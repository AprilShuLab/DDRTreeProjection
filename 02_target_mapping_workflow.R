############################################################
## 02. General target-cohort mapping workflow
##
## DDRTree is trained with adjusted reference phenotypes, but
## KNN mapping intentionally uses standardized RAW phenotypes.
## The reference-neighbor Y nodes provide tree-based subtype
## ratios and weights for projecting labels, Z coordinates,
## and K continuous subtype scores.
############################################################

############################################################
## 1. User configuration
############################################################

input_file <- file.path(
  "output",
  "Generic_phenotype_input.RData"
)

reference_model_file <- file.path(
  "output",
  "Generic_reference_DDRTree_model.RData"
)

mapping_output_file <- file.path(
  "output",
  "Generic_target_mapping_result.RData"
)

plot_dir <- file.path("output", "mapping_plots")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
## The number of neighbors.
mapping_k <- 3

## Optional stratified-bootstrap stability analysis.
# RUN_BOOTSTRAP <- FALSE
RUN_BOOTSTRAP <- TRUE
bootstrap_B <- 5
bootstrap_seed <- 123

############################################################
## 2. Packages and input
############################################################

library(FNN)
library(ggplot2)
library(igraph)

load(input_file)
load(reference_model_file)

id_col <- input_config$id_col
phenotype_cols <- input_config$phenotype_cols

subtype_levels <- reference_model$subtype_levels
K <- length(subtype_levels)

mapping_df_Reference <- reference_model$mapping_df_Reference
reference_subtype_scores <- reference_model$subtype_scores
root_to_Y_distance <- reference_model$root_to_Y_distance
subtype_scale <- reference_model$subtype_scale
subtype_colors <- reference_model$subtype_colors

reference_ids <- as.character(reference_data[[id_col]])
target_ids <- as.character(target_data[[id_col]])

if (mapping_k < 1 || mapping_k > nrow(reference_data)) {
  stop("mapping_k must be between 1 and the reference sample size.")
}

if (!all(reference_ids %in% rownames(reference_subtype_scores))) {
  stop("Reference subtype-score row names do not match reference IDs.")
}

############################################################
## 3. Standardize raw phenotype matrices using reference
##    center and scale
############################################################

reference_raw_matrix <- as.matrix(
  reference_data[, phenotype_cols, drop = FALSE]
)
rownames(reference_raw_matrix) <- reference_ids

target_raw_matrix <- as.matrix(
  target_data[, phenotype_cols, drop = FALSE]
)
rownames(target_raw_matrix) <- target_ids

reference_scaled <- scale(reference_raw_matrix)

raw_reference_center <- attr(
  reference_scaled,
  "scaled:center"
)

raw_reference_scale <- attr(
  reference_scaled,
  "scaled:scale"
)

target_scaled <- scale(
  target_raw_matrix,
  center = raw_reference_center,
  scale = raw_reference_scale
)

############################################################
## 4. General helper functions
############################################################

normalize_or_fallback <- function(
    branch_ratio,
    neighbor_subtypes,
    subtype_levels
) {
  if (
    all(is.finite(branch_ratio)) &&
    sum(branch_ratio) > 0
  ) {
    out <- branch_ratio / sum(branch_ratio)
  } else {
    label_count <- table(
      factor(
        neighbor_subtypes,
        levels = subtype_levels
      )
    )
    
    out <- as.numeric(label_count)
    names(out) <- subtype_levels
    
    if (sum(out) == 0) {
      stop("No valid subtype labels were found among the reference neighbors.")
    }
    
    out <- out / sum(out)
  }
  
  out[subtype_levels]
}


compute_branch_ratios <- function(
    neighbor_Y_indices,
    neighbor_subtypes,
    subtype_levels,
    root_to_Y_distance,
    subtype_scale
) {
  neighbor_Y_table <- as.data.frame(
    table(neighbor_Y_indices),
    stringsAsFactors = FALSE
  )
  
  colnames(neighbor_Y_table) <- c("Y_index", "count")
  neighbor_Y_table$Y_index <- as.integer(
    as.character(neighbor_Y_table$Y_index)
  )
  
  ## Every Y node has exactly one subtype in the reference model.
  Y_to_subtype <- setNames(
    reference_model$Y_label_table$Subtype,
    reference_model$Y_label_table$Y_index
  )
  
  neighbor_Y_table$Subtype <- Y_to_subtype[
    as.character(neighbor_Y_table$Y_index)
  ]
  
  neighbor_Y_table$root_distance <- root_to_Y_distance[
    paste0("Y_", neighbor_Y_table$Y_index)
  ]
  
  branch_ratio <- setNames(
    numeric(length(subtype_levels)),
    subtype_levels
  )
  
  for (g in subtype_levels) {
    df_g <- neighbor_Y_table[
      neighbor_Y_table$Subtype == g,
      ,
      drop = FALSE
    ]
    
    if (nrow(df_g) == 0) {
      branch_ratio[g] <- 0
    } else {
      weighted_graph_distance <- weighted.mean(
        x = df_g$root_distance,
        w = df_g$count
      )
      
      branch_ratio[g] <-
        weighted_graph_distance /
        subtype_scale[g]
    }
  }
  
  normalized_ratio <- normalize_or_fallback(
    branch_ratio = branch_ratio,
    neighbor_subtypes = neighbor_subtypes,
    subtype_levels = subtype_levels
  )
  
  list(
    raw_ratio = branch_ratio,
    normalized_ratio = normalized_ratio,
    neighbor_Y_table = neighbor_Y_table
  )
}


map_target_cohort <- function(
    reference_scaled,
    target_scaled,
    reference_Z,
    reference_ids,
    target_ids,
    mapping_df_Reference,
    reference_subtype_scores,
    subtype_levels,
    root_to_Y_distance,
    subtype_scale,
    mapping_k,
    show_progress = TRUE
) {
  n_target <- nrow(target_scaled)
  K_local <- length(subtype_levels)
  z_dimensions <- ncol(reference_Z)
  
  if (nrow(reference_Z) != nrow(reference_scaled)) {
    stop("reference_Z and reference_scaled must have the same row order and size.")
  }
  
  if (nrow(reference_subtype_scores) != nrow(reference_scaled)) {
    stop(
      "reference_subtype_scores and reference_scaled must have the same ",
      "row order and size."
    )
  }
  
  knn_result <- FNN::get.knnx(
    data = reference_scaled,
    query = target_scaled,
    k = mapping_k
  )
  
  target_label <- character(n_target)
  
  target_branch_ratios <- matrix(
    NA_real_,
    nrow = n_target,
    ncol = K_local,
    dimnames = list(target_ids, subtype_levels)
  )
  
  target_raw_branch_ratios <- matrix(
    NA_real_,
    nrow = n_target,
    ncol = K_local,
    dimnames = list(target_ids, subtype_levels)
  )
  
  target_subtype_scores <- matrix(
    NA_real_,
    nrow = n_target,
    ncol = K_local,
    dimnames = list(target_ids, subtype_levels)
  )
  
  target_Z <- matrix(
    NA_real_,
    nrow = n_target,
    ncol = z_dimensions,
    dimnames = list(
      target_ids,
      paste0("Dim", seq_len(z_dimensions))
    )
  )
  
  target_Z_covariance <- vector("list", n_target)
  names(target_Z_covariance) <- target_ids
  
  reference_Y_lookup <- mapping_df_Reference[
    match(reference_ids, mapping_df_Reference$individual),
    ,
    drop = FALSE
  ]
  
  if (anyNA(reference_Y_lookup$Y_index)) {
    stop("Some reference IDs have no Z-to-Y mapping.")
  }
  
  for (i in seq_len(n_target)) {
    if (show_progress && i %% 100 == 0) {
      message(
        "Progress: ", i, "/", n_target,
        " (", round(100 * i / n_target, 1), "%)"
      )
    }
    
    neighbor_indices <- knn_result$nn.index[i, ]
    neighbor_ids <- reference_ids[neighbor_indices]
    
    neighbor_Y_indices <- reference_Y_lookup$Y_index[
      neighbor_indices
    ]
    
    neighbor_subtypes <- reference_Y_lookup$Subtype[
      neighbor_indices
    ]
    
    ratio_result <- compute_branch_ratios(
      neighbor_Y_indices = neighbor_Y_indices,
      neighbor_subtypes = neighbor_subtypes,
      subtype_levels = subtype_levels,
      root_to_Y_distance = root_to_Y_distance,
      subtype_scale = subtype_scale
    )
    
    normalized_ratio <- ratio_result$normalized_ratio
    
    target_raw_branch_ratios[i, ] <- ratio_result$raw_ratio
    target_branch_ratios[i, ] <- normalized_ratio
    
    target_label[i] <- names(normalized_ratio)[
      which.max(normalized_ratio)
    ]
    
    ## Preserve the original mapping design: tree-ratio label
    ## weights are used without inverse raw-KNN distance weights.
    label_weight <- normalized_ratio[neighbor_subtypes]
    
    if (
      any(!is.finite(label_weight)) ||
      sum(label_weight) <= 0
    ) {
      label_weight <- rep(1, length(neighbor_indices))
    }
    
    knn_weight <- as.numeric(
      label_weight / sum(label_weight)
    )
    
    Z_neighbors <- reference_Z[
      neighbor_indices,
      ,
      drop = FALSE
    ]
    
    mapped_Z <- as.numeric(
      knn_weight %*% Z_neighbors
    )
    
    target_Z[i, ] <- mapped_Z
    
    if (nrow(Z_neighbors) == 1L) {
      target_Z_covariance[[i]] <- matrix(
        NA_real_,
        nrow = z_dimensions,
        ncol = z_dimensions
      )
    } else {
      Z_centered <- sweep(
        Z_neighbors,
        2,
        mapped_Z,
        "-"
      )
      
      target_Z_covariance[[i]] <-
        t(Z_centered) %*%
        diag(knn_weight) %*%
        Z_centered
    }
    
    neighbor_scores <- reference_subtype_scores[
      neighbor_indices,
      subtype_levels,
      drop = FALSE
    ]
    
    target_subtype_scores[i, ] <- as.numeric(
      knn_weight %*% neighbor_scores
    )
  }
  
  mapping_confidence <- apply(
    target_branch_ratios,
    1,
    max,
    na.rm = TRUE
  )
  
  list(
    target_label = factor(
      target_label,
      levels = subtype_levels
    ),
    mapping_confidence = mapping_confidence,
    target_raw_branch_ratios = target_raw_branch_ratios,
    target_branch_ratios = target_branch_ratios,
    target_subtype_scores = target_subtype_scores,
    target_Z = target_Z,
    target_Z_covariance = target_Z_covariance,
    knn_index = knn_result$nn.index,
    knn_distance = knn_result$nn.dist
  )
}

############################################################
## 5. Run target mapping
############################################################

mapping_result <- map_target_cohort(
  reference_scaled = reference_scaled,
  target_scaled = target_scaled,
  reference_Z = t(reference_model$ddr_results$Z),
  reference_ids = reference_ids,
  target_ids = target_ids,
  mapping_df_Reference = mapping_df_Reference,
  reference_subtype_scores = reference_subtype_scores,
  subtype_levels = subtype_levels,
  root_to_Y_distance = root_to_Y_distance,
  subtype_scale = subtype_scale,
  mapping_k = mapping_k
)

target_mapping_result <- data.frame(
  Individual = target_ids,
  Predicted_subtype = mapping_result$target_label,
  Mapping_confidence = mapping_result$mapping_confidence,
  mapping_result$target_Z,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ratio_output <- as.data.frame(
  mapping_result$target_branch_ratios,
  check.names = FALSE
)

colnames(ratio_output) <- paste0(
  colnames(ratio_output),
  "_ratio"
)

score_output <- as.data.frame(
  mapping_result$target_subtype_scores,
  check.names = FALSE
)

colnames(score_output) <- paste0(
  colnames(score_output),
  "_score"
)

target_mapping_result <- cbind(
  target_mapping_result,
  ratio_output,
  score_output
)

cat("Target subtype counts:\n")
print(table(target_mapping_result$Predicted_subtype))

############################################################
## 6.1 Plot mapped target coordinates and reference Y MST
############################################################

Y <- reference_model$ddr_results$Y
Y_edge_table <- reference_model$Y_edge_table

Y_edge_plot_data <- data.frame(
  x = Y[1, Y_edge_table$from],
  y = Y[2, Y_edge_table$from],
  xend = Y[1, Y_edge_table$to],
  yend = Y[2, Y_edge_table$to]
)

target_plot_data <- data.frame(
  Dim1 = mapping_result$target_Z[, 1],
  Dim2 = mapping_result$target_Z[, 2],
  Subtype = mapping_result$target_label
)

p_target <- ggplot() +
  geom_point(
    data = target_plot_data,
    aes(x = Dim1, y = Dim2, color = Subtype),
    size = 1.2,
    alpha = 0.65
  ) +
  geom_segment(
    data = Y_edge_plot_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black",
    linewidth = 0.8
  ) +
  geom_point(
    data = reference_model$Y_node_table,
    aes(x = Dim1, y = Dim2),
    color = "black",
    size = 1.5
  ) +
  scale_color_manual(values = subtype_colors) +
  coord_equal() +
  theme_classic() +
  theme(legend.position = "top") +
  labs(
    title = "Target cohort mapped to the reference DDRTree",
    x = "Dimension 1",
    y = "Dimension 2"
  )

ggsave(
  filename = file.path(plot_dir, "target_mapping_DDRTree.pdf"),
  plot = p_target,
  width = 9,
  height = 7
)

############################################################
## 6.2 Plot each projected target subtype score
############################################################

target_score_plot_list <- vector(
  "list",
  length(subtype_levels)
)

names(target_score_plot_list) <- subtype_levels

for (g in subtype_levels) {
  target_score_plot_data <- data.frame(
    Dim1 = mapping_result$target_Z[, 1],
    Dim2 = mapping_result$target_Z[, 2],
    Score = mapping_result$target_subtype_scores[, g],
    stringsAsFactors = FALSE
  )
  
  p_target_score <- ggplot(
    target_score_plot_data,
    aes(x = Dim1, y = Dim2, color = Score)
  ) +
    geom_point(size = 0.9, alpha = 0.65) +
    scale_color_viridis_c(
      option = "plasma",
      direction = -1,
      limits = c(0, 1)
    ) +
    coord_equal() +
    theme_classic() +
    theme(legend.position = "top") +
    labs(
      title = paste("Target subtype score:", g),
      x = "Dimension 1",
      y = "Dimension 2",
      color = "Score"
    )
  
  target_score_plot_list[[g]] <- p_target_score
  
  ggsave(
    filename = file.path(
      plot_dir,
      paste0("target_score_", make.names(g), ".pdf")
    ),
    plot = p_target_score,
    width = 8,
    height = 6
  )
}

############################################################
## 7. Optional general K-subtype bootstrap stability
############################################################

sample_idx_stratified_bootstrap <- function(y, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  y <- factor(y)
  index_by_subtype <- split(seq_along(y), y)
  
  as.integer(
    unlist(
      lapply(
        index_by_subtype,
        function(ii) sample(ii, size = length(ii), replace = TRUE)
      )
    )
  )
}


summarize_bootstrap_labels <- function(
    prediction_matrix,
    subtype_levels
) {
  assignment_counts <- sapply(
    subtype_levels,
    function(g) {
      rowSums(
        prediction_matrix == g,
        na.rm = TRUE
      )
    }
  )
  
  colnames(assignment_counts) <- paste0(
    "n_",
    make.names(subtype_levels)
  )
  
  max_count <- apply(assignment_counts, 1, max)
  final_index <- max.col(assignment_counts, ties.method = "first")
  
  data.frame(
    Individual = rownames(prediction_matrix),
    assignment_counts,
    Consistency = max_count / ncol(prediction_matrix),
    Final_subtype = subtype_levels[final_index],
    Combination = apply(
      assignment_counts,
      1,
      paste,
      collapse = "-"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


bootstrap_stability <- NULL

if (RUN_BOOTSTRAP) {
  prediction_matrix <- matrix(
    NA_character_,
    nrow = length(target_ids),
    ncol = bootstrap_B,
    dimnames = list(
      target_ids,
      paste0("Bootstrap_", seq_len(bootstrap_B))
    )
  )
  
  score_arrays <- setNames(
    lapply(
      subtype_levels,
      function(g) {
        matrix(
          NA_real_,
          nrow = length(target_ids),
          ncol = bootstrap_B,
          dimnames = list(
            target_ids,
            paste0("Bootstrap_", seq_len(bootstrap_B))
          )
        )
      }
    ),
    subtype_levels
  )
  
  for (b in seq_len(bootstrap_B)) {
    message("Bootstrap ", b, "/", bootstrap_B)
    
    boot_index <- sample_idx_stratified_bootstrap(
      y = reference_model$label_z,
      seed = bootstrap_seed + b
    )
    
    boot_reference_scaled <- reference_scaled[
      boot_index,
      ,
      drop = FALSE
    ]
    
    boot_reference_ids <- reference_ids[boot_index]
    
    boot_mapping_df <- mapping_df_Reference[
      match(boot_reference_ids, mapping_df_Reference$individual),
      ,
      drop = FALSE
    ]
    
    ## Unique row names are unnecessary here because lookup is by
    ## row position and individual ID; duplicated sampled IDs retain
    ## their original reference tree information.
    boot_scores <- reference_subtype_scores[
      boot_reference_ids,
      subtype_levels,
      drop = FALSE
    ]
    
    boot_result <- map_target_cohort(
      reference_scaled = boot_reference_scaled,
      target_scaled = target_scaled,
      reference_Z = t(reference_model$ddr_results$Z)[
        boot_index,
        ,
        drop = FALSE
      ],
      reference_ids = boot_reference_ids,
      target_ids = target_ids,
      mapping_df_Reference = boot_mapping_df,
      reference_subtype_scores = boot_scores,
      subtype_levels = subtype_levels,
      root_to_Y_distance = root_to_Y_distance,
      subtype_scale = subtype_scale,
      mapping_k = mapping_k,
      show_progress = FALSE
    )
    
    prediction_matrix[, b] <- as.character(
      boot_result$target_label
    )
    
    for (g in subtype_levels) {
      score_arrays[[g]][, b] <-
        boot_result$target_subtype_scores[, g]
    }
  }
  
  label_stability <- summarize_bootstrap_labels(
    prediction_matrix = prediction_matrix,
    subtype_levels = subtype_levels
  )
  
  score_stability <- do.call(
    cbind,
    lapply(
      subtype_levels,
      function(g) {
        out <- data.frame(
          rowMeans(score_arrays[[g]], na.rm = TRUE),
          apply(score_arrays[[g]], 1, sd, na.rm = TRUE)
        )
        
        colnames(out) <- c(
          paste0(g, "_mean"),
          paste0(g, "_sd")
        )
        
        out
      }
    )
  )
  
  rownames(score_stability) <- target_ids
  
  ##########################################################
  ## Mean individual-level bootstrap SD for each subtype
  ## score. This summarizes score variation across the B
  ## stratified bootstrap mappings.
  ##########################################################
  
  score_sd_summary <- data.frame(
    Subtype = subtype_levels,
    N_individuals = vapply(
      subtype_levels,
      function(g) {
        individual_sd <- apply(
          score_arrays[[g]],
          1,
          sd,
          na.rm = TRUE
        )
        
        sum(is.finite(individual_sd))
      },
      integer(1)
    ),
    Mean_bootstrap_SD = vapply(
      subtype_levels,
      function(g) {
        individual_sd <- apply(
          score_arrays[[g]],
          1,
          sd,
          na.rm = TRUE
        )
        
        mean(
          individual_sd[is.finite(individual_sd)],
          na.rm = TRUE
        )
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  combination_distribution <- as.data.frame(
    table(label_stability$Combination),
    stringsAsFactors = FALSE
  )
  
  colnames(combination_distribution) <- c(
    "Combination",
    "Count"
  )
  
  combination_distribution$Proportion <-
    combination_distribution$Count /
    nrow(label_stability)
  
  combination_distribution <- combination_distribution[
    order(-combination_distribution$Count),
    ,
    drop = FALSE
  ]
  
  bootstrap_stability <- list(
    prediction_matrix = prediction_matrix,
    score_arrays = score_arrays,
    label_stability = label_stability,
    score_stability = score_stability,
    score_sd_summary = score_sd_summary,
    combination_distribution = combination_distribution
  )
}

############################################################
## 8. Save mapping results
############################################################

save(
  target_mapping_result,
  mapping_result,
  bootstrap_stability,
  mapping_k,
  raw_reference_center,
  raw_reference_scale,
  file = mapping_output_file
)

cat("Saved:", mapping_output_file, "\n")
