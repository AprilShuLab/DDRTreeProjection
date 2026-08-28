############################################################
## 01. General reference DDRTree workflow
##
## Stage 1: fit DDRTree, construct the Y-node MST, establish
##          reference-sample Z-to-Y mapping, and inspect plots.
## Stage 2: define a complete Y-node subtype partition and
##          calculate reference labels and K subtype scores.
############################################################

############################################################
## 1. User configuration
############################################################

input_file <- file.path(
  "output",
  "Generic_phenotype_input.RData"
)

preliminary_model_file <- file.path(
  "output",
  "Generic_reference_DDRTree_preliminary.RData"
)

final_model_file <- file.path(
  "output",
  "Generic_reference_DDRTree_model.RData"
)

plot_dir <- file.path("output", "reference_plots")

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

## First run: keep FALSE. Inspect the numbered Y-node plot.
## Second run: enter root_node and subtype_Y_nodes below, then
## set DEFINE_SUBTYPES to TRUE.
DEFINE_SUBTYPES <- FALSE
# DEFINE_SUBTYPES <- TRUE

## This is only a template. Every Y node must appear exactly
## once across the lists when DEFINE_SUBTYPES is TRUE.
root_node <- NULL

subtype_Y_nodes <- NULL

## Example syntax only:
# root_node <- 11
# subtype_Y_nodes <- list(
#   Subtype_1 = c(26, 31, 9, 39, 33, 22, 17, 13, 46, 14, 45, 10, 24, 32, 40, 15, 12, 6, 35),
#   Subtype_2 = c(29, 43, 5, 1, 50, 44, 28, 23, 25, 42, 37, 38, 4, 36, 47, 2, 16, 20),
#   Subtype_3 = c(11, 7, 19, 21, 49,48, 18, 30, 3, 27, 34, 41, 8)
# )

## Path to the igraph-compatible orderCells source used in the
## original workflow. Set to NULL if no custom source is needed.
order_cells_source <- "order_cells_igraph214.R"

ddrtree_config <- list(
  dimensions = 2,
  maxIter = 30,
  sigma = 1e-3,
  lambda = NULL,
  ncenter = 50,
  gamma = 2,
  tol = 1e-3,
  seed = 123
)

############################################################
## 2. Packages and input
############################################################

library(DDRTree)
library(monocle)
library(igraph)
library(Matrix)
library(proxy)
library(ggplot2)

if (!is.null(order_cells_source)) {
  if (!file.exists(order_cells_source)) {
    stop("order_cells_source does not exist: ", order_cells_source)
  }
  source(order_cells_source)
}

load(input_file)

id_col <- input_config$id_col
phenotype_cols <- input_config$phenotype_cols
adjustment_formula <- input_config$adjustment_formula

reference_ids <- as.character(reference_data[[id_col]])

############################################################
## 3. Adjust reference phenotypes and fit DDRTree
############################################################

adjustment_result <- adjust_phenotypes(
  data = reference_data,
  id_col = id_col,
  phenotype_cols = phenotype_cols,
  adjustment_formula = adjustment_formula
)

reference_adjusted_matrix <- adjustment_result$adjusted_matrix
adjustment_models <- adjustment_result$fitted_models

X <- t(reference_adjusted_matrix)
colnames(X) <- reference_ids

set.seed(ddrtree_config$seed)

ddr_results <- DDRTree(
  X = X,
  dimensions = ddrtree_config$dimensions,
  maxIter = ddrtree_config$maxIter,
  sigma = ddrtree_config$sigma,
  lambda = ddrtree_config$lambda,
  ncenter = ddrtree_config$ncenter,
  param.gamma = ddrtree_config$gamma,
  tol = ddrtree_config$tol,
  verbose = TRUE
)

Z <- ddr_results$Z
Y <- ddr_results$Y
W <- ddr_results$W

colnames(Z) <- reference_ids

n_Y_nodes <- ncol(Y)
Y_names <- paste0("Y_", seq_len(n_Y_nodes))
colnames(Y) <- Y_names

############################################################
## 4. Helper: reconstruct the Monocle DDRTree object and
##    establish each reference sample's nearest Y node
############################################################

findNearestPointOnMST <- function(cds) {
  dp_mst <- minSpanningTree(cds)
  Z_local <- reducedDimS(cds)
  Y_local <- reducedDimK(cds)

  distances_Z_to_Y <- proxy::dist(
    t(Z_local),
    t(Y_local)
  )

  closest_vertex <- apply(
    distances_Z_to_Y,
    1,
    function(z) which(z == min(z))[1]
  )

  closest_vertex_df <- as.matrix(closest_vertex)
  rownames(closest_vertex_df) <- names(closest_vertex)

  cds@auxOrderingData[["DDRTree"]][[
    "pr_graph_cell_proj_closest_vertex"
  ]] <- closest_vertex_df

  cds
}


build_reference_cds <- function(
    X,
    reference_ids,
    ddr_results,
    Y_names
) {
  sample_metadata <- data.frame(
    row.names = reference_ids
  )

  feature_metadata <- data.frame(
    gene_short_name = rownames(X),
    row.names = rownames(X)
  )

  cds <- newCellDataSet(
    X,
    phenoData = new(
      "AnnotatedDataFrame",
      data = sample_metadata
    ),
    featureData = new(
      "AnnotatedDataFrame",
      data = feature_metadata
    ),
    lowerDetectionLimit = 0,
    expressionFamily = uninormal()
  )

  cds@dim_reduce_type <- "DDRTree"
  cds@reducedDimS <- ddr_results$Z
  colnames(cds@reducedDimS) <- reference_ids
  cds@reducedDimW <- ddr_results$W
  cds@reducedDimA <- matrix(nrow = 0, ncol = 0)
  cds@reducedDimK <- ddr_results$Y
  colnames(cds@reducedDimK) <- Y_names

  cds@auxOrderingData[["DDRTree"]]$objective_vals <-
    ddr_results$objective_vals

  adjusted_K <- Matrix::t(reducedDimK(cds))
  dp <- as.matrix(dist(adjusted_K))
  cellPairwiseDistances(cds) <- dp

  gp <- igraph::graph.adjacency(
    dp,
    mode = "undirected",
    weighted = TRUE
  )

  dp_mst <- igraph::minimum.spanning.tree(gp)
  minSpanningTree(cds) <- dp_mst
  cds@dim_reduce_type <- "DDRTree"
  cds <- findNearestPointOnMST(cds)

  cds
}


cds <- build_reference_cds(
  X = X,
  reference_ids = reference_ids,
  ddr_results = ddr_results,
  Y_names = Y_names
)

closest_vertex_df <- cds@auxOrderingData[["DDRTree"]][[
  "pr_graph_cell_proj_closest_vertex"
]]

ind_Y_index <- as.integer(closest_vertex_df[, 1])

mapping_df_Reference <- data.frame(
  individual = rownames(closest_vertex_df),
  Y_index = ind_Y_index,
  Y_name = Y_names[ind_Y_index],
  Y_dim_1 = Y[1, ind_Y_index],
  Y_dim_2 = Y[2, ind_Y_index],
  row.names = rownames(closest_vertex_df),
  stringsAsFactors = FALSE
)

############################################################
## 5. Extract the Y-node MST and graph-distance information
############################################################

Y_mst <- minSpanningTree(cds)

## Ensure vertex names are exactly Y_1, ..., Y_M.
if (is.null(igraph::V(Y_mst)$name)) {
  igraph::V(Y_mst)$name <- Y_names
}

Y_mst_vertex_order <- match(
  igraph::V(Y_mst)$name,
  Y_names
)

if (anyNA(Y_mst_vertex_order)) {
  stop("The reconstructed MST vertex names do not match the Y nodes.")
}

edge_ends <- igraph::ends(
  Y_mst,
  es = igraph::E(Y_mst),
  names = TRUE
)

Y_edge_table <- data.frame(
  from_name = edge_ends[, 1],
  to_name = edge_ends[, 2],
  stringsAsFactors = FALSE
)

Y_edge_table$from <- match(Y_edge_table$from_name, Y_names)
Y_edge_table$to <- match(Y_edge_table$to_name, Y_names)

Y_edge_table$edge_length <- sqrt(
  colSums(
    (
      Y[, Y_edge_table$from, drop = FALSE] -
        Y[, Y_edge_table$to, drop = FALSE]
    )^2
  )
)

igraph::E(Y_mst)$weight <- Y_edge_table$edge_length

Y_node_counts <- table(
  factor(
    mapping_df_Reference$Y_index,
    levels = seq_len(n_Y_nodes)
  )
)

Y_node_table <- data.frame(
  Y_index = seq_len(n_Y_nodes),
  Y_name = Y_names,
  Dim1 = Y[1, ],
  Dim2 = Y[2, ],
  n_reference = as.integer(Y_node_counts),
  stringsAsFactors = FALSE
)

############################################################
## 6. Reference diagnostic plots
############################################################

Y_edge_plot_data <- data.frame(
  x = Y[1, Y_edge_table$from],
  y = Y[2, Y_edge_table$from],
  xend = Y[1, Y_edge_table$to],
  yend = Y[2, Y_edge_table$to]
)

p_Y_numbered <- ggplot() +
  geom_segment(
    data = Y_edge_plot_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50",
    linewidth = 0.7
  ) +
  geom_point(
    data = Y_node_table,
    aes(x = Dim1, y = Dim2),
    color = "black",
    size = 2
  ) +
  geom_text(
    data = Y_node_table,
    aes(x = Dim1, y = Dim2, label = Y_index),
    nudge_y = 0.15,
    size = 3
  ) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "Numbered DDRTree Y-node MST",
    x = "Dimension 1",
    y = "Dimension 2"
  )

ggsave(
  filename = file.path(plot_dir, "01_numbered_Y_node_MST.pdf"),
  plot = p_Y_numbered,
  width = 9,
  height = 7
)

plot_data_Z <- data.frame(
  individual = reference_ids,
  Dim1 = Z[1, ],
  Dim2 = Z[2, ],
  stringsAsFactors = FALSE
)

p_Z_Y <- ggplot() +
  geom_point(
    data = plot_data_Z,
    aes(x = Dim1, y = Dim2),
    color = "#F8766D",
    alpha = 0.35,
    size = 0.8
  ) +
  geom_segment(
    data = Y_edge_plot_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "black",
    linewidth = 0.8
  ) +
  geom_point(
    data = Y_node_table,
    aes(x = Dim1, y = Dim2),
    color = "black",
    size = 1.8
  ) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "Reference Z coordinates and Y-node MST",
    x = "Dimension 1",
    y = "Dimension 2"
  )

ggsave(
  filename = file.path(plot_dir, "02_reference_Z_and_Y_MST.pdf"),
  plot = p_Z_Y,
  width = 9,
  height = 7
)

for (pheno in phenotype_cols) {
  pheno_plot_data <- plot_data_Z
  pheno_plot_data$phenotype_value <- reference_data[[pheno]]

  p_pheno <- ggplot(
    pheno_plot_data,
    aes(x = Dim1, y = Dim2, color = phenotype_value)
  ) +
    geom_point(size = 0.9, alpha = 0.65) +
    scale_color_viridis_c() +
    coord_equal() +
    theme_classic() +
    theme(legend.position = "top") +
    labs(
      title = paste("DDRTree gradient:", pheno),
      x = "Dimension 1",
      y = "Dimension 2",
      color = pheno
    )

  ggsave(
    filename = file.path(
      plot_dir,
      paste0("phenotype_gradient_", pheno, ".pdf")
    ),
    plot = p_pheno,
    width = 8,
    height = 6
  )
}

############################################################
## 7. Save preliminary model for first-stage inspection
############################################################

reference_preliminary_model <- list(
  input_config = input_config,
  ddrtree_config = ddrtree_config,
  adjustment_models = adjustment_models,
  reference_adjusted_matrix = reference_adjusted_matrix,
  ddr_results = ddr_results,
  cds = cds,
  Y_mst = Y_mst,
  Y_node_table = Y_node_table,
  Y_edge_table = Y_edge_table,
  mapping_df_Reference = mapping_df_Reference
)

save(
  reference_preliminary_model,
  file = preliminary_model_file
)

if (!DEFINE_SUBTYPES) {
  message(
    "Stage 1 completed. Inspect: ",
    file.path(plot_dir, "01_numbered_Y_node_MST.pdf"),
    "\nThen define root_node and subtype_Y_nodes, set ",
    "DEFINE_SUBTYPES <- TRUE, and rerun this script."
  )
} else {

  ############################################################
  ## 8. Validate the user-defined complete Y-node partition
  ############################################################

  validate_subtype_partition <- function(
      subtype_Y_nodes,
      root_node,
      n_Y_nodes,
      Y_mst,
      Y_names
  ) {
    if (is.null(root_node) || length(root_node) != 1) {
      stop("Define exactly one root_node.")
    }

    if (is.null(subtype_Y_nodes) || length(subtype_Y_nodes) < 2) {
      stop("Define subtype_Y_nodes as a named list with at least two subtypes.")
    }

    if (is.null(names(subtype_Y_nodes)) || any(names(subtype_Y_nodes) == "")) {
      stop("Every element of subtype_Y_nodes must have a subtype name.")
    }

    if (anyDuplicated(names(subtype_Y_nodes)) > 0) {
      stop("Subtype names must be unique.")
    }

    all_nodes <- as.integer(
      unlist(subtype_Y_nodes, use.names = FALSE)
    )

    valid_nodes <- seq_len(n_Y_nodes)
    duplicated_nodes <- unique(all_nodes[duplicated(all_nodes)])
    missing_nodes <- setdiff(valid_nodes, all_nodes)
    invalid_nodes <- setdiff(all_nodes, valid_nodes)

    if (length(duplicated_nodes) > 0) {
      stop(
        "Y nodes assigned to more than one subtype: ",
        paste(duplicated_nodes, collapse = ", ")
      )
    }

    if (length(missing_nodes) > 0) {
      stop(
        "Y nodes without a subtype assignment: ",
        paste(missing_nodes, collapse = ", ")
      )
    }

    if (length(invalid_nodes) > 0) {
      stop(
        "Invalid Y-node indices: ",
        paste(invalid_nodes, collapse = ", ")
      )
    }

    if (!(root_node %in% all_nodes)) {
      stop("root_node must belong to exactly one subtype.")
    }

    ## Connectivity is informative rather than mandatory.
    for (g in names(subtype_Y_nodes)) {
      node_names_g <- Y_names[subtype_Y_nodes[[g]]]
      graph_g <- igraph::induced_subgraph(Y_mst, vids = node_names_g)

      if (igraph::vcount(graph_g) > 1 && !igraph::is_connected(graph_g)) {
        warning(
          "Subtype '", g,
          "' contains disconnected Y-node regions."
        )
      }
    }

    invisible(TRUE)
  }

  validate_subtype_partition(
    subtype_Y_nodes = subtype_Y_nodes,
    root_node = root_node,
    n_Y_nodes = n_Y_nodes,
    Y_mst = Y_mst,
    Y_names = Y_names
  )

  subtype_levels <- names(subtype_Y_nodes)
  K <- length(subtype_levels)

  Y_label_table <- Y_node_table
  Y_label_table$Subtype <- NA_character_

  for (g in subtype_levels) {
    Y_label_table$Subtype[
      Y_label_table$Y_index %in% subtype_Y_nodes[[g]]
    ] <- g
  }

  mapping_df_Reference$Subtype <- Y_label_table$Subtype[
    match(mapping_df_Reference$Y_index, Y_label_table$Y_index)
  ]

  label_z <- factor(
    mapping_df_Reference$Subtype,
    levels = subtype_levels
  )
  names(label_z) <- mapping_df_Reference$individual

  ############################################################
  ## 9. Root-to-Y MST graph distances and subtype scales
  ############################################################

  root_name <- Y_names[root_node]

  root_distance_matrix <- igraph::distances(
    Y_mst,
    v = root_name,
    to = Y_names,
    weights = igraph::E(Y_mst)$weight
  )

  root_to_Y_distance <- as.numeric(root_distance_matrix[1, Y_names])
  names(root_to_Y_distance) <- Y_names

  subtype_terminal_nodes <- vapply(
    subtype_Y_nodes,
    function(nodes) {
      nodes <- as.integer(nodes)
      
      node_names <- Y_names[nodes]
      
      as.integer(
        nodes[
          which.max(
            root_to_Y_distance[node_names]
          )
        ]
      )
    },
    integer(1)
  )

  subtype_scale <- vapply(
    subtype_levels,
    function(g) {
      terminal_name <- Y_names[subtype_terminal_nodes[g]]
      root_to_Y_distance[terminal_name]
    },
    numeric(1)
  )

  if (any(!is.finite(subtype_scale)) || any(subtype_scale <= 0)) {
    stop(
      "Every subtype must contain at least one Y node with positive ",
      "graph distance from the root."
    )
  }

  ############################################################
  ## 10. Dynamically select one representative reference
  ##     individual at the farthest Y node of each subtype
  ############################################################

  representative_individuals <- setNames(
    character(K),
    subtype_levels
  )

  for (g in subtype_levels) {
    terminal_node <- subtype_terminal_nodes[g]

    candidate_ids <- mapping_df_Reference$individual[
      mapping_df_Reference$Y_index == terminal_node
    ]

    if (length(candidate_ids) > 0) {
      representative_individuals[g] <- candidate_ids[1]
    } else {
      terminal_coord <- Y[, terminal_node]
      distance_to_terminal <- sqrt(
        colSums((Z - terminal_coord)^2)
      )
      representative_individuals[g] <-
        reference_ids[which.min(distance_to_terminal)]
    }
  }

  if (anyDuplicated(representative_individuals) > 0) {
    stop(
      "At least two subtypes selected the same representative individual. ",
      "Inspect subtype_terminal_nodes and the Y-node partition."
    )
  }

  ############################################################
  ## 11. Calculate K reference subtype scores with orderCells
  ############################################################

  reference_subtype_scores <- matrix(
    NA_real_,
    nrow = length(reference_ids),
    ncol = K,
    dimnames = list(reference_ids, subtype_levels)
  )

  pseudotime_plot_list <- vector("list", K)
  names(pseudotime_plot_list) <- subtype_levels

  for (g_index in seq_along(subtype_levels)) {
    g <- subtype_levels[g_index]

    cds_g <- cds

    ## Preserve the original workflow's design while making the
    ## number of subtype-specific root states dynamic.
    pData(cds_g)$State <- K + 1L

    for (h_index in seq_along(subtype_levels)) {
      h <- subtype_levels[h_index]
      representative_row <- match(
        representative_individuals[h],
        rownames(pData(cds_g))
      )
      pData(cds_g)$State[representative_row] <- h_index
    }

    cds_g <- orderCells(
      cds_g,
      root_state = g_index
    )

    pseudotime <- pData(cds_g)$Pseudotime
    pseudotime_range <- range(pseudotime, na.rm = TRUE)

    if (
      any(!is.finite(pseudotime_range)) ||
      diff(pseudotime_range) == 0
    ) {
      stop("Invalid pseudotime range for subtype: ", g)
    }

    subtype_score <- 1 -
      (pseudotime - pseudotime_range[1]) /
      diff(pseudotime_range)

    reference_subtype_scores[rownames(pData(cds_g)), g] <-
      subtype_score

    pseudotime_plot_data <- data.frame(
      Dim1 = Z[1, ],
      Dim2 = Z[2, ],
      Score = reference_subtype_scores[reference_ids, g]
    )

    p_score <- ggplot(
      pseudotime_plot_data,
      aes(x = Dim1, y = Dim2, color = Score)
    ) +
      geom_point(size = 0.9, alpha = 0.65) +
      scale_color_viridis_c(option = "plasma", direction = -1) +
      coord_equal() +
      theme_classic() +
      theme(legend.position = "top") +
      labs(
        title = paste("Reference subtype score:", g),
        x = "Dimension 1",
        y = "Dimension 2"
      )

    pseudotime_plot_list[[g]] <- p_score

    ggsave(
      filename = file.path(
        plot_dir,
        paste0("reference_score_", make.names(g), ".pdf")
      ),
      plot = p_score,
      width = 8,
      height = 6
    )
  }

  ############################################################
  ## 12. Plot the final user-defined Y-node partition
  ############################################################

  subtype_colors <- setNames(
    scales::hue_pal()(K),
    subtype_levels
  )

  p_Y_subtypes <- ggplot() +
    geom_segment(
      data = Y_edge_plot_data,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      color = "grey50",
      linewidth = 0.7
    ) +
    geom_point(
      data = Y_label_table,
      aes(
        x = Dim1,
        y = Dim2,
        color = Subtype
      ),
      size = 2
    ) +
    geom_text(
      data = Y_label_table,
      aes(
        x = Dim1,
        y = Dim2,
        label = Y_index
      ),
      nudge_y = 0.15,
      size = 3,
      show.legend = FALSE
    ) +
    geom_point(
      data = Y_label_table[
        Y_label_table$Y_index == root_node,
      ],
      aes(
        x = Dim1,
        y = Dim2
      ),
      shape = 21,
      fill = NA,
      color = "black",
      stroke = 1.4,
      size = 6,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = subtype_colors
    ) +
    coord_equal() +
    theme_classic() +
    theme(
      legend.position = "top"
    ) +
    labs(
      title = "User-defined Y-node subtype partition",
      x = "Dimension 1",
      y = "Dimension 2",
      color = "Subtype"
    )

  ggsave(
    filename = file.path(plot_dir, "03_final_Y_node_subtypes.pdf"),
    plot = p_Y_subtypes,
    width = 9,
    height = 7
  )

  ############################################################
  ## 13. Save final reference model
  ############################################################

  reference_model <- list(
    input_config = input_config,
    ddrtree_config = ddrtree_config,
    adjustment_models = adjustment_models,
    reference_adjusted_matrix = reference_adjusted_matrix,
    ddr_results = ddr_results,
    cds = cds,
    Y_mst = Y_mst,
    Y_node_table = Y_node_table,
    Y_edge_table = Y_edge_table,
    Y_label_table = Y_label_table,
    mapping_df_Reference = mapping_df_Reference,
    root_node = root_node,
    subtype_Y_nodes = subtype_Y_nodes,
    subtype_levels = subtype_levels,
    subtype_terminal_nodes = subtype_terminal_nodes,
    representative_individuals = representative_individuals,
    root_to_Y_distance = root_to_Y_distance,
    subtype_scale = subtype_scale,
    subtype_colors = subtype_colors,
    label_z = label_z,
    subtype_scores = reference_subtype_scores
  )

  save(
    reference_model,
    file = final_model_file
  )

  cat("Reference subtype counts:\n")
  print(table(label_z))

  cat("Subtype farthest Y nodes:\n")
  print(subtype_terminal_nodes)

  cat("Subtype graph-distance scales:\n")
  print(subtype_scale)

  cat("Saved:", final_model_file, "\n")
}
