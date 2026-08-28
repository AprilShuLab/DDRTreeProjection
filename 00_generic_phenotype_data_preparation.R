############################################################
## 00. Generic phenotype data preparation
##
## This toy example creates a reference cohort and a target
## cohort with generic phenotype names (Pheno1, Pheno2, ...).
## Replace Section 3 with real data import when appropriate.
############################################################

set.seed(2026)

############################################################
## 1. User configuration
############################################################

id_col <- "ID"

phenotype_cols <- c(
  "Pheno1",
  "Pheno2"
)

## One-sided formula used only to construct the reference
## DDRTree input. Raw phenotypes are retained for KNN mapping.
adjustment_formula <- ~ Age + I(Age^2) + Sex

## Raw phenotype lower bound. Adjusted residuals are not
## restricted because residuals can legitimately be negative.
phenotype_lower_bound <- 0

## Output directory. Change this path for a real project.
output_dir <- "output"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

############################################################
## 2. General helper functions
############################################################

rtrunc_normal <- function(n, mean, sd, lower, upper) {
  out <- numeric(0)
  
  while (length(out) < n) {
    candidate <- rnorm(n, mean = mean, sd = sd)
    candidate <- candidate[
      candidate >= lower & candidate <= upper
    ]
    out <- c(out, candidate)
  }
  
  out[seq_len(n)]
}


## The latent components below are used only to create a
## branching toy distribution. They are not returned as known
## or true subtypes.
simulate_generic_cohort <- function(
    n,
    phenotype_cols,
    cohort_name,
    component_prob,
    component_directions,
    phenotype_center,
    phenotype_noise_sd,
    age_mean = 9,
    age_sd = 3.5,
    male_probability = 0.75,
    age_effect = NULL,
    sex_effect = NULL,
    phenotype_lower_bound = 0,
    max_generation_attempts = 10000L
) {
  p <- length(phenotype_cols)
  k_component <- length(component_prob)
  
  stopifnot(
    nrow(component_directions) == k_component,
    ncol(component_directions) == p,
    length(phenotype_center) == p,
    length(phenotype_noise_sd) == p
  )
  
  if (is.null(age_effect)) age_effect <- rep(0, p)
  if (is.null(sex_effect)) sex_effect <- rep(0, p)
  
  stopifnot(
    length(age_effect) == p,
    length(sex_effect) == p
  )
  
  if (length(phenotype_lower_bound) == 1L) {
    phenotype_lower_bound <- rep(phenotype_lower_bound, p)
  }
  
  stopifnot(
    length(phenotype_lower_bound) == p,
    all(is.finite(phenotype_lower_bound))
  )
  
  Sex <- rbinom(n, size = 1, prob = male_probability)
  
  Age_year <- rtrunc_normal(
    n = n,
    mean = age_mean,
    sd = age_sd,
    lower = 2,
    upper = 17.99
  )
  
  ## Age is stored in months in this toy example.
  Age <- Age_year * 12
  
  latent_component <- sample(
    seq_len(k_component),
    size = n,
    replace = TRUE,
    prob = component_prob
  )
  
  phenotype_matrix <- matrix(
    NA_real_,
    nrow = n,
    ncol = p,
    dimnames = list(NULL, phenotype_cols)
  )
  
  age_centered <- Age_year - mean(Age_year)
  
  for (i in seq_len(n)) {
    direction_i <- component_directions[latent_component[i], ]
    
    ########################################################
    ## Rejection sampling
    ##
    ## Draw a new latent radius and phenotype noise until the
    ## complete phenotype vector satisfies its lower bounds.
    ## Values are therefore generated within the valid range;
    ## they are not clipped to the boundary afterward.
    ########################################################
    
    accepted <- FALSE
    attempt <- 0L
    
    while (!accepted) {
      attempt <- attempt + 1L
      
      latent_radius_i <- rgamma(
        1,
        shape = 2.5,
        scale = 1.1
      )
      
      candidate_phenotypes <-
        phenotype_center +
        latent_radius_i * direction_i +
        rnorm(p, mean = 0, sd = phenotype_noise_sd) +
        age_centered[i] * age_effect +
        Sex[i] * sex_effect
      
      accepted <- all(
        is.finite(candidate_phenotypes) &
          candidate_phenotypes >= phenotype_lower_bound
      )
      
      if (accepted) {
        phenotype_matrix[i, ] <- candidate_phenotypes
      }
      
      if (!accepted && attempt >= max_generation_attempts) {
        stop(
          "Unable to generate a valid phenotype vector for individual ",
          i,
          " after ",
          max_generation_attempts,
          " attempts. Check phenotype_center, component_directions, ",
          "phenotype_noise_sd, and phenotype_lower_bound."
        )
      }
    }
  }
  
  out <- data.frame(
    ID = paste0(
      cohort_name,
      "_",
      sprintf("%04d", seq_len(n))
    ),
    Sex = Sex,
    Age = Age,
    phenotype_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(out)[1] <- "ID"
  out
}


adjust_phenotypes <- function(
    data,
    id_col,
    phenotype_cols,
    adjustment_formula
) {
  required_cols <- unique(
    c(id_col, phenotype_cols, all.vars(adjustment_formula))
  )
  
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns required for phenotype adjustment: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (any(!complete.cases(data[, required_cols, drop = FALSE]))) {
    stop(
      "Missing values were found in the ID, phenotype, or adjustment columns. ",
      "Impute or otherwise handle missing values before running this workflow."
    )
  }
  
  adjusted_matrix <- matrix(
    NA_real_,
    nrow = nrow(data),
    ncol = length(phenotype_cols),
    dimnames = list(
      as.character(data[[id_col]]),
      paste0(phenotype_cols, "_adjusted")
    )
  )
  
  fitted_models <- setNames(
    vector("list", length(phenotype_cols)),
    phenotype_cols
  )
  
  rhs_text <- paste(
    deparse(adjustment_formula[[2]]),
    collapse = ""
  )
  
  for (pheno in phenotype_cols) {
    model_formula <- as.formula(
      paste(pheno, "~", rhs_text)
    )
    
    fit <- lm(model_formula, data = data)
    
    adjusted_matrix[, paste0(pheno, "_adjusted")] <- residuals(fit)
    fitted_models[[pheno]] <- fit
  }
  
  list(
    adjusted_matrix = adjusted_matrix,
    fitted_models = fitted_models
  )
}


validate_input_data <- function(
    reference_data,
    target_data,
    id_col,
    phenotype_cols,
    adjustment_formula,
    phenotype_lower_bound = 0
) {
  reference_required <- unique(
    c(id_col, phenotype_cols, all.vars(adjustment_formula))
  )
  
  ## Target covariates are not required because target KNN uses
  ## raw phenotypes. They may still be present in target_data.
  target_required <- unique(c(id_col, phenotype_cols))
  
  missing_reference <- setdiff(
    reference_required,
    colnames(reference_data)
  )
  
  missing_target <- setdiff(
    target_required,
    colnames(target_data)
  )
  
  if (length(missing_reference) > 0) {
    stop(
      "Reference data are missing: ",
      paste(missing_reference, collapse = ", ")
    )
  }
  
  if (length(missing_target) > 0) {
    stop(
      "Target data are missing: ",
      paste(missing_target, collapse = ", ")
    )
  }
  
  if (anyDuplicated(reference_data[[id_col]]) > 0) {
    stop("Reference IDs must be unique.")
  }
  
  if (anyDuplicated(target_data[[id_col]]) > 0) {
    stop("Target IDs must be unique.")
  }
  
  reference_numeric <- vapply(
    reference_data[, phenotype_cols, drop = FALSE],
    is.numeric,
    logical(1)
  )
  
  target_numeric <- vapply(
    target_data[, phenotype_cols, drop = FALSE],
    is.numeric,
    logical(1)
  )
  
  if (!all(reference_numeric) || !all(target_numeric)) {
    stop("All phenotype columns must be numeric in both cohorts.")
  }
  
  if (length(phenotype_lower_bound) == 1L) {
    phenotype_lower_bound <- rep(
      phenotype_lower_bound,
      length(phenotype_cols)
    )
  }
  
  if (length(phenotype_lower_bound) != length(phenotype_cols)) {
    stop(
      "phenotype_lower_bound must have length 1 or match phenotype_cols."
    )
  }
  
  reference_below_bound <- vapply(
    seq_along(phenotype_cols),
    function(j) {
      any(
        reference_data[[phenotype_cols[j]]] <
          phenotype_lower_bound[j],
        na.rm = TRUE
      )
    },
    logical(1)
  )
  
  target_below_bound <- vapply(
    seq_along(phenotype_cols),
    function(j) {
      any(
        target_data[[phenotype_cols[j]]] <
          phenotype_lower_bound[j],
        na.rm = TRUE
      )
    },
    logical(1)
  )
  
  if (any(reference_below_bound) || any(target_below_bound)) {
    stop(
      "Raw phenotype values must be greater than or equal to ",
      "their configured phenotype_lower_bound."
    )
  }
  
  reference_sd <- vapply(
    reference_data[, phenotype_cols, drop = FALSE],
    sd,
    numeric(1),
    na.rm = TRUE
  )
  
  if (any(!is.finite(reference_sd)) || any(reference_sd == 0)) {
    stop("Reference phenotype columns must have finite, non-zero variance.")
  }
  
  if (any(!complete.cases(
    reference_data[, reference_required, drop = FALSE]
  ))) {
    stop("Reference data contain missing required values.")
  }
  
  if (any(!complete.cases(
    target_data[, target_required, drop = FALSE]
  ))) {
    stop("Target data contain missing ID or phenotype values.")
  }
  
  invisible(TRUE)
}

############################################################
## 3. Create generic toy data
##
## Replace this section with readRDS(), read.csv(), load(),
## or another data-import step for real datasets.
############################################################

## Rows are latent directions and columns correspond to
## phenotype_cols. These names have no biological meaning.
component_directions <- rbind(
  c(-4.2,  0.4),
  c( 2.8,  3.2),
  c( 2.4, -3.5)
)

reference_data <- simulate_generic_cohort(
  n = 6500,
  phenotype_cols = phenotype_cols,
  cohort_name = "REF",
  component_prob = c(0.45, 0.25, 0.30),
  component_directions = component_directions,
  phenotype_center = c(40, 40),
  phenotype_noise_sd = c(2.8, 2.8),
  age_mean = 9.18,
  age_sd = 3.43,
  male_probability = 0.75,
  age_effect = c(0.30, 0.20),
  sex_effect = c(1.20, 0.80),
  phenotype_lower_bound = phenotype_lower_bound
)

target_data <- simulate_generic_cohort(
  n = 2700,
  phenotype_cols = phenotype_cols,
  cohort_name = "TARGET",
  component_prob = c(0.55, 0.19, 0.26),
  component_directions = component_directions * 0.90,
  phenotype_center = c(39, 41),
  phenotype_noise_sd = c(3.0, 3.0),
  age_mean = 8.98,
  age_sd = 3.54,
  male_probability = 0.865,
  age_effect = c(0.30, 0.20),
  sex_effect = c(1.20, 0.80),
  phenotype_lower_bound = phenotype_lower_bound
)

############################################################
## 4. Validate and summarize
############################################################

validate_input_data(
  reference_data = reference_data,
  target_data = target_data,
  id_col = id_col,
  phenotype_cols = phenotype_cols,
  adjustment_formula = adjustment_formula,
  phenotype_lower_bound = phenotype_lower_bound
)

cat("Reference N:", nrow(reference_data), "\n")
cat("Target N:", nrow(target_data), "\n")

print(
  summary(reference_data[, phenotype_cols, drop = FALSE])
)

print(
  summary(target_data[, phenotype_cols, drop = FALSE])
)

############################################################
## 5. Save input data and configuration
############################################################

input_config <- list(
  id_col = id_col,
  phenotype_cols = phenotype_cols,
  adjustment_formula = adjustment_formula,
  phenotype_lower_bound = phenotype_lower_bound,
  output_dir = output_dir
)

save(
  reference_data,
  target_data,
  input_config,
  adjust_phenotypes,
  validate_input_data,
  file = file.path(output_dir, "Generic_phenotype_input.RData")
)

cat(
  "Saved:",
  file.path(output_dir, "Generic_phenotype_input.RData"),
  "\n"
)
