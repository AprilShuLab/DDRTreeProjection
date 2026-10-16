# Reference-Based DDRTree Subtype and Subtype Score Projection

This repository provides a general workflow for deriving categorical subtypes and continuous subtype scores in a reference cohort and projecting them to an independent target cohort. The framework combines a reference DDRTree, a user-defined partition of its principal Y-node tree, K-nearest-neighbor (KNN) mapping, tree-informed subtype similarities, and optional stratified-bootstrap stability analysis.

The workflow accompanies our study of autism spectrum disorder (ASD) support-need subtypes in the SPARK discovery cohort and their projection to the Simons Simplex Collection (SSC). Because individual-level SPARK and SSC data cannot be redistributed, the repository uses fully synthetic reference and target cohorts with generic phenotype names. The synthetic data demonstrate the computational workflow only and do not represent real participants or biological findings.

> **Manuscript:** Tan Q, et al. *Machine Learning-Derived Subtypes of Autism Reveal Distinct Genetic Contributions from Common and Rare Variants.* [Journal or preprint information to be added].

## Method overview

![Reference-based projection of categorical subtypes and continuous subtype scores](figures/method_overview.pdf)

The workflow consists of three scripts:

1. Synthetic reference and target cohorts are generated, validated, and saved with a shared input configuration.
2. Covariate-adjusted reference phenotypes are used to fit a DDRTree. After the root node and subtype-specific Y-node sets are defined, the workflow derives a categorical subtype label and one continuous score for each subtype for every reference individual.
3. Target individuals are mapped using K nearest reference neighbors identified in standardized raw phenotype space. The neighbors' reference-tree positions are used to project subtype similarities, categorical subtype assignments, DDRTree coordinates, and continuous subtype scores.

Importantly, the reference DDRTree is constructed from **covariate-adjusted phenotypes**, whereas KNN mapping is performed using **raw phenotypes standardized with the reference-cohort means and standard deviations**. The target cohort is not used to refit or optimize the reference DDRTree.

Although the accompanying study used three ASD support-need subtypes, the implementation supports a general partition of the reference tree into *K* subtypes.

## Repository contents

```text
.
├── 00_generic_phenotype_data_preparation.R
├── 01_reference_DDRTree_workflow.R
├── 02_target_mapping_workflow.R
├── order_cells_igraph214.R
├── figures/
│   └── method_overview.png
├── output/
├── README.md
└── LICENSE
```

| Script | Purpose |
|---|---|
| `00_generic_phenotype_data_preparation.R` | Generates and validates synthetic reference and target cohorts and saves the shared input configuration. Section 3 can be replaced with real-data import. |
| `01_reference_DDRTree_workflow.R` | Fits the reference DDRTree, reconstructs the Y-node minimum spanning tree (MST), supports manual root and subtype-branch definition, and calculates reference labels and continuous subtype scores. |
| `02_target_mapping_workflow.R` | Maps the target cohort to the fixed reference model, generates mapping and score plots, and optionally performs stratified-bootstrap stability analysis. |
| `order_cells_igraph214.R` | Provides the Monocle 2 `orderCells()` compatibility implementation used by the reference workflow with newer versions of `igraph`. |

## Software requirements

The workflow is written in R and uses the following software and packages:

| Software/package | Version |
|---|---:|
| `R` | 4.4.2 |
| `DDRTree` | 0.1.5 |
| `monocle` | 2.34.0 |
| `igraph` | 2.1.4 |
| `Matrix` | 1.7-2 |
| `proxy` | 0.4-27 |
| `ggplot2` | 3.5.2 |
| `FNN` | 1.1.4.1 |
| `scales` | 1.4.0 |

## Quick start

Clone the repository and run all scripts from the repository root.

```bash
git clone https://github.com/AprilShuLab/DDRTreeProjection.git
cd DDRTreeProjection
```

### Step 1: prepare the example data

```bash
Rscript 00_generic_phenotype_data_preparation.R
```

This script creates synthetic reference and target cohorts and saves:

```text
output/Generic_phenotype_input.RData
```

The saved file contains `reference_data`, `target_data`, `input_config`, and the input-preparation helper functions required by the downstream scripts.

### Step 2: construct and finalize the reference model

The reference workflow is intentionally run in two stages because the root node and complete subtype-specific Y-node partition must be defined after inspecting the fitted reference tree.

#### Stage 2A: fit and inspect the DDRTree

In `01_reference_DDRTree_workflow.R`, retain:

```r
DEFINE_SUBTYPES <- FALSE
root_node <- NULL
subtype_Y_nodes <- NULL
```

Then run:

```bash
Rscript 01_reference_DDRTree_workflow.R
```

The first-stage run saves a preliminary reference model and diagnostic plots, including:

```text
output/Generic_reference_DDRTree_preliminary.RData
output/reference_plots/01_numbered_Y_node_MST.pdf
output/reference_plots/02_reference_Z_and_Y_MST.pdf
```

Inspect `01_numbered_Y_node_MST.pdf` before defining the root and subtype-specific Y-node sets.

#### Stage 2B: define the root and subtype branches

Enter one root node and a named list containing the complete Y-node partition. For example:

```r
root_node <- 11

subtype_Y_nodes <- list(
  Subtype_1 = c(...),
  Subtype_2 = c(...),
  Subtype_3 = c(...)
)

DEFINE_SUBTYPES <- TRUE
```

Every Y node must appear **exactly once** across the subtype-specific lists. The script validates the root, the complete node partition, within-subtype connectivity, subtype terminal nodes, and the required graph-distance scales before calculating the final outputs.

Rerun:

```bash
Rscript 01_reference_DDRTree_workflow.R
```

The finalized reference model is saved as:

```text
output/Generic_reference_DDRTree_model.RData
```

It contains the fitted DDRTree, Y-node MST, reference-to-Y-node mapping, root node, subtype partition, subtype terminal nodes, root-to-Y-node graph distances, subtype-specific distance scales, reference subtype labels, and continuous reference subtype scores.

### Step 3: map the target cohort

Set the number of neighbors and bootstrap options in `02_target_mapping_workflow.R`:

```r
mapping_k <- 3

RUN_BOOTSTRAP <- TRUE
bootstrap_B <- 5
bootstrap_seed <- 123
```

Then run:

```bash
Rscript 02_target_mapping_workflow.R
```

The mapping results are saved as:

```text
output/Generic_target_mapping_result.RData
```

Target mapping and subtype-score plots are written to:

```text
output/mapping_plots/
```

## Input data

The synthetic example uses the following columns:

| Column | Cohort | Description |
|---|---|---|
| `ID` | Reference and target | Unique participant identifier. |
| `Age` | Reference and target | Age in months in the synthetic example. |
| `Sex` | Reference and target | Binary covariate in the synthetic example (`1` = male, `0` = female). |
| `Pheno1`, `Pheno2` | Reference and target | Generic phenotypes shared between the two cohorts. |

The default reference adjustment model is:

```r
adjustment_formula <- ~ Age + I(Age^2) + Sex
```

To use real data, replace Section 3 of `00_generic_phenotype_data_preparation.R` with the appropriate data-import code and update:

```r
id_col
phenotype_cols
adjustment_formula
phenotype_lower_bound
```

The following requirements apply:

- IDs must be unique within each cohort.
- Reference and target cohorts must contain the same phenotype variables with consistent definitions and ordering.
- Required IDs, phenotypes, and reference adjustment covariates must be complete before the workflow is run.
- Missing values should be imputed or otherwise handled before data preparation.
- Phenotype lower bounds apply to raw phenotypes only; adjusted residuals may legitimately be negative.
- Reference and target sample order should not be manually changed within saved workflow objects.

## Reference subtype labels and continuous scores

Each reference individual is assigned to its nearest Y node. The user-defined Y-node partition determines the categorical reference subtype label.

For every subtype, a representative individual at the distal end of the corresponding reference-tree branch is used to establish subtype-specific ordering with Monocle 2 `orderCells()`. The resulting pseudotime is reversed and rescaled to the interval from 0 to 1, producing one continuous score per subtype. Higher values indicate greater proximity along the learned reference trajectory to the corresponding subtype-specific phenotypic extreme.

These scores quantify relative position on the learned reference trajectory; they are not probabilities and are not required to sum to 1.

## Target-cohort projection

Raw reference phenotypes are standardized using their reference-cohort means and standard deviations. The same reference parameters are then applied to the target phenotypes. For each target individual, the workflow identifies `mapping_k` nearest reference neighbors in this standardized phenotype space.

The reference Y-node positions of the selected neighbors are used to calculate a tree-informed similarity for each subtype:

1. Neighbors are grouped by their mapped reference Y nodes.
2. For each represented subtype, the workflow calculates the count-weighted mean graph distance from the common root to the neighbors' Y nodes.
3. This distance is divided by the corresponding subtype-specific branch scale.
4. The resulting subtype ratios are normalized across all subtypes.

The categorical target subtype is the subtype with the largest normalized tree-informed ratio. `Mapping_confidence` is the largest normalized ratio for that target individual.

All selected reference neighbors are retained when projecting continuous subtype scores. Their contributions are weighted using the tree-informed subtype similarity associated with each neighbor's reference subtype. The target DDRTree coordinates are projected using the corresponding normalized neighbor weights.

## Main target outputs

`target_mapping_result` contains one row per target individual and includes:

| Output | Description |
|---|---|
| `Individual` | Target participant identifier. |
| `Predicted_subtype` | Projected categorical subtype. |
| `Mapping_confidence` | Maximum normalized tree-informed subtype ratio. |
| `Dim1`, `Dim2`, ... | Projected target coordinates in the reference DDRTree space. |
| `<Subtype>_ratio` | Normalized tree-informed similarity for each subtype. Ratios sum to 1 across subtypes. |
| `<Subtype>_score` | Projected continuous score for each subtype. Scores are bounded from 0 to 1 but are not required to sum to 1. |

The saved `mapping_result` object additionally contains the KNN indices and distances, projected score matrices, target branch ratios, target coordinates, and neighbor-weight details.

## Optional stratified-bootstrap stability analysis

When `RUN_BOOTSTRAP <- TRUE`, reference individuals are sampled with replacement independently within each reference subtype. This preserves the original subtype-specific sample size and subtype composition in every bootstrap replicate. The complete target-mapping procedure is then repeated using each bootstrap reference sample while keeping the learned reference-tree information fixed.

The saved `bootstrap_stability` object contains:

| Object | Description |
|---|---|
| `prediction_matrix` | Projected subtype label for every target individual in every bootstrap replicate. |
| `score_arrays` | Projected continuous subtype scores across bootstrap replicates. |
| `label_stability` | Subtype-specific assignment counts, maximum consistency, final subtype, and assignment-count combination for each target individual. |
| `score_stability` | Individual-level mean and standard deviation of each projected subtype score. |
| `score_sd_summary` | Mean individual-level bootstrap standard deviation for each subtype score. |
| `combination_distribution` | Frequency and proportion of each bootstrap label-count combination. |

The default `bootstrap_B = 5` matches the demonstration and the stability-analysis design used in the accompanying study. Larger values may be used in other applications when more precise bootstrap summaries are required.

## Monocle 2 and `igraph` compatibility

The reference workflow uses a locally adapted version of the Monocle 2 `orderCells()` implementation (`order_cells_igraph214.R`). This compatibility version was derived from the original Monocle 2 source code and updated to accommodate API changes introduced in newer versions of the `igraph` package. The modifications were limited to package-compatibility updates and did not intentionally alter the underlying pseudotime-ordering procedure.

Monocle 2 was developed by the Cole Trapnell laboratory and is distributed under the Artistic License 2.0. Users of this workflow should cite the original Monocle and DDRTree publications.

## Data availability

Individual-level SPARK and SSC data are not distributed in this repository because they are subject to controlled-access and data-use requirements. Researchers may request access through the corresponding SPARK and Simons Simplex Collection data-access procedures.

All participant-level data included in this repository are synthetic. They were created solely to demonstrate the computational workflow and must not be interpreted as representing real participants, cohort statistics, clinical characteristics, or biological findings.

## Citation

If you use this workflow, please cite the accompanying manuscript, Monocle 2, and DDRTree:

1. Tan Q, et al. *Machine Learning-Derived Subtypes of Autism Reveal Distinct Genetic Contributions from Common and Rare Variants.* [Journal, year, DOI to be added].
2. Qiu X, et al. Reversed graph embedding resolves complex single-cell trajectories. *Nature Methods.* 2017;14:979–982. https://doi.org/10.1038/nmeth.4402
3. Trapnell C, et al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. *Nature Biotechnology.* 2014;32:381–386. https://doi.org/10.1038/nbt.2859

## License

The original workflow code in this repository is released under the MIT License. The adapted Monocle 2 compatibility implementation in `order_cells_igraph214.R` remains subject to the Artistic License 2.0 and the attribution provided in that file.

## Contact

For questions about the workflow, please open a GitHub issue or contact:

**Qing Tan**  
Center for Genetic Epidemiology  
Department of Population and Public Health Sciences  
Keck School of Medicine of USC  
Email: qingtan@usc.edu
