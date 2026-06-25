# Differential Expression Analysis with Long Read RNA-Seq Data

Software accompanying the paper *“Differential Expression Analysis with Long Read RNA-Seq Data”* by Liu Z., Song C., Wang Z., Fan Q., Yu J., Hao N., Chen S., Sun X., Ding H. (2026)

---

## Overview

This R package implements a **hurdle negative binomial generalized linear model (hurdle-NB GLM)** for differential expression analysis of long-read RNA-Seq data. It provides:

- Data preparation and validation (`prepareDGE()`). Support for **matrix**, **data.frame**, **DGEList**, **DESeqDataSet**, or **SummarizedExperiment** inputs  
- Estimation of **size factors** for normalization (`sizeFactorsEst()`)  
- Estimation of **tag-wise dispersions** using emperical priors (`tagwiseEst()`)  
- Gene-level differential expression testing using:
  - **Likelihood Ratio Test** (`hurdle_LRT()`)  
  - **Wald Test** (`hurdle_Wald_Test()`)  

The package is optimized for long-read sequencing data with high dropout rates, where zero inflation and overdispersion are common.

---

## File Structure

| File / Folder                         | Description                                                                 |
|---------------------------------------|-----------------------------------------------------------------------------|
| `R/estimateSizeFactors.R`             | Functions for estimating normalization size factors for RNA-Seq counts.    |
| `R/estimateDispersion.R`              | Functions for prior estimation and tag-wise dispersion estimation.          |
| `R/hurdleNBTest.R`                    | Functions for hurdle-NB differential expression testing (Wald & LRT).      |
| `DESCRIPTION`                         | Package metadata, dependencies, and version information.                    |
| `NAMESPACE`                           | Exports, imports, and internal function definitions.                        |
| `README.md`                           | Package overview, usage, and installation instructions.                     |

---

## Installation

To install the stable **release** version of `LRDE` from Bioconductor:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LRDE")
```

To install the development version (currently required while under review):
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LRDE", version = "devel")
```

## Example usage
```r
# Load the package
library(LRDE)

# 1. Simulate a count matrix (negative binomial)
set.seed(123)
mat <- matrix(rnbinom(300, size = 5, mu = 5), nrow = 50)
grp <- factor(c("A", "A", "A", "B", "B", "B"))

# 2. Prepare the data object
y <- prepareDGE(mat, grp)

# 3. Estimate size factors for normalization
y <- sizeFactorsEst(y)

# 4. Estimate tag-wise gene dispersions
y <- tagwiseEst(y)

# 5a. Differential expression testing: Wald test
y <- hurdle.Wald.Test(y)

# 5b. Differential expression testing: Likelihood Ratio Test
y <- hurdle.LRT(y)

# 6. Access results
y$lrt_stats        # LRT statistics per gene
y$wald_stats       # Wald statistics per gene
y$p.values         # Corresponding p-values
y$tagwise.disp     # Estimated dispersions
y$zero_prob_matrix # Zero probabilities per group

```

---

## Contact
For questions, comments, or collaborations, please contact:

- **Hongxu Ding** — <hongxuding@arizona.edu>
- **Ziyang Liu** — <ziyang773@arizona.edu>

---
