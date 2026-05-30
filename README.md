# AAPA: Anchor-Assisted Pedigree Assignment

[![R Package](https://img.shields.io/badge/R%20package-v0.1.0-blue)](https://github.com/luansheng/AAPA)

High-throughput full-sib family assignment using known candidate parents and anchor individuals.

## Overview

AAPA combines **Mendelian compatibility checking** with **anchor-based kinship scoring** to assign test individuals to candidate families, with built-in confidence scoring and rejection mechanisms for unknown-family detection.

### Key Features

- **Mendelian conflict scoring**: Check offspring-parent genotype compatibility across all SNP loci
- **Anchor kinship scoring**: Leverage known family members (anchors) for IBS-based similarity scoring  
- **Composite scoring**: Weighted combination of conflict and kinship evidence
- **Top-k pruning**: Efficient candidate selection for large family sets
- **Rejection mechanism**: Multi-rule rejection for low-confidence and unknown-family individuals
- **Quality control**: Built-in SNP/sample filtering by missing rate and MAF

## Installation

```r
# Install from GitHub (development version)
# install.packages("devtools")
devtools::install_github("luansheng/AAPA")
```

## Quick Start

```r
library(aapa)

# Simulate example data
sim <- simulate_aapa_data(n_families = 10, n_snps = 500)

# Prepare inputs (see vignette for details)
# ... build parents and anchors objects ...

# Run assignment
result <- aapa_assign(genotype, parents, anchors, top_k = 5)

# View results
print(result)
summary(result)
```

See `vignette("getting-started")` for a complete walkthrough.

## Package Structure

```
aapa/
  R/
    data_io.R     # Data input: genotype, parents, anchors
    scoring.R     # Mendelian conflict, anchor kinship, composite score
    assign.R      # Main assignment pipeline with rejection logic
    qc.R          # Quality control filters
    plot.R        # Visualization functions
  src/            # C++ placeholders for future Rcpp acceleration
  tests/testthat/ # Unit tests
  vignettes/      # Getting started guide
```

## Research Plan

See [AAPA-research-plan.md](AAPA-research-plan.md) for the full research methodology, benchmark design, and development roadmap.

## License

MIT
