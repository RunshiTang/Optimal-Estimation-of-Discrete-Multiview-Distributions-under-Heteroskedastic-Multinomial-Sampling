# Simulation Code for Multiview Multinomial Density Estimation

This repository contains the R code used for the simulation study in the paper

**Optimal Estimation of Discrete Multiview Distributions under Heteroskedastic Multinomial Sampling**.
**Authors:** Runshi Tang, Julien Chhor, Olga Klopp, Alexandre B. Tsybakov, and Anru R. Zhang
The simulations compare several estimators for discrete multiview distributions under heteroskedastic multinomial sampling. The implemented methods include the pooled histogram estimator, the unscaled spectral estimator, oracle scaling, CP-based oracle-scaling approximation, oracle slice normalization, and estimated slice normalization.

## Repository structure

```text
.
├── multiview_simulation_functions_v5.R
├── Run_v5_compare_thinning.R
├── Plot_v5_compare_thinning.R
└── multiview_four_experiments_v5/
    ├── experiment1_unscaled_thinning_comparison/
    ├── experiment2_dense_heteroskedastic/
    ├── experiment3_vary_heteroskedastic/
    ├── experiment4_rank_dimension_scaling/
    └── figures_pdf/
```

The main files are:

- `multiview_simulation_functions_v5.R`: helper functions for generating multiview multinomial data, running the estimators, computing error metrics, and projecting estimates onto the probability simplex.
- `Run_v5_compare_thinning.R`: runs the four simulation experiments and saves raw and summarized CSV files.
- `Plot_v5_compare_thinning.R`: reads the saved CSV files and generates the PDF figures.
- `multiview_four_experiments_v5/`: contains saved simulation outputs and generated figures.

## Requirements

The code is written in R. Most estimation routines use only base R and `stats`. Some methods and plotting routines require additional packages:

```r
install.packages(c("ggplot2", "ggpubr", "dplyr", "rTensor"))
```

The package `rTensor` is only needed for the `oracle_cp` method. Parallel simulation on Linux or macOS uses `parallel::mclapply()`.

## Running the simulations

Before running the scripts, update the working directory at the top of the R files if needed:

```r
setwd("~/density/v5")
```

Then run the simulation script:

```r
source("Run_v5_compare_thinning.R")
```

By default, this runs all four experiments:

1. **Experiment 1:** comparison of the unscaled estimator with and without multinomial thinning.
2. **Experiment 2:** dense heteroskedastic multiview distributions with varying sample size.
3. **Experiment 3:** dense multiview distributions with varying heteroskedasticity strength.
4. **Experiment 4:** rank and dimension scaling experiment.

The default settings can be computationally expensive. To run a smaller test, reduce the number of Monte Carlo replications in `Run_v5_compare_thinning.R`, for example:

```r
NREP_EXP1 <- 2
NREP_EXP2 <- 2
NREP_EXP3 <- 2
NREP_EXP4 <- 2
```

The simulation results are saved as CSV files under `multiview_four_experiments_v5/`.

## Generating figures

After the simulation output has been created, run:

```r
source("Plot_v5_compare_thinning.R")
```

The figures are saved under:

```text
multiview_four_experiments_v5/figures_pdf/
```

## Output files

Each experiment folder contains some or all of the following files:

- `raw_results.csv`: raw Monte Carlo simulation results.
- `summary_l1.csv`: summary of entrywise $\ell_1$ errors.
- `summary_fro.csv`: summary of Frobenius errors.
- `summary_normalized_l1.csv`: summary of normalized $\ell_1$ errors.
- `summary_normalized_fro.csv`: summary of normalized Frobenius errors.

## Notes

- The script `Run_v5_compare_thinning.R` currently runs all four experiments when sourced.
- To run experiments separately, call the corresponding functions in `Run_v5_compare_thinning.R`, such as `run_experiment2_dense_heteroskedastic_simulation()`.
- The saved CSV and PDF files in this repository are the outputs from the current simulation configuration.

## Citation

If you use this code, please cite the accompanying paper:

> Optimal Estimation of Discrete Multiview Distributions under Heteroskedastic Multinomial Sampling.
