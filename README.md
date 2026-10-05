# Simulation code for censored multinomial adverse-event meta-analysis

This directory contains an illustrative implementation of the current simulation settings used to evaluate Bayesian meta-analytic models for censored, multi-severity adverse-event (AE) counts. It is intended to accompany the manuscript **“Bayesian Meta-analysis of Comparative Drug Safety”** (under review). 

The repository used the latest simulation configuration for illustration purposes. It is not an archive of every preliminary scenario considered during paper development.

## Overview

Each simulated dataset contains 85 AE types and 44 studies. Two studies are two-arm randomized controlled trials (RCTs), and the remaining studies are single-arm trials (SATs). Each study arm has a sample size of 100.

Two evidence configurations are considered:

- **Case 1:** 2 RCTs, 29 drug-A SATs, and 13 drug-B SATs.
- **Case 2:** 2 RCTs, 21 drug-A SATs, and 21 drug-B SATs.

The simulation includes four grade 1-5 odds ratios equal to 2, three grade 1-5 odds ratios equal to 1/2, and one grade 3-5 odds ratio equal to 3. The grade 3-5 signal is nested within an AE having a grade 1-5 odds ratio of 1/2, as specified in Web Appendix C. The AE labels assigned to these signals are shuffled within each replicate. Replicate-specific data generation is reproducible because the random-number seed is `2025 + 123 * replicate_id`.

## Files

### R simulation scripts

| Script | Purpose |
|---|---|
| `NEW_multiDG_multifit_HS_comp.R` | Complete-data multinomial model with horseshoe priors for AE-specific treatment effects. |
| `NEW_multiDG_multifit_FX_comp.R` | Complete-data multinomial model with fixed-effect priors for AE-specific treatment effects. |
| `NEW_multiDG_binfit_HS_comp.R` | Complete-data benchmark that fits grade 1-5 and grade 3-5 outcomes separately using binomial models. |
| `NEW_both_multiDG_multifit_HS_cens3_submit.R` | Censoring-aware multinomial model with 30% censoring. |
| `NEW_both_multiDG_multifit_HS_cens6_submit.R` | Censoring-aware multinomial model with 60% censoring. |
| `NEW_both_multiDG_multifit_HS_cens3_CCA.R` | Available-case benchmark for data generated with 30% censoring. |
| `NEW_both_multiDG_multifit_HS_cens6_CCA.R` | Available-case benchmark for data generated with 60% censoring. |

The filenames are retained to preserve correspondence with the analysis records used for the revision. In the files labeled `CCA`, observations with both severity outcomes censored are omitted, while an observed component is retained when only one outcome is censored; these scripts are therefore more precisely described as available-case benchmarks.

### JAGS model files

| Model file | Used by |
|---|---|
| `multinomf_ORsel_HS_comp.txt` | Complete-data multinomial horseshoe model. |
| `multinomf_ORsel_FX_comp.txt` | Complete-data multinomial fixed-effect model. |
| `binomf_ORsel_HS_comp.txt` | Separate binomial horseshoe benchmark. |
| `multinomf_ORsel_HS_cens.txt` | Censoring-aware multinomial horseshoe model. |
| `multinomf_ORsel_HS_cens_CCA.txt` | Available-case multinomial horseshoe benchmark. |

## Software requirements

- R
- JAGS 4.x
- R packages `R2jags` and `LaplacesDemon`

The required R packages can be installed with:

```r
install.packages(c("R2jags", "LaplacesDemon"))
```

JAGS must be installed separately and visible to `R2jags`.

## Running a simulation replicate

Each R script takes exactly two command-line arguments:

1. `case_id`: `1` or `2`;
2. `replicate_id`: any positive integer.

For example, from this directory:

```bash
Rscript NEW_multiDG_multifit_HS_comp.R 1 1
Rscript NEW_both_multiDG_multifit_HS_cens3_submit.R 2 1
```

The scripts locate their associated JAGS files relative to the script location, so they can be launched from another working directory. No user-specific path needs to be edited.

## Output

Each run creates one `.RData` file in `results/`. The saved object is named `res` and contains:

- the true grade 1-5 and grade 3-5 odds ratios;
- the true drug-specific AE probabilities;
- posterior summaries when JAGS fitting succeeds;
- signal-selection curves and areas under the corresponding curves.

The `results/` directory and generated `.RData` files are excluded from version control by `.gitignore`.

## Reproducibility notes

- The R data-generating mechanism and JAGS model are separated so that each model specification can be inspected independently.
- Input dimensions and the nesting constraint for AE probabilities are checked before model fitting.
- Censoring counts are computed explicitly with `floor(censoring_rate * n_AE)`.
- The supplied MCMC settings reproduce the analysis configuration and can require substantial computation. For large simulation batches, run different `(case_id, replicate_id)` combinations as independent jobs.

## Data and privacy

The simulations generate synthetic data and do not require patient-level or proprietary input files.
