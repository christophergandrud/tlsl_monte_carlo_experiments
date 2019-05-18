# Monte Carlo Experiments for Time-lagged Spatial-lagged Regressions

Cody A. Drolc, Christopher Gandrud, and Laron K. Williams

## About 

This repository contains all of the source code to replicate the Monte Carlo simulations presented in "Taking Time (and Space) Seriously:
How Scholars Falsely Infer Policy Diffusion from Model Misspecification".

## Reproduction instructions

### With GNU make

The simulations are run in R and coordinated with GNU make. If you have R and make installed, replicate the entire analysis by typing into your terminal:

```bash
cd THE_DIRECTORY_WITH_THE_REPOSITORY

make clean
make
```

### Manually without GNU make

You can also manually run each R file one at a time. Do this in the alphabetical order of the `.R` file names, e.g. start with `a_setup.R` and end with `z_mc_results_plots.R`.   

