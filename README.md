# Instructions
This repository contains the code required to reproduce the results from "Feature pre-selection for the development of epigenetic biomarkers".
These should be run in the order listed below to ensure each script has any dependencies generated from previous steps.

Each script can be run from within the repository using:
`Rscript path/to/script`

## Table of Contents
1. [Environment setup](#environment-setup)
2. [Configuration](#configuration)
3. [Preprocessing](#preprocessing)
4. [Fitting continuous trait models](#fitting-continuous-trait-models)
5. [Fitting EpiScore models](#fitting-episcore-models)
6. [Performance evaluation](#performance-evaluation)

## Environment setup (using Renv)
R version 4.3.0 is required. In addition, make sure renv 0.17.3 is installed.
The remaining dependencies can then be installed using renv by starting R in the top level directory of the repository and calling `renv::restore()`.

## Configuration
Edit the `config.yml` so the paths match those in your environment. These specify the location of:
- DNAm data files
- Phenotype target files
- MethylPipeR results output folder

## Preprocessing
Preprocessing scripts are located in `src/preprocessing/`.
`src/preprocessing/preprocessing.R` takes the EPIC array DNAm M-value files and target files for sets 1, 2, and 3 and produces preprocessed DNAm files and phenotype tables for application of pre-filtering methods.

## Fitting continuous trait models
Scripts for fitting the continuous trait models (part of the RTFS process) are located in `src/continuous_traits`.

## Fitting EpiScore models
Scripts for fitting the EpiScore models are located in `src/models`.

## Performance evaluation
Scripts for evaluating EpiScore model performance are located in `src/plots` and `src/tables`.