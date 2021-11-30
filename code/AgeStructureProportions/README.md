# Age Structure & Gonotrophic Cycles Model

The R script "AgeStructureGonotrophicCyclesModel.R" estimates mosquito age class proportions with
95% credible intervals. Details on estimation are available in the comments of the script. This
README file contains instructions on how to run the script.

## System requirements

The script runs on R software (https://cran.r-project.org/).
It has been tested on R version 4.1.1 for macOS X but should run on any version of R going back
at least 10 years and on any operating system, as it uses only base R functions and requires
no packages.

## Installation guide

- Place the script in a directory, with the data file if available. Open in any R editor/GUI/IDE.
- Submit the script to the R console.

## Instructions for use on real or demo data

- If the input data set is present, the script will load it to R and run the analysis.
- If the input data set is not present, the script will automatically run the analysis on demo data,
and give a warning to that effect.
- Outputs: estimates and 95% credible intervals for the age class proportions are printed to the
console and plotted as a bar plot. The script should run in less than 1 second.
