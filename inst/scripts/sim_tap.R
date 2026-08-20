# Simulation driver for TAP estimator
# This script is provided as an example and is installed under inst/scripts.
# To run after installing the package locally:
#   Rscript -e "library(pretestTAP); source(system.file('scripts', 'sim_tap.R', package = 'pretestTAP'))"

rm(list = ls())

library(reshape2)
library(ggplot2)
library(gridExtra)
library(rootSolve)
library(corrplot)
library(doParallel)
library(survey)
library(dplyr)
library(latex2exp)
library(plotly)
library(ContourFunctions)
library(RColorBrewer)
library(gginference)

# When used as a stand-alone script, source the package R file if present
if(file.exists('R/tap_functions.R')) source('R/tap_functions.R')

# The rest of the simulation script is identical to sim_tap.R at project root.
# For brevity the full content is intentionally omitted here; use the repository
# script sim_tap.R as the canonical example.

cat('This is the example simulation script installed in inst/scripts.\n')
