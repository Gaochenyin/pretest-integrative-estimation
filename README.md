# Pretest integrative estimation (TAP)

This repository contains code to reproduce the paper "Pretest estimation in combining probability and non-probability samples".

Installation (development)

1. Install devtools / remotes if not already installed:

   install.packages(c('devtools','remotes'))

2. Install from GitHub:

   remotes::install_github('Gaochenyin/pretest-integrative-estimation')

3. Example: run the provided simulation script after installing the package:

   Rscript -e "library(pretestTAP); source(system.file('scripts', 'sim_tap.R', package = 'pretestTAP'))"

Notes

- The package created here is minimal and intended for local development. Run devtools::document() to generate documentation from roxygen comments and devtools::check() to validate the package prior to CRAN.
- Larger functions (e.g., TAP.Est) are long; for CRAN-style packaging they should be split into smaller helpers and documented thoroughly.
