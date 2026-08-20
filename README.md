# Pretest integrative estimation (TAP)

An R package and set of scripts implementing the Test-and-Pool (TAP) integrative estimation method from the paper "Pretest estimation in combining probability and non-probability samples". The code provides tools to simulate synthetic populations, apply the TAP estimator to combine a probability sample (A) and a non-probability sample (B), and reproduce the figures and simulation results from the paper.

Key features
- Simulation utilities to generate synthetic populations and draw samples with different selection mechanisms.
- Analytic expressions and Monte Carlo approximations for TAP bias, variance, and MSE.
- A main TAP estimator (TAP.Est) that runs propensity estimation, point estimation, bootstrap inference, and adaptive confidence interval construction.
- Example simulation driver that produces the figures used in the manuscript.

Installation (development)
1. Install prerequisites (if needed):

   install.packages(c("devtools", "remotes"))

2. Install the package from GitHub (development version):

   remotes::install_github('Gaochenyin/pretest-integrative-estimation')

3. During development, load the package directly from the project root:

   # from an R session in the project directory
   devtools::load_all('.')

Quick start examples
- Run the provided example simulation (installed under inst/scripts):

   Rscript -e "library(pretestTAP); source(system.file('scripts', 'sim_tap.R', package = 'pretestTAP'))"

- Or run the repository's driver script directly (produces the MSE surface and adaptive-interval figures and performs Monte Carlo simulations):

   Rscript sim_tap.R

Why the R/ directory and inst/scripts/ are both present
- R/: contains the package's R source files that are built and installed as the package API. Place core functions you want available via library(pretestTAP) here (for example, GenerateSimuDta and the analytic helpers).
- inst/scripts/: contains example scripts that are shipped with the installed package. Files here are copied into the installed package under <installed-package>/scripts/ and can be found at runtime with system.file('scripts', 'your_script.R', package = 'pretestTAP').

Practical note on the two sim_tap.R files in this repo
- Canonical (full) driver: the top-level sim_tap.R in the repository is the full simulation driver used during development and to reproduce the paper figures.
- Installed wrapper: inst/scripts/sim_tap.R is a lightweight wrapper (installed with the package). It intentionally avoids duplicating the full driver so the two copies don't drift; when you clone the repository you can run the full sim_tap.R directly, and the installed wrapper will instruct users how to run the canonical script or run a short demo when appropriate.

Repository layout
```text
R/                 # Core implementation (tap_functions.R) with GenerateSimuDta, TAP.Est (stub for installed package), and analytic helpers
inst/scripts/      # Installed wrapper scripts (inst/scripts/sim_tap.R) that point to the canonical driver or provide brief demos
sim_tap.R          # Canonical simulation driver used to reproduce figures and experiments (development)
DESCRIPTION         # R package metadata
NAMESPACE           # R package namespace
README.md           # This file
```

Important files and functions
- R/tap_functions.R: main implementation. Key exported functions include:
  - GenerateSimuDta(...): create a synthetic population and draw samples A (probability) and B (non-probability).
  - TAP.Est(...): main estimator and inference routine (propensity estimation, point estimation, bootstrap, adaptive CI) — note: the installed package may contain a lightweight stub; to use the full implementation, run the repository version or load the package from the project root during development.
  - LambdaCgammaMSE, LambdaCgammaMSE.MC, OptimLambdaCgamma.MC: analytic and Monte Carlo routines for evaluating/optimizing lambda and c_gamma tuning parameters.
  - NormalTruncatedFirstMom / NormalTruncatedSecondMom: helpers used in analytic bias/variance expressions.

How to run the canonical simulation driver (from a clone)
```bash
# clone the repo
git clone https://github.com/Gaochenyin/pretest-integrative-estimation.git
cd pretest-integrative-estimation
# run the full driver (this may take a long time depending on B/K and niter)
Rscript sim_tap.R
```

How to run the installed wrapper (after remotes::install_github)
```r
# run the lightweight wrapper included in the installed package
Rscript -e "library(pretestTAP); source(system.file('scripts', 'sim_tap.R', package = 'pretestTAP'))"
```

Development notes
- This package is intentionally minimal for development. To generate Rd documentation from roxygen comments run:

   devtools::document()

- To run package checks (CRAN style):

   devtools::check()

- Large functions are separated in this repository for readability during development; for CRAN submission consider splitting them into smaller helpers and adding more granular unit tests.

License and citation
- License: MIT (see LICENSE file)
- If you use this code in published work, please cite the accompanying paper "Pretest estimation in combining probability and non-probability samples" and include a reference to this GitHub repository.

Contact
- Author: Gaochenyin Gao (781747089@qq.com)
