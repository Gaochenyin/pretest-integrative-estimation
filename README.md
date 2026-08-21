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

License and citation
- License: MIT (see LICENSE file)
- If you use this code in published work, please cite the accompanying paper "Pretest estimation in combining probability and non-probability samples" and include a reference to this GitHub repository.

Contact
- Author: Gaochenyin Gao (781747089@qq.com)
