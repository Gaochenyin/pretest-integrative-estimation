# Lightweight wrapper installed with the package.
# It prefers to source the canonical top-level sim_tap.R when available
# (developer clone). When installed from GitHub, the package ships this
# lightweight wrapper so users can run a short demo or get instructions.

# If running from the package installation and the full driver was shipped,
# attempt to source it from the installed scripts directory.
pkg_script <- system.file("scripts", "sim_tap.R", package = "pretestTAP")
if (nzchar(pkg_script) && file.exists(pkg_script) && !identical(pkg_script, "")) {
  # Avoid re-sourcing this wrapper itself if it ended up installed in scripts/
  # Check for a marker (we expect the canonical driver to be at the project root
  # during development; installed packages normally won't include the full driver).
  if (basename(pkg_script) != "sim_tap.R" || file.size(pkg_script) > 200) {
    try(source(pkg_script), silent = FALSE)
  } else {
    cat("This is the installed lightweight wrapper for the canonical sim_tap.R driver.\n")
    cat("When developing locally, run the canonical driver from the repository root:\n")
    cat("  Rscript sim_tap.R\n\n")
    cat("After installing the package, run:\n")
    cat("  Rscript -e \"library(pretestTAP); source(system.file('scripts', 'sim_tap.R', package = 'pretestTAP'))\"\n")
  }
} else if (file.exists(file.path('.', 'sim_tap.R'))) {
  # Developer environment: prefer the project's top-level driver
  source(file.path('.', 'sim_tap.R'))
} else {
  cat("This is a lightweight wrapper for the canonical sim_tap.R driver.\n")
  cat("During development, run the canonical driver from the repository root:\n")
  cat("  Rscript sim_tap.R\n\n")
  cat("After installing the package, run the installed wrapper with:\n")
  cat("  Rscript -e \"library(pretestTAP); source(system.file('scripts', 'sim_tap.R', package = 'pretestTAP'))\"\n")
}
