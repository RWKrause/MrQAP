# Regenerate the golden-master fixture.
#
#   Rscript data-raw/make_golden.R
#
# ONLY run this deliberately, when a numeric change is intended and
# understood. A failing test-golden.R means a refactor changed results --
# diagnose it, do not re-baseline it.

devtools::load_all(".", quiet = TRUE)
source("tests/testthat/helper-golden.R")

dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
golden <- golden_run_all()
saveRDS(golden, "tests/testthat/fixtures/golden.rds", version = 3)

cat("wrote tests/testthat/fixtures/golden.rds\n")
cat("cases:", length(golden), "\n")
cat("elapsed:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n\n")

for (nm in names(golden)) {
  g <- golden[[nm]]
  cf <- if (is.list(g$coefficients)) g$coefficients[[1]] else g$coefficients
  cat(sprintf("  %-24s %2d coef  b1=%9.5f\n", nm, length(cf),
              if (length(cf) > 1) cf[[2]] else NA_real_))
}
