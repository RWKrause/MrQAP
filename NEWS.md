# MrQAP 0.3.0

## Breaking changes

* **`QAPglm()` and `QAPcss()` have been replaced by a single `QAP()`.**
  Whether the data are network matrices or CSS arrays is a property of the
  data, not a reason for two functions, and is now detected from the number
  of dimensions (override with `css =`). Existing scripts need
  `QAPglm(...)` / `QAPcss(...)` changed to `QAP(...)`; no other argument
  changed name.

* **`nullhyp` now defaults to `"qapspp"` everywhere.** `QAPglm()` already
  defaulted to semi-partialling plus; `QAPcss()` defaulted to `"qapy"`.
  CSS models called without an explicit `nullhyp` will now use the stronger
  (and more expensive) test. Pass `nullhyp = "qapy"` to keep the old CSS
  behaviour.

* **The returned object has the same shape for every input shape.** CSS fits
  previously exposed results only under `$base`; they now also carry
  `$coefficients`, `$t` and the family extras at the top level, as matrix
  fits always did.

## Results that change

Two fixed bugs alter numbers produced by earlier versions. **Re-run any
analysis affected by these.**

* **`mode = "undirected"` was silently ignored for network matrices.** The
  `mode` argument was accepted and never read, so every dyad entered the
  model twice: n was doubled, standard errors were roughly a factor of
  sqrt(2) too small, and all p-values were too optimistic. CSS data were not
  affected.

* **Robust standard errors were wrong for every non-gaussian family.**
  `use_robust_errors = TRUE` fed *deviance* residuals and unweighted
  leverages into an OLS HC3 formula. On a binomial fit this inflated
  t-values by a factor of 2 to 3. GLMs now use `sandwich::vcovHC()`, fixest
  models use their own heteroskedasticity-robust covariance, and the
  closed-form HC3 is used only for linear models, where it is exact.
  `sandwich` is a new (suggested) dependency.

* **Multinomial p-values were wrong wherever they appeared at all.** A
  multinomial fit's coefficients are a `(ncat - 1) x k` matrix, so the
  permutation result matrices stack a block of coefficient comparisons
  above a block of t-comparisons. `print()` read row 2 and labelled it
  `Pr(<=t)` -- for three categories that is the *second category's
  coefficient* row -- and recycled a length-`k` vector into a
  `(ncat-1)*k` table. `summary()` and `as.data.frame()` reported `NA`
  throughout. Both now unpack the block. Under the default
  `nullhyp = "qapspp"` a multinomial fit did not run at all; see below.

* **GMM estimates were sensitive to their starting values.** `gmm()` was
  started from `rnorm()`, redrawn on every permutation. On a well-scaled
  design this was harmless, but where the moment conditions have several
  local optima it was not: 25 random starts on one Poisson design gave `x1`
  estimates from -80.6 to +184.8 (sd 24.6) where the truth was 0.6, and
  that spread went straight into the null distribution. Starting values now
  come from the corresponding GLM fit, which is deterministic and
  data-driven; the same design gives 0.528. **Re-run any GMM analysis.**

## New features

* **`predict()`**, including `type = "matrix"`, which puts the predictions
  back into the shape of the network with `NA` in every structurally
  excluded cell -- the diagonal, the redundant triangle in undirected mode,
  anything dropped for missingness. Accepts `newdata` as a named list of
  matrices shaped like the fitted data.

* **`plot()`** draws the permutation distribution of each tested
  coefficient with the observed statistic marked: the picture the p-value
  summarises.

* **`weights` and `offset`** arguments, supported for every family and for
  both estimators. For GMM they enter the moment conditions themselves;
  for Poisson and logit, where the moments are the score equations and the
  system is just-identified, the result matches the weighted, offset MLE
  exactly. Weights and offsets are properties of the cell, not of the
  actors, so they stay in place while the outcome (or the tested
  predictor) is permuted -- see `?QAP`.

* **`data` accepts more shapes.** A dyadic or triadic data frame is passed
  to `df_to_mat()` for you (`sender=`, `receiver=`, `perceiver=`,
  `split_by=`), and `igraph` or `network` objects anywhere in the list are
  converted to adjacency matrices (`edge_attr=` to read a weight).

* **`tidy()` and `glance()`** register with **broom** when it is installed.
  `tidy()` deliberately omits `std.error`: a QAP fit has no single standard
  error that the permutation p-value corresponds to.

* **Standard extractor methods**: `coef()`, `summary()`, `confint()`,
  `vcov()`, `nobs()`, `fitted()`, `residuals()`, `logLik()` (hence `AIC()`
  and `BIC()`) and `as.data.frame()`. Results no longer have to be read off
  the printed table.

  `summary()` prints the same report as the fit, with significance codes,
  and `summary(fit)$coefficients` gives a data frame of estimates,
  statistics and the three p-values. All of them handle `comparison=`
  models, returning one result per comparison.

* **`confint()` uses the permutation distribution.** The permuted
  coefficients are now retained (`fit$null_dist`, a few tens of kilobytes;
  suppressed by `less_mem`), and intervals are formed as
  `estimate ± z * sd(permutation draws)`. Because permutation preserves the
  network's dependence structure, this reflects dyadic dependence that a
  model-based standard error ignores.

  Note this is *not* a percentile interval of the permutation distribution:
  that distribution is centred on zero under the null, so its percentiles
  bracket zero rather than the estimate. Only its spread is used, under a
  normal approximation.

* **Permutation p-values print honestly.** A p-value of 0 means "no
  permutation was this extreme", i.e. below `1/reps` — it is now shown as
  e.g. `<0.005` rather than the `<2e-16` that the conventional coefficient
  formatter produces, and the resolution is stated beneath the table.

* Fitted objects gain a shared `"QAP"` parent class, which the methods
  above dispatch on. The specific classes (`QAPRegression`, `QAPGLM`,
  `QAPCSS`) are unchanged and still drive `print()`.

* **Multi-mode data are supported**: each array dimension may be its own
  node set, and each is then permuted independently. Controlled by
  `multi_mode`, which infers `TRUE` when the dimensions differ in size.
  A square matrix is genuinely ambiguous, so it is assumed to be one-mode
  unless you say otherwise. Self-ties and undirected mode are undefined for
  multi-mode data and are now rejected rather than silently misapplied.

* `QAP()` validates `family`, `mode`, `nullhyp` and `estimator` by name, so
  a typo is an immediate error rather than a confusing downstream failure.

* Predictors are checked for conformability with the outcome. A
  wrong-sized predictor used to recycle a mismatched index and produce
  silently incorrect data.

* Perfectly collinear predictors now warn once and report `NA`, instead of
  desynchronising the coefficient and t-value vectors.

## Other bug fixes

* `sample()` was called on permutation groups of size one, where
  `sample(7)` returns a permutation of `1:7` rather than `7`. Any grouped
  design with a singleton group was affected.

* `RMPerm()` dropped `CSS = TRUE` when recursing over a list of arrays.

* `combine_qap_estimates()` never actually pooled comparison models: it
  indexed `res[[comparison]]$lower` where results are stored as
  `res$lower[[comparison]]`, so it incremented `$reps` while leaving the
  p-values untouched.

* Per-network permutation groups were handed to every network instead of
  the one they belonged to.

* `QAPcss()` built its model formula before the checks that switch random
  intercepts off, so the formula kept terms the warnings said were dropped.

* A CSS permutation that could not find a usable draw aborted the entire
  run from inside a parallel worker; it is now dropped and reported like
  any other failed permutation.

* `array_to_vector()` looped over dimension 1 while slicing dimension 3.
  Cubic input hid this completely; on a 2x3x4 array it silently returned 12
  of 24 values.

* CSS fits always reported "permutations were performed within groups only",
  because `groups` was internally defaulted to a single all-in-one group.

* The rows of `$lower`, `$larger` and `$abs` are now labelled `"b"` and
  `"t"` on every path. Previously the outcome-permutation path leaked
  `rbind()`'s argument names (`"perm_coefs"`, `"perm_t"`) into the result,
  the GPU path had no row labels at all, and the semi-partialling path had
  none either. Values are unchanged; code indexing these matrices by row
  *number* is unaffected.

* **`QAP(family = "multinom")` errored on its own default settings.**
  `qap_init_pmats()` sized the result matrices with `length(coefficients)`
  -- `(ncat-1)*k` for a matrix -- instead of `ncol()`, and took `names()`
  of a matrix, which is `NULL`. The per-predictor assignment then failed,
  so multinomial fits were unusable under `nullhyp = "qapspp"`.

* **`less_mem = TRUE` did not save any memory.** It skipped only the
  top-level copy of the fitted model while the same object survived in
  `fit$base$base_model`: at n = 120 the result was 4.68 MB of 4.69 MB,
  even as `fitted()`, `vcov()` and the rest refused to run and advised
  refitting without it. It now drops the model for real (9.37 MB to
  0.01 MB), and `coef()`, `nobs()` and `summary()` still work.

* **`nobs()` returned `NA` under `less_mem`.** The count is now recorded
  when the data are vectorised, so it no longer depends on the model
  object.

* **`combine_qap_estimates()` would pool unrelated models.** It checked
  only the class, `nullhyp`, comparison names and dimnames, so two fits to
  *different networks* pooled without complaint. It now also requires the
  family, estimator, mode, diagonal, data-shape flags and reference to
  match, and compares the coefficients -- which differ if and only if the
  data did. It also pools `$null_dist`, so `confint()` on a combined
  object uses every run's draws rather than the first run's while
  reporting the summed `reps`; and it weights by the permutations that
  actually contributed (`fit$n_valid`, new) rather than by the number
  requested, which over-weighted a run that lost permutations to
  non-convergence. A single fit passed alone is itself a list and used to
  be read as a list of models.

* **Mismatched node orderings are rejected.** `validate_qap_input()`
  compared `dim()` but never `dimnames()`, so matrices built from separate
  sources in different actor orders produced a silently meaningless model.
  Where both sides label a dimension, the labels must now agree; unlabelled
  data is unaffected.

* `confint()` and `plot()` explain why they are unavailable for
  multinomial fits instead of blaming `less_mem`.

## Performance

The permutation loop is roughly an order of magnitude faster:

| case | before | after |
|---|---|---|
| network matrix, gaussian, qapy | 4.00 ms/rep | 0.33 ms/rep |
| network matrix, gaussian, qapspp | 6.10 ms/rep | 0.90 ms/rep |
| CSS, gaussian, qapy | 11.40 ms/rep | 0.80 ms/rep |

Linear models with complete data now take a direct-solve path that
precomputes the invariant structure once instead of rebuilding the model
frame and refitting through the formula interface every permutation.
Set `options(MrQAP.disable_fast_ols = TRUE)` to force the general path;
the test suite uses this to check the two agree to 1e-12.

Semi-partialling gained a second speed-up on top of that. Only one column
of the design changes per permutation, and only that column's statistics
are ever read, so the fit is now obtained by Frisch-Waugh-Lovell from a
decomposition of the *fixed* columns formed once per predictor, rather
than by re-solving the whole design each time. Robust standard errors use
the matching leverage decomposition. Both are exact algebra, held to
`lm()`, to the full solver and to `HC3()` at 1e-10 in the test suite.

| case (n = 60, 3 predictors, 400 reps) | before | after |
|---|---|---|
| qapspp | 1.882 ms/rep | 0.480 ms/rep |
| qapspp, robust SEs | 2.928 ms/rep | 0.524 ms/rep |

which puts semi-partialling at parity with outcome permutation per
replication, having been three to six times dearer.

`summary()` was being computed three times per gaussian fit, and R-squared
on every permutation although only the baseline's is ever reported.
`array_to_vector()` and `make_css_data()` no longer rebuild invariant index
structures per call.

The GPU path (`use_gpu = TRUE`) now batches the semi-partialling case
instead of issuing one small GPU call per permutation, and ships only the
column that the permutation actually changes, expanding the rest of the
design on the device. Previously it assembled and transferred the whole
design matrix once per permutation, which made the GPU path *slower* than
the CPU one for semi-partialling.

Measured on an RTX 3070 against the CPU fast path, 500 permutations,
3 predictors:

| network size | qapy | qapspp |
|---|---|---|
| n = 40 (1.5k dyads) | 1.8x slower | 4.7x faster |
| n = 100 (9.9k dyads) | 3.4x faster | 4.4x faster |
| n = 200 (39.8k dyads) | 1.3x faster | 2.9x faster |

The GPU is not worth it on small networks, where the transfer overhead
exceeds the work. With CPU-only `torch` it is slower than the plain R fast
path throughout, since there is no device to amortise against.

Batch size adapts to the problem instead of sitting at a fixed 500, and is
tunable with `options(MrQAP.gpu_batch_mb = ...)` (default 512 MB). The
previous 64 MB internal default throttled throughput badly on a large card.

The tensor algebra is verified against the CPU solver in the test suite
whenever `torch` is installed, on both the CPU and CUDA devices, and each
run re-checks its first batch at runtime and stops if the two disagree.

## Housekeeping

* The package license is now GPL (>= 3). `DESCRIPTION` previously declared
  `MIT + file LICENSE` while the `LICENSE` file contained the GPL-3 text.
* `QAP()` and `probabilistic_confusion_matrix()` restore the caller's RNG
  state on exit instead of leaving `set.seed()` applied globally.
* `QAP()` only changes the active `future` plan when `ncores` is given;
  previously it forced `sequential` and discarded a plan you had set.
* **Progress reporting is now opt-in**, as **progressr** intends. The
  package emits a progressor; you switch reporting on by wrapping the call:
  `progressr::with_progress(QAP(...))`. It previously called
  `with_progress()` internally, which made progress compulsory and nested
  badly inside a user's own. Signals are batched at about 1% of the run
  rather than one condition throw per permutation.
* Continuous integration on Windows, macOS and Linux (release, devel and
  oldrel), plus test coverage reporting and a **pkgdown** site.
* `R CMD check --as-cran` is clean: 0 errors, 0 warnings, 0 notes.


# MrQAP 0.2.0

* Initial development version.
