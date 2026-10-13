## Submission

This is a new submission of MrQAP 0.3.0.

The package fits multiple regression quadratic assignment procedure (MRQAP)
models on network matrices, cognitive social structure arrays and multi-mode
data, with permutation-based inference.

## Test environments

* local Windows 11, R 4.6.1 -- 0 errors, 0 warnings, 1 note
* GitHub Actions: ubuntu-latest (R-devel, release, oldrel-1),
  macos-latest (release), windows-latest (release)
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Robert W. Krause <robert.w.krause@uky.edu>'

New submission
```

This note is the expected one for a first submission.

## Notes for the reviewer

* The two DOIs in the Description field were checked against the published
  records and resolve correctly. They are given in the
  `authors (year) <doi:...>` form.

* Long-running tests are wrapped in `testthat::skip_on_cran()`, so the
  check runs well inside the time limit on CRAN's machines. The full suite
  runs on every push in GitHub Actions, where `NOT_CRAN` is set.

* `torch` is suggested and used only for optional GPU-accelerated batch
  estimation. Every code path that touches it is guarded by
  `requireNamespace()` plus a check that a usable LibTorch back end is
  actually present, and the corresponding tests skip when it is not. The
  package is fully functional without it.

* All other suggested packages (`fixest`, `gmm`, `glmmTMB`, `lme4`, `nnet`,
  `pscl`, `sandwich`, `igraph`, `network`, `broom`) are used conditionally
  and their tests skip when they are unavailable.

## Downstream dependencies

There are none; this is a new package.
