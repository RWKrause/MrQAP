# MrQAP

Multiple Regression Quadratic Assignment Procedure (MRQAP) for network data
in R.

Dyadic observations are not independent: the tie from *i* to *j* shares a
sender with the tie from *i* to *k*, and standard regression standard errors
are badly wrong as a result. MRQAP handles this by permuting the network —
relabelling nodes and refitting many times — and reading the p-value off the
resulting null distribution rather than off a t-table.

`MrQAP` provides one function, `QAP()`, for

- **network matrices** — ordinary one-mode socio-matrices,
- **cognitive social structures** — 3D arrays of `[sender, receiver, perceiver]`,
- **multi-mode data** — where each dimension is its own node set,

with gaussian, binomial, Poisson, negative binomial, zero-inflated Poisson
and multinomial families, optional fixed or random effects, prior weights
and offsets, and heteroskedasticity-consistent standard errors.

## Installation

```r
# install.packages("remotes")
remotes::install_github("RWKrause/MrQAP")
```

## Usage

```r
library(MrQAP)

set.seed(1)
n  <- 20
x1 <- matrix(rnorm(n^2), n, n)          # some dyadic covariate
x2 <- matrix(rnorm(n^2), n, n)          # another
y  <- 1 + 0.8 * x1 - 0.4 * x2 + matrix(rnorm(n^2), n, n)

fit <- QAP(y ~ x1 + x2,
           data = list(y = y, x1 = x1, x2 = x2),
           reps = 1000)
fit
```

The kind of data is detected from its shape — pass a 3D array and you get a
CSS model, pass a rectangular matrix and you get a two-mode model, with no
change to the call:

```r
# cognitive social structure: [sender, receiver, perceiver]
QAP(y ~ x1, data = list(y = css_array, x1 = pred_array), reps = 1000)

# two-mode: 40 people x 12 organisations
QAP(y ~ x1, data = list(y = affil, x1 = pred), reps = 1000)
```

A dyadic data frame can go straight in, as can `igraph` and `network`
objects:

```r
d <- data.frame(sender = ..., receiver = ..., weight = ..., shared = ...)
QAP(weight ~ shared, data = d, sender = "sender", receiver = "receiver")

QAP(y ~ x1, data = list(y = some_igraph, x1 = another_igraph))
```

## Working with a fit

```r
summary(fit)                # the report, with significance codes
summary(fit)$coefficients   # a data frame, one row per term
confint(fit)                # intervals from the permutation distribution
plot(fit)                   # the null distribution the p-value summarises
predict(fit, type = "matrix")   # fitted values, back in network shape
```

`coef()`, `vcov()`, `nobs()`, `fitted()`, `residuals()`, `logLik()`,
`AIC()`, `BIC()` and `as.data.frame()` all work as usual, and `tidy()` /
`glance()` are registered with **broom** if you have it.

See `vignette("mrqap")` for a full walkthrough.

## Which null hypothesis?

`nullhyp = "qapspp"` (the default) is Dekker's semi-partialling plus
procedure: each predictor is residualised on the others before being
permuted. `nullhyp = "qapy"` permutes the outcome instead. Semi-partialling
plus is the more reliable test when predictors are correlated with one
another — the usual situation — and is what you should normally use. It
costs `reps` permutations *per predictor*.

## Performance

The permutation loop is the whole cost of MRQAP, so it is optimised
accordingly. Linear models with complete data take a fast path that
precomputes everything invariant across permutations and solves the normal
equations directly, roughly 12x faster than the general path on network
matrices and 14x on CSS arrays. Results are identical either way — the test
suite checks the two paths against each other, and against a stored fixture,
at a tolerance of 1e-12.

Semi-partialling costs `reps` permutations per predictor, so it gets a
further saving: only one column of the design changes per permutation, and
only that column's statistics are ever read, so each fit is obtained by
Frisch-Waugh-Lovell from a decomposition of the columns that do not change.
That brings it to roughly the same cost per permutation as outcome
permutation, from three to six times dearer.

For larger jobs, `ncores` parallelises over permutations via the
[future](https://future.futureverse.org/) framework. If you leave `ncores`
unset, whatever `future` plan you have set yourself is used unchanged. To
watch a long run, wrap it: `progressr::with_progress(QAP(...))`.

`use_gpu = TRUE` offloads batched linear models to a GPU through
[torch](https://torch.mlverse.org/), which pays off from a few thousand
dyads upwards. Each run re-checks its first batch against the CPU solver
and stops if they disagree.

## Citation

```r
citation("MrQAP")
```

## License

GPL (>= 3)
