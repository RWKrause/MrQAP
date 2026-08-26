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
and multinomial families, optional fixed or random effects, and
heteroskedasticity-consistent standard errors.

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

Building the input from a dyadic data frame:

```r
d <- data.frame(sender = ..., receiver = ..., weight = ..., shared = ...)
mats <- df_to_mat(d, sender = "sender", receiver = "receiver")
QAP(weight ~ shared, data = mats, reps = 1000)
```

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

For larger jobs, `ncores` parallelises over permutations via the
[future](https://future.futureverse.org/) framework. If you leave `ncores`
unset, whatever `future` plan you have set yourself is used unchanged.

## Citation

```r
citation("MrQAP")
```

## License

GPL (>= 3)
