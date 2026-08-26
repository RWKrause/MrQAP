# ============================================================
# GPU-accelerated batch OLS.
#
# The torch algebra lives in gpu_solve_fixed_x() and gpu_solve_varying_x(),
# deliberately kept as small functions over plain R arrays so they can be
# checked directly against the CPU solver -- see tests/testthat/test-gpu.R,
# which does exactly that whenever torch is installed.
#
# gpu_batch_ols() additionally re-checks its own first batch at runtime and
# stops loudly on disagreement, so a torch or driver change that breaks the
# translation surfaces immediately rather than as quietly wrong p-values.
# ============================================================

#' GPU-accelerated batch OLS for QAP permutations
#'
#' Solves many permuted OLS problems at once on a GPU via \pkg{torch}.  Only
#' applies to \code{family = "gaussian"} with no random effects, no fixed
#' effects and no comparisons.
#'
#' Under \code{qapy} the design matrix is fixed, so its projection is formed
#' once and each batch is a single matrix product.  Under \code{qapspp} one
#' column of the design changes every permutation, so a batch of designs is
#' stacked into a 3D tensor and solved with batched operations rather than
#' one small solve per permutation.
#'
#' @param tmpl Output of \code{qap_ols_template()}.
#' @param data Named list of matrices or arrays.
#' @param dep Character; dependent variable name.
#' @param main Character; predictor names.
#' @param groups Permutation grouping.
#' @param reps Integer; number of permutations.
#' @param baseline_fit Baseline fit from \code{fit_qap_model()}.
#' @param perm_var Character or NULL; the variable to permute (qapspp).  When
#'   NULL the outcome is permuted (qapy).
#' @param css,multi_mode Structure flags, passed to \code{perm_networks()}.
#' @param use_robust_errors Logical; must match the baseline fit.
#' @param batch_size Integer or NULL; permutations per batch.  \code{NULL}
#'   sizes the batch from the problem so one batch stays near 64 MB.
#' @param device Character; \code{"cuda"} or \code{"cpu"}.
#'
#' @return A list with \code{lower}, \code{larger}, \code{abs}.
#' @keywords internal

gpu_batch_ols <- function(tmpl, data, dep, main, groups, reps, baseline_fit,
                          perm_var = NULL, css = FALSE, multi_mode = FALSE,
                          use_robust_errors = FALSE,
                          batch_size = NULL, device = "cuda") {

  if (!requireNamespace("torch", quietly = TRUE))
    stop("The 'torch' package is required for GPU acceleration. ",
         "Install it with: install.packages('torch')")

  if (use_robust_errors)
    stop("use_gpu = TRUE does not support use_robust_errors; ",
         "run without the GPU for HC3 standard errors.")

  if (device == "cuda" && !torch::cuda_is_available()) {
    message("CUDA not available. Falling back to CPU torch.")
    device <- "cpu"
  }

  X       <- tmpl$X
  n_obs   <- tmpl$n_obs
  p       <- tmpl$p
  extract <- tmpl$extract

  if (is.null(batch_size)) batch_size <- gpu_batch_size(n_obs, p)

  base_coefs <- baseline_fit$coefficients
  base_t     <- baseline_fit$t
  bres       <- rbind(base_coefs, base_t)

  lower_sum <- larger_sum <- abs_sum <- rep(0, length(base_coefs) * 2)
  n_used    <- 0L
  checked   <- FALSE
  # Permuted statistics, retained so confint() can use the same draws the
  # p-values come from.
  draw_b <- matrix(NA_real_, nrow = reps, ncol = length(base_coefs),
                   dimnames = list(NULL, names(base_coefs)))
  draw_t <- draw_b

  target <- if (is.null(perm_var)) dep else perm_var
  obj    <- data[[target]]
  col    <- if (is.null(perm_var)) NA_integer_ else match(perm_var, colnames(X))

  done <- 0L
  while (done < reps) {
    nb <- min(batch_size, reps - done)
    Y <- NULL; Pc <- NULL

    if (is.null(perm_var)) {
      # ---- batch of permuted responses, design fixed ----
      Y <- vapply(seq_len(nb),
                  function(j) extract(perm_networks(obj, groups, CSS = css,
                                                    multi_mode = multi_mode)),
                  numeric(n_obs))
      sol <- gpu_solve_fixed_x(X, Y, device = device)

    } else {
      # ---- batch of permuted designs, response fixed ----
      Pc  <- gpu_perm_columns(nb, obj, extract, groups, css, multi_mode,
                              n_obs = n_obs)
      sol <- gpu_solve_varying_x(X, col, Pc, tmpl$y, device = device)
    }
    Bc <- sol$coefficients
    Tc <- sol$t

    # Cheap insurance: one batch per run is re-solved on the CPU and
    # compared, so a torch or driver change that breaks the translation
    # fails loudly instead of producing quietly wrong p-values.
    if (!checked) {
      ref <- if (is.null(perm_var)) {
        qap_ols_solve(X, Y[, 1])
      } else {
        Xp <- X; Xp[, col] <- Pc[, 1]
        qap_ols_solve(Xp, tmpl$y)
      }
      bad <- suppressWarnings(max(abs(Bc[, 1] - ref$coefficients),
                                  abs(Tc[, 1] - ref$t), na.rm = TRUE))
      if (!is.finite(bad) || bad > 1e-6)
        stop("GPU batch OLS disagrees with the CPU solver (max difference ",
             format(bad), "). This is a bug: please report it, and use ",
             "use_gpu = FALSE in the meantime.")
      checked <- TRUE
    }

    for (j in seq_len(nb)) {
      pres <- rbind(Bc[, j], Tc[, j])
      if (!all(is.finite(pres))) next
      n_used     <- n_used + 1L
      lower_sum  <- lower_sum  + as.vector(pres <= bres)
      larger_sum <- larger_sum + as.vector(pres >= bres)
      abs_sum    <- abs_sum    + as.vector(abs(pres) >= abs(bres))
      draw_b[n_used, ] <- pres[1, ]
      draw_t[n_used, ] <- pres[2, ]
    }
    done <- done + nb
  }

  if (n_used == 0L) stop("All permutations failed to converge.")
  if (n_used < reps)
    warning(reps - n_used, " of ", reps,
            " permutations failed and were excluded.")

  mk <- function(v) matrix(v / n_used, nrow = 2, ncol = length(base_coefs),
                           dimnames = list(qap_stat_rows(),
                                           names(base_coefs)))
  keep <- seq_len(n_used)
  list(lower = mk(lower_sum), larger = mk(larger_sum), abs = mk(abs_sum),
       draws = list(b = draw_b[keep, , drop = FALSE],
                    t = draw_t[keep, , drop = FALSE]),
       n_valid = n_used)
}


#' Batched OLS with a fixed design matrix
#'
#' Solves \code{ncol(Y)} regressions of the columns of \code{Y} on a single
#' \code{X}.  Because the design does not change, its projection is formed
#' once and the whole batch is one matrix product.
#'
#' Kept separate from \code{gpu_batch_ols()} so the torch algebra can be
#' checked directly against \code{qap_ols_solve()} -- see
#' \code{tests/testthat/test-gpu.R}.
#'
#' @param X Numeric matrix (n x p), including the intercept column.
#' @param Y Numeric matrix (n x B) of responses, one per permutation.
#' @param device Character; \code{"cuda"} or \code{"cpu"}.
#'
#' @return A list with \code{coefficients} and \code{t}, each \code{p x B}.
#' @keywords internal

gpu_solve_fixed_x <- function(X, Y, device = "cpu") {
  n_obs <- nrow(X)
  p     <- ncol(X)
  tt <- function(m) torch::torch_tensor(m, dtype = torch::torch_float64(),
                                        device = device)

  X_t         <- tt(X)
  XtX         <- torch::torch_mm(torch::torch_t(X_t), X_t)
  XtXinv      <- torch::torch_inverse(XtX)
  M           <- torch::torch_mm(XtXinv, torch::torch_t(X_t))
  XtXinv_diag <- torch::torch_diag(XtXinv)

  Y_t <- tt(Y)
  B   <- torch::torch_mm(M, Y_t)
  E   <- Y_t - torch::torch_mm(X_t, B)
  MSE <- torch::torch_sum(E^2, dim = 1) / (n_obs - p)
  SE  <- torch::torch_sqrt(torch::torch_ger(XtXinv_diag, MSE))

  list(coefficients = as.matrix(B$cpu()),
       t            = as.matrix((B / SE)$cpu()))
}


#' Batched OLS with a design matrix that changes per permutation
#'
#' Solves \code{dim(Xb)[1]} regressions of a single \code{y} on the stacked
#' designs in \code{Xb}, using batched tensor operations rather than one
#' small solve per permutation.
#'
#' Kept separate from \code{gpu_batch_ols()} so the torch algebra can be
#' checked directly against \code{qap_ols_solve()}.
#'
#' @param X Numeric matrix (n x p); the design, whose column \code{col} is
#'   replaced per permutation.
#' @param col Integer; index of the column that varies.
#' @param Pcols Numeric matrix (n x B); the permuted values of that column,
#'   one column per permutation.
#' @param y Numeric response vector of length n.
#' @param device Character; \code{"cuda"} or \code{"cpu"}.
#'
#' @return A list with \code{coefficients} and \code{t}, each \code{p x B}.
#' @keywords internal

gpu_solve_varying_x <- function(X, col, Pcols, y, device = "cpu") {
  nb    <- ncol(Pcols)
  n_obs <- nrow(X)
  p     <- ncol(X)
  tt <- function(m) torch::torch_tensor(m, dtype = torch::torch_float64(),
                                        device = device)

  # Only the varying column crosses the bus. Building the full B x n x p
  # stack on the CPU would copy (and transfer) p times more data than the
  # permutation actually changes.
  Xb_t <- tt(X)$unsqueeze(1L)$expand(c(nb, n_obs, p))$clone()
  Xb_t[, , col] <- tt(t(Pcols))
  y_t  <- tt(matrix(y, ncol = 1))

  Xt   <- Xb_t$transpose(2L, 3L)
  XtX  <- torch::torch_bmm(Xt, Xb_t)
  yb   <- y_t$unsqueeze(1L)$expand(c(nb, n_obs, 1L))
  Xty  <- torch::torch_bmm(Xt, yb)
  Bt   <- torch::linalg_solve(XtX, Xty)
  Et   <- yb - torch::torch_bmm(Xb_t, Bt)
  s2   <- torch::torch_sum(Et^2, dim = 2) / (n_obs - p)
  dinv <- torch::torch_diagonal(torch::linalg_inv(XtX), dim1 = 2L, dim2 = 3L)
  SE   <- torch::torch_sqrt(dinv * s2)

  list(coefficients = t(as.matrix(Bt$squeeze(3L)$cpu())),
       t            = t(as.matrix((Bt$squeeze(3L) / SE)$cpu())))
}


#' Build a batch of permuted predictor columns
#'
#' Only the permuted column is materialised: the rest of the design is
#' identical across permutations and is expanded on the device instead.
#' Pure R, so it is testable without \pkg{torch}.
#'
#' @param nb Integer; batch size.
#' @param obj The variable being permuted (matrix, array, or list thereof).
#' @param extract Vectoriser from \code{qap_ols_template()}.
#' @param groups,css,multi_mode Passed to \code{perm_networks()}.
#' @param n_obs Integer; length of the vectorised column.
#'
#' @return An \code{n_obs x nb} matrix.
#' @keywords internal

gpu_perm_columns <- function(nb, obj, extract, groups,
                             css = FALSE, multi_mode = FALSE, n_obs) {
  vapply(seq_len(nb),
         function(j) extract(perm_networks(obj, groups, CSS = css,
                                           multi_mode = multi_mode)),
         numeric(n_obs))
}


#' Choose a batch size that keeps one batch near a memory budget
#'
#' The qapspp path stacks \code{batch} copies of the design matrix, so batch
#' size trades GPU memory against throughput.  Too small and the transfer
#' overhead dominates: on an 8 GB card, a 200-node network ran 3.7x slower
#' at a 64 MB budget than at 1 GB.  Past roughly 1000 permutations per batch
#' the returns disappear, because the batch is then assembled on the CPU
#' faster than the GPU can be kept busy.
#'
#' Override with \code{options(MrQAP.gpu_batch_mb = ...)} if your card has
#' less memory to spare, or more.
#'
#' @param n_obs Integer; rows in the design matrix.
#' @param p Integer; columns in the design matrix.
#' @param budget_mb Numeric; target size of one batch, in megabytes.
#'
#' @return Integer batch size, at least 1 and at most 1000.
#' @keywords internal

gpu_batch_size <- function(n_obs, p,
                           budget_mb = getOption("MrQAP.gpu_batch_mb", 512)) {
  bytes_per_rep <- 8 * n_obs * p          # float64
  max(1L, min(1000L, as.integer(budget_mb * 1e6 / bytes_per_rep)))
}


#' Check GPU availability
#'
#' @return Logical; \code{TRUE} if both \pkg{torch} and CUDA are available.
#' @export
#'
#' @examples
#' gpu_available()

gpu_available <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) return(FALSE)
  isTRUE(try(torch::cuda_is_available(), silent = TRUE))
}
