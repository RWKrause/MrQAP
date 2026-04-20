#' GPU-accelerated batch OLS for QAP permutations
#'
#' Uses the \pkg{torch} package to compute OLS regressions for many
#' permutations simultaneously on a GPU (or CPU via torch).  Only applicable
#' when \code{family = "gaussian"}, no random effects, no fixed effects, and
#' no comparisons.
#'
#' The approach pre-computes \eqn{M = (X'X)^{-1} X'}, then for \code{reps}
#' permutations constructs a matrix \eqn{Y_{perm}} (qapy) or swaps a single
#' predictor column (qapspp) and recomputes \eqn{M} per batch.  All betas and
#' t-values are obtained via batch matrix multiplications.
#'
#' @param data Named list of matrices (same as in \code{QAPglm}).
#' @param parsed Output of \code{parse_qap_formula()}.
#' @param mode Character; "digraph" or "graph".
#' @param diag Logical; include diagonal?
#' @param groups Vector or NULL; permutation groups.
#' @param reps Integer; number of permutations.
#' @param baseline_fit Baseline fit from \code{fit_qap_model()}.
#' @param perm_var Character or NULL; variable to permute for qapspp.
#'   When NULL, Y is permuted (qapy).
#' @param batch_size Integer; permutations per GPU batch (default 500).
#' @param device Character; "cuda" for GPU, "cpu" for CPU torch.
#'
#' @return A list with \code{lower}, \code{larger}, \code{abs} proportion
#'   matrices (same format as QAPglm output).
#' @keywords internal

gpu_batch_ols <- function(data, parsed, mode, diag, groups, reps,
                          baseline_fit, perm_var = NULL,
                          batch_size = 500, device = "cuda") {

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The 'torch' package is required for GPU acceleration. ",
         "Install it with: install.packages('torch')")
  }

  if (device == "cuda" && !torch::cuda_is_available()) {
    message("CUDA not available. Falling back to CPU torch.")
    device <- "cpu"
  }

  dep  <- parsed$dependent
  main <- parsed$main

  # Create baseline data frame once
  pred0 <- make_qap_data(y    = data[[dep]],
                          x    = data[main],
                          g    = groups,
                          diag = diag,
                          mode = mode,
                          net  = 1,
                          perm = FALSE,
                          xi   = NULL)

  y_vec <- pred0$yv
  X_mat <- cbind(1, as.matrix(pred0[, main, drop = FALSE]))
  n_obs <- nrow(X_mat)
  p     <- ncol(X_mat)

  base_coefs <- baseline_fit$coefficients
  base_t     <- baseline_fit$t

  # Accumulators
  lower_sum  <- rep(0, length(base_coefs) * 2)
  larger_sum <- rep(0, length(base_coefs) * 2)
  abs_sum    <- rep(0, length(base_coefs) * 2)

  dim_out <- c(2, length(base_coefs))

  if (is.null(perm_var)) {
    # ---- qapy: permute Y, X stays fixed ----
    # Pre-compute projection matrix once
    X_t      <- torch::torch_tensor(X_mat, dtype = torch::torch_float64(),
                                    device = device)
    XtX      <- torch::torch_mm(torch::torch_t(X_t), X_t)
    XtXinv   <- torch::torch_inverse(XtX)
    M        <- torch::torch_mm(XtXinv, torch::torch_t(X_t))
    XtXinv_diag <- torch::torch_diag(XtXinv)

    reps_done <- 0
    while (reps_done < reps) {
      current_batch <- min(batch_size, reps - reps_done)

      Y_batch <- matrix(NA_real_, nrow = n_obs, ncol = current_batch)
      for (j in seq_len(current_batch)) {
        y_perm <- RMPerm(data[[dep]], groups)
        perm_pred <- make_qap_data(y    = y_perm,
                                    x    = data[main],
                                    g    = groups,
                                    diag = diag,
                                    mode = mode,
                                    net  = 1,
                                    perm = FALSE,
                                    xi   = NULL)
        Y_batch[, j] <- perm_pred$yv
      }

      Y_t <- torch::torch_tensor(Y_batch, dtype = torch::torch_float64(),
                                  device = device)

      B   <- torch::torch_mm(M, Y_t)
      E   <- Y_t - torch::torch_mm(X_t, B)
      MSE <- torch::torch_sum(E^2, dim = 1) / (n_obs - p)
      SE  <- torch::torch_sqrt(torch::torch_ger(XtXinv_diag, MSE))
      T_vals <- B / SE

      B_cpu <- as.matrix(B$cpu())
      T_cpu <- as.matrix(T_vals$cpu())

      for (j in seq_len(current_batch)) {
        pres <- rbind(B_cpu[, j], T_cpu[, j])
        bres <- rbind(base_coefs, base_t)
        lower_sum  <- lower_sum  + as.vector(pres <= bres)
        larger_sum <- larger_sum + as.vector(pres >= bres)
        abs_sum    <- abs_sum    + as.vector(abs(pres) >= abs(bres))
      }

      reps_done <- reps_done + current_batch
    }

  } else {
    # ---- qapspp: permute one predictor, rebuild X each time ----
    # Y is fixed on the device
    y_t <- torch::torch_tensor(matrix(y_vec, ncol = 1),
                               dtype = torch::torch_float64(),
                               device = device)

    reps_done <- 0
    while (reps_done < reps) {
      current_batch <- min(batch_size, reps - reps_done)

      # Collect permuted betas and t-values
      B_batch <- matrix(NA_real_, nrow = p, ncol = current_batch)
      T_batch <- matrix(NA_real_, nrow = p, ncol = current_batch)

      for (j in seq_len(current_batch)) {
        d_perm <- data
        d_perm[[perm_var]] <- RMPerm(d_perm[[perm_var]], groups)

        perm_pred <- make_qap_data(y    = d_perm[[dep]],
                                    x    = d_perm[main],
                                    g    = groups,
                                    diag = diag,
                                    mode = mode,
                                    net  = 1,
                                    perm = FALSE,
                                    xi   = NULL)
        X_perm <- cbind(1, as.matrix(perm_pred[, main, drop = FALSE]))
        Xp_t   <- torch::torch_tensor(X_perm, dtype = torch::torch_float64(),
                                       device = device)

        XpXp      <- torch::torch_mm(torch::torch_t(Xp_t), Xp_t)
        XpXp_inv  <- torch::torch_inverse(XpXp)
        b_perm    <- torch::torch_mm(XpXp_inv,
                                     torch::torch_mm(torch::torch_t(Xp_t), y_t))
        e_perm    <- y_t - torch::torch_mm(Xp_t, b_perm)
        mse_perm  <- (torch::torch_sum(e_perm^2) / (n_obs - p))$item()
        XpXp_diag <- as.numeric(torch::torch_diag(XpXp_inv)$cpu())
        se_perm   <- sqrt(XpXp_diag * mse_perm)

        b_cpu <- as.numeric(b_perm$cpu())
        t_cpu <- b_cpu / se_perm

        B_batch[, j] <- b_cpu
        T_batch[, j] <- t_cpu
      }

      for (j in seq_len(current_batch)) {
        pres <- rbind(B_batch[, j], T_batch[, j])
        bres <- rbind(base_coefs, base_t)
        lower_sum  <- lower_sum  + as.vector(pres <= bres)
        larger_sum <- larger_sum + as.vector(pres >= bres)
        abs_sum    <- abs_sum    + as.vector(abs(pres) >= abs(bres))
      }

      reps_done <- reps_done + current_batch
    }
  }

  list(
    lower  = matrix(lower_sum / reps,  nrow = dim_out[1], ncol = dim_out[2],
                    dimnames = list(NULL, names(base_coefs))),
    larger = matrix(larger_sum / reps, nrow = dim_out[1], ncol = dim_out[2],
                    dimnames = list(NULL, names(base_coefs))),
    abs    = matrix(abs_sum / reps,    nrow = dim_out[1], ncol = dim_out[2],
                    dimnames = list(NULL, names(base_coefs)))
  )
}


#' GPU-accelerated batch OLS for CSS permutations
#'
#' CSS variant of \code{gpu_batch_ols}.  Handles 3D arrays and CSS-specific
#' permutation (\code{RMPerm(..., CSS = TRUE)}).
#'
#' @param data Named list of 3D arrays (same as in \code{QAPcss}).
#' @param parsed Output of \code{parse_qap_formula()}.
#' @param mode Character; "directed" or "undirected".
#' @param diag Logical; include diagonal?
#' @param groups Vector or NULL; permutation groups.
#' @param reps Integer; number of permutations.
#' @param baseline_fit Baseline fit from \code{fit_qap_model()}.
#' @param perm_var Character or NULL; variable to permute for qapspp.
#' @param batch_size Integer; permutations per GPU batch (default 500).
#' @param device Character; "cuda" for GPU, "cpu" for CPU torch.
#'
#' @return A list with \code{lower}, \code{larger}, \code{abs} proportion
#'   matrices.
#' @keywords internal

gpu_batch_ols_css <- function(data, parsed, mode, diag, groups, reps,
                              baseline_fit, perm_var = NULL,
                              batch_size = 500, device = "cuda") {

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The 'torch' package is required for GPU acceleration. ",
         "Install it with: install.packages('torch')")
  }

  if (device == "cuda" && !torch::cuda_is_available()) {
    message("CUDA not available. Falling back to CPU torch.")
    device <- "cpu"
  }

  dep  <- parsed$dependent
  main <- parsed$main
  data_vars <- parsed$all_data_vars

  # Create baseline data frame once
  x_list <- lapply(data_vars, function(v) data[[v]])
  names(x_list) <- data_vars
  cssd  <- make_css_data(y = data[[dep]], x = x_list,
                         nets = 1, diag = diag, mode = mode)
  pred0 <- cssd$pred

  y_vec <- pred0$yv
  X_mat <- cbind(1, as.matrix(pred0[, main, drop = FALSE]))
  n_obs <- nrow(X_mat)
  p     <- ncol(X_mat)

  base_coefs <- baseline_fit$coefficients
  base_t     <- baseline_fit$t

  # Accumulators
  lower_sum  <- rep(0, length(base_coefs) * 2)
  larger_sum <- rep(0, length(base_coefs) * 2)
  abs_sum    <- rep(0, length(base_coefs) * 2)

  dim_out <- c(2, length(base_coefs))

  # Helper: build vectorised data from (possibly permuted) CSS data
  build_css_pred <- function(d) {
    xl <- lapply(data_vars, function(v) d[[v]])
    names(xl) <- data_vars
    make_css_data(y = d[[dep]], x = xl,
                  nets = 1, diag = diag, mode = mode)$pred
  }

  if (is.null(perm_var)) {
    # ---- qapy: permute Y, X stays fixed ----
    X_t      <- torch::torch_tensor(X_mat, dtype = torch::torch_float64(),
                                    device = device)
    XtX      <- torch::torch_mm(torch::torch_t(X_t), X_t)
    XtXinv   <- torch::torch_inverse(XtX)
    M        <- torch::torch_mm(XtXinv, torch::torch_t(X_t))
    XtXinv_diag <- torch::torch_diag(XtXinv)

    reps_done <- 0
    while (reps_done < reps) {
      current_batch <- min(batch_size, reps - reps_done)

      Y_batch <- matrix(NA_real_, nrow = n_obs, ncol = current_batch)
      for (j in seq_len(current_batch)) {
        d_perm <- data
        d_perm[[dep]] <- RMPerm(d_perm[[dep]], groups, CSS = TRUE)
        perm_pred <- build_css_pred(d_perm)
        Y_batch[, j] <- perm_pred$yv
      }

      Y_t <- torch::torch_tensor(Y_batch, dtype = torch::torch_float64(),
                                  device = device)

      B   <- torch::torch_mm(M, Y_t)
      E   <- Y_t - torch::torch_mm(X_t, B)
      MSE <- torch::torch_sum(E^2, dim = 1) / (n_obs - p)
      SE  <- torch::torch_sqrt(torch::torch_ger(XtXinv_diag, MSE))
      T_vals <- B / SE

      B_cpu <- as.matrix(B$cpu())
      T_cpu <- as.matrix(T_vals$cpu())

      for (j in seq_len(current_batch)) {
        pres <- rbind(B_cpu[, j], T_cpu[, j])
        bres <- rbind(base_coefs, base_t)
        lower_sum  <- lower_sum  + as.vector(pres <= bres)
        larger_sum <- larger_sum + as.vector(pres >= bres)
        abs_sum    <- abs_sum    + as.vector(abs(pres) >= abs(bres))
      }

      reps_done <- reps_done + current_batch
    }

  } else {
    # ---- qapspp: permute one predictor, rebuild X each time ----
    y_t <- torch::torch_tensor(matrix(y_vec, ncol = 1),
                               dtype = torch::torch_float64(),
                               device = device)

    reps_done <- 0
    while (reps_done < reps) {
      current_batch <- min(batch_size, reps - reps_done)

      B_batch <- matrix(NA_real_, nrow = p, ncol = current_batch)
      T_batch <- matrix(NA_real_, nrow = p, ncol = current_batch)

      for (j in seq_len(current_batch)) {
        d_perm <- data
        d_perm[[perm_var]] <- RMPerm(d_perm[[perm_var]], groups, CSS = TRUE)

        perm_pred <- build_css_pred(d_perm)
        X_perm <- cbind(1, as.matrix(perm_pred[, main, drop = FALSE]))
        Xp_t   <- torch::torch_tensor(X_perm, dtype = torch::torch_float64(),
                                       device = device)

        XpXp      <- torch::torch_mm(torch::torch_t(Xp_t), Xp_t)
        XpXp_inv  <- torch::torch_inverse(XpXp)
        b_perm    <- torch::torch_mm(XpXp_inv,
                                     torch::torch_mm(torch::torch_t(Xp_t), y_t))
        e_perm    <- y_t - torch::torch_mm(Xp_t, b_perm)
        mse_perm  <- (torch::torch_sum(e_perm^2) / (n_obs - p))$item()
        XpXp_diag <- as.numeric(torch::torch_diag(XpXp_inv)$cpu())
        se_perm   <- sqrt(XpXp_diag * mse_perm)

        b_cpu <- as.numeric(b_perm$cpu())
        t_cpu <- b_cpu / se_perm

        B_batch[, j] <- b_cpu
        T_batch[, j] <- t_cpu
      }

      for (j in seq_len(current_batch)) {
        pres <- rbind(B_batch[, j], T_batch[, j])
        bres <- rbind(base_coefs, base_t)
        lower_sum  <- lower_sum  + as.vector(pres <= bres)
        larger_sum <- larger_sum + as.vector(pres >= bres)
        abs_sum    <- abs_sum    + as.vector(abs(pres) >= abs(bres))
      }

      reps_done <- reps_done + current_batch
    }
  }

  list(
    lower  = matrix(lower_sum / reps,  nrow = dim_out[1], ncol = dim_out[2],
                    dimnames = list(NULL, names(base_coefs))),
    larger = matrix(larger_sum / reps, nrow = dim_out[1], ncol = dim_out[2],
                    dimnames = list(NULL, names(base_coefs))),
    abs    = matrix(abs_sum / reps,    nrow = dim_out[1], ncol = dim_out[2],
                    dimnames = list(NULL, names(base_coefs)))
  )
}


#' Check GPU availability
#'
#' @return Logical; TRUE if torch and CUDA are available.
#' @export

gpu_available <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) return(FALSE)
  torch::cuda_is_available()
}
