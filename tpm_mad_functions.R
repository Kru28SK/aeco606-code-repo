# =============================================================================
# Transition Probability Matrix Estimation: Minimum Absolute Deviations (MAD)
# Method: Lee, Judge & Zellner (1970)
# Approach: Iteratively Reweighted Least Squares (IRLS) with simplex projection
# =============================================================================
# DATA FORMAT:
#   data.frame / matrix  |  rows = years  |  cols = crops  |  values = areas
# =============================================================================

# ------------------------------------------------------------------------------
# HELPER 1: Simplex projection (Duchi et al. 2008)
# ------------------------------------------------------------------------------
project_simplex <- function(v) {
  n   <- length(v)
  u   <- sort(v, decreasing = TRUE)
  css <- cumsum(u)
  rho <- max(which(u > (css - 1) / seq_len(n)))
  lam <- (1 - css[rho]) / rho
  pmax(v + lam, 0)
}

# ------------------------------------------------------------------------------
# HELPER 2: Sum of Absolute Deviations
# ------------------------------------------------------------------------------
compute_sad <- function(X, Y, P) sum(abs(Y - X %*% P))

# ------------------------------------------------------------------------------
# HELPER 3: Stationary distribution + full convergence profile
#   Returns k values at three practical thresholds AND the full decay table.
# ------------------------------------------------------------------------------
compute_stationary <- function(P, max_k = 200) {
  
  n <- nrow(P)
  
  # Power iteration to find pi
  v <- rep(1 / n, n)
  for (iter in seq_len(2000)) {
    v_new <- as.vector(v %*% P)
    v_new <- pmax(v_new, 0)
    v_new <- v_new / sum(v_new)
    if (max(abs(v_new - v)) < 1e-14) break
    v <- v_new
  }
  pi <- v_new
  
  # Verification
  pi_check  <- as.vector(matrix(pi, nrow = 1) %*% P)
  max_error <- max(abs(pi_check - pi))
  
  # Build full convergence profile: max row-deviation of P^k from pi
  thresholds  <- c(1e-2, 1e-3, 1e-4)          # practical tolerances
  k_converge  <- setNames(rep(NA_integer_, 3), c("1e-2", "1e-3", "1e-4"))
  profile_k   <- integer(0)
  profile_dev <- numeric(0)
  
  Pk <- P
  for (k in seq_len(max_k)) {
    dev <- max(apply(Pk, 1, function(row) max(abs(row - pi))))
    profile_k   <- c(profile_k,   k)
    profile_dev <- c(profile_dev, dev)
    
    for (nm in names(k_converge))
      if (is.na(k_converge[nm]) && dev < as.numeric(nm))
        k_converge[nm] <- k
    
    if (!is.na(k_converge["1e-4"])) break   # stop once tightest threshold met
    Pk <- Pk %*% P
  }
  
  list(pi          = pi,
       pi_check    = pi_check,
       max_error   = max_error,
       k_converge  = k_converge,          # named: "1e-2", "1e-3", "1e-4"
       profile     = data.frame(k = profile_k, max_dev = profile_dev))
}

# ------------------------------------------------------------------------------
# MAIN FUNCTION: estimate_tpm_mad()
# ------------------------------------------------------------------------------
estimate_tpm_mad <- function(data,
                             max_outer = 30,
                             max_inner = 500,
                             tol_outer = 1e-6,
                             tol_inner = 1e-8,
                             eps       = 1e-4,
                             max_k     = 200,
                             verbose   = TRUE) {
  
  mat        <- as.matrix(data)
  crop_names <- if (!is.null(colnames(mat))) colnames(mat) else paste0("Crop", seq_len(ncol(mat)))
  n_years    <- nrow(mat)
  n          <- ncol(mat)
  nT         <- n_years - 1
  
  if (nT < 2 || n < 2)
    stop("Need at least 3 years (2 transitions) and 2 crops.")
  
  S <- mat / rowSums(mat)
  X <- S[1:nT,       , drop = FALSE]
  Y <- S[2:(nT + 1), , drop = FALSE]
  
  P <- matrix(0.2 / (n - 1), nrow = n, ncol = n)
  diag(P) <- 0.8
  
  outer_iter <- 0
  
  for (outer in seq_len(max_outer)) {
    outer_iter <- outer
    P_old <- P
    
    W   <- 1 / pmax(abs(Y - X %*% P), eps)
    trQ <- sum(rowMeans(W) * rowSums(X^2))
    if (trQ == 0) trQ <- 1
    lr  <- 1.5 / trQ
    
    for (inner in seq_len(max_inner)) {
      grad     <- t(X) %*% (W * (X %*% P - Y))
      P_new    <- P - lr * grad
      max_diff <- 0
      for (i in seq_len(n)) {
        p_proj   <- project_simplex(P_new[i, ])
        max_diff <- max(max_diff, max(abs(p_proj - P[i, ])))
        P[i, ]   <- p_proj
      }
      if (max_diff < tol_inner) break
    }
    
    max_diff_outer <- max(abs(P - P_old))
    if (verbose)
      cat(sprintf("  Outer iter %2d | max|Delta P| = %.2e | SAD = %.6f\n",
                  outer, max_diff_outer, compute_sad(X, Y, P)))
    
    if (max_diff_outer < tol_outer) break
  }
  
  rownames(P) <- crop_names
  colnames(P) <- crop_names
  
  stat                 <- compute_stationary(P, max_k = max_k)
  names(stat$pi)       <- crop_names
  names(stat$pi_check) <- crop_names
  
  list(P          = P,
       SAD        = compute_sad(X, Y, P),
       outer_iter = outer_iter,
       stationary = stat,
       crop_names = crop_names)
}

# ------------------------------------------------------------------------------
# DISPLAY FUNCTION: print_tpm_results()
# ------------------------------------------------------------------------------
print_tpm_results <- function(res, digits_P = 4, digits_pi = 6) {
  
  sep  <- strrep("=", 65)
  sep2 <- strrep("-", 65)
  
  cat("\n", sep, "\n", sep = "")
  cat("  TRANSITION PROBABILITY MATRIX  (MAD Estimator)\n")
  cat("  Method: Lee, Judge & Zellner (1970)\n")
  cat(sep, "\n\n", sep = "")
  
  # 1. TPM
  cat("1. ESTIMATED TPM\n")
  cat("   P(i,j) = prob. of area moving from row-crop to col-crop\n")
  cat(sep2, "\n", sep = "")
  print(round(res$P, digits_P))
  cat("\n   Row sums (must = 1):\n")
  print(round(rowSums(res$P), 10))
  
  # 2. Fit summary
  cat("\n", sep2, "\n", sep = "")
  cat("2. FIT SUMMARY\n")
  cat(sep2, "\n", sep = "")
  cat(sprintf("   Sum of Absolute Deviations (SAD)  : %.10f\n", res$SAD))
  cat(sprintf("   IRLS outer iterations completed   : %d\n",    res$outer_iter))
  
  # 3. Stationary distribution
  cat("\n", sep2, "\n", sep = "")
  cat("3. STATIONARY (EQUILIBRIUM) DISTRIBUTION\n")
  cat("   Solves pi*P = pi, sum(pi) = 1  (power iteration)\n")
  cat("   Long-run equilibrium crop area shares\n")
  cat(sep2, "\n", sep = "")
  
  pi_df <- data.frame(
    Crop       = res$crop_names,
    Stationary = round(res$stationary$pi,                                  digits_pi),
    pi_times_P = round(res$stationary$pi_check,                            digits_pi),
    Difference = formatC(abs(res$stationary$pi - res$stationary$pi_check),
                         format = "e", digits = 3),
    check.names = FALSE
  )
  print(pi_df, row.names = FALSE)
  
  cat(sprintf("\n   Sum of pi               : %.10f\n", sum(res$stationary$pi)))
  max_e <- res$stationary$max_error
  cat(sprintf("   Max |pi*P - pi|         : %.3e  %s\n",
              max_e, if (is.finite(max_e) && max_e < 1e-8) "[PASS]" else "[CHECK]"))
  
  # 4. Convergence profile
  cat("\n", sep2, "\n", sep = "")
  cat("4. CONVERGENCE TO STATIONARY DISTRIBUTION\n")
  cat("   Criterion: max over all rows of  max|P^k[i,] - pi| < tol\n")
  cat(sep2, "\n", sep = "")
  
  kv <- res$stationary$k_converge
  cat(sprintf("   %-30s  k = %s\n", "Practical  (tol = 0.01)  :",
              if (!is.na(kv["1e-2"])) kv["1e-2"] else "> max_k"))
  cat(sprintf("   %-30s  k = %s\n", "Applied    (tol = 0.001) :",
              if (!is.na(kv["1e-3"])) kv["1e-3"] else "> max_k"))
  cat(sprintf("   %-30s  k = %s\n", "Precise    (tol = 0.0001):",
              if (!is.na(kv["1e-4"])) kv["1e-4"] else "> max_k"))
  
  cat("\n   Full decay profile (max row-deviation from pi):\n\n")
  prof <- res$stationary$profile
  # Print a compact table: k and deviation, flag threshold crossings
  cat(sprintf("   %4s   %12s   %s\n", "k", "Max deviation", ""))
  cat(sprintf("   %s\n", strrep("-", 35)))
  crossed <- c(`1e-2` = FALSE, `1e-3` = FALSE, `1e-4` = FALSE)
  for (i in seq_len(nrow(prof))) {
    flag <- ""
    if (!crossed["1e-2"] && prof$max_dev[i] < 1e-2) {
      flag <- "<- tol = 0.01  (practical)";  crossed["1e-2"] <- TRUE
    } else if (!crossed["1e-3"] && prof$max_dev[i] < 1e-3) {
      flag <- "<- tol = 0.001 (applied)";    crossed["1e-3"] <- TRUE
    } else if (!crossed["1e-4"] && prof$max_dev[i] < 1e-4) {
      flag <- "<- tol = 0.0001 (precise)";   crossed["1e-4"] <- TRUE
    }
    cat(sprintf("   %4d   %12.6f   %s\n", prof$k[i], prof$max_dev[i], flag))
  }
  
  cat("\n", sep, "\n\n", sep = "")
  invisible(res)
}