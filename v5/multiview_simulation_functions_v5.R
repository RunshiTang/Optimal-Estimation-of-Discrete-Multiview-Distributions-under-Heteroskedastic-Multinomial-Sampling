
# ============================================================
# Multiview density estimation simulation helpers (v5)
#
# This script implements the helper functions needed for the
# simulation study in the multiview density estimation paper.
#
# Methods included:
#   1. histogram      : pooled histogram baseline
#   2. unscaled       : Algorithm 1 with M = 1
#   3. oracle         : Algorithm 1 with true oracle scaling M
#   4. oracle_cp      : Algorithm 1 with M estimated by CP decomposition
#                       of a pilot / pooled histogram using rTensor::cp()
#   5. slice_oracle   : Algorithm 1 with true slice normalization
#   6. slice_est      : Algorithm 1 with estimated slice normalization
#
# Dependency note:
#   - All methods except "oracle_cp" use only base R / stats.
#   - Method "oracle_cp" requires the CRAN package rTensor.
#
# Suggested usage:
#   source("multiview_simulation_functions_v5.R")
#
# v5 update:
#   The six estimation methods can now be run with either
#   algorithm = "no_thinning" or algorithm = "thinning".
# ============================================================

# -----------------------------
# Basic utilities
# -----------------------------

# Function: %||%
# Purpose : Return y if x is NULL; otherwise return x.
# Inputs  : x = any R object; y = fallback value.
# Output  : x if x is not NULL, else y.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Function: assert_array
# Purpose : Stop if the input is not an array.
# Inputs  : X = object to check; name = label used in the error message.
# Output  : Invisibly returns TRUE if X is an array.
assert_array <- function(X, name = deparse(substitute(X))) {
  if (is.null(dim(X))) stop(name, " must be an array.")
  invisible(TRUE)
}

# Function: frobenius_norm
# Purpose : Compute the Frobenius norm of an array.
# Inputs  : X = numeric array.
# Output  : Nonnegative scalar equal to sqrt(sum(X^2)).
frobenius_norm <- function(X) {
  sqrt(sum(X^2))
}

# Function: tensor_l1_norm
# Purpose : Compute the entrywise l1 norm of an array.
# Inputs  : X = numeric array.
# Output  : Nonnegative scalar equal to sum(abs(X)).
tensor_l1_norm <- function(X) {
  sum(abs(X))
}

# Function: tensor_l1_error
# Purpose : Compute entrywise l1 error between two arrays.
# Inputs  : X, Y = numeric arrays with the same dimensions.
# Output  : Nonnegative scalar equal to sum(abs(X - Y)).
tensor_l1_error <- function(X, Y) {
  sum(abs(X - Y))
}

# Function: tensor_fro_error
# Purpose : Compute Frobenius error between two arrays.
# Inputs  : X, Y = numeric arrays with the same dimensions.
# Output  : Nonnegative scalar equal to sqrt(sum((X - Y)^2)).
tensor_fro_error <- function(X, Y) {
  sqrt(sum((X - Y)^2))
}

# Function: safe_divide
# Purpose : Divide A by B with denominator lower-bounded by eps.
# Inputs  : A, B = numeric arrays / scalars; eps = small positive constant.
# Output  : Array of the same shape as A / B with denominator clipped below by eps.
safe_divide <- function(A, B, eps = 1e-15) {
  A / pmax(B, eps)
}

# Function: project_to_simplex
# Purpose : Project a vector onto the probability simplex.
# Inputs  : x = numeric vector.
# Output  : Numeric vector with nonnegative entries summing to 1.
project_to_simplex <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  u <- sort(x, decreasing = TRUE)
  cssv <- cumsum(u) - 1
  rho <- max(which(u - cssv / seq_len(n) > 0))
  theta <- cssv[rho] / rho
  pmax(x - theta, 0)
}

# Function: project_tensor_to_probability
# Purpose : Project all entries of a tensor onto the global probability simplex.
# Inputs  : X = numeric array.
# Output  : Numeric array with the same dimensions as X, nonnegative and summing to 1.
project_tensor_to_probability <- function(X) {
  out <- project_to_simplex(as.numeric(X))
  array(out, dim = dim(X))
}

# Function: normalize_probability_tensor
# Purpose : Normalize a nonnegative array so its entries sum to 1.
# Inputs  : X = numeric array; eps = threshold used to avoid division by zero.
# Output  : Probability tensor with the same dimensions as X.
normalize_probability_tensor <- function(X, eps = 1e-15) {
  X <- pmax(X, 0)
  s <- sum(X)
  if (s <= eps) {
    X[] <- 1 / length(X)
    return(X)
  }
  X / s
}

# Function: normalize_probability_vector
# Purpose : Normalize a vector into a probability vector using a positive-part / abs fallback.
# Inputs  : x = numeric vector; eps = threshold used to avoid division by zero.
# Output  : Numeric vector with nonnegative entries summing to 1.
normalize_probability_vector <- function(x, eps = 1e-15) {
  x <- as.numeric(x)
  y <- pmax(x, 0)
  if (sum(y) <= eps) y <- abs(x)
  if (sum(y) <= eps) y <- rep(1, length(x))
  y / sum(y)
}

# -----------------------------
# Tensor algebra
# -----------------------------

# Function: unfold_tensor
# Purpose : Matricize a tensor along a given mode.
# Inputs  : X = numeric array; mode = integer in {1, ..., length(dim(X))}.
# Output  : Matrix whose rows correspond to the specified mode.
unfold_tensor <- function(X, mode) {
  assert_array(X)
  dims <- dim(X)
  d <- length(dims)
  if (mode < 1 || mode > d) stop("Invalid mode.")
  perm <- c(mode, setdiff(seq_len(d), mode))
  Xp <- aperm(X, perm)
  matrix(Xp, nrow = dims[mode], ncol = prod(dims[-mode]))
}

# Function: fold_tensor
# Purpose : Refold a mode-unfolded matrix back into a tensor.
# Inputs  : M = matrix; dims = target tensor dimensions; mode = unfolding mode.
# Output  : Numeric array with dimensions given by dims.
fold_tensor <- function(M, dims, mode) {
  d <- length(dims)
  if (mode < 1 || mode > d) stop("Invalid mode.")
  perm <- c(mode, setdiff(seq_len(d), mode))
  tmp_dims <- dims[perm]
  Xp <- array(M, dim = tmp_dims)
  aperm(Xp, order(perm))
}

# Function: tensor_times_matrix
# Purpose : Multiply a tensor by a matrix along one mode.
# Inputs  : X = tensor; A = matrix of size m x dim(X)[mode]; mode = multiplication mode.
# Output  : Tensor with the same modes as X except mode replaced by m.
tensor_times_matrix <- function(X, A, mode) {
  assert_array(X)
  dims <- dim(X)
  if (ncol(A) != dims[mode]) {
    stop("Dimension mismatch in tensor_times_matrix().")
  }
  M <- unfold_tensor(X, mode)
  AM <- A %*% M
  new_dims <- dims
  new_dims[mode] <- nrow(A)
  fold_tensor(AM, new_dims, mode)
}

# Function: rank1_tensor
# Purpose : Build an outer-product tensor from a list of vectors.
# Inputs  : v_list = list of numeric vectors, one per mode.
# Output  : Rank-1 tensor whose (i1,...,id) entry is product_k v_list[[k]][ik].
rank1_tensor <- function(v_list) {
  d <- length(v_list)
  if (d < 1) stop("v_list must be non-empty.")

  dims <- vapply(v_list, length, integer(1))
  out <- as.numeric(v_list[[1]])

  if (d >= 2) {
    for (k in 2:d) {
      out <- as.vector(outer(out, as.numeric(v_list[[k]])))
    }
  }

  array(out, dim = dims)
}

# -----------------------------
# Tensor norms / marginals
# -----------------------------

# Function: slice_marginals
# Purpose : Compute the slice sums along one mode.
# Inputs  : X = numeric tensor; mode = target mode.
# Output  : Numeric vector of length dim(X)[mode].
slice_marginals <- function(X, mode) {
  rowSums(unfold_tensor(X, mode))
}

# Function: fiber_l1_by_mode
# Purpose : Compute l1 norms of fibers for a fixed mode.
# Inputs  : X = numeric tensor; mode = target mode.
# Output  : Numeric vector containing l1 norms of all mode-specific fibers.
fiber_l1_by_mode <- function(X, mode) {
  colSums(abs(unfold_tensor(X, mode)))
}

# Function: slice_l1_by_mode
# Purpose : Compute l1 norms of slices for a fixed mode.
# Inputs  : X = numeric tensor; mode = target mode.
# Output  : Numeric vector of slice l1 norms.
slice_l1_by_mode <- function(X, mode) {
  rowSums(abs(unfold_tensor(X, mode)))
}

# Function: fiber_l1_max
# Purpose : Compute the maximum fiber l1 norm over all modes.
# Inputs  : X = numeric tensor.
# Output  : Nonnegative scalar.
fiber_l1_max <- function(X) {
  max(vapply(seq_along(dim(X)), function(k) max(fiber_l1_by_mode(X, k)), numeric(1)))
}

# Function: slice_l1_max
# Purpose : Compute the maximum slice l1 norm over all modes.
# Inputs  : X = numeric tensor.
# Output  : Nonnegative scalar.
slice_l1_max <- function(X) {
  max(vapply(seq_along(dim(X)), function(k) max(slice_l1_by_mode(X, k)), numeric(1)))
}

# Function: summarize_tensor_signal
# Purpose : Summarize basic signal-size quantities of a tensor.
# Inputs  : X = numeric tensor.
# Output  : List with dimensions, l1 norm, Frobenius norm, max entry,
#           maximum fiber l1 norm, and maximum slice l1 norm.
summarize_tensor_signal <- function(X) {
  list(
    dims = dim(X),
    l1 = tensor_l1_norm(X),
    fro = frobenius_norm(X),
    max_entry = max(X),
    fiber_l1_max = fiber_l1_max(X),
    slice_l1_max = slice_l1_max(X)
  )
}

# -----------------------------
# Data generation
# -----------------------------

# Function: rdirichlet1
# Purpose : Draw one Dirichlet random vector.
# Inputs  : alpha = positive concentration vector.
# Output  : Probability vector of the same length as alpha.
rdirichlet1 <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  if (sum(x) <= 0) x[] <- 1
  x / sum(x)
}

# Function: generate_weights
# Purpose : Create a vector of CP mixture weights.
# Inputs  : R = rank; type = "balanced", "geometric", or "custom";
#           values = user-supplied weights if type = "custom";
#           ratio = common ratio if type = "geometric".
# Output  : Probability vector of length R.
generate_weights <- function(R,
                             type = c("balanced", "geometric", "custom"),
                             values = NULL,
                             ratio = 0.5) {
  type <- match.arg(type)
  if (type == "balanced") {
    w <- rep(1 / R, R)
  } else if (type == "geometric") {
    w <- ratio^(seq_len(R) - 1)
    w <- w / sum(w)
  } else {
    if (is.null(values)) stop("For custom weights, provide 'values'.")
    if (length(values) != R) stop("'values' must have length R.")
    if (any(values < 0)) stop("Weights must be nonnegative.")
    w <- values / sum(values)
  }
  w
}

# Function: generate_dense_factors
# Purpose : Generate dense probability factors with Dirichlet columns.
# Inputs  : p_vec = vector of mode dimensions; R = rank; alpha = Dirichlet parameter.
# Output  : List of factor matrices; the k-th matrix has dimension p_vec[k] x R.
generate_dense_factors <- function(p_vec, R, alpha = 1) {
  lapply(p_vec, function(p) {
    out <- sapply(seq_len(R), function(r) rdirichlet1(rep(alpha, p)))
    if (R == 1) out <- matrix(out, ncol = 1)
    out
  })
}

# Function: generate_sparse_factor
# Purpose : Generate one sparse probability vector.
# Inputs  : p = vector length; s = support size; alpha = Dirichlet parameter on support.
# Output  : Probability vector of length p with at most s nonzero entries.
generate_sparse_factor <- function(p, s, alpha = 0.2) {
  if (s > p) stop("s must be <= p.")
  idx <- sort(sample.int(p, s))
  v <- numeric(p)
  v[idx] <- rdirichlet1(rep(alpha, s))
  v
}

# Function: generate_sparse_factors
# Purpose : Generate sparse factor matrices for all modes.
# Inputs  : p_vec = vector of mode dimensions; R = rank; s_vec = support sizes;
#           alpha = Dirichlet parameter on each support.
# Output  : List of sparse factor matrices; the k-th matrix has dimension p_vec[k] x R.
generate_sparse_factors <- function(p_vec, R, s_vec = NULL, alpha = 0.2) {
  s_vec <- s_vec %||% pmax(2, ceiling(p_vec / 5))
  if (length(s_vec) == 1) s_vec <- rep(s_vec, length(p_vec))
  Map(function(p, s) {
    out <- sapply(seq_len(R), function(r) generate_sparse_factor(p, s, alpha))
    if (R == 1) out <- matrix(out, ncol = 1)
    out
  }, p = p_vec, s = s_vec)
}

# Function: build_cp_tensor
# Purpose : Construct a CP probability tensor from weights and factor matrices.
# Inputs  : weights = probability vector of length R;
#           factors = list of mode-specific factor matrices with common column count R.
# Output  : Probability tensor with dimensions equal to the row sizes of factors.
build_cp_tensor <- function(weights, factors) {
  d <- length(factors)
  R <- ncol(factors[[1]])
  dims <- vapply(factors, nrow, integer(1))
  if (length(weights) != R) stop("weights length must equal number of columns in factors.")

  P <- array(0, dim = dims)
  for (r in seq_len(R)) {
    vecs <- lapply(factors, function(A) A[, r])
    P <- P + weights[r] * rank1_tensor(vecs)
  }
  normalize_probability_tensor(P)
}

# Function: generate_multiview_tensor
# Purpose : Generate one low-rank multiview probability tensor and its factors.
# Inputs  : p_vec = vector of dimensions;
#           R = CP rank;
#           factor_type = "dense" or "sparse";
#           alpha = Dirichlet parameter;
#           s_vec = support sizes for sparse factors;
#           weight_type = "balanced", "geometric", or "custom";
#           weight_values = custom weights if needed;
#           geometric_ratio = ratio for geometric weights.
# Output  : List containing P, weights, factors, and factor_type.
generate_multiview_tensor <- function(p_vec,
                                      R,
                                      factor_type = c("dense", "sparse"),
                                      alpha = 1,
                                      s_vec = NULL,
                                      weight_type = c("balanced", "geometric", "custom"),
                                      weight_values = NULL,
                                      geometric_ratio = 0.5) {
  factor_type <- match.arg(factor_type)
  weight_type <- match.arg(weight_type)

  weights <- generate_weights(
    R = R,
    type = weight_type,
    values = weight_values,
    ratio = geometric_ratio
  )

  factors <- if (factor_type == "dense") {
    generate_dense_factors(p_vec, R, alpha = alpha)
  } else {
    generate_sparse_factors(p_vec, R, s_vec = s_vec, alpha = alpha)
  }

  P <- build_cp_tensor(weights, factors)

  list(
    P = P,
    weights = weights,
    factors = factors,
    factor_type = factor_type
  )
}

# Function: rmultinomial_tensor
# Purpose : Draw one multinomial count tensor with cell probabilities P.
# Inputs  : n = total sample size; P = probability tensor.
# Output  : Count tensor with the same dimensions as P and total count n.
rmultinomial_tensor <- function(n, P) {
  assert_array(P)
  counts <- rmultinom(1, size = n, prob = c(P))
  array(as.numeric(counts[, 1]), dim = dim(P))
}

# Function: simulate_histograms
# Purpose : Generate multiple independent multinomial histograms from P.
# Inputs  : P = probability tensor; n = sample size in each split; nsplit = number of histograms.
# Output  : List of length nsplit containing count tensors.
simulate_histograms <- function(P, n, nsplit = 3) {
  replicate(nsplit, rmultinomial_tensor(n, P), simplify = FALSE)
}

# Function: simulate_one_replication
# Purpose : Generate one full simulation replication.
# Inputs  : p_vec, R, n, factor_type, alpha, s_vec, weight_type,
#           weight_values, geometric_ratio = data-generating parameters;
#           nsplit = legacy argument; ignored by the no-thinning generator;
#           pilot_n = optional pilot histogram sample size for estimating M.
# Output  : List containing the true tensor P, weights, factors, Y/Y_list,
#           pilot_Y (possibly NULL), and n.
simulate_one_replication <- function(p_vec,
                                     R,
                                     n,
                                     factor_type = c("dense", "sparse"),
                                     alpha = 1,
                                     s_vec = NULL,
                                     weight_type = c("balanced", "geometric", "custom"),
                                     weight_values = NULL,
                                     geometric_ratio = 0.5,
                                     nsplit = 1,
                                     pilot_n = NULL) {
  model <- generate_multiview_tensor(
    p_vec = p_vec,
    R = R,
    factor_type = factor_type,
    alpha = alpha,
    s_vec = s_vec,
    weight_type = weight_type,
    weight_values = weight_values,
    geometric_ratio = geometric_ratio
  )

  # No multinomial thinning: draw one histogram tensor and use it in all
  # estimation stages.  The legacy name Y_list is retained only so older
  # wrappers keep working; it now stores the same single tensor Y.
  Y <- rmultinomial_tensor(n, model$P)
  pilot_Y <- if (!is.null(pilot_n) && pilot_n > 0) rmultinomial_tensor(pilot_n, model$P) else NULL

  c(model, list(Y = Y, Y_list = Y, pilot_Y = pilot_Y, n = n))
}

# -----------------------------
# Scaling tensors
# -----------------------------

# Function: make_identity_scaling
# Purpose : Create the all-ones scaling tensor.
# Inputs  : dims = tensor dimensions.
# Output  : Array of ones with dimensions dims.
make_identity_scaling <- function(dims) {
  array(1, dim = dims)
}

# Function: make_oracle_scaling_from_factors
# Purpose : Construct the oracle scaling tensor M from factor matrices.
# Inputs  : factors = list of factor matrices;
#           gamma = optional list of per-mode weights used in sigma_k;
#           eps = small positive constant for numerical stability.
# Output  : Scaling tensor M with the same dimensions as the target tensor.
make_oracle_scaling_from_factors <- function(factors, gamma = NULL, eps = 1e-15) {
  d <- length(factors)
  p_vec <- vapply(factors, nrow, integer(1))

  if (is.null(gamma)) {
    gamma <- lapply(p_vec, function(p) rep(1, p))
  }
  if (length(gamma) != d) stop("gamma must be a list with one vector per mode.")

  b_list <- vector("list", d)
  for (k in seq_len(d)) {
    sigma_k <- gamma[[k]] * apply(factors[[k]], 1, max)
    sigma_k <- pmax(sigma_k, eps)
    b_list[[k]] <- sqrt(sum(sigma_k) / sigma_k)
  }

  rank1_tensor(b_list)
}

# Function: make_slice_scaling_from_tensor
# Purpose : Construct the slice-normalization scaling tensor from the true tensor.
# Inputs  : P = probability tensor; eps = small positive constant.
# Output  : Scaling tensor with the same dimensions as P.
make_slice_scaling_from_tensor <- function(P, eps = 1e-15) {
  dims <- dim(P)
  d <- length(dims)
  b_list <- vector("list", d)
  for (k in seq_len(d)) {
    sl <- slice_marginals(P, k)
    b_list[[k]] <- 1 / sqrt(pmax(sl, 1 / dims[k], eps))
  }
  rank1_tensor(b_list)
}

# Function: make_slice_scaling_from_histogram
# Purpose : Construct the slice-normalization scaling tensor from a histogram.
# Inputs  : Y = count tensor or nonnegative tensor estimate; eps = small positive constant.
# Output  : Scaling tensor with the same dimensions as Y.
make_slice_scaling_from_histogram <- function(Y, eps = 1e-15) {
  make_slice_scaling_from_tensor(normalize_probability_tensor(Y), eps = eps)
}

# -----------------------------
# CP-decomposition-based estimated oracle scaling
# -----------------------------

# Function: require_rTensor
# Purpose : Stop with an informative message if the rTensor package is unavailable.
# Inputs  : None.
# Output  : Invisibly returns TRUE if rTensor is installed.
require_rTensor <- function() {
  if (!requireNamespace("rTensor", quietly = TRUE)) {
    stop(
      paste(
        "Method 'oracle_cp' requires the CRAN package 'rTensor'.",
        "Please run install.packages('rTensor') and try again."
      )
    )
  }
  invisible(TRUE)
}

# Function: cp_extract_probability_factors
# Purpose : Convert rTensor::cp() output into per-mode probability factor matrices.
# Inputs  : cp_fit = object returned by rTensor::cp(); eps = small positive constant.
# Output  : List of factor matrices with nonnegative columns summing to 1.
cp_extract_probability_factors <- function(cp_fit, eps = 1e-15) {
  U_list <- cp_fit$U
  lapply(U_list, function(Uk) {
    out <- apply(Uk, 2, function(v) normalize_probability_vector(v, eps = eps))
    if (is.null(dim(out))) out <- matrix(out, ncol = 1)
    out
  })
}

# Function: cp_estimated_tensor_to_array
# Purpose : Safely extract the estimated tensor array from rTensor::cp() output.
# Inputs  : cp_fit = object returned by rTensor::cp().
# Output  : Numeric array with the same dimensions as the decomposed tensor.
cp_estimated_tensor_to_array <- function(cp_fit) {
  est <- cp_fit$est
  if (methods::is(est, "Tensor")) {
    return(est@data)
  }
  if (is.array(est)) return(est)
  stop("Could not extract the estimated tensor from cp_fit$est.")
}

# Function: fit_cp_with_rtensor
# Purpose : Run CP decomposition on a tensor using rTensor::cp().
# Inputs  : X = nonnegative tensor (counts or probabilities);
#           rank = target CP rank;
#           cp_max_iter = maximum ALS iterations;
#           cp_tol = relative Frobenius tolerance.
# Output  : List with cp_fit, factors_est, and P_est_cp.
fit_cp_with_rtensor <- function(X,
                                rank,
                                cp_max_iter = 50,
                                cp_tol = 1e-5) {
  require_rTensor()
  X_prob <- normalize_probability_tensor(X)
  tnsr <- rTensor::as.tensor(X_prob)
  cp_fit <- rTensor::cp(
    tnsr,
    num_components = rank,
    max_iter = cp_max_iter,
    tol = cp_tol
  )
  factors_est <- cp_extract_probability_factors(cp_fit)
  P_est_cp <- normalize_probability_tensor(cp_estimated_tensor_to_array(cp_fit))
  list(
    cp_fit = cp_fit,
    factors_est = factors_est,
    P_est_cp = P_est_cp
  )
}

# Function: make_oracle_scaling_from_cp_tensor
# Purpose : Estimate oracle scaling M by CP decomposition of a tensor estimate.
# Inputs  : X = count tensor or nonnegative pilot tensor;
#           rank = target CP rank;
#           cp_max_iter = maximum ALS iterations;
#           cp_tol = relative Frobenius tolerance;
#           eps = numerical-stability constant passed to make_oracle_scaling_from_factors().
# Output  : List with M, cp_fit, factors_est, and P_est_cp.
make_oracle_scaling_from_cp_tensor <- function(X,
                                               rank,
                                               cp_max_iter = 50,
                                               cp_tol = 1e-5,
                                               eps = 1e-15) {
  cp_obj <- fit_cp_with_rtensor(
    X = X,
    rank = rank,
    cp_max_iter = cp_max_iter,
    cp_tol = cp_tol
  )
  M <- make_oracle_scaling_from_factors(
    factors = cp_obj$factors_est,
    eps = eps
  )
  c(list(M = M), cp_obj)
}

# -----------------------------
# Initializers for Algorithm 1
# -----------------------------

# Function: .make_diag_matrix
# Purpose : Build a diagonal matrix, preserving the 1x1 case.
# Inputs  : x = numeric scalar or vector.
# Output  : Square diagonal matrix.
.make_diag_matrix <- function(x) {
  if (length(x) == 1L) {
    matrix(x, nrow = 1, ncol = 1)
  } else {
    diag(x, nrow = length(x), ncol = length(x))
  }
}

# Function: top_eigenvectors
# Purpose : Compute the top r eigenvectors of a symmetric matrix.
# Inputs  : G = square matrix; r = requested rank.
# Output  : Matrix with orthonormal columns.
top_eigenvectors <- function(G, r) {
  ee <- eigen((G + t(G)) / 2, symmetric = TRUE)
  r_eff <- min(r, ncol(ee$vectors))
  U <- ee$vectors[, seq_len(r_eff), drop = FALSE]
  qr.Q(qr(U))
}

# Function: heteropca
# Purpose : Simple iterative HeteroPCA-style initializer.
# Inputs  : G = square Gram matrix; r = target rank; niter = number of iterations.
# Output  : Matrix with r orthonormal columns.
heteropca <- function(G, r, niter = 15) {
  G <- (G + t(G)) / 2
  diag(G) <- 0
  G_imp <- G

  for (it in seq_len(niter)) {
    ee <- eigen((G_imp + t(G_imp)) / 2, symmetric = TRUE)
    r_eff <- min(r, length(ee$values))
    U <- ee$vectors[, seq_len(r_eff), drop = FALSE]
    lam <- pmax(ee$values[seq_len(r_eff)], 0)
    G_low <- U %*% .make_diag_matrix(lam) %*% t(U)
    diag(G_imp) <- diag(G_low)
  }

  qr.Q(qr(U))
}

# Function: deflated_heteropca
# Purpose : Practical deflation-based variant of the HeteroPCA initializer.
# Inputs  : G = square Gram matrix; r = target rank; niter = number of iterations.
# Output  : Matrix with r orthonormal columns.
deflated_heteropca <- function(G, r, niter = 15) {
  G_resid <- (G + t(G)) / 2
  diag(G_resid) <- 0
  p <- nrow(G_resid)
  U <- matrix(0, nrow = p, ncol = r)

  for (j in seq_len(r)) {
    u_j <- heteropca(G_resid, r = 1, niter = niter)
    if (j > 1) {
      U_tmp <- qr.Q(qr(cbind(U[, seq_len(j - 1), drop = FALSE], u_j)))
      u_j <- U_tmp[, j, drop = FALSE]
    }
    lam_j <- as.numeric(t(u_j) %*% G_resid %*% u_j)
    U[, j] <- as.numeric(u_j)
    G_resid <- G_resid - lam_j * tcrossprod(u_j)
    diag(G_resid) <- 0
  }

  qr.Q(qr(U))
}

# Function: top_left_singular_vectors
# Purpose : Compute the top left singular vectors of a matrix.
# Inputs  : M = matrix; r = requested rank.
# Output  : Matrix with r orthonormal columns.
top_left_singular_vectors <- function(M, r) {
  r_eff <- min(r, nrow(M), ncol(M))
  if (r_eff < 1) stop("Requested rank is invalid.")
  sv <- svd(M, nu = r_eff, nv = 0)
  U <- sv$u[, seq_len(r_eff), drop = FALSE]
  qr.Q(qr(U))
}

# -----------------------------
# Algorithm 1 estimator
# -----------------------------

# Function: .as_single_histogram_tensor
# Purpose : Coerce an input into one count tensor for the no-thinning algorithm.
# Inputs  : Y = either one count tensor, a length-one list containing a tensor,
#           or a legacy list of split histograms.
# Output  : A single count tensor.  Legacy split lists are pooled.
.as_single_histogram_tensor <- function(Y) {
  if (is.list(Y)) {
    if (length(Y) < 1) stop("Y must contain at least one tensor.")
    if (length(Y) == 1) {
      Y <- Y[[1]]
    } else {
      # Backward compatibility: if old code passes split histograms, use the
      # pooled histogram rather than separate samples for separate stages.
      Y <- Reduce(`+`, Y)
    }
  }
  assert_array(Y, "Y")
  Y
}

# Function: multinomial_thin_tensor
# Purpose : Split one multinomial histogram tensor into three conditionally
#           independent multinomial-thinned tensors. The implementation uses
#           the equivalent sequential-binomial representation and is vectorized
#           over tensor cells.
# Inputs  : Y = one count tensor; probs = three thinning probabilities.
# Output  : A length-three list of count tensors with the same dimensions as Y.
multinomial_thin_tensor <- function(Y, probs = rep(1 / 3, 3)) {
  Y <- .as_single_histogram_tensor(Y)
  if (length(probs) != 3L) stop("probs must have length 3.")
  if (any(!is.finite(probs)) || any(probs < 0) || sum(probs) <= 0) {
    stop("probs must be nonnegative and have positive sum.")
  }
  probs <- probs / sum(probs)

  y <- as.numeric(Y)
  dims <- dim(Y)

  Y1 <- stats::rbinom(n = length(y), size = y, prob = probs[1])
  rem1 <- y - Y1

  p2_cond <- if (probs[2] + probs[3] > 0) probs[2] / (probs[2] + probs[3]) else 0
  Y2 <- stats::rbinom(n = length(y), size = rem1, prob = p2_cond)
  Y3 <- rem1 - Y2

  list(
    array(as.numeric(Y1), dim = dims),
    array(as.numeric(Y2), dim = dims),
    array(as.numeric(Y3), dim = dims)
  )
}

# Function: make_algorithm_histograms
# Purpose : Create the stage-specific histograms required by either algorithm.
# Inputs  : Y = one count tensor or legacy list of count tensors;
#           algorithm = "no_thinning" or "thinning";
#           thinning_probs = probabilities used for multinomial thinning.
# Output  : List with Y1, Y2, Y3 and their stage sample sizes.
make_algorithm_histograms <- function(Y,
                                      algorithm = c("no_thinning", "thinning"),
                                      thinning_probs = rep(1 / 3, 3)) {
  algorithm <- match.arg(algorithm)
  Y <- .as_single_histogram_tensor(Y)

  parts <- if (algorithm == "thinning") {
    multinomial_thin_tensor(Y, probs = thinning_probs)
  } else {
    list(Y, Y, Y)
  }

  names(parts) <- c("Y1", "Y2", "Y3")
  list(
    Y1 = parts[[1]],
    Y2 = parts[[2]],
    Y3 = parts[[3]],
    n1 = sum(parts[[1]]),
    n2 = sum(parts[[2]]),
    n3 = sum(parts[[3]]),
    n_total = sum(Y),
    algorithm = algorithm
  )
}


# Function: algorithm1_recover_split
# Purpose : Run tensor estimation using three stage-specific histograms.
#           Setting Y1 = Y2 = Y3 recovers the no-thinning algorithm, while
#           using multinomial-thinned Y1, Y2, Y3 recovers the sample-splitting
#           / thinning algorithm.
# Inputs  : Y1 = histogram for the off-diagonal Gram initialization;
#           Y2 = histogram for the refinement SVD step;
#           Y3 = histogram for the final projection step;
#           rank_vec = scalar rank or vector of mode-specific ranks;
#           init = "deflated" or "hetero";
#           niter_init = number of initializer iterations.
# Output  : List with Xhat (recovered scaled tensor), U0, and U.
algorithm1_recover_split <- function(Y1,
                                     Y2,
                                     Y3,
                                     rank_vec,
                                     init = c("deflated", "hetero"),
                                     niter_init = 15) {
  init <- match.arg(init)
  Y1 <- .as_single_histogram_tensor(Y1)
  Y2 <- .as_single_histogram_tensor(Y2)
  Y3 <- .as_single_histogram_tensor(Y3)

  dims <- dim(Y1)
  if (!identical(dim(Y2), dims) || !identical(dim(Y3), dims)) {
    stop("Y1, Y2, and Y3 must have the same dimensions.")
  }

  d <- length(dims)
  if (length(rank_vec) == 1) rank_vec <- rep(rank_vec, d)
  if (length(rank_vec) != d) stop("rank_vec must be length 1 or length d.")

  # Initial subspace estimates from P_off_diag(M_k(Y1) M_k(Y1)^T).
  U0 <- vector("list", d)
  for (k in seq_len(d)) {
    MkY <- unfold_tensor(Y1, k)
    G0 <- MkY %*% t(MkY)
    diag(G0) <- 0
    U0[[k]] <- if (init == "deflated") {
      deflated_heteropca(G0, r = rank_vec[k], niter = niter_init)
    } else {
      heteropca(G0, r = rank_vec[k], niter = niter_init)
    }
  }

  # Refined mode-k subspace estimates using Y2.
  U <- vector("list", d)
  for (k in seq_len(d)) {
    Zk <- Y2
    for (h in setdiff(seq_len(d), k)) {
      Zk <- tensor_times_matrix(Zk, t(U0[[h]]), mode = h)
    }
    MkZ <- unfold_tensor(Zk, k)
    U[[k]] <- top_left_singular_vectors(MkZ, r = rank_vec[k])
  }

  # Final projection using Y3.
  Xhat <- Y3
  for (k in seq_len(d)) {
    Pk <- U[[k]] %*% t(U[[k]])
    Xhat <- tensor_times_matrix(Xhat, Pk, mode = k)
  }

  list(Xhat = Xhat, U0 = U0, U = U)
}

# Function: algorithm1_recover
# Purpose : Backward-compatible wrapper for tensor estimation without
#           multinomial thinning on one histogram.
algorithm1_recover <- function(Y,
                               rank_vec,
                               init = c("deflated", "hetero"),
                               niter_init = 15) {
  init <- match.arg(init)
  Y <- .as_single_histogram_tensor(Y)
  algorithm1_recover_split(
    Y1 = Y,
    Y2 = Y,
    Y3 = Y,
    rank_vec = rank_vec,
    init = init,
    niter_init = niter_init
  )
}

# -----------------------------
# Estimators used in the simulation
# -----------------------------

# Function: estimate_histogram_baseline
# Purpose : Compute the histogram baseline estimator.
# Inputs  : Y_list = one count tensor, or a legacy list of count tensors;
#           project = whether to project to the simplex;
#           algorithm = "no_thinning" or "thinning";
#           thinning_probs = probabilities used if thinning is requested;
#           algorithm_histograms = optional precomputed output from
#                                  make_algorithm_histograms().
# Output  : List with Phat, method, algorithm, and stage sample sizes.
estimate_histogram_baseline <- function(Y_list,
                                        project = TRUE,
                                        algorithm = c("no_thinning", "thinning"),
                                        thinning_probs = rep(1 / 3, 3),
                                        algorithm_histograms = NULL) {
  algorithm <- match.arg(algorithm)
  stage <- algorithm_histograms %||% make_algorithm_histograms(
    Y = Y_list,
    algorithm = algorithm,
    thinning_probs = thinning_probs
  )
  if (stage$n3 <= 0) stop("The final-stage histogram has zero total count. Increase n.")
  Phat <- stage$Y3 / stage$n3
  if (project) Phat <- project_tensor_to_probability(Phat)
  list(
    Phat = Phat,
    method = "histogram",
    algorithm = algorithm,
    n_stage = c(n1 = stage$n1, n2 = stage$n2, n3 = stage$n3, n_total = stage$n_total)
  )
}

# Function: estimate_scaled_algorithm1
# Purpose : Apply either the no-thinning or multinomial-thinning tensor estimator
#           after scaling the histogram(s) by M.
# Inputs  : Y_list = one count tensor, or a legacy list of count tensors;
#           rank = scalar or vector rank;
#           M = scaling tensor;
#           algorithm = "no_thinning" or "thinning";
#           thinning_probs = thinning probabilities for Y1, Y2, Y3;
#           algorithm_histograms = optional precomputed output from
#                                  make_algorithm_histograms();
#           init = initializer type;
#           niter_init = number of initializer iterations;
#           project = whether to project the final estimate to the simplex.
# Output  : List with Phat, Qhat, scaling, fit, algorithm, and stage sample sizes.
estimate_scaled_algorithm1 <- function(Y_list,
                                       rank,
                                       M,
                                       algorithm = c("no_thinning", "thinning"),
                                       thinning_probs = rep(1 / 3, 3),
                                       algorithm_histograms = NULL,
                                       init = c("deflated", "hetero"),
                                       niter_init = 15,
                                       project = TRUE) {
  algorithm <- match.arg(algorithm)
  init <- match.arg(init)
  Y <- .as_single_histogram_tensor(Y_list)
  if (!identical(dim(Y), dim(M))) stop("Y and M must have the same dimensions.")

  stage <- algorithm_histograms %||% make_algorithm_histograms(
    Y = Y,
    algorithm = algorithm,
    thinning_probs = thinning_probs
  )
  if (stage$n3 <= 0) stop("The final-stage histogram has zero total count. Increase n.")

  fit <- algorithm1_recover_split(
    Y1 = stage$Y1 * M,
    Y2 = stage$Y2 * M,
    Y3 = stage$Y3 * M,
    rank_vec = rank,
    init = init,
    niter_init = niter_init
  )
  Qhat <- fit$Xhat / stage$n3
  Phat <- safe_divide(Qhat, M)
  if (project) Phat <- project_tensor_to_probability(Phat)
  list(
    Phat = Phat,
    Qhat = Qhat,
    scaling = M,
    fit = fit,
    algorithm = algorithm,
    n_stage = c(n1 = stage$n1, n2 = stage$n2, n3 = stage$n3, n_total = stage$n_total)
  )
}

# Function: estimate_multiview
# Purpose : Dispatch to one of the simulation estimators.
# Inputs  : Y_list = count tensor or legacy list of count tensors;
#           rank = scalar or vector rank;
#           method = one of "histogram", "unscaled", "oracle", "oracle_cp",
#                    "slice_oracle", "slice_est";
#           P_true = true tensor, needed for slice_oracle;
#           factors = true factors, needed for oracle;
#           pilot_Y = pilot histogram, optional for oracle_cp and slice_est;
#           init = initializer type;
#           niter_init = number of initializer iterations;
#           project = whether to project final estimate to the simplex;
#           cp_max_iter = maximum CP-ALS iterations for oracle_cp;
#           cp_tol = CP-ALS stopping tolerance.
# Output  : Method-specific fit object, always containing Phat and method.
estimate_multiview <- function(Y_list,
                               rank,
                               method = c("histogram", "unscaled", "oracle", "oracle_cp", "slice_oracle", "slice_est"),
                               P_true = NULL,
                               factors = NULL,
                               pilot_Y = NULL,
                               algorithm = c("no_thinning", "thinning"),
                               thinning_probs = rep(1 / 3, 3),
                               algorithm_histograms = NULL,
                               init = c("deflated", "hetero"),
                               niter_init = 15,
                               project = TRUE,
                               cp_max_iter = 50,
                               cp_tol = 1e-5) {
  method <- match.arg(method)
  algorithm <- match.arg(algorithm)
  init <- match.arg(init)
  Y <- .as_single_histogram_tensor(Y_list)
  dims <- dim(Y)

  if (method == "histogram") {
    return(estimate_histogram_baseline(
      Y,
      project = project,
      algorithm = algorithm,
      thinning_probs = thinning_probs,
      algorithm_histograms = algorithm_histograms
    ))
  }

  cp_obj <- NULL

  if (method == "unscaled") {
    M <- make_identity_scaling(dims)
  } else if (method == "oracle") {
    if (is.null(factors)) stop("factors must be provided for method = 'oracle'.")
    M <- make_oracle_scaling_from_factors(factors)
  } else if (method == "oracle_cp") {
    pilot_for_cp <- pilot_Y %||% Y
    cp_obj <- make_oracle_scaling_from_cp_tensor(
      X = pilot_for_cp,
      rank = if (length(rank) == 1) rank else max(rank),
      cp_max_iter = cp_max_iter,
      cp_tol = cp_tol
    )
    M <- cp_obj$M
  } else if (method == "slice_oracle") {
    if (is.null(P_true)) stop("P_true must be provided for method = 'slice_oracle'.")
    M <- make_slice_scaling_from_tensor(P_true)
  } else if (method == "slice_est") {
    pilot_for_slice <- pilot_Y %||% Y
    M <- make_slice_scaling_from_histogram(pilot_for_slice)
  } else {
    stop("Unknown method.")
  }

  out <- estimate_scaled_algorithm1(
    Y_list = Y,
    rank = rank,
    M = M,
    algorithm = algorithm,
    thinning_probs = thinning_probs,
    algorithm_histograms = algorithm_histograms,
    init = init,
    niter_init = niter_init,
    project = project
  )
  out$method <- method
  out$algorithm <- algorithm
  if (!is.null(cp_obj)) {
    out$cp_fit <- cp_obj$cp_fit
    out$factors_est <- cp_obj$factors_est
    out$P_est_cp <- cp_obj$P_est_cp
  }
  out
}

# -----------------------------
# Evaluation and simulation wrappers
# -----------------------------

# Function: evaluate_estimate
# Purpose : Compute accuracy metrics for one estimator.
# Inputs  : Phat = estimated probability tensor; P_true = true probability tensor.
# Output  : One-row data frame with l1, Frobenius, and max-absolute errors.
evaluate_estimate <- function(Phat, P_true) {
  data.frame(
    l1_error = tensor_l1_error(Phat, P_true),
    fro_error = tensor_fro_error(Phat, P_true),
    max_abs_error = max(abs(Phat - P_true)),
    stringsAsFactors = FALSE
  )
}

# Function: compare_methods_one_rep
# Purpose : Run several methods on one replication and collect errors.
# Inputs  : P_true = true probability tensor;
#           Y_list = count tensor or legacy list of count tensors;
#           rank = scalar or vector rank;
#           factors = true factors for oracle;
#           pilot_Y = optional pilot histogram;
#           methods = character vector of method names;
#           init = initializer type;
#           niter_init = number of initializer iterations;
#           project = whether to project the final estimate to the simplex;
#           cp_max_iter = maximum CP-ALS iterations for oracle_cp;
#           cp_tol = CP-ALS stopping tolerance.
# Output  : Data frame with one row per method.

compare_methods_one_rep <- function(P_true,
                                    Y_list,
                                    rank,
                                    factors = NULL,
                                    pilot_Y = NULL,
                                    methods = c("histogram", "unscaled", "oracle", "oracle_cp", "slice_oracle", "slice_est"),
                                    algorithm = "no_thinning",
                                    thinning_probs = rep(1 / 3, 3),
                                    init = c("deflated", "hetero"),
                                    niter_init = 15,
                                    project = TRUE,
                                    cp_max_iter = 50,
                                    cp_tol = 1e-5) {
  algorithm <- match.arg(
    algorithm,
    choices = c("no_thinning", "thinning"),
    several.ok = TRUE
  )
  init <- match.arg(init)

  res_list <- unlist(lapply(algorithm, function(alg) {
    stage <- make_algorithm_histograms(
      Y = Y_list,
      algorithm = alg,
      thinning_probs = thinning_probs
    )
    lapply(methods, function(meth) {
      fit <- estimate_multiview(
        Y_list = Y_list,
        rank = rank,
        method = meth,
        P_true = P_true,
        factors = factors,
        pilot_Y = pilot_Y,
        algorithm = alg,
        thinning_probs = thinning_probs,
        algorithm_histograms = stage,
        init = init,
        niter_init = niter_init,
        project = project,
        cp_max_iter = cp_max_iter,
        cp_tol = cp_tol
      )
      cbind(
        data.frame(method = meth, algorithm = alg, stringsAsFactors = FALSE),
        evaluate_estimate(fit$Phat, P_true)
      )
    })
  }), recursive = FALSE)

  do.call(rbind, res_list)
}

# Function: run_simulation
# Purpose : Run a Monte Carlo simulation across replications.
# Inputs  : nrep = number of replications;
#           p_vec, R, n, factor_type, alpha, s_vec, weight_type,
#           weight_values, geometric_ratio = data-generating parameters;
#           methods = character vector of method names;
#           init = initializer type;
#           niter_init = number of initializer iterations;
#           project = whether to project final estimates;
#           pilot_n = optional pilot histogram sample size;
#           seed = random seed;
#           verbose = whether to print progress;
#           cp_max_iter = maximum CP-ALS iterations for oracle_cp;
#           cp_tol = CP-ALS stopping tolerance.
# Output  : Data frame stacking all methods and all replications.
run_simulation <- function(nrep,
                           p_vec,
                           R,
                           n,
                           factor_type = c("dense", "sparse"),
                           alpha = 1,
                           s_vec = NULL,
                           weight_type = c("balanced", "geometric", "custom"),
                           weight_values = NULL,
                           geometric_ratio = 0.5,
                           methods = c("histogram", "unscaled", "oracle", "oracle_cp", "slice_oracle", "slice_est"),
                           algorithm = "no_thinning",
                           thinning_probs = rep(1 / 3, 3),
                           init = c("deflated", "hetero"),
                           niter_init = 15,
                           project = TRUE,
                           pilot_n = NULL,
                           seed = NULL,
                           verbose = TRUE,
                           cp_max_iter = 50,
                           cp_tol = 1e-5) {
  factor_type <- match.arg(factor_type)
  weight_type <- match.arg(weight_type)
  algorithm <- match.arg(algorithm, several.ok = TRUE)
  init <- match.arg(init)

  if (!is.null(seed)) set.seed(seed)

  out <- vector("list", nrep)
  for (b in seq_len(nrep)) {
    if (verbose && (b %% max(1, floor(nrep / 10)) == 0 || b == 1)) {
      message(sprintf("Replication %d / %d", b, nrep))
    }

    sim <- simulate_one_replication(
      p_vec = p_vec,
      R = R,
      n = n,
      factor_type = factor_type,
      alpha = alpha,
      s_vec = s_vec,
      weight_type = weight_type,
      weight_values = weight_values,
      geometric_ratio = geometric_ratio,
      nsplit = 1,
      pilot_n = pilot_n
    )

    tmp <- compare_methods_one_rep(
      P_true = sim$P,
      Y_list = sim$Y %||% sim$Y_list,
      rank = R,
      factors = sim$factors,
      pilot_Y = sim$pilot_Y,
      methods = methods,
      algorithm = algorithm,
      thinning_probs = thinning_probs,
      init = init,
      niter_init = niter_init,
      project = project,
      cp_max_iter = cp_max_iter,
      cp_tol = cp_tol
    )

    tmp$rep <- b
    tmp$n <- n
    tmp$R <- R
    tmp$factor_type <- factor_type
    tmp$weight_type <- weight_type
    tmp$p <- paste(p_vec, collapse = "x")
    out[[b]] <- tmp
  }

  do.call(rbind, out)
}

# Function: summarize_simulation
# Purpose : Summarize Monte Carlo errors by method.
# Inputs  : sim_df = output data frame from run_simulation().
# Output  : Data frame with method-wise means and standard errors.
summarize_simulation <- function(sim_df) {
  split_key <- if ("algorithm" %in% names(sim_df)) {
    interaction(sim_df$method, sim_df$algorithm, drop = TRUE, lex.order = TRUE)
  } else {
    sim_df$method
  }
  split_df <- split(sim_df, split_key)
  out <- lapply(split_df, function(df) {
    ans <- data.frame(
      method = unique(df$method),
      mean_l1 = mean(df$l1_error),
      se_l1 = sd(df$l1_error) / sqrt(nrow(df)),
      mean_fro = mean(df$fro_error),
      se_fro = sd(df$fro_error) / sqrt(nrow(df)),
      mean_max_abs = mean(df$max_abs_error),
      se_max_abs = sd(df$max_abs_error) / sqrt(nrow(df))
    )
    if ("algorithm" %in% names(df)) ans$algorithm <- unique(df$algorithm)
    ans
  })
  do.call(rbind, out)
}

# -----------------------------
# Example settings
# -----------------------------

# Function: example_settings
# Purpose : Return a few canned simulation settings.
# Inputs  : None.
# Output  : Named list of setting lists.
example_settings <- function() {
  list(
    balanced_dense = list(
      p_vec = c(40, 40, 40),
      R = 4,
      n = 40 * 40 * 4,
      factor_type = "dense",
      alpha = 1,
      weight_type = "balanced"
    ),
    sparse_balanced = list(
      p_vec = c(40, 40, 40),
      R = 4,
      n = 40 * 40 * 4,
      factor_type = "sparse",
      s_vec = c(8, 8, 8),
      alpha = 0.2,
      weight_type = "balanced"
    ),
    imbalanced_dense = list(
      p_vec = c(40, 40, 40),
      R = 4,
      n = 40 * 40 * 4,
      factor_type = "dense",
      alpha = 1,
      weight_type = "custom",
      weight_values = c(0.85, 0.10, 0.04, 0.01)
    )
  )
}
