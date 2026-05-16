# ============================================================
# Run the first four simulation experiments (v5) for the multiview
# density estimation paper.
#
# This script does simulation ONLY.
# It saves raw and summarized CSV files for later plotting.
#
# Linux / macOS parallelism is supported through mclapply().
# ============================================================

# -----------------------------
# User configuration
# -----------------------------
setwd("~/density/v5")
source("multiview_simulation_functions_v5.R")

METHODS <- c("histogram", "unscaled", "oracle", "oracle_cp", "slice_oracle", "slice_est")
ALGORITHM_OPTIONS <- c("no_thinning", "thinning")
INIT_METHOD <- "deflated"
NITER_INIT <- 15
PROJECT_TO_SIMPLEX <- TRUE
CP_MAX_ITER <- 50
CP_TOL <- 1e-5
MC_CORES <- max(1L, parallel::detectCores() - 1L)

# These defaults are still moderately expensive. Increase or decrease as needed.
NREP_EXP1 <- 30
NREP_EXP2 <- 50
NREP_EXP3 <- 50
NREP_EXP4 <- 30

# Experiment 1 compares the unscaled estimator with and without thinning.
# It uses smaller p and n than Experiment 2 to keep the extra comparison cheap.
DEFAULT_N_GRID_EXP1 <- c(10000L, 30000L, 100000L, 300000L, 1000000L)

# Experiment 2 keeps the original larger sample-size grid.
DEFAULT_N_GRID_EXP2 <- c(300000L, 1000000L, 3000000L, 10000000L, 30000000L, 100000000L)

OUTPUT_DIR <- "multiview_four_experiments_v5"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Small file / summary helpers
# -----------------------------

# Function: write_results_csv
# Purpose : Save a data frame as a CSV file without row names.
# Inputs  : df = data frame to save;
#           filename = output CSV path.
# Output  : Invisibly returns filename.
write_results_csv <- function(df, filename) {
  utils::write.csv(df, file = filename, row.names = FALSE)
  invisible(filename)
}

# Function: summarize_numeric_by
# Purpose : Summarize a numeric vector by selected grouping variables.
# Inputs  : values = numeric vector of length nrow(df);
#           df = data frame containing grouping columns;
#           group_vars = character vector of grouping-column names;
#           metric_name = label stored in the metric column.
# Output  : Data frame with grouping columns plus n_rep, metric, mean, sd, and se.
summarize_numeric_by <- function(values,
                                 df,
                                 group_vars,
                                 metric_name = "custom_metric") {
  if (length(values) != nrow(df)) {
    stop("values must have length nrow(df).")
  }
  if (!all(group_vars %in% names(df))) {
    stop("Some group_vars are not columns in df.")
  }

  tmp_df <- df[, group_vars, drop = FALSE]
  tmp_df$.value_to_summarize <- values
  split_key <- do.call(interaction, c(tmp_df[group_vars], list(drop = TRUE, lex.order = TRUE)))
  split_list <- split(tmp_df, split_key)

  out <- lapply(split_list, function(dat) {
    one <- dat[1, group_vars, drop = FALSE]
    one$n_rep <- nrow(dat)
    one$metric <- metric_name
    one$mean <- mean(dat$.value_to_summarize)
    one$sd <- stats::sd(dat$.value_to_summarize)
    one$se <- one$sd / sqrt(one$n_rep)
    one
  })

  out_df <- do.call(rbind, out)
  rownames(out_df) <- NULL
  out_df
}

# Function: summarize_metric_by
# Purpose : Summarize one simulation error metric by arbitrary grouping variables.
# Inputs  : sim_df = simulation results data frame;
#           group_vars = grouping columns;
#           metric = one of "l1_error", "fro_error", or "max_abs_error".
# Output  : Data frame with grouping columns plus n_rep, metric, mean, sd, and se.
summarize_metric_by <- function(sim_df,
                                group_vars = c("method"),
                                metric = c("l1_error", "fro_error", "max_abs_error")) {
  metric <- match.arg(metric)
  if (!all(group_vars %in% names(sim_df))) {
    stop("Some group_vars are not columns in sim_df.")
  }

  split_key <- do.call(interaction, c(sim_df[group_vars], list(drop = TRUE, lex.order = TRUE)))
  split_list <- split(sim_df, split_key)

  out <- lapply(split_list, function(df) {
    one <- df[1, group_vars, drop = FALSE]
    one$n_rep <- nrow(df)
    one$metric <- metric
    one$mean <- mean(df[[metric]])
    one$sd <- stats::sd(df[[metric]])
    one$se <- one$sd / sqrt(one$n_rep)
    one
  })

  out_df <- do.call(rbind, out)
  rownames(out_df) <- NULL
  out_df
}

# Function: add_normalized_error_columns
# Purpose : Add normalized l1 and Frobenius errors used in scaling plots.
# Inputs  : sim_df = simulation data frame;
#           p_col = column containing common mode dimension p;
#           R_col = column containing rank R;
#           n_col = column containing sample size n.
# Output  : The same data frame with columns normalized_l1 and normalized_fro added.
add_normalized_error_columns <- function(sim_df,
                                         p_col = "p_scalar",
                                         R_col = "R",
                                         n_col = "n") {
  if (!all(c(p_col, R_col, n_col) %in% names(sim_df))) {
    stop("p_col, R_col, and n_col must all be present in sim_df.")
  }
  sim_df$normalized_fiber_l1_max = sim_df$fiber_l1_max * sim_df[[p_col]]^2
  sim_df$normalized_slice_l1_max = sim_df$slice_l1_max * sim_df[[p_col]]
  sim_df$normalized_l1 <- sim_df$l1_error * sqrt(sim_df[[n_col]] / (sim_df[[p_col]] * sim_df[[R_col]]))
  sim_df$normalized_fro <- sim_df$fro_error * sqrt(sim_df[[n_col]] / (sim_df$fiber_l1_max * sim_df[[R_col]]))
  sim_df
}

# Function: compute_ratio_to_reference
# Purpose : Compute summary-value ratios relative to a reference method.
# Inputs  : summary_df = summarized data frame;
#           reference_method = method used in the denominator;
#           group_vars = grouping columns excluding method;
#           method_col = method-column name;
#           value_col = summarized-value column.
# Output  : Data frame with added columns ref_value and ratio_to_ref.
compute_ratio_to_reference <- function(summary_df,
                                       reference_method = "oracle",
                                       group_vars,
                                       method_col = "method",
                                       value_col = "mean") {
  if (!all(c(group_vars, method_col, value_col) %in% names(summary_df))) {
    stop("group_vars, method_col, and value_col must all be present in summary_df.")
  }

  ref_df <- summary_df[summary_df[[method_col]] == reference_method, c(group_vars, value_col), drop = FALSE]
  names(ref_df)[names(ref_df) == value_col] <- "ref_value"

  out <- merge(summary_df, ref_df, by = group_vars, all.x = TRUE, sort = FALSE)
  out$ratio_to_ref <- out[[value_col]] / out$ref_value
  out
}

# -----------------------------
# Parallel simulation helpers
# -----------------------------

# Function: run_and_tag_standard_simulation
# Purpose : Run one design point using the standard generator in
#           multiview_simulation_functions_wo_thinning.R and append metadata columns.
# Inputs  : sim_args = named list passed to run_custom_simulation_parallel();
#           metadata = named list of extra columns to append;
#           mc.cores = number of worker processes.
# Output  : Data frame of raw simulation results for this design point.
# run_and_tag_standard_simulation <- function(sim_args,
#                                             metadata = list(),
#                                             mc.cores = MC_CORES) {
#   sim_df <- do.call(
#     run_custom_simulation_parallel,
#     c(sim_args, list(mc.cores = mc.cores))
#   )
# 
#   for (nm in names(metadata)) {
#     sim_df[[nm]] <- metadata[[nm]]
#   }
#   sim_df
# }

# Function: generate_dense_heteroskedastic_factors
# Purpose : Generate dense probability-factor matrices with strongly varying,
#           but still fully dense, coordinate profiles.
# Inputs  : p_vec = vector of mode dimensions;
#           R = CP rank;
#           alpha = global scale for Dirichlet concentrations;
#           hetero_strength = max/min scale ratio of the mean profile;
#           random_permutation = whether to randomly permute the profile for
#                                each component and mode.
# Output  : List of factor matrices; each matrix has dimension p_vec[k] x R.
generate_dense_heteroskedastic_factors <- function(p_vec,
                                                   R,
                                                   alpha = 0.8,
                                                   hetero_strength = 50,
                                                   random_permutation = TRUE) {
  if (hetero_strength < 1) stop("hetero_strength must be >= 1.")

  base_profile <- exp(seq(0, log(hetero_strength), length.out = max(p_vec)))
  base_profile <- base_profile / mean(base_profile)

  lapply(p_vec, function(p) {
    out <- sapply(seq_len(R), function(r) {
      prof <- base_profile[seq_len(p)]
      if (random_permutation) {
        prof <- prof[sample.int(p)]
      }
      alpha_vec <- alpha * prof
      rdirichlet1(alpha_vec)
    })
    if (R == 1) out <- matrix(out, ncol = 1)
    out
  })
}

# Function: generate_multiview_tensor_dense_heteroskedastic
# Purpose : Generate one dense but heteroskedastic low-rank probability tensor.
# Inputs  : p_vec = vector of mode dimensions;
#           R = CP rank;
#           alpha = Dirichlet scale parameter;
#           hetero_strength = max/min scale ratio in each factor profile;
#           weight_type = "balanced", "geometric", or "custom";
#           weight_values = custom weights if needed;
#           geometric_ratio = common ratio for geometric weights.
# Output  : List containing P, weights, factors, and factor_type.
generate_multiview_tensor_dense_heteroskedastic <- function(p_vec,
                                                            R,
                                                            alpha = 0.8,
                                                            hetero_strength = 50,
                                                            weight_type = c("balanced", "geometric", "custom"),
                                                            weight_values = NULL,
                                                            geometric_ratio = 0.5) {
  weight_type <- match.arg(weight_type)

  weights <- generate_weights(
    R = R,
    type = weight_type,
    values = weight_values,
    ratio = geometric_ratio
  )

  factors <- generate_dense_heteroskedastic_factors(
    p_vec = p_vec,
    R = R,
    alpha = alpha,
    hetero_strength = hetero_strength,
    random_permutation = TRUE
  )

  P <- build_cp_tensor(weights, factors)

  list(
    P = P,
    weights = weights,
    factors = factors,
    factor_type = "dense_heteroskedastic"
  )
}

# Function: simulate_one_replication_dense_heteroskedastic
# Purpose : Generate one full replication under the dense heteroskedastic model.
# Inputs  : p_vec = vector of mode dimensions;
#           R = CP rank;
#           n = sample size for Y_list;
#           alpha = Dirichlet scale parameter;
#           hetero_strength = profile heterogeneity strength;
#           weight_type = weight regime;
#           weight_values = custom weight vector if needed;
#           geometric_ratio = common ratio for geometric weights;
#           nsplit = legacy argument; ignored by the no-thinning generator;
#           pilot_n = pilot-sample size for M estimation.
# Output  : List containing P, weights, factors, Y/Y_list, pilot_Y, n,
#           factor_type, and weight_type.
simulate_one_replication_dense_heteroskedastic <- function(p_vec,
                                                           R,
                                                           n,
                                                           alpha = 0.8,
                                                           hetero_strength = 50,
                                                           weight_type = c("balanced", "geometric", "custom"),
                                                           weight_values = NULL,
                                                           geometric_ratio = 0.5,
                                                           nsplit = 1,
                                                           pilot_n = NULL) {
  weight_type <- match.arg(weight_type)

  model <- generate_multiview_tensor_dense_heteroskedastic(
    p_vec = p_vec,
    R = R,
    alpha = alpha,
    hetero_strength = hetero_strength,
    weight_type = weight_type,
    weight_values = weight_values,
    geometric_ratio = geometric_ratio
  )

  # No multinomial thinning: draw one histogram tensor and use it in all
  # estimation stages. The legacy name Y_list is retained only for wrappers.
  Y <- rmultinomial_tensor(n, model$P)
  pilot_Y <- if (!is.null(pilot_n) && pilot_n > 0) rmultinomial_tensor(pilot_n, model$P) else NULL

  c(model, list(
    Y = Y,
    Y_list = Y,
    pilot_Y = pilot_Y,
    n = n,
    weight_type = weight_type,
    hetero_strength = hetero_strength
  ))
}

# Function: run_custom_simulation_parallel
# Purpose : Run Monte Carlo replications using a user-supplied simulation generator.
# Inputs  : nrep = number of replications;
#           generator_fun = function returning a list with P, factors, Y/Y_list,
#                           pilot_Y, n, and optionally weight_type / factor_type;
#           generator_args = named list of arguments passed to generator_fun;
#           methods = estimation methods to compare;
#           algorithm = one or more of "no_thinning" and "thinning";
#           thinning_probs = probabilities used by the thinning algorithm;
#           init = initializer type for Algorithm 1;
#           niter_init = iterations for the initializer;
#           project = whether to project the final estimate to the simplex;
#           seed = base seed; replication b uses seed + b - 1;
#           verbose = whether to print progress messages;
#           cp_max_iter = maximum ALS iterations for oracle_cp;
#           cp_tol = stopping tolerance for oracle_cp;
#           mc.cores = number of worker processes;
#           mc.preschedule = passed to mclapply().
# Output  : Data frame stacking all methods and replications.
run_custom_simulation_parallel <- function(nrep,
                                           generator_fun,
                                           generator_args,
                                           methods = METHODS,
                                           algorithm = "no_thinning",
                                           thinning_probs = rep(1 / 3, 3),
                                           init = c("deflated", "hetero"),
                                           niter_init = NITER_INIT,
                                           project = PROJECT_TO_SIMPLEX,
                                           seed = NULL,
                                           verbose = TRUE,
                                           cp_max_iter = CP_MAX_ITER,
                                           cp_tol = CP_TOL,
                                           mc.cores = MC_CORES,
                                           mc.preschedule = TRUE) {
  algorithm <- match.arg(algorithm, choices = ALGORITHM_OPTIONS, several.ok = TRUE)
  init <- match.arg(init)

  worker_fun <- function(b) {
    if (!is.null(seed)) {
      set.seed(seed + b - 1L)
    }

    sim <- do.call(generator_fun, generator_args)

    tmp <- compare_methods_one_rep(
      P_true = sim$P,
      Y_list = sim$Y %||% sim$Y_list,
      rank = generator_args$R,
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
    tmp$n <- sim$n
    tmp$R <- generator_args$R
    tmp$factor_type <- sim$factor_type %||% "custom_generator"
    tmp$weight_type <- sim$weight_type %||% NA_character_
    tmp$p <- paste(generator_args$p_vec, collapse = "x")
    tmp$fiber_l1_max = fiber_l1_max(sim$P)
    tmp$slice_l1_max = slice_l1_max(sim$P)
    
    tmp
  }

  use_parallel <- (.Platform$OS.type != "windows") && (mc.cores > 1L)

  if (verbose) {
    if (use_parallel) {
      message(sprintf("Running %d replications in parallel with mc.cores = %d", nrep, mc.cores))
    } else {
      message(sprintf("Running %d replications sequentially", nrep))
    }
  }

  if (use_parallel) {
    out <- parallel::mclapply(
      X = seq_len(nrep),
      FUN = worker_fun,
      mc.cores = mc.cores,
      mc.preschedule = mc.preschedule,
      mc.set.seed = FALSE
    )
  } else {
    out <- lapply(seq_len(nrep), function(b) {
      if (verbose && (b %% max(1, floor(nrep / 10)) == 0 || b == 1)) {
        message(sprintf("Replication %d / %d", b, nrep))
      }
      worker_fun(b)
    })
  }

  do.call(rbind, out)
}


# -----------------------------
# Experiment 1
# -----------------------------

# Function: run_experiment1_unscaled_thinning_comparison_simulation
# Purpose : Compare the unscaled estimator under the no-thinning and
#           multinomial-thinning algorithms.  The data-generating model is
#           the same dense heteroskedastic regime as Experiment 2, but with
#           smaller p and n for a lighter diagnostic experiment.
# Inputs  : out_dir = root output directory;
#           nrep = number of Monte Carlo replications per sample size;
#           p = common mode dimension;
#           R = CP rank;
#           n_grid = vector of sample sizes;
#           alpha = Dirichlet scale parameter for dense heteroskedastic factors;
#           hetero_strength = max/min scale ratio in each factor profile.
# Output  : Invisibly returns a list with raw and summarized data frames.
run_experiment1_unscaled_thinning_comparison_simulation <- function(out_dir,
                                                                    nrep = NREP_EXP1,
                                                                    p = 30,
                                                                    R = 4,
                                                                    n_grid = DEFAULT_N_GRID_EXP1,
                                                                    alpha = 0.8,
                                                                    hetero_strength = 50) {
  exp_dir <- file.path(out_dir, "experiment1_unscaled_thinning_comparison")
  dir.create(exp_dir, recursive = TRUE, showWarnings = FALSE)

  raw_list <- vector("list", length(n_grid))
  for (i in seq_along(n_grid)) {
    n_cur <- n_grid[i]
    message(sprintf("Experiment 1 | n = %d | algorithms = %s",
                    n_cur, paste(ALGORITHM_OPTIONS, collapse = ", ")))

    raw_list[[i]] <- run_custom_simulation_parallel(
      nrep = nrep,
      generator_fun = simulate_one_replication_dense_heteroskedastic,
      generator_args = list(
        p_vec = c(p, p, p),
        R = R,
        n = n_cur,
        alpha = alpha,
        hetero_strength = hetero_strength,
        weight_type = "balanced",
        nsplit = 1,
        pilot_n = NULL
      ),
      methods = "unscaled",
      algorithm = ALGORITHM_OPTIONS,
      thinning_probs = rep(1 / 3, 3),
      init = INIT_METHOD,
      niter_init = NITER_INIT,
      project = PROJECT_TO_SIMPLEX,
      seed = 10000 + 1000 * i,
      verbose = TRUE,
      cp_max_iter = CP_MAX_ITER,
      cp_tol = CP_TOL,
      mc.cores = MC_CORES
    )

    raw_list[[i]]$experiment <- "unscaled_thinning_comparison"
    raw_list[[i]]$p_scalar <- p
    raw_list[[i]]$sample_multiplier <- n_cur / (p * R)
    raw_list[[i]]$hetero_strength <- hetero_strength
    raw_list[[i]]$alpha_profile <- alpha
  }

  raw_df <- do.call(rbind, raw_list)
  raw_df <- add_normalized_error_columns(raw_df)

  group_vars <- c("method", "algorithm", "n", "sample_multiplier", "hetero_strength")
  summary_l1 <- summarize_metric_by(raw_df, group_vars = group_vars, metric = "l1_error")
  summary_fro <- summarize_metric_by(raw_df, group_vars = group_vars, metric = "fro_error")
  summary_normalized_l1 <- summarize_numeric_by(
    values = raw_df$normalized_l1,
    df = raw_df,
    group_vars = group_vars,
    metric_name = "normalized_l1"
  )
  summary_normalized_fro <- summarize_numeric_by(
    values = raw_df$normalized_fro,
    df = raw_df,
    group_vars = group_vars,
    metric_name = "normalized_fro"
  )

  write_results_csv(raw_df, file.path(exp_dir, "raw_results.csv"))
  write_results_csv(summary_l1, file.path(exp_dir, "summary_l1.csv"))
  write_results_csv(summary_fro, file.path(exp_dir, "summary_fro.csv"))
  write_results_csv(summary_normalized_l1, file.path(exp_dir, "summary_normalized_l1.csv"))
  write_results_csv(summary_normalized_fro, file.path(exp_dir, "summary_normalized_fro.csv"))

  invisible(list(
    raw = raw_df,
    summary_l1 = summary_l1,
    summary_fro = summary_fro,
    summary_normalized_l1 = summary_normalized_l1,
    summary_normalized_fro = summary_normalized_fro
  ))
}

# Backward-compatible alias: Experiment 1 is now the unscaled thinning comparison.
run_experiment1_balanced_dense_simulation <- run_experiment1_unscaled_thinning_comparison_simulation

# -----------------------------
# Experiment 2
# -----------------------------

# Function: run_experiment2_dense_heteroskedastic_simulation
# Purpose : Run Experiment 2 using a dense but heteroskedastic tensor and save CSV files.
# Inputs  : out_dir = root output directory;
#           nrep = number of Monte Carlo replications per sample size;
#           p = common mode dimension;
#           R = CP rank;
#           n_grid = vector of sample sizes;
#           alpha = Dirichlet scale parameter for dense heteroskedastic factors;
#           hetero_strength = max/min scale ratio in each factor profile.
# Output  : Invisibly returns a list with raw and summarized data frames.
run_experiment2_dense_heteroskedastic_simulation <- function(out_dir,
                                                             nrep = NREP_EXP2,
                                                             p = 50,
                                                             R = 4,
                                                             n_grid = DEFAULT_N_GRID_EXP2,
                                                             alpha = 0.8,
                                                             hetero_strength = 100) {
  exp_dir <- file.path(out_dir, "experiment2_dense_heteroskedastic")
  dir.create(exp_dir, recursive = TRUE, showWarnings = FALSE)

  raw_list <- vector("list", length(n_grid))
  for (i in seq_along(n_grid)) {
    n_cur <- n_grid[i]
    message(sprintf("Experiment 2 | n = %d", n_cur))

    raw_list[[i]] <- run_custom_simulation_parallel(
      nrep = nrep,
      generator_fun = simulate_one_replication_dense_heteroskedastic,
      generator_args = list(
        p_vec = c(p, p, p),
        R = R,
        n = n_cur,
        alpha = alpha,
        hetero_strength = hetero_strength,
        weight_type = "balanced",
        nsplit = 1,
        pilot_n = NULL
      ),
      methods = METHODS,
      init = INIT_METHOD,
      niter_init = NITER_INIT,
      project = PROJECT_TO_SIMPLEX,
      seed = 2000 + 1000 * i,
      verbose = TRUE,
      cp_max_iter = CP_MAX_ITER,
      cp_tol = CP_TOL,
      mc.cores = MC_CORES
    )

    raw_list[[i]]$experiment <- "dense_heteroskedastic"
    raw_list[[i]]$p_scalar <- p
    raw_list[[i]]$sample_multiplier <- n_cur / (p * R)
    raw_list[[i]]$hetero_strength <- hetero_strength
    raw_list[[i]]$alpha_profile <- alpha
  }

  raw_df <- do.call(rbind, raw_list)
  raw_df <- add_normalized_error_columns(raw_df)

  summary_l1 <- summarize_metric_by(raw_df, group_vars = c("method", "n", "sample_multiplier", "hetero_strength"), metric = "l1_error")
  summary_fro <- summarize_metric_by(raw_df, group_vars = c("method", "n", "sample_multiplier", "hetero_strength"), metric = "fro_error")
  summary_normalized_l1 <- summarize_numeric_by(
    values = raw_df$normalized_l1,
    df = raw_df,
    group_vars = c("method", "n", "sample_multiplier", "hetero_strength"),
    metric_name = "normalized_l1"
  )
  summary_normalized_fro <- summarize_numeric_by(
    values = raw_df$normalized_fro,
    df = raw_df,
    group_vars = c("method", "n", "sample_multiplier", "hetero_strength"),
    metric_name = "normalized_fro"
  )

  write_results_csv(raw_df, file.path(exp_dir, "raw_results.csv"))
  write_results_csv(summary_l1, file.path(exp_dir, "summary_l1.csv"))
  write_results_csv(summary_fro, file.path(exp_dir, "summary_fro.csv"))
  write_results_csv(summary_normalized_l1, file.path(exp_dir, "summary_normalized_l1.csv"))
  write_results_csv(summary_normalized_fro, file.path(exp_dir, "summary_normalized_fro.csv"))

  invisible(list(
    raw = raw_df,
    summary_l1 = summary_l1,
    summary_fro = summary_fro,
    summary_normalized_l1 = summary_normalized_l1,
    summary_normalized_fro = summary_normalized_fro
  ))
}

# -----------------------------
# Experiment 3
# -----------------------------

# Function: run_experiment3_vary_heteroskedastic_simulation
# Purpose : Run Experiment 3 across several heteroskedastic regimes and save CSV files.
# Inputs  : out_dir = root output directory;
#           nrep = number of Monte Carlo replications per heteroskedastic setting;
#           p = common mode dimension;
#           R = CP rank;
#           n = sample size;
#           weight_list = named list of heteroskedastic vectors.
# Output  : Invisibly returns a list with raw and summarized data frames.
run_experiment3_vary_heteroskedastic_simulation <- function(out_dir,
                                                        nrep = NREP_EXP3,
                                                        p = 50,
                                                        R = 4,
                                                        n = 200000,
                                                        alpha = 0.8,
                                                        #hetero_strength_vec = seq(1, 100, length.out = 6)
                                                        #hetero_strength_vec = c(1, 10, 100, 1000, 10000, 100000)
                                                        hetero_strength_vec = c(1, 3, 10, 30, 100, 300, 1000)
                                                        ) {
  exp_dir <- file.path(out_dir, "experiment3_vary_heteroskedastic")
  dir.create(exp_dir, recursive = TRUE, showWarnings = FALSE)

  raw_list <- vector("list", length(hetero_strength_vec))

  for (i in seq_along(hetero_strength_vec)) {
    hetero_strength <- hetero_strength_vec[i]
    message(sprintf("Experiment 3 | hetero_strength = %s", hetero_strength))

    raw_list[[i]] <- run_custom_simulation_parallel(
      nrep = nrep,
      generator_fun = simulate_one_replication_dense_heteroskedastic,
      generator_args = list(
        p_vec = c(p, p, p),
        R = R,
        n = n,
        alpha = alpha,
        hetero_strength = hetero_strength,
        weight_type = "balanced",
        nsplit = 1,
        pilot_n = NULL
      ),
      methods = METHODS,
      init = INIT_METHOD,
      niter_init = NITER_INIT,
      project = PROJECT_TO_SIMPLEX,
      seed = 20000 + 1000 * i,
      verbose = TRUE,
      cp_max_iter = CP_MAX_ITER,
      cp_tol = CP_TOL,
      mc.cores = MC_CORES
    )
    
    raw_list[[i]]$experiment <- "dense_heteroskedastic"
    raw_list[[i]]$p_scalar <- p
    raw_list[[i]]$sample_multiplier <- n / (p * R)
    raw_list[[i]]$hetero_strength <- hetero_strength
    raw_list[[i]]$alpha_profile <- alpha
  }

  raw_df <- do.call(rbind, raw_list)
  raw_df <- add_normalized_error_columns(raw_df)

  summary_l1 <- summarize_metric_by(raw_df, group_vars = c("method", "n", "sample_multiplier", "hetero_strength"), metric = "l1_error")
  summary_fro <- summarize_metric_by(raw_df, group_vars = c("method", "n", "sample_multiplier", "hetero_strength"), metric = "fro_error")
  summary_normalized_l1 <- summarize_numeric_by(
    values = raw_df$normalized_l1,
    df = raw_df,
    group_vars = c("method", "n", "sample_multiplier", "hetero_strength"),
    metric_name = "normalized_l1"
  )
  summary_normalized_fro <- summarize_numeric_by(
    values = raw_df$normalized_fro,
    df = raw_df,
    group_vars = c("method", "n", "sample_multiplier", "hetero_strength"),
    metric_name = "normalized_fro"
  )

  write_results_csv(raw_df, file.path(exp_dir, "raw_results.csv"))
  write_results_csv(summary_l1, file.path(exp_dir, "summary_l1.csv"))
  write_results_csv(summary_fro, file.path(exp_dir, "summary_fro.csv"))
  write_results_csv(summary_normalized_l1, file.path(exp_dir, "summary_normalized_l1.csv"))
  write_results_csv(summary_normalized_fro, file.path(exp_dir, "summary_normalized_fro.csv"))
  
  invisible(list(
    raw = raw_df,
    summary_l1 = summary_l1,
    summary_fro = summary_fro,
    summary_normalized_l1 = summary_normalized_l1,
    summary_normalized_fro = summary_normalized_fro
  ))
}

# -----------------------------
# Experiment 4
# -----------------------------

# Function: run_experiment4_rank_dimension_scaling_simulation
# Purpose : Run Experiment 4 varying p and R with n proportional to pR, then save CSV files.
# Inputs  : out_dir = root output directory;
#           nrep = number of Monte Carlo replications per grid point;
#           p_grid = vector of common mode dimensions;
#           R_grid = vector of CP ranks;
#           sample_multiplier = constant c in n = c * p * R.
# Output  : Invisibly returns a list with raw and summarized data frames.
run_experiment4_rank_dimension_scaling_simulation <- function(out_dir,
                                                              nrep = NREP_EXP4,
                                                              p_grid = c(60, 80, 100, 120),
                                                              R_grid = c(2, 4, 6, 8),
                                                              sample_multiplier = 120,
                                                              alpha = 0.8) {
  exp_dir <- file.path(out_dir, "experiment4_rank_dimension_scaling")
  dir.create(exp_dir, recursive = TRUE, showWarnings = FALSE)

  raw_list <- vector("list", length(p_grid) * length(R_grid))
  idx <- 1L

  for (R in R_grid) {
    for (p in p_grid) {
      n_cur <- sample_multiplier * p * R
      message(sprintf("Experiment 4 | p = %d, R = %d, n = %d", p, R, n_cur))
      
      raw_list[[idx]] <- run_custom_simulation_parallel(
        nrep = nrep,
        generator_fun = simulate_one_replication_dense_heteroskedastic,
        generator_args = list(
          p_vec = c(p, p, p),
          R = R,
          n = n_cur,
          alpha = alpha,
          hetero_strength = 1,
          weight_type = "balanced",
          nsplit = 1,
          pilot_n = NULL
        ),
        methods = METHODS,
        init = INIT_METHOD,
        niter_init = NITER_INIT,
        project = PROJECT_TO_SIMPLEX,
        seed = 20000 + 1000 * idx,
        verbose = TRUE,
        cp_max_iter = CP_MAX_ITER,
        cp_tol = CP_TOL,
        mc.cores = MC_CORES
      )
      
      raw_list[[idx]]$experiment <- "dense_heteroskedastic_RP"
      raw_list[[idx]]$p_scalar <- p
      raw_list[[idx]]$sample_multiplier <- sample_multiplier
      raw_list[[idx]]$hetero_strength <- 1
      raw_list[[idx]]$alpha_profile <- alpha
      
      idx <- idx + 1L
    }
  }

  raw_df <- do.call(rbind, raw_list)
  raw_df <- add_normalized_error_columns(raw_df)

  summary_normalized_l1 <- summarize_numeric_by(
    values = raw_df$normalized_l1,
    df = raw_df,
    group_vars = c("method", "p_scalar", "R"),
    metric_name = "normalized_l1"
  )
  summary_normalized_fro <- summarize_numeric_by(
    values = raw_df$normalized_fro,
    df = raw_df,
    group_vars = c("method", "p_scalar", "R"),
    metric_name = "normalized_fro"
  )

  write_results_csv(raw_df, file.path(exp_dir, "raw_results.csv"))
  write_results_csv(summary_normalized_l1, file.path(exp_dir, "summary_normalized_l1.csv"))
  write_results_csv(summary_normalized_fro, file.path(exp_dir, "summary_normalized_fro.csv"))

  invisible(list(
    raw = raw_df,
    summary_normalized_l1 = summary_normalized_l1,
    summary_normalized_fro = summary_normalized_fro
  ))
}

# -----------------------------
# Master runner
# -----------------------------

# Function: run_all_four_experiments_simulation
# Purpose : Run all four experiments and save CSV outputs under out_dir.
# Inputs  : out_dir = root output directory.
# Output  : Invisibly returns a named list of experiment outputs.
run_all_four_experiments_simulation <- function(out_dir = OUTPUT_DIR) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  exp1 <- run_experiment1_unscaled_thinning_comparison_simulation(out_dir = out_dir)
  exp2 <- run_experiment2_dense_heteroskedastic_simulation(out_dir = out_dir)
  exp3 <- run_experiment3_vary_heteroskedastic_simulation(out_dir = out_dir)
  exp4 <- run_experiment4_rank_dimension_scaling_simulation(out_dir = out_dir)

  invisible(list(
    experiment1 = exp1,
    experiment2 = exp2,
    experiment3 = exp3,
    experiment4 = exp4
  ))
}

# Uncomment the next line to run everything immediately after sourcing:
# all_results <- run_all_four_experiments_simulation()

all_results <- run_all_four_experiments_simulation()

