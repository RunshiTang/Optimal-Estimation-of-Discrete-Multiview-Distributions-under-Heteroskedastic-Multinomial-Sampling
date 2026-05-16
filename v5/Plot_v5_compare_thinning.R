# ============================================================
# Plot the first four simulation experiments (v5) for the multiview
# density estimation paper.
#
# This script does plotting ONLY.
# It reads CSV files created by run_four_experiments_simulation_v5.R
# and saves PDF figures using ggplot2.
#
# v5 fix:
#   Save all figures under root_dir/figures_pdf/<experiment_name>/
#   so figures are separated from the CSV output folders.
#
# v3 fixes:
#   1. Robust type conversion after read.csv().
#   2. Explicit ordering of rows by x-variable within each method.
#   3. Correct handling of numeric x-axes on log scale.
#   4. Error-bar width bug fixed (width = 0 for numeric x).
#   5. Clear tick marks for all requested x values.
# ============================================================
setwd("~/density/v5")
set.seed(0)
# -----------------------------
# User configuration
# -----------------------------

# Install once if needed:
# install.packages("ggplot2")

SIM_OUTPUT_DIR <- "multiview_four_experiments_v5"
FIGURE_ROOT_DIR_NAME <- "figures_pdf"

# Function: get_figure_dir
# Purpose : Build a figure-output directory that is separate from the CSV folders.
# Inputs  : root_dir = root simulation-output directory;
#           experiment_name = experiment subdirectory name.
# Output  : Character scalar giving the figure directory path.
get_figure_dir <- function(root_dir, experiment_name) {
  file.path(root_dir, FIGURE_ROOT_DIR_NAME, experiment_name)
}
METHOD_LEVELS <- c("histogram", "unscaled_no_thinning", "unscaled_thinning", "unscaled", "oracle", "oracle_cp", "slice_oracle", "slice_est")

# -----------------------------
# Package helper
# -----------------------------

# Function: require_ggplot2
# Purpose : Stop with an informative message if ggplot2 is unavailable.
# Inputs  : None.
# Output  : Invisibly returns TRUE when ggplot2 is installed.
require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please run install.packages('ggplot2').")
  }
  invisible(TRUE)
}

require_ggpubr <- function() {
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required. Please run install.packages('ggpubr').")
  }
  invisible(TRUE)
}

require_dplyr <- function() {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Please run install.packages('dplyr').")
  }
  suppressPackageStartupMessages(library(dplyr))
  invisible(TRUE)
}

# -----------------------------
# File / factor helpers
# -----------------------------

# Function: format_large_number
# Purpose : Pretty labels for large numeric axes.
# Inputs  : x = numeric vector.
# Output  : Character vector with commas and no scientific notation.
format_large_number <- function(x) {
  format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
}

# Function: read_experiment_csv
# Purpose : Read one CSV file from an experiment subdirectory and
#           robustly convert columns to their natural types.
# Inputs  : root_dir = root simulation-output directory;
#           experiment_dir = experiment subdirectory name;
#           filename = CSV filename.
# Output  : Data frame read from the CSV file.
read_experiment_csv <- function(root_dir, experiment_dir, filename) {
  path <- file.path(root_dir, experiment_dir, filename)
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path))
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  df[] <- lapply(df, function(z) utils::type.convert(z, as.is = TRUE))
  df
}

# Function: ensure_method_factor
# Purpose : Convert the method column to a factor with a fixed plotting order.
# Inputs  : df = data frame;
#           method_col = name of the method column.
# Output  : The same data frame with the method column converted to a factor.
ensure_method_factor <- function(df, method_col = "method") {
  if (method_col %in% names(df)) {
    df[[method_col]] <- factor(df[[method_col]], levels = METHOD_LEVELS)
  }
  df
}

# Function: add_algorithm_method_label
# Purpose : For Experiment 1, convert method + algorithm into a single plotting
#           method so the two unscaled algorithm options appear as separate curves.
# Inputs  : df = data frame with method and algorithm columns.
# Output  : Updated data frame with a relabeled method column.
add_algorithm_method_label <- function(df) {
  if (all(c("method", "algorithm") %in% names(df))) {
    df$method <- ifelse(
      df$method == "unscaled" & df$algorithm == "no_thinning",
      "unscaled_no_thinning",
      ifelse(
        df$method == "unscaled" & df$algorithm == "thinning",
        "unscaled_thinning",
        df$method
      )
    )
  }
  ensure_method_factor(df)
}

# Function: coerce_numeric_columns
# Purpose : Force selected columns to numeric when present.
# Inputs  : df = data frame;
#           cols = character vector of column names.
# Output  : Updated data frame.
coerce_numeric_columns <- function(df, cols) {
  for (nm in cols) {
    if (nm %in% names(df)) {
      df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    }
  }
  df
}

# Function: order_for_plotting
# Purpose : Sort rows by method, optional facet variable, and x variable.
# Inputs  : df = data frame;
#           x_var = x-axis column name;
#           facet_var = optional facet column name.
# Output  : Sorted data frame.
order_for_plotting <- function(df, x_var, facet_var = NULL) {
  ord_cols <- c("method", facet_var, x_var)
  ord_cols <- ord_cols[!is.na(ord_cols) & ord_cols %in% names(df)]
  if (length(ord_cols) > 0) {
    # IMPORTANT: strip names before do.call(order, ...).
    # Otherwise a column named "method" is passed as order(method = ...),
    # which collides with base::order(method = c("auto", "shell", "radix"))
    # and triggers: Error in match.arg(method).
    ord_args <- unname(as.list(df[ord_cols]))
    ord <- do.call(base::order, ord_args)
    df <- df[ord, , drop = FALSE]
  }
  rownames(df) <- NULL
  df
}

# Function: add_deterministic_method_offset
# Purpose : Add a deterministic method-specific x-offset for numeric curves.
#           This replaces random jitter. All geoms for the same method use the
#           same horizontal offset, so lines, points, and error bars align.
# Inputs  : df = data frame;
#           x_var = original x-axis column;
#           method_col = method column;
#           use_log10_x = whether the x-axis will be shown on a log10 scale;
#           offset_frac = additive offset as a fraction of the minimum x-gap;
#           log10_offset = multiplicative offset size in log10 units.
# Output  : Updated data frame with .x_plot column.
add_deterministic_method_offset <- function(df,
                                            x_var,
                                            method_col = "method",
                                            use_log10_x = FALSE,
                                            offset_frac = 0.08,
                                            log10_offset = 0.06) {
  if (!(x_var %in% names(df))) stop("x_var must be a column of df.")
  if (!(method_col %in% names(df))) stop("method_col must be a column of df.")

  method_factor <- droplevels(as.factor(df[[method_col]]))
  method_levels <- levels(method_factor)
  n_methods <- length(method_levels)

  if (n_methods <= 1L) {
    df$.x_plot <- df[[x_var]]
    return(df)
  }

  method_index <- match(as.character(df[[method_col]]), method_levels)
  offset_index <- method_index - (n_methods + 1) / 2
  offset_index <- offset_index / max(abs(offset_index))

  x <- as.numeric(df[[x_var]])
  if (use_log10_x && all(x > 0, na.rm = TRUE)) {
    # On a log10 axis, use multiplicative offsets. This gives visually
    # comparable separation at small and large sample sizes.
    df$.x_plot <- 10^(log10(x) + offset_index * log10_offset)
  } else {
    ux <- sort(unique(x[is.finite(x)]))
    gap <- if (length(ux) >= 2L) min(diff(ux)) else NA_real_
    if (!is.finite(gap) || gap <= 0) {
      gap <- max(abs(x), na.rm = TRUE)
    }
    if (!is.finite(gap) || gap <= 0) gap <- 1
    df$.x_plot <- x + offset_index * offset_frac * gap
  }

  df
}


# Function: save_plot_pdf
# Purpose : Save a ggplot object as a PDF file.
# Inputs  : plot_obj = ggplot object;
#           filename = output path ending in .pdf;
#           width, height = figure size in inches.
# Output  : Invisibly returns filename.
save_plot_pdf <- function(plot_obj, filename, width = 8, height = 5) {
  require_ggplot2()
  ggplot2::ggsave(
    filename = filename,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf
  )
  invisible(filename)
}

# -----------------------------
# Generic plotting helpers
# -----------------------------

# Function: plot_summary_curve
# Purpose : Create a ggplot line plot with points and standard-error bars.
# Inputs  : summary_df = summary data frame containing x_var, method, mean, and se;
#           x_var = x-axis column name;
#           title = plot title;
#           xlab = x-axis label;
#           ylab = y-axis label;
#           facet_var = optional faceting variable;
#           use_log10_x = whether to apply a log10 x-axis.
# Output  : ggplot object.
plot_summary_curve <- function(summary_df,
                               x_var,
                               title = NULL,
                               xlab = NULL,
                               ylab = NULL,
                               facet_var = NULL,
                               use_log10_x = FALSE) {
  require_ggplot2()
  
  summary_df <- ensure_method_factor(summary_df)
  summary_df$ymin <- pmax(summary_df$mean - summary_df$sd, 0)
  summary_df$ymax <- summary_df$mean + summary_df$sd

  numeric_x <- is.numeric(summary_df[[x_var]]) || is.integer(summary_df[[x_var]])
  if (numeric_x) {
    summary_df[[x_var]] <- as.numeric(summary_df[[x_var]])
  }

  summary_df <- order_for_plotting(summary_df, x_var = x_var, facet_var = facet_var)

  if (numeric_x) {
    # Use deterministic method-specific offsets instead of random jitter.
    # This keeps all geoms of the same method/color aligned and reproducible.
    summary_df <- add_deterministic_method_offset(
      summary_df,
      x_var = x_var,
      method_col = "method",
      use_log10_x = use_log10_x
    )

    p <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = .data[[".x_plot"]], y = .data[["mean"]], color = .data[["method"]], group = .data[["method"]])
    ) +
      ggplot2::geom_line(linewidth = 0.7, alpha = 0.9) +
      ggplot2::geom_point(size = 2.2) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data[["ymin"]], ymax = .data[["ymax"]]),
        linewidth = 0.4,
        alpha = 0.8
      )

    if (use_log10_x) {
      p <- p + ggplot2::scale_x_continuous(
        trans = "log10",
        breaks = sort(unique(summary_df[[x_var]])),
        labels = format_large_number
      )
    } else {
      p <- p + ggplot2::scale_x_continuous(
        breaks = sort(unique(summary_df[[x_var]])),
        labels = format_large_number
      )
    }
  } else {
    dodge <- ggplot2::position_dodge(width = 0.35)
    p <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = .data[[x_var]], y = .data[["mean"]], 
                   color = .data[["method"]], group = .data[["method"]])
    ) +
      ggplot2::geom_line(linewidth = 0.7, alpha = 0.9, position = dodge) +
      ggplot2::geom_point(size = 2.2, position = dodge) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[["ymin"]], ymax = .data[["ymax"]]),
        width = 0.15,
        linewidth = 0.4,
        alpha = 0.8,
        position = dodge
      )
  }

  method_colors <- c(
    "histogram" = "#F8766D",
    "unscaled_no_thinning" = "#F564E3",
    "unscaled_thinning" = "#00BFC4",
    "unscaled" = "#F564E3",
    "oracle" = "#B79F00",
    "oracle_cp" = "#00BA38",
    "slice_est" = "#00BFC4",
    "slice_oracle" = "#619CFF"
  )
  method_labels <- c(
    "histogram" = "Histogram",
    "unscaled_no_thinning" = "Unscaled (no thinning)",
    "unscaled_thinning" = "Unscaled (thinning)",
    "unscaled" = "Unscaled",
    "oracle" = "Oracle",
    "oracle_cp" = "Oracle CP",
    "slice_est" = "Slice est.",
    "slice_oracle" = "Slice oracle"
  )
  
  p <- p +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = "Method") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = if (numeric_x && use_log10_x) 45 else 0, hjust = 1)
    )+ 
    ggplot2::scale_color_manual(values = method_colors, labels = method_labels)

  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)))
  }

  p
}

# Function: plot_ratio_curve
# Purpose : Plot ratios relative to the oracle method.
# Inputs  : ratio_df = data frame containing method, x_var, ratio_to_ref, and ref_value;
#           x_var = x-axis column name;
#           title = plot title;
#           xlab = x-axis label;
#           ylab = y-axis label.
# Output  : ggplot object.
plot_ratio_curve <- function(ratio_df,
                             x_var,
                             title = NULL,
                             xlab = NULL,
                             ylab = NULL,
                             use_log10_x = FALSE) {
  require_ggplot2()

  ratio_df <- ratio_df[ratio_df$method != "oracle", , drop = FALSE]
  ratio_df <- ensure_method_factor(ratio_df)

  numeric_x <- is.numeric(ratio_df[[x_var]]) || is.integer(ratio_df[[x_var]])
  if (numeric_x) {
    ratio_df[[x_var]] <- as.numeric(ratio_df[[x_var]])
  }
  ratio_df <- order_for_plotting(ratio_df, x_var = x_var)

  if (numeric_x) {
    p <- ggplot2::ggplot(
      ratio_df,
      ggplot2::aes(x = .data[[x_var]], y = .data[["ratio_to_ref"]], color = .data[["method"]], group = .data[["method"]])
    ) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 2.2)

    if (use_log10_x) {
      p <- p + ggplot2::scale_x_continuous(
        trans = "log10",
        breaks = sort(unique(ratio_df[[x_var]])),
        labels = format_large_number
      )
    } else {
      p <- p + ggplot2::scale_x_continuous(
        breaks = sort(unique(ratio_df[[x_var]])),
        labels = format_large_number
      )
    }
  } else {
    dodge <- ggplot2::position_dodge(width = 0.35)
    p <- ggplot2::ggplot(
      ratio_df,
      ggplot2::aes(x = .data[[x_var]], y = .data[["ratio_to_ref"]], color = .data[["method"]], group = .data[["method"]])
    ) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2) +
      ggplot2::geom_line(linewidth = 0.7, position = dodge) +
      ggplot2::geom_point(size = 2.2, position = dodge)
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = "Method") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = if (numeric_x && use_log10_x) 45 else 0, hjust = 1)
    )
}

# -----------------------------
# Experiment-specific plotters
# -----------------------------


# Function: plot_experiment1
# Purpose : Plot Experiment 1, comparing the unscaled estimator under the
#           no-thinning and multinomial-thinning algorithm options.
# Inputs  : root_dir = root directory containing experiment CSV files.
# Output  : Invisibly returns a named list of ggplot objects.
plot_experiment1 <- function(root_dir = SIM_OUTPUT_DIR) {
  require_dplyr()
  exp_name <- "experiment1_unscaled_thinning_comparison"
  fig_dir <- get_figure_dir(root_dir, exp_name)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  sum_l1 <- read_experiment_csv(root_dir, exp_name, "summary_l1.csv")
  sum_fro <- read_experiment_csv(root_dir, exp_name, "summary_fro.csv")
  sum_nl1 <- read_experiment_csv(root_dir, exp_name, "summary_normalized_l1.csv")
  sum_nfro <- read_experiment_csv(root_dir, exp_name, "summary_normalized_fro.csv")

  for (obj_name in c("sum_l1", "sum_fro", "sum_nl1", "sum_nfro")) {
    obj <- get(obj_name)
    obj <- coerce_numeric_columns(obj, c("n", "sample_multiplier", "hetero_strength", "mean", "sd", "se"))
    obj <- add_algorithm_method_label(obj)
    assign(obj_name, obj)
  }

  require_ggpubr()

  p_l1 <- plot_summary_curve(
    sum_l1,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = TRUE
  )
  p_fro <- plot_summary_curve(
    sum_fro,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean Frobenius error",
    use_log10_x = TRUE
  )
  p_nl1 <- plot_summary_curve(
    sum_nl1,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = TRUE
  )
  p_nfro <- plot_summary_curve(
    sum_nfro,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = TRUE
  )

  save_plot_pdf(p_l1, file.path(fig_dir, "experiment1_unscaled_thinning_l1_vs_n.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_fro, file.path(fig_dir, "experiment1_unscaled_thinning_fro_vs_n.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_nl1, file.path(fig_dir, "experiment1_unscaled_thinning_normalized_l1_vs_n.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_nfro, file.path(fig_dir, "experiment1_unscaled_thinning_normalized_fro_vs_n.pdf"), width = 8.8, height = 5.4)

  p_f <- ggpubr::ggarrange(p_fro, p_nfro, nrow = 1, common.legend = TRUE)
  p_l <- ggpubr::ggarrange(p_l1, p_nl1, nrow = 1, common.legend = TRUE)

  save_plot_pdf(p_f, file.path(fig_dir, "to_use_experiment1_unscaled_thinning_F_error.pdf"), width = 9, height = 4)
  save_plot_pdf(p_l, file.path(fig_dir, "to_use_experiment1_unscaled_thinning_l1_error.pdf"), width = 9, height = 4)

  invisible(list(l1 = p_l1, fro = p_fro, normalized_l1 = p_nl1, normalized_fro = p_nfro))
}

# Function: plot_experiment2
# Purpose : Plot Experiment 2 from saved CSV files and write PDF figures.
# Inputs  : root_dir = root directory containing experiment CSV files.
# Output  : Invisibly returns a named list of ggplot objects.
plot_experiment2 <- function(root_dir = SIM_OUTPUT_DIR) {
  require_dplyr()
  exp_dir <- file.path(root_dir, "experiment2_dense_heteroskedastic")
  fig_dir <- get_figure_dir(root_dir, "experiment2_dense_heteroskedastic")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  sum_l1 <- read_experiment_csv(root_dir, "experiment2_dense_heteroskedastic", "summary_l1.csv")
  sum_fro <- read_experiment_csv(root_dir, "experiment2_dense_heteroskedastic", "summary_fro.csv")
  sum_nl1 <- read_experiment_csv(root_dir, "experiment2_dense_heteroskedastic", "summary_normalized_l1.csv")
  sum_nfro <- read_experiment_csv(root_dir, "experiment2_dense_heteroskedastic", "summary_normalized_fro.csv")
  
  sum_l1$method = factor(sum_l1$method, levels = c("oracle", 
                                                   "slice_est", 
                                                   "slice_oracle",
                                                   "histogram",
                                                   "oracle_cp",
                                                   "unscaled"))
  sum_fro$method = factor(sum_fro$method, levels = c("oracle", 
                                                   "slice_est", 
                                                   "slice_oracle",
                                                   "histogram",
                                                   "oracle_cp",
                                                   "unscaled"))
  sum_nl1$method = factor(sum_nl1$method, levels = c("oracle", 
                                                   "slice_est", 
                                                   "slice_oracle",
                                                   "histogram",
                                                   "oracle_cp",
                                                   "unscaled"))
  sum_nfro$method = factor(sum_nfro$method, levels = c("oracle", 
                                                   "slice_est", 
                                                   "slice_oracle",
                                                   "histogram",
                                                   "oracle_cp",
                                                   "unscaled"))
  
  for (obj_name in c("sum_l1", "sum_fro", "sum_nl1", "sum_nfro")) {
    obj <- get(obj_name)
    obj <- coerce_numeric_columns(obj, c("n", "sample_multiplier", "mean", "sd", "se"))
    assign(obj_name, obj)
  }
  
  require_ggpubr()
  
  p_fro <- plot_summary_curve(
    sum_fro%>% filter(method %in% c("unscaled", "histogram")),
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean Frobenius error",
    use_log10_x = TRUE
  )
  p_nfro <- plot_summary_curve(
    sum_nfro %>% filter(method %in% c("unscaled", "histogram")),
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = TRUE
  )
  p = ggpubr::ggarrange(p_fro, p_nfro, nrow = 1, common.legend = T)
  
  save_plot_pdf(p, file.path(fig_dir, "to_use_experiment2_F_error.pdf"), width = 9, height = 4)
  
  
  p_l1 <- plot_summary_curve(
    sum_l1%>% filter(!(method %in% c("unscaled","oracle_cp")), n != 1e+08),
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = TRUE
  )
  
  p_l11 <- plot_summary_curve(
    sum_l1%>% filter(!(method %in% c("unscaled","oracle_cp", "histogram")), n != 1e+08),
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = TRUE
  )
  p_nl11 <- plot_summary_curve(
    sum_nl1 %>% filter(!(method %in% c("unscaled","oracle_cp", "histogram")), n != 1e+08),
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = TRUE
  )
  
  p = ggpubr::ggarrange(p_l1 , p_l11, p_nl11, nrow = 1, common.legend = T)
  
  save_plot_pdf(p, file.path(fig_dir, "to_use_experiment2_l1_error.pdf"), width = 12, height = 4)
  
  
  p_l1 <- plot_summary_curve(
    sum_l1,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = TRUE
  )
  p_fro <- plot_summary_curve(
    sum_fro,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean Frobenius error",
    use_log10_x = TRUE
  )
  p_nl1 <- plot_summary_curve(
    sum_nl1,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = TRUE
  )
  p_nfro <- plot_summary_curve(
    sum_nfro,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = TRUE
  )

  save_plot_pdf(p_l1, file.path(fig_dir, "experiment2_l1_vs_n.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_fro, file.path(fig_dir, "experiment2_fro_vs_n.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_nl1, file.path(fig_dir, "experiment2_normalized_l1_vs_n.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_nfro, file.path(fig_dir, "experiment2_normalized_fro_vs_n.pdf"), width = 8.8, height = 5.4)
  
  
  sum_l1 = sum_l1 %>% filter(method != "oracle_cp")
  sum_fro = sum_fro %>% filter(method != "oracle_cp")
  sum_nl1 = sum_nl1 %>% filter(method != "oracle_cp")
  sum_nfro = sum_nfro %>% filter(method != "oracle_cp")
  
  p_l1 <- plot_summary_curve(
    sum_l1,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = TRUE
  )
  p_fro <- plot_summary_curve(
    sum_fro,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean Frobenius error",
    use_log10_x = TRUE
  )
  p_nl1 <- plot_summary_curve(
    sum_nl1,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = TRUE
  )
  p_nfro <- plot_summary_curve(
    sum_nfro,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = TRUE
  )
  
  save_plot_pdf(p_l1, file.path(fig_dir, "experiment2_l1_vs_n_wo_cp.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_fro, file.path(fig_dir, "experiment2_fro_vs_n_wo_cp.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_nl1, file.path(fig_dir, "experiment2_normalized_l1_vs_n_wo_cp.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_nfro, file.path(fig_dir, "experiment2_normalized_fro_vs_n_wo_cp.pdf"), width = 8.8, height = 5.4)
  
  sum_l1 = sum_l1 %>% filter(method != "histogram")
  sum_fro = sum_fro %>% filter(method != "histogram")
  sum_nl1 = sum_nl1 %>% filter(method != "histogram")
  sum_nfro = sum_nfro %>% filter(method != "histogram")
  
  p_l1 <- plot_summary_curve(
    sum_l1,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = TRUE
  )
  p_fro <- plot_summary_curve(
    sum_fro,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean Frobenius error",
    use_log10_x = TRUE
  )
  p_nl1 <- plot_summary_curve(
    sum_nl1,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = TRUE
  )
  p_nfro <- plot_summary_curve(
    sum_nfro,
    x_var = "n",
    title = NULL,
    xlab = "Sample size n",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = TRUE
  )
  
  save_plot_pdf(p_l1, file.path(fig_dir, "experiment2_l1_vs_n_wo_cp_hist.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_fro, file.path(fig_dir, "experiment2_fro_vs_n_wo_cp_hist.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_nl1, file.path(fig_dir, "experiment2_normalized_l1_vs_n_wo_cp_hist.pdf"), width = 8.8, height = 5.4)
  save_plot_pdf(p_nfro, file.path(fig_dir, "experiment2_normalized_fro_vs_n_wo_cp_hist.pdf"), width = 8.8, height = 5.4)
  
  
  invisible(list(l1 = p_l1, fro = p_fro, normalized_l1 = p_nl1, normalized_fro = p_nfro))
}

# Function: plot_experiment3
# Purpose : Plot Experiment 3 from saved CSV files and write PDF figures.
# Inputs  : root_dir = root directory containing experiment CSV files.
# Output  : Invisibly returns a named list of ggplot objects.
plot_experiment3 <- function(root_dir = SIM_OUTPUT_DIR) {
  require_dplyr()
  exp_dir <- file.path(root_dir, "experiment3_vary_heteroskedastic")
  fig_dir <- get_figure_dir(root_dir, "experiment3_vary_heteroskedastic")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  sum_l1 <- read_experiment_csv(root_dir, "experiment3_vary_heteroskedastic", "summary_l1.csv")
  sum_fro <- read_experiment_csv(root_dir, "experiment3_vary_heteroskedastic", "summary_fro.csv")
  sum_nl1 <- read_experiment_csv(root_dir, "experiment3_vary_heteroskedastic", "summary_normalized_l1.csv")
  sum_nfro <- read_experiment_csv(root_dir, "experiment3_vary_heteroskedastic", "summary_normalized_fro.csv")

  for (obj_name in c("sum_l1", "sum_fro", "sum_nl1", "sum_nfro")) {
    obj <- get(obj_name)
    obj <- coerce_numeric_columns(obj, c("n", "sample_multiplier", "mean", "sd", "se"))
    assign(obj_name, obj)
  }
  
  
  p_fro <- plot_summary_curve(
    sum_fro%>% filter(method %in% c("unscaled", "histogram")),
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = "Mean Frobenius error",
    use_log10_x = TRUE
  )
  p_nfro <- plot_summary_curve(
    sum_nfro %>% filter(method %in% c("unscaled", "histogram")),
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = TRUE
  )+ expand_limits(y = 0)
  p = ggpubr::ggarrange(p_fro, p_nfro, nrow = 1, common.legend = T)
  
  save_plot_pdf(p, file.path(fig_dir, "to_use_experiment3_F_error.pdf"), width = 9, height = 4)
  
  
  p_l1 <- plot_summary_curve(
    sum_l1 %>% filter(!(method %in% c("unscaled"))),
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = T
  )
  p_l11 <- plot_summary_curve(
    sum_nl1 %>% filter(!(method %in% c("unscaled", "oracle_cp"))),
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = T
  )+ expand_limits(y = 0)
  # p_l111 <- plot_summary_curve(
  #   sum_l1 %>% filter(!(method %in% c("unscaled", "histogram", "oracle_cp"))),
  #   x_var = "hetero_strength",
  #   title = NULL,
  #   xlab = "Hetero strength",
  #   ylab = latex2exp::TeX("Mean $l_1$ error"), 
  #   use_log10_x = T
  # )
  p = ggpubr::ggarrange(p_l1, p_l11, nrow = 1, common.legend = T)
  
  save_plot_pdf(p, file.path(fig_dir, "to_use_experiment3_l1_error.pdf"), width = 9, height = 4)
  
  

  p_l1 <- plot_summary_curve(
    sum_l1,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = latex2exp::TeX("Mean $l_1$ error"),  # "Mean l1 error",
    use_log10_x = T
  )
  p_fro <- plot_summary_curve(
    sum_fro,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = "Mean Frobenius error",
    use_log10_x = T
  )
  p_ratio_l1 <- plot_summary_curve(
    sum_nl1,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = T
  )
  p_ratio_fro <- plot_summary_curve(
    sum_nfro,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = T
  )

  save_plot_pdf(p_l1, file.path(fig_dir, "experiment3_l1_vs_hetero_strength.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_fro, file.path(fig_dir, "experiment3_fro_vs_hetero_strength.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_ratio_l1, file.path(fig_dir, "experiment3_normalized_l1_vs_hetero_strength.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_ratio_fro, file.path(fig_dir, "experiment3_normalized_fro_vs_hetero_strength.pdf"), width = 8.2, height = 5.2)
  
  sum_l1 = sum_l1 %>% filter(method != "oracle_cp")
  sum_fro = sum_fro %>% filter(method != "oracle_cp")
  sum_nl1 = sum_nl1 %>% filter(method != "oracle_cp")
  sum_nfro = sum_nfro %>% filter(method != "oracle_cp")
  
  p_l1 <- plot_summary_curve(
    sum_l1,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = T
  )
  p_fro <- plot_summary_curve(
    sum_fro,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = "Mean Frobenius error",
    use_log10_x = T
  )
  p_ratio_l1 <- plot_summary_curve(
    sum_nl1,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = T
  )
  p_ratio_fro <- plot_summary_curve(
    sum_nfro,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = T
  )
  
  save_plot_pdf(p_l1, file.path(fig_dir, "experiment3_l1_vs_hetero_strength_wo_cp.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_fro, file.path(fig_dir, "experiment3_fro_vs_hetero_strength_wo_cp.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_ratio_l1, file.path(fig_dir, "experiment3_normalized_l1_vs_hetero_strength_wo_cp.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_ratio_fro, file.path(fig_dir, "experiment3_normalized_fro_vs_hetero_strength_wo_cp.pdf"), width = 8.2, height = 5.2)
  
  sum_l1 = sum_l1 %>% filter(method != "histogram")
  sum_fro = sum_fro %>% filter(method != "histogram")
  sum_nl1 = sum_nl1 %>% filter(method != "histogram")
  sum_nfro = sum_nfro %>% filter(method != "histogram")
  
  p_l1 <- plot_summary_curve(
    sum_l1,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = latex2exp::TeX("Mean $l_1$ error"), 
    use_log10_x = T
  )
  p_fro <- plot_summary_curve(
    sum_fro,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = "Mean Frobenius error",
    use_log10_x = T
  )
  p_ratio_l1 <- plot_summary_curve(
    sum_nl1,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    use_log10_x = T
  )
  p_ratio_fro <- plot_summary_curve(
    sum_nfro,
    x_var = "hetero_strength",
    title = NULL,
    xlab = "Hetero strength",
    ylab = "Mean normalized Frobenius error",
    use_log10_x = T
  )
  
  save_plot_pdf(p_l1, file.path(fig_dir, "experiment3_l1_vs_hetero_strength_wo_cp_hist.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_fro, file.path(fig_dir, "experiment3_fro_vs_hetero_strength_wo_cp_hist.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_ratio_l1, file.path(fig_dir, "experiment3_normalized_l1_vs_hetero_strength_wo_cp_hist.pdf"), width = 8.2, height = 5.2)
  save_plot_pdf(p_ratio_fro, file.path(fig_dir, "experiment3_normalized_fro_vs_hetero_strength_wo_cp_hist.pdf"), width = 8.2, height = 5.2)
  
  invisible(list(l1 = p_l1, fro = p_fro, ratio_l1 = p_ratio_l1, ratio_fro = p_ratio_fro))
}

# Function: plot_experiment4
# Purpose : Plot Experiment 4 from saved CSV files and write PDF figures.
# Inputs  : root_dir = root directory containing experiment CSV files.
# Output  : Invisibly returns a named list of ggplot objects.
plot_experiment4 <- function(root_dir = SIM_OUTPUT_DIR) {
  require_dplyr()
  exp_dir <- file.path(root_dir, "experiment4_rank_dimension_scaling")
  fig_dir <- get_figure_dir(root_dir, "experiment4_rank_dimension_scaling")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  sum_nl1 <- read_experiment_csv(root_dir, "experiment4_rank_dimension_scaling", "summary_normalized_l1.csv")
  sum_nfro <- read_experiment_csv(root_dir, "experiment4_rank_dimension_scaling", "summary_normalized_fro.csv")

  for (obj_name in c("sum_nl1", "sum_nfro")) {
    obj <- get(obj_name)
    obj <- coerce_numeric_columns(obj, c("p_scalar", "R", "mean", "sd", "se"))
    obj$R <- factor(obj$R, levels = sort(unique(obj$R)))
    assign(obj_name, obj)
  }
  
  sum_nl1$R = paste0("R = ", sum_nl1$R)
  sum_nfro$R = paste0("R = ", sum_nfro$R)
  
  p_nl1 <- plot_summary_curve(
    sum_nl1%>% filter(method != "unscaled"),
    x_var = "p_scalar",
    title = NULL,
    xlab = "Common mode dimension p",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    facet_var = "R"
  ) + ggplot2::scale_x_continuous(breaks = unique(sum_nl1$p_scalar))

  p_nfro <- plot_summary_curve(
    sum_nfro,
    x_var = "p_scalar",
    title = NULL,
    xlab = "Common mode dimension p",
    ylab = "Mean normalized Frobenius error",
    facet_var = "R"
  ) + ggplot2::scale_x_continuous(breaks = unique(sum_nfro$p_scalar))

  save_plot_pdf(p_nl1, file.path(fig_dir, "experiment4_normalized_l1_vs_p_facet_R.pdf"), width = 10, height = 6)
  save_plot_pdf(p_nfro, file.path(fig_dir, "experiment4_normalized_fro_vs_p_facet_R.pdf"), width = 10, height = 6)
  
  sum_nl1 = sum_nl1 %>% filter(method != "oracle_cp")
  sum_nfro = sum_nfro %>% filter(method != "oracle_cp")
  
  p_nl1 <- plot_summary_curve(
    sum_nl1%>% filter(method != "unscaled"),
    x_var = "p_scalar",
    title = NULL,
    xlab = "Common mode dimension p",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    facet_var = "R"
  ) + ggplot2::scale_x_continuous(breaks = unique(sum_nl1$p_scalar))
  
  p_nfro <- plot_summary_curve(
    sum_nfro,
    x_var = "p_scalar",
    title = NULL,
    xlab = "Common mode dimension p",
    ylab = "Mean normalized Frobenius error",
    facet_var = "R"
  ) + ggplot2::scale_x_continuous(breaks = unique(sum_nfro$p_scalar))
  
  save_plot_pdf(p_nl1, file.path(fig_dir, "experiment4_normalized_l1_vs_p_facet_R_wo_cp.pdf"), width = 10, height = 6)
  save_plot_pdf(p_nfro, file.path(fig_dir, "experiment4_normalized_fro_vs_p_facet_R_wo_cp.pdf"), width = 10, height = 6)
  
  sum_nl1 = sum_nl1 %>% filter(method != "histogram")
  sum_nfro = sum_nfro %>% filter(method != "histogram")
  
  p_nl1 <- plot_summary_curve(
    sum_nl1%>% filter(method != "unscaled"),
    x_var = "p_scalar",
    title = NULL,
    xlab = "Common mode dimension p",
    ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
    facet_var = "R"
  ) + ggplot2::scale_x_continuous(breaks = unique(sum_nl1$p_scalar))
  
  p_nfro <- plot_summary_curve(
    sum_nfro,
    x_var = "p_scalar",
    title = NULL,
    xlab = "Common mode dimension p",
    ylab = "Mean normalized Frobenius error",
    facet_var = "R"
  ) + ggplot2::scale_x_continuous(breaks = unique(sum_nfro$p_scalar))
  
  save_plot_pdf(p_nl1, file.path(fig_dir, "experiment4_normalized_l1_vs_p_facet_R_wo_cp_hist.pdf"), width = 10, height = 6)
  save_plot_pdf(p_nfro, file.path(fig_dir, "experiment4_normalized_fro_vs_p_facet_R_wo_cp_hist.pdf"), width = 10, height = 6)
  
  
  
  # p_nl1 <- plot_summary_curve(
  #   sum_nl1,
  #   x_var = "p_scalar",
  #   title = "Experiment 4: normalized l1 error vs p",
  #   xlab = "Common mode dimension p",
  #   ylab = latex2exp::TeX("Mean normalized $l_1$ error"),
  #   facet_var = "R"
  # ) + ggplot2::scale_x_continuous(breaks = c(20, 40, 60, 80)) + 
  #   lims(y = c(0,2.5))
  # 
  # save_plot_pdf(p_nl1, file.path(fig_dir, "experiment4_normalized_l1_vs_p_facet_R_2.pdf"), width = 10, height = 6)
  # 
  invisible(list(normalized_l1 = p_nl1, normalized_fro = p_nfro))
}

# -----------------------------
# Master plotter
# -----------------------------

# Function: plot_all_four_experiments
# Purpose : Read all CSV outputs and save all PDF figures.
# Inputs  : root_dir = root simulation-output directory.
# Output  : Invisibly returns a named list of ggplot objects.
plot_all_four_experiments <- function(root_dir = SIM_OUTPUT_DIR) {
  require_ggplot2()
  require_dplyr()

  exp1 <- plot_experiment1(root_dir = root_dir)
  exp2 <- plot_experiment2(root_dir = root_dir)
  exp3 <- plot_experiment3(root_dir = root_dir)
  exp4 <- plot_experiment4(root_dir = root_dir)

  invisible(list(
    experiment1 = exp1,
    experiment2 = exp2,
    experiment3 = exp3,
    experiment4 = exp4
  ))
}

# Uncomment the next line to create all PDF figures after sourcing:
all_plots <- plot_all_four_experiments()
