# sccomp_outlier_plot_helpers.R
# Plotting helpers for sccomp results with censored outlier sample × cell-group pairs shown as red squares.

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(forcats)
library(scales)
library(magrittr)
library(sccomp)

# Attach sccomp-censored outlier sample × cell-group pairs to sccomp model object.
attach_sccomp_outliers <- function(fit_sccomp, outliers_tbl) {
  sample_col <- as.character(attr(fit_sccomp, ".sample"))
  cell_col   <- as.character(attr(fit_sccomp, ".cell_group"))

  out_tbl <- outliers_tbl %>%
    as.data.frame() %>%
    dplyr::rename(
      sample = dplyr::any_of(c("sample_id", "unique_sample_ID")),
      cell_group = dplyr::any_of(c("cell_type", "celltype", "Celltype"))
    ) %>%
    dplyr::transmute(
      !!sample_col := as.character(.data[["sample"]]),
      !!cell_col   := as.character(.data[["cell_group"]]),
      outlier = TRUE
    ) %>%
    dplyr::distinct()

  attr(fit_sccomp, "outliers") <- out_tbl
  fit_sccomp
}


# Get observed proportions for each cell-group in each sample.
get_sccomp_data_proportion <- function(x) {
  count_data <- attr(x, "count_data")

  .cell_group <- attr(x, ".cell_group")
  .count      <- attr(x, ".count")
  .sample     <- attr(x, ".sample")

  data_proportion <- x %>%
    dplyr::select(-`factor`) %>%
    tidyr::pivot_wider(
      names_from = parameter,
      values_from = c(dplyr::contains("c_"), dplyr::contains("v_"))
    ) %>%
    dplyr::left_join(count_data, by = dplyr::join_by(!!.cell_group)) %>%
    dplyr::with_groups(
      !!.sample,
      ~ dplyr::mutate(.x, proportion = (!!.count) / sum(!!.count))
    )

  if (!is.null(attr(x, "outliers"))) {
    data_proportion <- data_proportion %>%
      dplyr::left_join(attr(x, "outliers"), by = dplyr::join_by(!!.sample, !!.cell_group)) %>%
      dplyr::mutate(outlier = tidyr::replace_na(outlier, FALSE))
  } else {
    data_proportion <- data_proportion %>%
      dplyr::mutate(outlier = FALSE)
  }

  data_proportion
}


# Scatterplot of observed proportions vs. factor of interest, with outliers displayed as red squares.
plot_scatterplot_patched <- function(
  .data,
  data_proportion,
  factor_of_interest,
  .cell_group,
  .sample,
  significance_threshold = 0.05,
  my_theme = sccomp_theme(),
  marginal_pred_tbl = NULL,
  marginal_x_col = factor_of_interest
) {
  proportion <- outlier <- NULL

  S_sqrt <- function(x) sign(x) * sqrt(abs(x))
  IS_sqrt <- function(x) x^2 * sign(x)
  S_sqrt_trans <- function() scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

  .cell_group <- rlang::enquo(.cell_group)
  .sample <- rlang::enquo(.sample)
  cell_group_col <- rlang::as_name(.cell_group)

  p <- ggplot2::ggplot()

  # Optional marginal sccomp prediction line/ribbon.
  # This replaces the previous geom_smooth/sccomp_replicate visual.
  if (!is.null(marginal_pred_tbl)) {
    pred_df <- marginal_pred_tbl %>%
      as.data.frame()

    # Harmonize x column to factor_of_interest if needed
    if (marginal_x_col != factor_of_interest && marginal_x_col %in% colnames(pred_df)) {
      pred_df <- pred_df %>%
        dplyr::rename(!!factor_of_interest := dplyr::all_of(marginal_x_col))
    }

    stopifnot(
      all(c(
        cell_group_col,
        factor_of_interest,
        "proportion_mean",
        "proportion_lower",
        "proportion_upper"
      ) %in% colnames(pred_df))
    )

    pred_sum <- pred_df %>%
      dplyr::group_by(
        .data[[cell_group_col]],
        .data[[factor_of_interest]]
      ) %>%
      dplyr::summarise(
        proportion_mean  = mean(proportion_mean, na.rm = TRUE),
        proportion_lower = mean(proportion_lower, na.rm = TRUE),
        proportion_upper = mean(proportion_upper, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data[[cell_group_col]], .data[[factor_of_interest]])

    p <- p +
      ggplot2::geom_ribbon(
        data = pred_sum,
        ggplot2::aes(
          x = .data[[factor_of_interest]],
          ymin = proportion_lower,
          ymax = proportion_upper
        ),
        fill = "grey70",
        alpha = 0.45,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(
        data = pred_sum,
        ggplot2::aes(
          x = .data[[factor_of_interest]],
          y = proportion_mean
        ),
        color = "black",
        linewidth = 0.45,
        inherit.aes = FALSE
      )
  }

  p +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data[[factor_of_interest]],
        y = proportion,
        shape = outlier,
        color = outlier
      ),
      data = data_proportion,
      position = ggplot2::position_jitter(width = 0.15, height = 0),
      size = 0.65,
      alpha = 0.75
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(!!.cell_group),
      scales = "free_y",
      nrow = 4
    ) +

    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "#e11f28")) +
    ggplot2::scale_y_continuous(trans = S_sqrt_trans(), labels = scales::label_percent(accuracy = 0.1)) +
    ggplot2::xlab(factor_of_interest) +
    ggplot2::ylab("Cell-group proportion") +
    ggplot2::labs(
      color = "sccomp outlier",
      shape = "sccomp outlier"
    ) +
    my_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
      title = ggplot2::element_text(size = 3),
      legend.position = "bottom",
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))
    )
}


# Boxplot of observed proportions vs. factor of interest, with outliers displayed as red squares.
# Significance threshold is displayed as a blue fill.
plot_cat_sccomp_patched <- function(
  .data,
  data_proportion,
  factor_of_interest,
  .cell_group,
  significance_threshold = 0.05,
  significance_statistic = c("pH0", "FDR"),
  my_theme = sccomp_theme()
) {
  significance_statistic <- match.arg(significance_statistic)

  .cell_group <- rlang::enquo(.cell_group)

  S_sqrt <- function(x) sign(x) * sqrt(abs(x))
  IS_sqrt <- function(x) x^2 * sign(x)
  S_sqrt_trans <- function() scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

  # significance label per cell group for this factor
  sig_col <- if (significance_statistic == "FDR") "c_FDR" else "c_pH0"

  # This significance-fill mapping is tailored to Para_bin-style contrasts with "1-3" and "4+" levels.
  sig_df <- .data %>%
    dplyr::filter(
      factor == factor_of_interest,
      !is.na(.data[[sig_col]])
    ) %>%
    dplyr::mutate(
      sig_level = dplyr::case_when(
        parameter == paste0(factor_of_interest, "1-3") ~ "1-3",
        parameter == paste0(factor_of_interest, "4+")  ~ "4+",
        TRUE ~ NA_character_
      ),
      significant = .data[[sig_col]] < significance_threshold
    ) %>%
    dplyr::filter(!is.na(sig_level)) %>%
    dplyr::group_by(!!.cell_group, sig_level) %>%
    dplyr::summarise(
      significant = any(significant, na.rm = TRUE),
      .groups = "drop"
    )

  plot_df <- data_proportion %>%
    dplyr::mutate(
      sig_level = as.character(.data[[factor_of_interest]])
    ) %>%
    dplyr::left_join(
      sig_df,
      by = c(rlang::as_name(.cell_group), "sig_level")
    ) %>%
    dplyr::mutate(
      significant = tidyr::replace_na(significant, FALSE),
      outlier = tidyr::replace_na(outlier, FALSE)
    )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[factor_of_interest]], y = proportion)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = significant),
      outlier.shape = NA,
      linewidth = 0.25,
      width = 0.65,
      alpha = 0.45
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = outlier, color = outlier),
      position = ggplot2::position_jitter(width = 0.15, height = 0),
      size = 0.65,
      alpha = 0.75
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(!!.cell_group),
      scales = "free_y",
      nrow = 4
    ) +
    ggplot2::scale_y_continuous(trans = S_sqrt_trans(), labels = scales::label_percent(accuracy = 0.1)) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#9ecae1")) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "#e11f28")) +
    ggplot2::labs(
      x = factor_of_interest,
      y = "Cell-group proportion",
      fill = paste0(significance_statistic, " < ", significance_threshold),
      color = "sccomp outlier",
      shape = "sccomp outlier"
    ) +
    my_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
      title = ggplot2::element_text(size = 3),
      legend.position = "bottom",
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))
    )
}


# Main plotting function for sccomp models with outliers displayed as red squares.
plot_sccomp_patched <- function(
  x,
  significance_threshold = 0.05,
  significance_statistic = c("pH0", "FDR"),
  show_fdr_message = FALSE,
  marginal_predictions = NULL,
  marginal_x_cols = NULL,
  ...
) {
  significance_statistic <- match.arg(significance_statistic)

  .cell_group <- attr(x, ".cell_group")
  .sample <- attr(x, ".sample")

  if (x %>% dplyr::select(dplyr::ends_with("FDR")) %>% ncol() == 0) {
    stop("sccomp says: to produce plots, run sccomp_test() first.")
  }

  data_proportion <- get_sccomp_data_proportion(x)

  my_factors <- x %>%
    dplyr::filter(!is.na(`factor`)) %>%
    dplyr::distinct(`factor`) %>%
    dplyr::pull(`factor`)

  plots <- list()

  plots$boxplot <- purrr::map(
    my_factors,
    function(fac) {
      marginal_pred_tbl <- NULL
      marginal_x_col <- fac

      if (!is.null(marginal_predictions) && fac %in% names(marginal_predictions)) {
        marginal_pred_tbl <- marginal_predictions[[fac]]
      }

      if (!is.null(marginal_x_cols) && fac %in% names(marginal_x_cols)) {
        marginal_x_col <- marginal_x_cols[[fac]]
      }

      if (data_proportion %>%
          dplyr::select(dplyr::all_of(fac)) %>%
          dplyr::pull(1) %>%
          is.numeric()) {

        p <- plot_scatterplot_patched(
          .data = x,
          data_proportion = data_proportion,
          factor_of_interest = fac,
          .cell_group = !!rlang::sym(.cell_group),
          .sample = !!rlang::sym(.sample),
          my_theme = sccomp_theme(),
          significance_threshold = significance_threshold,
          marginal_pred_tbl = marginal_pred_tbl,
          marginal_x_col = marginal_x_col
        )

      } else {
        p <- plot_cat_sccomp_patched(
          .data = x,
          data_proportion = data_proportion,
          factor_of_interest = fac,
          .cell_group = !!rlang::sym(.cell_group),
          significance_threshold = significance_threshold,
          significance_statistic = significance_statistic,
          my_theme = sccomp_theme()
        )
      }

      # p + ggplot2::ggtitle(
      #   sprintf(
      #     "Grouped by %s (for multi-factor models, associations may be hard to see with unidimensional stratification)",
      #     fac
      #   )
      # )
    }
  )

  names(plots$boxplot) <- my_factors

  plots$credible_intervals_1D <- plot_1D_intervals(
    .data = x,
    significance_threshold = significance_threshold,
    significance_statistic = significance_statistic,
    show_fdr_message = show_fdr_message
  )

  if ("v_effect" %in% colnames(x) && nrow(x %>% dplyr::filter(!is.na(v_effect))) > 0) {
    plots$credible_intervals_2D <- plot_2D_intervals(
      .data = x,
      significance_threshold = significance_threshold,
      significance_statistic = significance_statistic,
      show_fdr_message = show_fdr_message
    )
  }

  plots
}
