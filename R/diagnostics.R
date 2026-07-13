## pepVet Score Diagnostics
##
## Functions for quantifying multicollinearity, dimensionality, and
## component contributions in pepVet scoring models. Used to validate
## that the five (or six) component scores measure orthogonal
## digest-quality dimensions.

## Component color palette shared between diagnostics and visualization

.diagnostic_known_components <- function() {
  unique(c(
    names(.default_scoring_weights$protein_only),
    names(.default_scoring_weights$proteome_aware)
  ))
}

.validate_diagnostic_batch_result <- function(batch_result) {
  if (!inherits(batch_result, "data.frame")) {
    .abort(
      "{.arg batch_result} must be a tibble from {.fn batch_evaluate}.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  .validate_unique_columns(
    batch_result,
    "batch_result",
    class = "pepvet_error_invalid_diagnostics_input"
  )

  if (nrow(batch_result) < 2L) {
    .abort(
      "{.arg batch_result} must contain at least two proteins.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  comp_cols <- grep("^S_", names(batch_result), value = TRUE)

  if (length(comp_cols) < 2L) {
    .abort(
      paste0(
        "{.arg batch_result} must contain at least two ",
        "component score columns (S_*)."
      ),
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  if (anyDuplicated(comp_cols) > 0L) {
    .abort(
      "Component score columns in {.arg batch_result} must be unique.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  unknown_components <- setdiff(
    comp_cols, .diagnostic_known_components()
  )

  if (length(unknown_components) > 0L) {
    .abort(
      c(
        "{.arg batch_result} contains unknown component score columns.",
        "i" = "Unknown: {.val {unknown_components}}"
      ),
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  if (!all(vapply(batch_result[comp_cols], is.numeric, logical(1)))) {
    .abort(
      "Component score columns in {.arg batch_result} must be numeric.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  if (any(vapply(
    batch_result[comp_cols],
    function(values) {
      anyNA(values) ||
        any(!is.finite(values)) ||
        any(values < 0) ||
        any(values > 1)
    },
    logical(1)
  ))) {
    .abort(
      "Component scores in {.arg batch_result} must be finite values between 0 and 1.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  if (any(vapply(
    batch_result[comp_cols],
    function(values) length(unique(values)) < 2L,
    logical(1)
  ))) {
    .abort(
      "Component score columns in {.arg batch_result} must vary across proteins.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  if ("composite_score" %in% names(batch_result)) {
    composite <- batch_result$composite_score
    if (
      !is.numeric(composite) ||
        anyNA(composite) ||
        any(!is.finite(composite)) ||
        any(composite < 0) ||
        any(composite > 1)
    ) {
      .abort(
        "{.arg batch_result} contains invalid composite scores.",
        class = "pepvet_error_invalid_diagnostics_input"
      )
    }
  }

  if ("verdict" %in% names(batch_result)) {
    verdict <- batch_result$verdict
    if (
      !is.character(verdict) ||
        anyNA(verdict) ||
        any(!verdict %in% c("Good", "Moderate", "Poor"))
    ) {
      .abort(
        "{.arg batch_result} contains invalid verdict labels.",
        class = "pepvet_error_invalid_diagnostics_input"
      )
    }
  }

  comp_cols
}

.validate_diagnostic_weights <- function(weights, component_names) {
  if (
    !is.numeric(weights) ||
      length(weights) != length(component_names) ||
      is.null(names(weights)) ||
      anyNA(weights) ||
      any(!is.finite(weights))
  ) {
    .abort(
      "{.arg weights} must be a named finite numeric vector.",
      class = "pepvet_error_invalid_weights"
    )
  }

  weight_names <- names(weights)
  if (
    anyNA(weight_names) ||
      any(!nzchar(trimws(weight_names))) ||
      anyDuplicated(weight_names) > 0L ||
      !setequal(weight_names, component_names)
  ) {
    .abort(
      "Named {.arg weights} must match the component score columns.",
      class = "pepvet_error_invalid_weights"
    )
  }

  normalized_weights <- as.numeric(weights[match(component_names, weight_names)])
  if (any(normalized_weights < 0) || sum(normalized_weights) <= 0) {
    .abort(
      "{.arg weights} must contain non-negative values with a positive sum.",
      class = "pepvet_error_invalid_weights"
    )
  }

  names(normalized_weights) <- component_names
  normalized_weights / sum(normalized_weights)
}

#' Score diagnostics for pepVet scoring models
#'
#' `score_diagnostics()` runs three analyses on a [batch_evaluate()] result
#' to quantify multicollinearity, dimensionality, and component
#' contributions in the scoring model. Use it to validate that the
#' five (or six) component scores are measuring orthogonal dimensions
#' and that each contributes meaningfully to the composite.
#'
#' @param batch_result A tibble returned by [batch_evaluate()]. Must contain
#'   at least two rows and two component-score columns (prefixed `S_`). If
#'   `NULL`, too few rows, or missing required columns, raises an error.
#' @param weights Optional named numeric vector of component weights.
#'   When `NULL` (default), weights are inferred from the detected
#'   component columns using pepVet's default scoring weights
#'   (protein-only or proteome-aware, depending on presence of
#'   `S_unique`).
#'
#' @details Three analyses are performed:
#'
#' **VIF (Variance Inflation Factor):** For each component, a linear
#' regression is fit using all other components as predictors. VIF =
#' 1 / (1 - R^2). VIF < 5 is acceptable. VIF > 10 indicates problematic
#' collinearity that may warrant component removal or merger. VIF
#' requires more than `n_components + 1` observations; when this condition is
#' not met for an input with at least two rows, NA values are returned with a
#' warning. Perfect collinearity is
#' represented by `Inf` rather than an unclassed regression warning.
#'
#' **PCA (Principal Component Analysis):** Centered and scaled PCA on
#' the component-score matrix. The proportion of variance explained
#' by each principal component reveals the effective dimensionality
#' of the scoring model. If PC1 + PC2 explain > 80% of variance, the
#' model has low effective dimensionality (expected for a well-designed
#' multi-criteria score).
#'
#' **Ablation:** Each component is set to 0 (its minimum possible value)
#' for all proteins while others remain at their actual values, and the
#' composite score is recomputed. The mean and maximum drop in composite
#' score, plus the number of verdicts that flip, quantify each
#' component's marginal contribution. Components with mean drops
#' near zero are candidates for removal or down-weighting.
#'
#' @section Limitations:
#' Diagnostics are meaningful only on batch results with enough proteins
#' (at least 10-20 recommended; VIF requires more than `n_components + 1`).
#' Results are specific to the enzyme, missed-cleavage setting, and
#' protein set used. Running diagnostics on a different batch may
#' give different VIF and ablation values.
#'
#' @return A named list with six elements:
#' \describe{
#'   \item{\code{vif}}{A named numeric vector of VIF values, one per
#'     component. `NA` when too few proteins for reliable estimation and
#'     `Inf` for perfect collinearity.}
#'   \item{\code{pca}}{A list with \code{var_explained} (numeric vector),
#'     \code{loadings} (rotation matrix), \code{sdev} (standard deviations),
#'     and \code{x} (PCA scores matrix).}
#'   \item{\code{ablation}}{A data.frame with columns \code{component},
#'     \code{weight}, \code{mean_drop}, \code{sd_drop}, \code{max_drop},
#'     and \code{n_verdict_flipped}.}
#'   \item{\code{n_proteins}}{Number of proteins in the input.}
#'   \item{\code{n_components}}{Number of component scores detected.}
#'   \item{\code{weights}}{The resolved weight vector used for ablation.}
#' }
#'
#' @seealso [batch_evaluate()] for upstream batch evaluation,
#'   [plot_score_diagnostics()] for visualization.
#'
#' @family diagnostics
#'
#' @examples
#' small_fasta <- system.file(
#'   "extdata", "small_proteome_50_proteins.fasta",
#'   package = "pepVet"
#' )
#' batch <- batch_evaluate(small_fasta, enzyme = "trypsin")
#' diag <- score_diagnostics(batch)
#' diag$vif
#' diag$ablation
#'
#' # Proteome-aware mode
#' prot_digest <- digest_protein(small_fasta, enzyme = "trypsin")
#' batch_pa <- batch_evaluate(small_fasta, enzyme = "trypsin",
#'   proteome = prot_digest)
#' diag_pa <- score_diagnostics(batch_pa)
#' diag_pa$vif
#' @export
score_diagnostics <- function(batch_result, weights = NULL) {
  comp_cols <- .validate_diagnostic_batch_result(batch_result)

  if (is.null(weights)) {
    w0 <- if ("S_unique" %in% comp_cols) {
      .default_scoring_weights$proteome_aware
    } else {
      .default_scoring_weights$protein_only
    }
    w_sub <- w0[comp_cols]
  } else {
    w_sub <- .validate_diagnostic_weights(weights, comp_cols)
  }
  w_sub <- w_sub / sum(w_sub)

  component_matrix <- as.matrix(batch_result[, comp_cols])
  n_prot <- nrow(batch_result)
  n_comp <- length(comp_cols)

  ## VIF
  vif_vals <- setNames(numeric(n_comp), comp_cols)

  if (n_prot <= n_comp + 1L) {
    vif_vals[] <- NA_real_
    cli::cli_warn(
      paste0(
        "Too few proteins ({n_prot}) for reliable VIF estimation ",
        "(need > {n_comp + 1L}). Returning NA."
      ),
      class = "pepvet_warning_diagnostics_vif"
    )
  } else {
    for (i in seq_along(comp_cols)) {
      response <- comp_cols[i]
      predictors <- comp_cols[-i]
      formula_str <- paste(response, "~", paste(predictors, collapse = " + "))
      fit <- suppressWarnings(
        stats::lm(stats::as.formula(formula_str), data = batch_result)
      )
      r_squared <- suppressWarnings(summary(fit)$r.squared)
      vif_vals[i] <- if (!is.finite(r_squared)) {
        NA_real_
      } else if (r_squared >= 1 - .Machine$double.eps) {
        Inf
      } else {
        1 / (1 - max(0, r_squared))
      }
    }
  }

  ## PCA
  pca_result <- tryCatch(
    stats::prcomp(component_matrix, center = TRUE, scale. = TRUE),
    error = function(error) {
      .abort(
        "{.arg batch_result} could not be decomposed by PCA.",
        class = "pepvet_error_invalid_diagnostics_input"
      )
    }
  )
  var_exp <- pca_result$sdev^2 / sum(pca_result$sdev^2)

  ## Ablation: set each component to 0 and measure composite drop
  true_composite <- drop(component_matrix %*% w_sub)

  if ("verdict" %in% names(batch_result)) {
    true_verdict <- batch_result$verdict
  } else {
    true_verdict <- ifelse(
      true_composite >= .get_param("verdict_good"), "Good",
      ifelse(
        true_composite >= .get_param("verdict_moderate"),
        "Moderate", "Poor"
      )
    )
  }

  ablation_tbl <- data.frame(
    component          = comp_cols,
    weight             = unname(w_sub[comp_cols]),
    mean_drop          = NA_real_,
    sd_drop            = NA_real_,
    max_drop           = NA_real_,
    n_verdict_flipped  = NA_integer_,
    stringsAsFactors   = FALSE
  )

  for (i in seq_along(comp_cols)) {
    perturbed <- component_matrix
    perturbed[, i] <- 0
    perturbed_composite <- drop(perturbed %*% w_sub)
    drops <- true_composite - perturbed_composite

    ablation_tbl$mean_drop[i] <- mean(drops)
    ablation_tbl$sd_drop[i] <- stats::sd(drops)
    ablation_tbl$max_drop[i] <- max(drops)

    perturbed_verdict <- ifelse(
      perturbed_composite >= .get_param("verdict_good"), "Good",
      ifelse(
        perturbed_composite >= .get_param("verdict_moderate"),
        "Moderate", "Poor"
      )
    )
    ablation_tbl$n_verdict_flipped[i] <- sum(perturbed_verdict != true_verdict)
  }

  list(
    vif           = vif_vals,
    pca           = list(
      var_explained = var_exp,
      loadings      = pca_result$rotation,
      sdev          = pca_result$sdev,
      x             = pca_result$x
    ),
    ablation      = ablation_tbl,
    n_proteins    = n_prot,
    n_components  = n_comp,
    weights       = w_sub
  )
}


.validate_diagnostic_title <- function(title) {
  if (
    !is.null(title) &&
      (!is.character(title) || length(title) != 1L || is.na(title))
  ) {
    .abort(
      "{.arg title} must be a single character string or {.val NULL}.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  invisible(title)
}

.validate_diagnostic_result <- function(x) {
  required <- c(
    "vif", "pca", "ablation", "n_proteins", "n_components", "weights"
  )

  if (!is.list(x) || !all(required %in% names(x))) {
    .abort(
      "{.arg x} must be a list returned by {.fn score_diagnostics}.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  if (
      !is.numeric(x$n_proteins) ||
      length(x$n_proteins) != 1L ||
      is.na(x$n_proteins) ||
      !is.finite(x$n_proteins) ||
      x$n_proteins < 2 ||
      x$n_proteins > .Machine$integer.max ||
      x$n_proteins != floor(x$n_proteins) ||
      !is.numeric(x$n_components) ||
      length(x$n_components) != 1L ||
      is.na(x$n_components) ||
      !is.finite(x$n_components) ||
      x$n_components < 2 ||
      x$n_components > .Machine$integer.max ||
      x$n_components != floor(x$n_components)
  ) {
    .abort(
      "The diagnostics result contains invalid dimensions.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  n_prot <- as.integer(x$n_proteins)
  n_comp <- as.integer(x$n_components)

  if (
      !is.numeric(x$vif) ||
      length(x$vif) != n_comp ||
      is.null(names(x$vif)) ||
      anyNA(names(x$vif)) ||
      any(!nzchar(names(x$vif))) ||
      anyDuplicated(names(x$vif)) > 0L ||
      any(!names(x$vif) %in% .diagnostic_known_components()) ||
      any(is.nan(x$vif)) ||
      any(!is.na(x$vif) & x$vif < 0)
  ) {
    .abort(
      "The diagnostics result contains invalid VIF values.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  if (
    !is.numeric(x$weights) ||
      length(x$weights) != n_comp ||
      is.null(names(x$weights)) ||
      !identical(names(x$weights), names(x$vif)) ||
      anyNA(x$weights) ||
      any(!is.finite(x$weights)) ||
      any(x$weights < 0) ||
      !isTRUE(all.equal(sum(x$weights), 1, tolerance = 1e-10))
  ) {
    .abort(
      "The diagnostics result contains invalid resolved weights.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  if (!is.list(x$pca) ||
      !all(c("var_explained", "loadings", "sdev", "x") %in% names(x$pca))) {
    .abort(
      "The diagnostics result contains an invalid PCA result.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  var_explained <- x$pca$var_explained
  loadings <- x$pca$loadings
  sdev <- x$pca$sdev
  scores <- x$pca$x

  if (
    !is.numeric(var_explained) ||
      length(var_explained) < 2L ||
      anyNA(var_explained) ||
      any(!is.finite(var_explained)) ||
      any(var_explained < 0) ||
      !isTRUE(all.equal(sum(var_explained), 1, tolerance = 1e-10)) ||
      !is.matrix(loadings) ||
      nrow(loadings) != n_comp ||
      ncol(loadings) != length(var_explained) ||
      !is.numeric(loadings) ||
      anyNA(loadings) ||
      any(!is.finite(loadings)) ||
      !is.numeric(sdev) ||
      length(sdev) != length(var_explained) ||
      anyNA(sdev) ||
      any(!is.finite(sdev)) ||
      any(sdev < 0) ||
      !is.matrix(scores) ||
      nrow(scores) != n_prot ||
      ncol(scores) != length(var_explained) ||
      !is.numeric(scores) ||
      anyNA(scores) ||
      any(!is.finite(scores))
  ) {
    .abort(
      "The diagnostics result contains invalid PCA dimensions or values.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  required_ablation <- c(
    "component", "weight", "mean_drop", "sd_drop", "max_drop",
    "n_verdict_flipped"
  )

  if (
    !inherits(x$ablation, "data.frame") ||
      !all(required_ablation %in% names(x$ablation)) ||
      nrow(x$ablation) != n_comp ||
      !is.character(x$ablation$component) ||
      !identical(x$ablation$component, names(x$vif)) ||
      !is.numeric(x$ablation$weight) ||
      !is.numeric(x$ablation$mean_drop) ||
      !is.numeric(x$ablation$sd_drop) ||
      !is.numeric(x$ablation$max_drop) ||
      !is.numeric(x$ablation$n_verdict_flipped) ||
      anyNA(x$ablation[c(
        "weight", "mean_drop", "sd_drop", "max_drop", "n_verdict_flipped"
      )]) ||
      any(!is.finite(x$ablation$weight)) ||
      any(!is.finite(x$ablation$mean_drop)) ||
      any(!is.finite(x$ablation$sd_drop)) ||
      any(!is.finite(x$ablation$max_drop)) ||
      any(!is.finite(x$ablation$n_verdict_flipped)) ||
      any(x$ablation$weight < 0) ||
      any(x$ablation$mean_drop < 0) ||
      any(x$ablation$sd_drop < 0) ||
      any(x$ablation$max_drop < 0) ||
      any(x$ablation$n_verdict_flipped < 0) ||
      any(x$ablation$n_verdict_flipped > n_prot) ||
      any(x$ablation$n_verdict_flipped != floor(x$ablation$n_verdict_flipped)) ||
      !isTRUE(all.equal(
        unname(x$ablation$weight), unname(x$weights), tolerance = 1e-10
      ))
  ) {
    .abort(
      "The diagnostics result contains invalid ablation values.",
      class = "pepvet_error_invalid_diagnostics_input"
    )
  }

  invisible(x)
}


#' Plot score diagnostics
#'
#' `plot_score_diagnostics()` visualises the result of
#' [score_diagnostics()] as a three-panel figure: (A) VIF bar chart
#' with threshold lines at 5 and 10, (B) PCA scree plot with cumulative
#' variance line and an 80% reference, (C) ablation waterfall showing
#' the mean composite drop with error bars when each component is set
#' to zero, annotated with the number of verdicts that flip.
#'
#' @param x A list returned by [score_diagnostics()]. If `NULL` or
#'   missing required elements, raises an error.
#' @param title Optional character string for an overall figure title.
#'   When `NULL` (default), generates "Score Diagnostics" with the
#'   protein and component counts.
#'
#' @return A `patchwork` object with three panels: VIF (A), PCA (B),
#'   ablation (C).
#'
#' @section Limitations:
#' Shows what score_diagnostics() returns. Interpretation needs domain
#' knowledge.
#'
#' @seealso [score_diagnostics()] for the upstream analysis.
#'
#' @family diagnostics
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   small_fasta <- system.file(
#'     "extdata", "small_proteome_50_proteins.fasta",
#'     package = "pepVet"
#'   )
#'   batch <- batch_evaluate(small_fasta, enzyme = "trypsin")
#'   diag <- score_diagnostics(batch)
#'   p <- plot_score_diagnostics(diag)
#'   print(p)
#' }
#' @export
plot_score_diagnostics <- function(x, title = NULL) {
  .validate_diagnostic_title(title)
  .validate_diagnostic_result(x)

  rlang::check_installed("ggplot2",
    reason = "to produce score diagnostics plots"
  )
  rlang::check_installed("patchwork",
    reason = "to assemble score diagnostics panels"
  )

  n_prot <- x$n_proteins
  n_comp <- x$n_components
  figure_title <- title %||% paste0(
    "Score Diagnostics  |  ", n_prot, " proteins  |  ",
    n_comp, " components"
  )

  ## Panel A: VIF bar chart with severity coloring
  finite_vif <- x$vif[is.finite(x$vif)]
  vif_cap <- if (length(finite_vif) == 0L) {
    10
  } else {
    max(10, finite_vif)
  }
  plot_vif <- x$vif
  plot_vif[is.na(plot_vif)] <- 0
  plot_vif[is.infinite(plot_vif)] <- vif_cap * 1.1
  vif_label <- ifelse(
    is.na(x$vif), "N/A",
    ifelse(is.infinite(x$vif), "Inf", sprintf("%.1f", x$vif))
  )

  vif_df <- data.frame(
    component = factor(names(x$vif), levels = rev(names(x$vif))),
    vif       = x$vif,
    plot_vif  = plot_vif,
    vif_label = vif_label,
    severity  = ifelse(is.na(x$vif), "unknown",
                  ifelse(x$vif >= 10, "high",
                    ifelse(x$vif >= 5, "moderate", "low"))),
    stringsAsFactors = FALSE
  )
  vif_df$severity <- factor(vif_df$severity,
    levels = c("low", "moderate", "high", "unknown"))

  severity_colors <- c(
    low      = .pepvet_pal$good,
    moderate = .pepvet_pal$moderate,
    high     = .pepvet_pal$poor,
    unknown  = .pepvet_pal$na_gray
  )

  pa <- ggplot2::ggplot(vif_df,
      ggplot2::aes(y = .data$component, x = .data$plot_vif)) +
    ggplot2::geom_vline(xintercept = 5,
      linetype = "dashed", color = severity_colors[["moderate"]],
      linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 10,
      linetype = "dashed", color = severity_colors[["high"]],
      linewidth = 0.4) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$severity),
      alpha = 0.85, width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$vif_label),
      hjust = -0.1, size = 2.5) +
    ggplot2::scale_fill_manual(
      values = severity_colors,
      labels = c("low"      = "< 5  (acceptable)",
                 "moderate" = "5-10  (concerning)",
                 "high"     = "> 10  (problematic)",
                 "unknown"  = "N/A"),
      name = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(
      title = "A  Multicollinearity (VIF)",
      x = "Variance Inflation Factor", y = NULL) +
    .pepvet_theme() +
    ggplot2::theme(
      legend.position = .get_plot_theme_value("legend.position", "bottom"),
      legend.key.size = ggplot2::unit(8, "pt"),
      legend.text = ggplot2::element_text(size = 6),
      legend.margin = ggplot2::margin(t = -4, b = 0),
      legend.spacing = ggplot2::unit(0, "pt"))

  ## Panel B: PCA scree + cumulative with PC1+PC2 annotation
  cum2 <- sum(x$pca$var_explained[1:2])
  var_df <- data.frame(
    pc = factor(paste0("PC", seq_along(x$pca$var_explained)),
      levels = paste0("PC", seq_along(x$pca$var_explained))),
    var_exp = x$pca$var_explained,
    cum_var = cumsum(x$pca$var_explained),
    stringsAsFactors = FALSE
  )

  pb <- ggplot2::ggplot(var_df, ggplot2::aes(x = .data$pc)) +
    ggplot2::geom_col(ggplot2::aes(y = .data$var_exp,
      fill = .data$pc), alpha = 0.85, width = 0.6) +
    ggplot2::geom_point(ggplot2::aes(y = .data$cum_var),
      color = .pepvet_pal$brand, size = 1.5) +
    ggplot2::geom_line(ggplot2::aes(y = .data$cum_var, group = 1),
      color = .pepvet_pal$brand, linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 0.8,
      linetype = "dashed", color = "darkgrey", linewidth = 0.3) +
    ggplot2::annotate("text", x = 1, y = 0.8,
      label = "80%", hjust = -0.2, vjust = -0.5,
      size = 2.5, color = "darkgrey") +
    ggplot2::scale_fill_manual(values = rep(.pepvet_pal$brand_light,
      length(x$pca$var_explained)), guide = "none") +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(x * 100), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(
      title = "B  PCA variance explained",
      subtitle = paste0("PC1+PC2 = ", round(cum2 * 100), "%"),
      x = NULL, y = "Variance explained") +
    .pepvet_theme()

  ## Panel C: Ablation waterfall with error bars and flip counts
  abl_df <- x$ablation
  abl_df$component_f <- factor(abl_df$component,
    levels = rev(abl_df$component))

  pc <- ggplot2::ggplot(abl_df, ggplot2::aes(y = .data$component_f,
    x = .data$mean_drop)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$component),
      width = 0.6, alpha = 0.85) +
    ggplot2::geom_errorbar(ggplot2::aes(
      xmin = .data$mean_drop - .data$sd_drop,
      xmax = .data$mean_drop + .data$sd_drop),
      width = 0.1, linewidth = 0.3,
      color = .pepvet_pal$text_subtitle) +
    ggplot2::geom_text(ggplot2::aes(
      label = sprintf("drop %.3f\n(%d/%d flips)",
        .data$mean_drop, .data$n_verdict_flipped, n_prot)),
      hjust = -0.15, size = 2.5,
      color = .pepvet_pal$text_dark, lineheight = 0.9) +
    ggplot2::scale_fill_manual(values = .pepvet_pal$component,
      guide = "none") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.55))) +
    ggplot2::labs(
      title = "C  Ablation (component = 0)",
      x = "Mean composite drop", y = NULL) +
    .pepvet_theme()

  combined <- pa + pb + pc +
    patchwork::plot_layout(ncol = 3, widths = c(1.2, 1, 1.2))

  combined <- combined + patchwork::plot_annotation(
    title = figure_title,
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 10, face = "bold",
        color = .pepvet_pal$brand_dark,
        margin = ggplot2::margin(b = 4))
    )
  )

  combined
}
