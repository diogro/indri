#' Plot an association scan
#'
#' A Manhattan plot (`-log10(p)` along the genome) or a QQ plot (observed
#' against expected `-log10(p)` under the null), drawn with ggplot2 and
#' returned as a `ggplot` object to print, save or extend.
#'
#' @param x An `indri_gwas` scan.
#' @param type `"manhattan"` or `"qq"`.
#' @param which The p-values to draw: `"joint"` (`pvalues`, the default),
#'   `"dominance"` (`pvalues_dominance`), `"gxe"` (`pvalues_gxe`) or
#'   `"empirical"` (`pvalues_empirical`).
#' @param trait For single-trait scans, the trait to draw; `NULL` (the
#'   default) draws each SNP's smallest p-value over traits, adjusted when
#'   the scan was (`pvalues_adjust`). The QQ plot of a single-trait scan
#'   uses every SNP and trait when `trait` is `NULL`.
#' @param threshold Significance line on the Manhattan plot, as a p-value;
#'   defaults to Bonferroni over the SNPs drawn (`0.05 / M`). `NA` for none.
#' @param highlight Names of SNPs to mark and label.
#' @param ... Unused.
#' @return A `ggplot` object.
#' @details The Manhattan plot places SNPs by `chrom` and `pos` from
#'   `snp_info` when both are there, chromosomes end to end; otherwise in
#'   column order. SNPs with `NA` p-values (monomorphic, or without a
#'   dominance estimate) are left out.
#' @seealso [indri_gwas()]
#' @export
plot.indri_gwas <- function(x, type = c("manhattan", "qq"),
                            which = c("joint", "dominance", "gxe",
                                      "empirical"),
                            trait = NULL, threshold = NULL, highlight = NULL,
                            ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Plotting an indri_gwas scan needs ggplot2: ",
         "install.packages(\"ggplot2\").", call. = FALSE)
  type <- match.arg(type)
  which <- match.arg(which)
  p <- .indri_gwas_plot_p(x, which, trait, all_traits = type == "qq")
  label <- .indri_gwas_plot_label(x, which, trait)
  if (type == "qq") return(.indri_gwas_qq(p, label))

  snps <- .indri_gwas_snp_names(x)
  data <- .indri_gwas_positions(x$snp_info, snps)
  data$logp <- -log10(unname(p))
  data$highlight <- FALSE
  if (!is.null(highlight)) {
    unknown <- setdiff(highlight, snps)
    if (length(unknown))
      stop("`highlight` names SNPs not in the scan: ",
           paste(utils::head(unknown, 5L), collapse = ", "), call. = FALSE)
    data$highlight <- data$snp %in% highlight
  }
  data <- data[!is.na(data$logp), , drop = FALSE]
  if (is.null(threshold)) threshold <- 0.05 / nrow(data)

  aes <- ggplot2::aes
  plot <- ggplot2::ggplot(data, aes(x = .data$x, y = .data$logp)) +
    ggplot2::geom_point(aes(colour = .data$band), size = 0.9,
                        show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(odd = "grey25", even = "grey60")) +
    ggplot2::labs(x = NULL, y = expression(-log[10](p)), title = label) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
  if (!all(is.na(data$chrom))) {
    ticks <- tapply(data$x, data$chrom, mean)
    plot <- plot + ggplot2::scale_x_continuous(breaks = unname(ticks),
                                               labels = names(ticks))
  }
  if (!is.na(threshold))
    plot <- plot + ggplot2::geom_hline(
      data = data.frame(yintercept = -log10(threshold)),
      ggplot2::aes(yintercept = .data$yintercept),
      linetype = "dashed", colour = "firebrick")
  if (any(data$highlight)) {
    marked <- data[data$highlight, , drop = FALSE]
    plot <- plot +
      ggplot2::geom_point(data = marked, colour = "firebrick", size = 2) +
      ggplot2::geom_text(data = marked, ggplot2::aes(label = .data$snp),
                         vjust = -0.8, size = 3)
  }
  plot
}

#' Summarise an association scan
#'
#' The SNPs with the smallest p-values and the genomic inflation factor
#' `lambda_gc` (the median test statistic over its expectation under the
#' null; near 1 for a calibrated scan of mostly null SNPs).
#'
#' @param object An `indri_gwas` scan.
#' @param n Number of top SNPs (or SNP-trait pairs) to list.
#' @param ... Unused.
#' @return An object of class `summary.indri_gwas` with `top` (a data frame),
#'   `lambda_gc` (one value, or one per trait for single-trait scans),
#'   `n_tested` and the scan's description.
#' @export
summary.indri_gwas <- function(object, n = 10L, ...) {
  df <- as.data.frame(object)
  df <- df[!is.na(df$pvalue), , drop = FALSE]
  # keep the table readable: per-trait effect columns stay in the scan
  df <- df[, !grepl("^(a|d)_", names(df)), drop = FALSE]
  top <- utils::head(df[order(df$pvalue), , drop = FALSE], n)
  rownames(top) <- NULL
  df_test <- if (object$mode == "multi_trait") length(object$traits) else 1L
  df_test <- df_test * (1L + isTRUE(object$dominance) +
                          length(object$gxe_terms))
  lambda <- function(p) {
    p <- p[!is.na(p)]
    stats::qchisq(stats::median(p), df_test, lower.tail = FALSE) /
      stats::qchisq(0.5, df_test, lower.tail = FALSE)
  }
  lambda_gc <- if (object$mode == "multi_trait") lambda(object$pvalues)
               else apply(object$pvalues, 2L, lambda)
  structure(list(top = top, lambda_gc = lambda_gc,
                 n_tested = sum(!is.na(as.matrix(object$pvalues)[, 1L])),
                 mode = object$mode, model = object$model,
                 dominance = isTRUE(object$dominance), delta = object$delta,
                 fdr_threshold = object$fdr_threshold, fdr = object$fdr),
            class = "summary.indri_gwas")
}

#' @export
print.summary.indri_gwas <- function(x, ...) {
  what <- if (x$mode == "multi_trait") "multi-trait" else "single-trait"
  cat("indri association scan:", what, toupper(x$model),
      if (x$dominance) "with dominance", "\n")
  cat(sprintf("  %d SNPs tested\n", x$n_tested))
  if (length(x$lambda_gc) == 1L) {
    cat(sprintf("  genomic inflation (lambda_gc): %.3f\n", x$lambda_gc))
  } else {
    cat(sprintf("  genomic inflation (lambda_gc) over traits: %.3f to %.3f\n",
                min(x$lambda_gc), max(x$lambda_gc)))
  }
  if (!is.null(x$fdr)) cat(.indri_gwas_fdr_line(x$fdr, x$fdr_threshold))
  cat("\nsmallest p-values:\n")
  print(x$top, row.names = FALSE, digits = 4)
  invisible(x)
}

# ---------------------------------------------------------------------------

.indri_gwas_snp_names <- function(x) {
  if (!is.null(x$snp_info$snp)) as.character(x$snp_info$snp) else rownames(x$a)
}

.indri_gwas_plot_p <- function(x, which, trait, all_traits) {
  p <- switch(which,
    joint = x$pvalues,
    dominance = {
      if (!isTRUE(x$dominance))
        stop("This scan has no dominance tests; rerun indri_gwas() with ",
             "`dominance = TRUE`.", call. = FALSE)
      x$pvalues_dominance
    },
    gxe = {
      if (is.null(x$gxe_terms))
        stop("This scan has no genotype-by-environment tests; rerun ",
             "indri_gwas() with `gxe`.", call. = FALSE)
      x$pvalues_gxe
    },
    empirical = {
      if (is.null(x$pvalues_empirical))
        stop("This scan has no empirical p-values; rerun indri_gwas() with ",
             "`n_permutations`.", call. = FALSE)
      x$pvalues_empirical
    })
  if (x$mode == "multi_trait") return(p)
  if (!is.null(trait)) {
    if (!trait %in% colnames(p))
      stop("`trait` must be one of the scan's traits: ",
           paste(utils::head(colnames(p), 5L), collapse = ", "), call. = FALSE)
    return(p[, trait])
  }
  if (all_traits) return(c(p))
  if (which == "joint" && !is.null(x$pvalues_adjust)) p <- x$pvalues_adjust
  apply(p, 1L, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
}

.indri_gwas_plot_label <- function(x, which, trait) {
  test <- switch(which, joint = if (isTRUE(x$dominance) || !is.null(x$gxe_terms))
                   "joint test" else "association",
                 dominance = "dominance given additive",
                 gxe = paste0("interaction with ",
                              paste(x$gxe_terms, collapse = ", ")),
                 empirical = "empirical p-values")
  scope <- if (x$mode == "multi_trait") {
    sprintf("multi-trait %s, %d traits", toupper(x$model), length(x$traits))
  } else if (!is.null(trait)) {
    sprintf("single-trait %s, %s", toupper(x$model), trait)
  } else {
    sprintf("single-trait %s, best trait per SNP%s", toupper(x$model),
            if (which == "joint" && !is.null(x$pvalues_adjust))
              " (adjusted)" else "")
  }
  paste0(scope, ": ", test)
}

# x positions: chromosomes end to end when snp_info has chrom and pos,
# column order otherwise. `band` alternates by chromosome for the colours.
.indri_gwas_positions <- function(info, snps) {
  M <- length(snps)
  if (!is.null(info) && all(c("chrom", "pos") %in% names(info))) {
    chrom <- info$chrom
    levels <- unique(chrom[order(suppressWarnings(as.numeric(as.character(chrom))),
                                 as.character(chrom))])
    pos <- as.numeric(info$pos)
    span <- tapply(pos, factor(chrom, levels = levels), max)
    gap <- 0.02 * sum(span)
    offset <- c(0, cumsum(span + gap))[seq_along(levels)]
    names(offset) <- levels
    x <- pos + offset[as.character(chrom)]
    band <- ifelse(match(as.character(chrom), levels) %% 2L == 1L, "odd", "even")
    return(data.frame(snp = snps, chrom = chrom, x = unname(x), band = band,
                      stringsAsFactors = FALSE))
  }
  data.frame(snp = snps, chrom = NA, x = seq_len(M), band = "odd",
             stringsAsFactors = FALSE)
}

.indri_gwas_qq <- function(p, label) {
  p <- sort(unname(p[!is.na(p)]))
  data <- data.frame(expected = -log10(stats::ppoints(length(p))),
                     observed = -log10(p))
  ggplot2::ggplot(data, ggplot2::aes(x = .data$expected, y = .data$observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::labs(x = expression(expected ~ -log[10](p)),
                  y = expression(observed ~ -log[10](p)),
                  title = paste("QQ,", label)) +
    ggplot2::theme_minimal()
}
