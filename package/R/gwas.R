#' Genome-wide association scan with the assembled covariances
#'
#' Tests each SNP for association with the traits in `Y`, either jointly
#' (one multi-trait "any effect" test per SNP) or trait by trait, in a linear
#' mixed model with relationship matrix `A` or in a linear model. The
#' multi-trait mixed model uses the genetic and residual covariances `G` and
#' `E`, usually from [fit_G()], which is what makes it feasible for many
#' traits.
#'
#' @section Models:
#' For `N` individuals, `P` traits and one SNP `g`,
#' `Y = X B + g beta' + U + Eps`. Every test profiles an overall scale
#' `sigma^2` and is a likelihood-ratio test:
#'
#' * multi-trait LMM: `vec(U + Eps) ~ N(0, sigma^2 (G (x) A + delta E (x) I))`,
#'   `P` degrees of freedom. `G` and `E` set the shape; the ratio `delta` is
#'   fitted once by REML under the null (no SNP) and held fixed for every SNP.
#'   With `G` and `E` fitted to the same data it should be near 1, and is
#'   returned so that it can be checked.
#' * multi-trait LM: `sigma^2 E (x) I`, `P` degrees of freedom, with `E` by
#'   default the covariance of the residuals of `Y` on `X`.
#' * single-trait LMM: per trait, `sigma^2 (A + delta I)`, one degree of
#'   freedom, with `delta` fitted by REML under the null
#'   (`single_trait_variance = "reml"`) or taken as `E[k, k] / G[k, k]`
#'   (`"fit"`).
#' * single-trait LM: ordinary least squares.
#'
#' @section Dominance:
#' With `dominance = TRUE` each SNP enters as its additive dosage (0/1/2) and a
#' heterozygote indicator (1 for genotype 1, 0 for both homozygotes).
#' `pvalues` is then the joint test of both terms (2 degrees of freedom per
#' trait), `pvalues_dominance` the test of the dominance term with the
#' additive one in the null (1 per trait), and `a` and `d` the additive and
#' dominance effects of the joint model. A SNP observed in only two genotype
#' classes has a dominance term collinear with its dosage: it keeps the
#' additive-only test in `pvalues`, gets `NA` dominance results, and is
#' `FALSE` in `dominance_estimable`. Missing genotypes are imputed by the
#' column mean in each coding. The variance ratio is fitted without any SNP,
#' as without dominance.
#'
#' @section Genotype-by-environment:
#' With `gxe`, each SNP's dosage `g` also enters multiplied by each named
#' column `e` of `X` (`g * e`), next to the environments' main effects, which
#' are in `X` already. `pvalues` is then the joint test of the SNP's main
#' effect and its interactions (and dominance), with
#' `P * (1 + n_gxe [+ 1])` degrees of freedom for `P` traits, and
#' `pvalues_gxe` tests the interactions with the main effect (and dominance)
#' in the null. `a` is the main effect (the SNP's effect where every named
#' environment is 0) and `gxe` holds the interaction effects (SNPs by traits
#' by environments). With `dominance = TRUE` the dominance term is a main
#' effect only and `pvalues_dominance` tests it given the rest. A SNP whose
#' interactions cannot be estimated (for instance every carrier in one
#' environment) keeps the test without them, gets `NA` interaction results,
#' and is `FALSE` in `gxe_estimable`. The covariance is the same as without
#' `gxe`: GxE in the background (a polygenic G x E, or a residual variance
#' that differs between environments) is not modelled.
#'
#' @section Environment-specific background:
#' With `fit` from [fit_G_env()], the polygenic background has its own
#' genetic covariance in each group of a discrete environment and between
#' groups, and its own residual covariance in each: the covariance of the
#' traits of individuals `i` and `j` is `A_ij G[(env_i), (env_j)]` plus, for
#' `i = j`, `delta E[(env_i)]`. This is the background to use with `gxe` on
#' that environment: unmodelled polygenic GxE inflates the interaction
#' tests. The covariance is no Kronecker product, so the scan rotates by a
#' generalized eigendecomposition of its genetic and residual parts, once
#' per connected component of `A`: fast for families (a block-diagonal `A`),
#' cubic in the number of individuals times traits for a dense `A`. `delta`
#' is fitted by REML as usual; `single_trait_variance` does not apply.
#'
#' @section Permutations:
#' With `n_permutations > 0` the genotype rows are permuted across
#' individuals and the scan rerun, keeping the null variance ratio of the
#' observed scan. `pvalues_empirical` is then
#' `(1 + #{permuted p <= observed p}) / (n_permutations + 1)`, per SNP (and
#' trait); with `dominance = TRUE` it refers to the joint test.
#'
#' The permuted scans also estimate the false discovery rate, per SNP (for
#' single-trait scans, from each SNP's smallest adjusted p-value, or raw
#' without `adjust`): at a threshold `t`,
#' `FDR(t) = (mean # permuted p <= t per permutation) / (# observed p <= t)`.
#' `qvalues_empirical` is, for each SNP, the smallest `FDR(t)` over
#' thresholds at or above its own p-value, and with `fdr`, `fdr_threshold` is
#' the largest observed p-value with `FDR(t) <= fdr` (`NA` if none): the SNPs
#' at or below it are the discoveries at that rate. (LiMMBo's "empirical
#' FDR" was the `fdr`-quantile of the permuted p-values, which is close to
#' `fdr` itself and so an uncorrected threshold.)
#' Permutations are drawn in the main process (`set.seed(seed)`, then
#' `n_permutations` calls to `sample.int(N)`), so results do not depend on
#' `cores`.
#'
#' **Relatedness.** Permuting genotypes across individuals treats them as
#' exchangeable, which related individuals are not: the permuted null breaks
#' the link between genotype and relatedness that the mixed model accounts
#' for. The empirical p-values follow LiMMBo and are a check on the analytic
#' ones, not a replacement for them in structured samples.
#'
#' The covariance is diagonalised once (eigenvectors of `A`, then of `G`
#' relative to `E`), so the scan costs one `N x N` rotation of the genotypes
#' plus `O(N P)` per SNP.
#'
#' @param Y Numeric matrix, individuals by traits, with no missing values.
#' @param genotypes Numeric matrix, individuals by SNPs (typically 0/1/2).
#'   Missing values are replaced by the SNP's mean; monomorphic SNPs get `NA`.
#' @param fit An `indri_fit` from [fit_G()] supplying `G` and `E`, or an
#'   `indri_fit_env` from [fit_G_env()] for an environment-specific
#'   background; see that section.
#' @param G,E Genetic and residual covariance matrices of the traits, in the
#'   units of `Y`, instead of `fit`. Matched to the columns of `Y` by name
#'   when both are named. For the multi-trait LM, `E` replaces the residual
#'   covariance.
#' @param A Relationship matrix of the individuals; required for
#'   `model = "lmm"`.
#' @param X Fixed-effects design with full column rank that can represent a
#'   constant, as for [fit_G()]; defaults to an intercept. Principal
#'   components go here.
#' @param snp_info Optional data frame with one row per column of
#'   `genotypes` (for example `snp`, `chrom`, `pos`), carried into the
#'   results.
#' @param mode `"multi_trait"` for one joint test per SNP, `"single_trait"`
#'   for one test per SNP and trait.
#' @param model `"lmm"` (with `A`) or `"lm"`.
#' @param single_trait_variance For the single-trait LMM, `"reml"` fits each
#'   trait's variance ratio, `"fit"` takes it from `G` and `E`.
#' @param delta Optional fixed variance ratio for the LMMs, instead of the
#'   fitted one: a number for the multi-trait test, one value per trait for
#'   the single-trait test.
#' @param dominance If `TRUE`, also fit a dominance (heterozygote) effect for
#'   each SNP; see the Dominance section. Needs genotypes coded 0, 1 and 2.
#' @param gxe Names of columns of `X` (environments, such as sex, diet or a
#'   covariate) whose interaction with each SNP's additive effect is also
#'   fitted; see the Genotype-by-environment section. A factor enters as the
#'   indicator columns `X` has for it.
#' @param adjust Single-trait p-value adjustment for the number of traits:
#'   `"bonferroni"` multiplies by `P`, `"effective"` by
#'   [indri_effective_tests()] of the residuals after `X`; both cap at 1.
#' @param n_permutations Number of genotype permutations for empirical
#'   p-values (`0`, the default, for none); see the Permutations section.
#' @param fdr Optional false discovery rate for a p-value threshold estimated
#'   from the permutations; without `n_permutations` it runs 10
#'   permutations (LiMMBo's number).
#' @param seed Optional seed for the permutations; the caller's random stream
#'   is restored afterwards.
#' @param cores Number of processes for the permuted scans (forked; one on
#'   Windows).
#' @param chunk_size Number of SNPs rotated at a time.
#'
#' @return An object of class `indri_gwas`: `pvalues` (a vector over SNPs,
#'   or SNPs by traits), `a` (the additive effects, SNPs by traits, in the
#'   units of `Y` per allele; not to be confused with the relationship
#'   matrix `A`), `delta`, `pvalues_adjust` and `n_tests` when
#'   adjusted, `snp_info`, `mode`, `model` and `call`. With
#'   `dominance = TRUE` also `pvalues_dominance`, `d` (the dominance effects: the
#'   heterozygote's deviation from the homozygotes' midpoint) and
#'   `dominance_estimable`; `pvalues_adjust` then adjusts the joint test.
#'   With `gxe`, also `pvalues_gxe`, `gxe` (SNPs by traits by environments),
#'   `gxe_estimable` and `gxe_terms`; `pvalues_adjust` then adjusts the joint
#'   test.
#'   With permutations, `pvalues_empirical`, `qvalues_empirical` and
#'   `n_permutations`, and with `fdr`, `fdr` and `fdr_threshold`.
#' @seealso [fit_G()], [indri_effective_tests()]
#' @references
#' Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
#' a simple, scalable approach for linear mixed models in high-dimensional
#' genetic association studies. bioRxiv, preprint.
#' \doi{10.1101/255497}
#' @export
indri_gwas <- function(Y, genotypes, fit = NULL, G = NULL, E = NULL, A = NULL,
                       X = NULL, snp_info = NULL,
                       mode = c("multi_trait", "single_trait"),
                       model = c("lmm", "lm"),
                       single_trait_variance = c("reml", "fit"),
                       delta = NULL, dominance = FALSE, gxe = NULL,
                       adjust = c("none", "bonferroni", "effective"),
                       n_permutations = 0L, fdr = NULL, seed = NULL,
                       cores = 1L, chunk_size = 1000L) {
  call <- match.call()
  mode <- match.arg(mode)
  model <- match.arg(model)
  single_trait_variance <- match.arg(single_trait_variance)
  adjust <- match.arg(adjust)
  chunk_size <- .indri_count(chunk_size, "chunk_size")
  if (!is.numeric(n_permutations) || length(n_permutations) != 1L ||
      is.na(n_permutations) || n_permutations < 0 ||
      n_permutations != round(n_permutations))
    stop("`n_permutations` must be a whole number >= 0.", call. = FALSE)
  n_permutations <- as.integer(n_permutations)
  if (!is.null(fdr)) {
    if (!is.numeric(fdr) || length(fdr) != 1L || is.na(fdr) || fdr <= 0 ||
        fdr >= 1)
      stop("`fdr` must be a number between 0 and 1.", call. = FALSE)
    if (n_permutations == 0L) n_permutations <- 10L   # LiMMBo's default
  }
  cores <- .indri_count(cores, "cores")

  Y <- .indri_gwas_matrix(Y, "Y")
  if (anyNA(Y))
    stop("`Y` has missing values; indri_gwas() needs complete phenotypes ",
         "(drop or impute the incomplete individuals first).", call. = FALSE)
  N <- nrow(Y)
  P <- ncol(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("trait", seq_len(P))
  traits <- colnames(Y)

  genotypes <- .indri_gwas_matrix(genotypes, "genotypes")
  if (is.null(colnames(genotypes)))
    colnames(genotypes) <- paste0("snp", seq_len(ncol(genotypes)))
  genotypes <- .indri_gwas_rows(genotypes, rownames(Y), N, "genotypes")
  M <- ncol(genotypes)

  if (is.null(X)) {
    X <- matrix(1, N, 1L, dimnames = list(rownames(Y), "(Intercept)"))
  } else {
    X <- .indri_gwas_rows(.indri_gwas_matrix(X, "X"), rownames(Y), N, "X")
  }
  .indri_check_design(X)

  if (model == "lmm") {
    if (is.null(A))
      stop("`A` is required for model = \"lmm\".", call. = FALSE)
    A <- .indri_gwas_matrix(A, "A")
    if (is.null(rownames(A)) && !is.null(colnames(A))) rownames(A) <- colnames(A)
    if (is.null(colnames(A))) colnames(A) <- rownames(A)
    A <- .indri_gwas_rows(A, rownames(Y), N, "A")
    A <- t(.indri_gwas_rows(t(A), rownames(Y), N, "A"))
  }

  if (!is.null(snp_info)) {
    snp_info <- as.data.frame(snp_info)
    if (nrow(snp_info) != M)
      stop("`snp_info` must have one row per column of `genotypes` (", M,
           "), not ", nrow(snp_info), ".", call. = FALSE)
  } else {
    snp_info <- data.frame(snp = colnames(genotypes), stringsAsFactors = FALSE)
  }

  env_fit <- NULL
  if (inherits(fit, "indri_fit_env")) {
    if (!is.null(G) || !is.null(E))
      stop("Pass either `fit` or `G` and `E`, not both.", call. = FALSE)
    if (model != "lmm")
      stop("An environment-specific fit (fit_G_env()) needs ",
           "model = \"lmm\".", call. = FALSE)
    if (!setequal(fit$traits, traits))
      stop("`fit` describes the traits ",
           paste(utils::head(fit$traits, 5L), collapse = ", "),
           ", not the traits of `Y`.", call. = FALSE)
    fe <- fit$env
    if (!is.null(names(fe)) && !is.null(rownames(Y))) {
      missing <- setdiff(rownames(Y), names(fe))
      if (length(missing))
        stop("`fit` has no environment for ", length(missing),
             " individual(s) of `Y`: ",
             paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
      fe <- fe[rownames(Y)]
    } else if (length(fe) != N) {
      stop("`fit` has an environment for ", length(fe), " individuals, `Y` ",
           "has ", N, ".", call. = FALSE)
    }
    env_fit <- list(Gs = fit$G, Es = fit$E, env = as.character(fe))
    fit <- NULL
  }
  if (!is.null(fit)) {
    if (!is.null(G) || !is.null(E))
      stop("Pass either `fit` or `G` and `E`, not both.", call. = FALSE)
    G <- fit$G
    E <- fit$E
    if (is.null(E) && (model == "lmm" && (mode == "multi_trait" ||
                                          single_trait_variance == "fit")))
      stop("`fit` has no `E`: it was made before indri assembled the ",
           "residual covariance. Refit it with fit_G().", call. = FALSE)
  }
  needs_GE <- model == "lmm" && is.null(env_fit) &&
    (mode == "multi_trait" || single_trait_variance == "fit")
  if (needs_GE && (is.null(G) || is.null(E)))
    stop(if (mode == "multi_trait") "The multi-trait LMM" else
           "`single_trait_variance = \"fit\"`",
         " needs `G` and `E` (or an indri `fit`).", call. = FALSE)
  if (!is.null(G)) G <- .indri_gwas_trait_matrix(G, traits, "G")
  if (!is.null(E)) E <- .indri_gwas_trait_matrix(E, traits, "E")

  if (adjust != "none" && mode == "multi_trait")
    stop("`adjust` corrects single-trait p-values for the number of traits; ",
         "use it with mode = \"single_trait\".", call. = FALSE)
  if (!is.null(delta)) {
    if (model != "lmm")
      stop("`delta` applies to model = \"lmm\" only.", call. = FALSE)
    delta <- .indri_gwas_delta(delta, mode, traits)
  }

  if (!is.logical(dominance) || length(dominance) != 1L || is.na(dominance))
    stop("`dominance` must be TRUE or FALSE.", call. = FALSE)
  dom <- NULL
  estimable <- NULL
  if (dominance) {
    seen <- genotypes[!is.na(genotypes)]
    if (any(abs(seen - round(seen)) > 1e-8 | seen < 0 | seen > 2))
      stop("`dominance = TRUE` needs genotypes coded 0, 1 and 2 (the ",
           "heterozygote is genotype 1); these are dosages.", call. = FALSE)
    estimable <- apply(genotypes, 2L,
                       function(g) length(unique(g[!is.na(g)])) == 3L)
    names(estimable) <- colnames(genotypes)
    dom <- .indri_gwas_impute((genotypes == 1) * 1)
  }

  env <- NULL
  if (!is.null(gxe)) {
    if (!is.character(gxe) || !length(gxe) || anyNA(gxe))
      stop("`gxe` must be the names of columns of `X`.", call. = FALSE)
    if (is.null(colnames(X)))
      stop("`gxe` names columns of `X`, but `X` has no column names.",
           call. = FALSE)
    unknown <- setdiff(gxe, colnames(X))
    if (length(unknown))
      stop("`gxe` names columns that `X` does not have: ",
           paste(unknown, collapse = ", "), ". The environments must be in ",
           "`X`, so that their main effects are in the model.", call. = FALSE)
    env <- X[, gxe, drop = FALSE]
    constant <- apply(env, 2L, function(e) max(e) - min(e) < 1e-12)
    if (any(constant))
      stop("`gxe` names a constant column of `X` (",
           paste(gxe[constant], collapse = ", "), "): an interaction with it ",
           "is the SNP's main effect.", call. = FALSE)
  }

  genotypes <- .indri_gwas_impute(genotypes)
  monomorphic <- apply(genotypes, 2L, function(g) max(g) - min(g) < 1e-12)
  if (any(monomorphic))
    warning(sum(monomorphic), " monomorphic SNP(s) get NA results: ",
            paste(utils::head(colnames(genotypes)[monomorphic], 5L),
                  collapse = ", "),
            if (sum(monomorphic) > 5L) ", ...", call. = FALSE)

  run_scan <- function(genotypes, dom, delta) {
    if (!is.null(env_fit)) {
      .indri_gwas_env(Y, X, genotypes, env_fit, A, mode, delta, chunk_size,
                      monomorphic, dom, estimable, env)
    } else if (mode == "multi_trait") {
      .indri_gwas_multi(Y, X, genotypes, G, E, A, model, delta, chunk_size,
                        monomorphic, dom, estimable, env)
    } else {
      .indri_gwas_single(Y, X, genotypes, G, E, A, model, delta,
                         single_trait_variance, chunk_size, monomorphic, dom,
                         estimable, env)
    }
  }
  scan <- run_scan(genotypes, dom, delta)

  out <- list(pvalues = scan$pvalues, a = scan$betas, delta = scan$delta,
              pvalues_adjust = NULL, n_tests = NULL, snp_info = snp_info,
              dominance = dominance,
              pvalues_dominance = scan$pvalues_dominance,
              d = scan$betas_dominance,
              dominance_estimable = if (is.null(scan$dominance_estimable))
                estimable else scan$dominance_estimable,
              gxe_terms = gxe,
              pvalues_gxe = scan$pvalues_gxe,
              gxe = scan$betas_gxe,
              gxe_estimable = scan$gxe_estimable,
              background = if (is.null(env_fit)) "kronecker" else "environment",
              n_components = scan$n_components,
              mode = mode, model = model, n_individuals = N, traits = traits,
              call = call)
  if (adjust != "none") {
    out$n_tests <- if (adjust == "bonferroni") P else indri_effective_tests(Y, X)
    out$pvalues_adjust <- pmin(scan$pvalues * out$n_tests, 1)
  }
  if (n_permutations > 0L) {
    # What each permuted scan contributes to the FDR pool: the p-value per
    # SNP, or for single-trait scans each SNP's smallest (adjusted) one.
    per_snp <- function(p) {
      if (mode == "multi_trait") return(p)
      if (!is.null(out$n_tests)) p <- pmin(p * out$n_tests, 1)
      apply(p, 1L, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
    }
    restore <- .indri_use_seed(seed)
    perms <- vapply(seq_len(n_permutations), function(i) sample.int(N),
                    integer(N))
    restore()
    one <- function(i) {
      o <- perms[, i]
      run_scan(genotypes[o, , drop = FALSE],
               if (is.null(dom)) NULL else dom[o, , drop = FALSE],
               scan$delta)$pvalues
    }
    null_p <- if (cores > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(seq_len(n_permutations), one, mc.cores = cores)
    } else {
      lapply(seq_len(n_permutations), one)
    }
    failed <- vapply(null_p, inherits, logical(1), "try-error")
    if (any(failed))
      stop("A permuted scan failed: ", as.character(null_p[[which(failed)[1]]]),
           call. = FALSE)
    count <- Reduce(`+`, lapply(null_p, function(p) p <= scan$pvalues))
    out$pvalues_empirical <- (count + 1) / (n_permutations + 1)
    out$n_permutations <- n_permutations
    est <- .indri_gwas_permutation_fdr(per_snp(scan$pvalues),
                                       lapply(null_p, per_snp), fdr)
    out$qvalues_empirical <- stats::setNames(est$q, colnames(genotypes))
    if (!is.null(fdr)) {
      out$fdr <- fdr
      out$fdr_threshold <- est$threshold
    }
  }
  structure(out, class = "indri_gwas")
}

#' Effective number of independent traits
#'
#' Galwey's (2009) effective number of tests, from the eigenvalues
#' `lambda` of the Spearman correlation matrix of the traits:
#' `(sum(sqrt(lambda)))^2 / sum(lambda)`. It is 1 for perfectly correlated
#' traits and `ncol(Y)` for independent ones, and is the factor
#' `adjust = "effective"` multiplies single-trait p-values by in
#' [indri_gwas()].
#'
#' With a design `X`, the correlation is that of the residuals of `Y` on `X`:
#' the single-trait tests see the traits after the fixed effects, and
#' covariates shared by the traits (sex, age, line) would otherwise make
#' them look more correlated than they are to the tests.
#'
#' @param Y Numeric matrix, individuals by traits.
#' @param X Optional fixed-effects design, as in [indri_gwas()].
#' @return A number between 1 and `ncol(Y)`.
#' @references Galwey, N. W. (2009). A new measure of the effective number of
#'   tests, a practical tool for comparing families of non-independent
#'   significance tests. Genetic Epidemiology 33, 559-568.
#'   \doi{10.1002/gepi.20408}
#' @export
indri_effective_tests <- function(Y, X = NULL) {
  Y <- as.matrix(Y)
  if (!is.null(X)) Y <- stats::lm.fit(as.matrix(X), Y)$residuals
  R <- stats::cor(Y, method = "spearman", use = "pairwise.complete.obs")
  lambda <- pmax(eigen(R, symmetric = TRUE, only.values = TRUE)$values, 0)
  sum(sqrt(lambda))^2 / sum(lambda)
}

# ---------------------------------------------------------------------------
# The scans. Both reduce to independent weighted regressions, one per column
# of a rotated response, with per-entry variances V[i, c] (see
# .indri_gwas_columns()).

.indri_gwas_multi <- function(Y, X, genotypes, G, E, A, model, delta,
                              chunk_size, monomorphic, dom = NULL,
                              estimable = NULL, env = NULL) {
  N <- nrow(Y); P <- ncol(Y)
  if (model == "lm" && is.null(E)) {
    resid <- stats::lm.fit(X, Y)$residuals
    E <- crossprod(resid) / (N - ncol(X))
  }
  L <- .indri_gwas_chol(E, "E")
  if (model == "lmm") {
    # E^-1/2 G E^-T/2 = W diag(d) W'; T = L^-T W whitens E and diagonalises G.
    Gw <- backsolve(L, t(backsolve(L, G, upper.tri = FALSE)), upper.tri = FALSE)
    eg <- eigen((Gw + t(Gw)) / 2, symmetric = TRUE)
    d <- pmax(eg$values, 0)
    Tm <- backsolve(t(L), eg$vectors)
    rot <- .indri_gwas_rotation(A)
  } else {
    d <- rep(0, P)
    Tm <- backsolve(t(L), diag(P))
    rot <- list(Q = NULL, lambda = rep(0, N))
  }
  Ys <- .indri_gwas_rotate(rot$Q, Y %*% Tm)
  Xs <- .indri_gwas_rotate(rot$Q, X)
  V_of <- function(dl) outer(rot$lambda, d) + if (model == "lmm") dl else 1
  if (model == "lmm" && is.null(delta))
    delta <- .indri_gwas_reml_delta(Ys, Xs, V_of)
  V <- V_of(if (model == "lmm") delta else 1)

  M <- ncol(genotypes)
  null <- .indri_gwas_null(Ys, Xs, V)
  if (!is.null(env)) {
    out <- .indri_gwas_gxe(null, Xs, V, rot, genotypes, dom, env, estimable,
                           monomorphic, chunk_size, multi = TRUE, N = N,
                           Ti = solve(Tm), snps = colnames(genotypes),
                           traits = colnames(Y))
    out$delta <- if (model == "lmm") unname(delta) else NULL
    return(out)
  }
  stat <- stat_dom <- rep(NA_real_, M)
  df <- rep(P, M)
  beta_tilde <- beta_tilde_dom <- matrix(NA_real_, M, P)
  for (idx in .indri_gwas_chunks(M, chunk_size)) {
    Gs <- .indri_gwas_rotate(rot$Q, genotypes[, idx, drop = FALSE])
    Ds <- if (!is.null(dom)) .indri_gwas_rotate(rot$Q, dom[, idx, drop = FALSE])
    alt <- .indri_gwas_columns(null, Xs, Gs, V, Ds, estimable[idx])
    q_a <- colSums(alt$q1)
    if (is.null(dom)) {
      stat[idx] <- N * P * log(sum(null$q0) / q_a)
      beta_tilde[idx, ] <- t(alt$beta)
    } else {
      # the joint model on every trait, or the additive one where the
      # dominance term is not estimable in some column
      ok <- colSums(!alt$ok) == 0L
      q_ad <- ifelse(ok, colSums(alt$q2), q_a)
      stat[idx] <- N * P * log(sum(null$q0) / q_ad)
      stat_dom[idx] <- ifelse(ok, N * P * log(q_a / q_ad), NA_real_)
      df[idx] <- ifelse(ok, 2L * P, P)
      b_a <- t(alt$beta_joint)
      b_a[!ok, ] <- t(alt$beta)[!ok, ]
      b_d <- t(alt$beta_dom)
      b_d[!ok, ] <- NA_real_
      beta_tilde[idx, ] <- b_a
      beta_tilde_dom[idx, ] <- b_d
    }
  }
  stat[monomorphic] <- stat_dom[monomorphic] <- NA_real_
  beta_tilde[monomorphic, ] <- beta_tilde_dom[monomorphic, ] <- NA_real_
  # beta' T = beta_tilde', so beta' = beta_tilde' T^-1
  Ti <- solve(Tm)
  named <- function(b) {
    dimnames(b) <- list(colnames(genotypes), colnames(Y))
    b
  }
  pvalues <- stats::setNames(stats::pchisq(stat, df, lower.tail = FALSE),
                             colnames(genotypes))
  out <- list(pvalues = pvalues, betas = named(beta_tilde %*% Ti),
              delta = if (model == "lmm") unname(delta) else NULL)
  if (!is.null(dom)) {
    out$pvalues_dominance <- stats::setNames(
      stats::pchisq(stat_dom, P, lower.tail = FALSE), colnames(genotypes))
    out$betas_dominance <- named(beta_tilde_dom %*% Ti)
  }
  out
}

.indri_gwas_single <- function(Y, X, genotypes, G, E, A, model, delta,
                               single_trait_variance, chunk_size,
                               monomorphic, dom = NULL, estimable = NULL,
                               env = NULL) {
  N <- nrow(Y); P <- ncol(Y)
  rot <- if (model == "lmm") .indri_gwas_rotation(A)
         else list(Q = NULL, lambda = rep(0, N))
  Ys <- .indri_gwas_rotate(rot$Q, Y)
  Xs <- .indri_gwas_rotate(rot$Q, X)
  if (model == "lmm" && is.null(delta)) {
    delta <- if (single_trait_variance == "fit") {
      g <- diag(G)
      if (any(g <= 0))
        stop("`single_trait_variance = \"fit\"` needs a positive genetic ",
             "variance for every trait.", call. = FALSE)
      diag(E) / g
    } else {
      vapply(seq_len(P), function(k)
        .indri_gwas_reml_delta(Ys[, k, drop = FALSE], Xs,
                               function(dl) matrix(rot$lambda + dl, ncol = 1L)),
        numeric(1))
    }
    names(delta) <- colnames(Y)
  }
  V <- if (model == "lmm") outer(rot$lambda, rep(1, P)) +
         matrix(delta, N, P, byrow = TRUE)
       else matrix(1, N, P)

  M <- ncol(genotypes)
  empty <- matrix(NA_real_, M, P,
                  dimnames = list(colnames(genotypes), colnames(Y)))
  pvalues <- betas <- p_dom <- b_dom <- empty
  null <- .indri_gwas_null(Ys, Xs, V)
  if (!is.null(env)) {
    out <- .indri_gwas_gxe(null, Xs, V, rot, genotypes, dom, env, estimable,
                           monomorphic, chunk_size, multi = FALSE, N = N,
                           Ti = NULL, snps = colnames(genotypes),
                           traits = colnames(Y))
    out$delta <- if (model == "lmm") delta else NULL
    return(out)
  }
  for (idx in .indri_gwas_chunks(M, chunk_size)) {
    Gs <- .indri_gwas_rotate(rot$Q, genotypes[, idx, drop = FALSE])
    Ds <- if (!is.null(dom)) .indri_gwas_rotate(rot$Q, dom[, idx, drop = FALSE])
    alt <- .indri_gwas_columns(null, Xs, Gs, V, Ds, estimable[idx])
    if (is.null(dom)) {
      stat <- N * log(null$q0 / alt$q1)
      pvalues[idx, ] <- t(stats::pchisq(stat, 1, lower.tail = FALSE))
      betas[idx, ] <- t(alt$beta)
    } else {
      # one rule per SNP for every trait: joint where estimable everywhere
      ok <- matrix(colSums(!alt$ok) == 0L, P, length(idx), byrow = TRUE)
      q_ad <- ifelse(ok, alt$q2, alt$q1)
      joint <- N * log(null$q0 / q_ad)
      pvalues[idx, ] <- t(stats::pchisq(joint, ifelse(ok, 2, 1),
                                        lower.tail = FALSE))
      p_dom[idx, ] <- t(ifelse(ok, stats::pchisq(N * log(alt$q1 / q_ad), 1,
                                                 lower.tail = FALSE), NA_real_))
      betas[idx, ] <- t(ifelse(ok, alt$beta_joint, alt$beta))
      b_dom[idx, ] <- t(ifelse(ok, alt$beta_dom, NA_real_))
    }
  }
  pvalues[monomorphic, ] <- betas[monomorphic, ] <- NA_real_
  p_dom[monomorphic, ] <- b_dom[monomorphic, ] <- NA_real_
  out <- list(pvalues = pvalues, betas = betas,
              delta = if (model == "lmm") delta else NULL)
  if (!is.null(dom)) {
    out$pvalues_dominance <- p_dom
    out$betas_dominance <- b_dom
  }
  out
}

# Weighted least squares of every column of Ys on Xs under the null, kept for
# the SNP tests: per column c, weights 1 / V[, c], the residuals r0 and their
# weighted sum of squares q0.
.indri_gwas_null <- function(Ys, Xs, V) {
  C <- ncol(Ys)
  W <- 1 / V
  r0 <- matrix(0, nrow(Ys), C)
  XtWX_inv <- vector("list", C)
  for (c in seq_len(C)) {
    A <- crossprod(Xs, Xs * W[, c])
    XtWX_inv[[c]] <- chol2inv(chol(A))
    b <- XtWX_inv[[c]] %*% crossprod(Xs, Ys[, c] * W[, c])
    r0[, c] <- Ys[, c] - Xs %*% b
  }
  list(W = W, r0 = r0, q0 = colSums(W * r0^2), XtWX_inv = XtWX_inv)
}

# Add each SNP to every column's regression (Frisch-Waugh: residualise the
# rotated SNP on Xs with the column's weights). Returns, per column (rows)
# and SNP (columns), the additive-only effect and residual sum of squares
# (`beta`, `q1`). With a dominance coding Ds, also the joint additive +
# dominance fit (`beta_joint`, `beta_dom`, `q2`) and `ok`, FALSE where the
# dominance term is not estimable (fewer than three genotype classes, or
# numerically collinear with the dosage).
.indri_gwas_columns <- function(null, Xs, Gs, V, Ds = NULL, estimable = NULL) {
  C <- ncol(null$r0); m <- ncol(Gs)
  beta <- q1 <- matrix(NA_real_, C, m)
  if (!is.null(Ds)) {
    beta_joint <- beta_dom <- q2 <- matrix(NA_real_, C, m)
    ok <- matrix(TRUE, C, m)
  }
  for (c in seq_len(C)) {
    w <- null$W[, c]
    wr <- w * null$r0[, c]
    resid_on_X <- function(M)
      M - Xs %*% (null$XtWX_inv[[c]] %*% crossprod(Xs, M * w))
    Ga <- resid_on_X(Gs)
    ra <- colSums(Ga * wr)
    Saa <- colSums(Ga^2 * w)
    beta[c, ] <- ra / Saa
    q1[c, ] <- null$q0[c] - ra^2 / Saa
    if (!is.null(Ds)) {
      Gd <- resid_on_X(Ds)
      rd <- colSums(Gd * wr)
      Sdd <- colSums(Gd^2 * w)
      Sad <- colSums(Ga * Gd * w)
      det <- Saa * Sdd - Sad^2
      good <- estimable & Sdd > 0 & det > 1e-10 * Saa * Sdd
      good[is.na(good)] <- FALSE
      ba <- (Sdd * ra - Sad * rd) / det
      bd <- (Saa * rd - Sad * ra) / det
      ok[c, ] <- good
      beta_joint[c, ] <- ifelse(good, ba, NA_real_)
      beta_dom[c, ] <- ifelse(good, bd, NA_real_)
      q2[c, ] <- ifelse(good, null$q0[c] - (ba * ra + bd * rd), NA_real_)
    }
  }
  out <- list(beta = beta, q1 = q1)
  if (!is.null(Ds))
    out <- c(out, list(beta_joint = beta_joint, beta_dom = beta_dom, q2 = q2,
                       ok = ok))
  out
}

# Permutation estimate of the false discovery rate. At a threshold t,
# FDR(t) = (mean number of permuted values <= t per permutation) /
#          (number of observed values <= t), capped at 1;
# the q-value of an observed value is the smallest FDR(t) over t at or above
# it, and the threshold the largest observed value with FDR(t) <= fdr (NA if
# none). `observed` holds one value per SNP (NA for untested SNPs), and
# `null_list` the same for each permuted scan.
.indri_gwas_permutation_fdr <- function(observed, null_list, fdr = NULL) {
  q <- rep(NA_real_, length(observed))
  ok <- !is.na(observed)
  t <- observed[ok]
  pooled <- sort(unlist(null_list), na.last = NA)
  n <- length(null_list)
  at <- pmin(1, findInterval(t, pooled) / n /
               findInterval(t, sort(t)))
  o <- order(t, decreasing = TRUE)
  q_ok <- numeric(length(t))
  q_ok[o] <- cummin(at[o])
  q[ok] <- q_ok
  threshold <- NULL
  if (!is.null(fdr)) {
    passing <- t[at <= fdr]
    threshold <- if (length(passing)) max(passing) else NA_real_
  }
  list(q = q, threshold = threshold)
}

# Scans with genotype-by-environment terms. Each SNP brings the terms a
# (dosage), d (heterozygote indicator, with dominance) and one g * e per
# environment; the nested models a, a + d, a + gxe and a + d + gxe are
# fitted per rotated column, and each SNP uses the largest one whose terms
# are estimable in every column (dominance needs three genotype classes).
# Multi-trait scans sum the columns; single-trait scans keep them apart.
.indri_gwas_gxe <- function(null, Xs, V, rot, genotypes, dom, env, estimable,
                            monomorphic, chunk_size, multi, N, Ti, snps,
                            traits) {
  M <- ncol(genotypes)
  C <- ncol(null$r0)
  m_env <- ncol(env)
  has_d <- !is.null(dom)
  i_g <- (2L + has_d):(1L + has_d + m_env)
  subsets <- list(a = 1L, ag = c(1L, i_g))
  if (has_d) subsets <- c(subsets, list(ad = 1:2, adg = c(1:2, i_g)))
  store <- lapply(subsets, function(s)
    list(q = matrix(NA_real_, C, M), ok = matrix(FALSE, C, M),
         beta = array(NA_real_, c(C, M, length(s)))))
  for (idx in .indri_gwas_chunks(M, chunk_size)) {
    g <- genotypes[, idx, drop = FALSE]
    terms <- c(list(.indri_gwas_rotate(rot$Q, g)),
               if (has_d) list(.indri_gwas_rotate(rot$Q, dom[, idx, drop = FALSE])),
               lapply(seq_len(m_env), function(k)
                 .indri_gwas_rotate(rot$Q, g * env[, k])))
    fits <- .indri_gwas_nested(null, Xs, terms, V, subsets)
    for (s in names(subsets)) {
      store[[s]]$q[, idx] <- fits[[s]]$q
      store[[s]]$ok[, idx] <- fits[[s]]$ok
      store[[s]]$beta[, idx, ] <- fits[[s]]$beta
    }
  }
  # one model per SNP, for every column
  all_ok <- function(s) colSums(!store[[s]]$ok) == 0L
  d_ok <- if (has_d) estimable & all_ok("ad") else rep(FALSE, M)
  g_ok <- ifelse(d_ok, if (has_d) all_ok("adg") else FALSE, all_ok("ag"))
  fit_of <- ifelse(d_ok, ifelse(g_ok, "adg", "ad"), ifelse(g_ok, "ag", "a"))
  no_g_of <- ifelse(d_ok, "ad", "a")
  no_d_of <- ifelse(g_ok, "ag", "a")
  pick <- function(which, field = "q") {
    vapply(seq_len(M), function(j) store[[which[j]]][[field]][, j],
           numeric(C))
  }
  q_fit <- matrix(pick(fit_of), C)
  q_no_g <- matrix(pick(no_g_of), C)
  q_no_d <- matrix(pick(no_d_of), C)
  coef <- function(j, pos) store[[fit_of[j]]]$beta[, j, pos]
  pos_d <- 2L
  pos_g <- function(j) if (d_ok[j]) 2L + seq_len(m_env) else 1L + seq_len(m_env)
  b_a <- vapply(seq_len(M), function(j) coef(j, 1L), numeric(C))
  b_d <- vapply(seq_len(M), function(j)
    if (d_ok[j]) coef(j, pos_d) else rep(NA_real_, C), numeric(C))
  b_g <- array(NA_real_, c(C, M, m_env))
  for (j in which(g_ok)) b_g[, j, ] <- store[[fit_of[j]]]$beta[, j, pos_g(j)]
  b_a <- matrix(b_a, C); b_d <- matrix(b_d, C)
  na <- monomorphic

  n_terms <- 1 + d_ok + m_env * g_ok
  if (multi) {
    P <- C
    stat <- N * P * log(sum(null$q0) / colSums(q_fit))
    stat_g <- ifelse(g_ok, N * P * log(colSums(q_no_g) / colSums(q_fit)), NA)
    stat_d <- ifelse(d_ok, N * P * log(colSums(q_no_d) / colSums(q_fit)), NA)
    p <- stats::pchisq(stat, P * n_terms, lower.tail = FALSE)
    p_g <- stats::pchisq(stat_g, P * m_env, lower.tail = FALSE)
    p_d <- stats::pchisq(stat_d, P, lower.tail = FALSE)
    p[na] <- p_g[na] <- p_d[na] <- NA_real_
    back <- function(b) {
      out <- t(b) %*% Ti
      out[na, ] <- NA_real_
      dimnames(out) <- list(snps, traits)
      out
    }
    gxe <- array(NA_real_, c(M, P, m_env),
                 dimnames = list(snps, traits, colnames(env)))
    for (k in seq_len(m_env)) gxe[, , k] <- back(b_g[, , k])
    out <- list(pvalues = stats::setNames(p, snps), betas = back(b_a),
                pvalues_gxe = stats::setNames(p_g, snps), betas_gxe = gxe,
                gxe_estimable = stats::setNames(g_ok, snps))
    if (has_d) {
      out$pvalues_dominance <- stats::setNames(p_d, snps)
      out$betas_dominance <- back(b_d)
      out$dominance_estimable <- stats::setNames(d_ok, snps)
    }
    return(out)
  }
  per <- function(x) {
    x <- t(x)
    x[na, ] <- NA_real_
    dimnames(x) <- list(snps, traits)
    x
  }
  rep_cols <- function(v) matrix(v, C, M, byrow = TRUE)
  stat <- N * log(null$q0 / q_fit)
  stat_g <- ifelse(rep_cols(g_ok), N * log(q_no_g / q_fit), NA)
  stat_d <- ifelse(rep_cols(d_ok), N * log(q_no_d / q_fit), NA)
  p <- stats::pchisq(stat, rep_cols(n_terms), lower.tail = FALSE)
  gxe <- array(NA_real_, c(M, C, m_env),
               dimnames = list(snps, traits, colnames(env)))
  for (k in seq_len(m_env)) gxe[, , k] <- per(b_g[, , k])
  out <- list(pvalues = per(p), betas = per(b_a),
              pvalues_gxe = per(stats::pchisq(stat_g, m_env,
                                              lower.tail = FALSE)),
              betas_gxe = gxe, gxe_estimable = stats::setNames(g_ok, snps))
  if (has_d) {
    out$pvalues_dominance <- per(stats::pchisq(stat_d, 1, lower.tail = FALSE))
    out$betas_dominance <- per(b_d)
    out$dominance_estimable <- stats::setNames(d_ok, snps)
  }
  out
}

# Scans with an environment-specific background (fit_G_env()). The
# covariance of vec(Y) (trait-major) is V(delta) = Vg + delta Ve, with
# Vg[(k,i),(l,j)] = A_ij Gs[(k, env_i), (l, env_j)] and Ve block diagonal by
# individual. That is no Kronecker product, but it is a pencil of two
# matrices: with Ve = L L' and L^-1 Vg L^-T = U diag(lambda) U', the rotation
# R = U' L^-1 gives R V(delta) R' = diag(lambda + delta) for every delta. A
# block-diagonal A (families) makes V block diagonal, so the decomposition
# runs per connected component of A. After the rotation each SNP is a
# weighted regression of the rotated response on rotated terms, one column
# per trait and SNP term, as in .indri_gwas_nested().
.indri_gwas_env <- function(Y, X, genotypes, env_fit, A, mode, delta,
                            chunk_size, monomorphic, dom, estimable, env) {
  N <- nrow(Y)
  traits <- colnames(Y)
  comps <- .indri_components(A)
  problems <- if (mode == "multi_trait") list(traits) else as.list(traits)
  has_d <- !is.null(dom)
  m_env <- if (is.null(env)) 0L else ncol(env)
  M <- ncol(genotypes)
  snps <- colnames(genotypes)
  results <- lapply(problems, function(tr) {
    p <- length(tr)
    rot <- lapply(comps, function(idx)
      .indri_env_rotation(A[idx, idx, drop = FALSE], env_fit$env[idx],
                          env_fit$Gs, env_fit$Es, tr))
    rotate_terms <- function(T) {
      # one rotated column per trait for each column of T (individuals x m)
      lapply(seq_len(p), function(k) do.call(rbind, lapply(seq_along(comps),
        function(c) {
          idx <- comps[[c]]
          block <- (k - 1L) * length(idx) + seq_along(idx)
          rot[[c]]$R[, block, drop = FALSE] %*% T[idx, , drop = FALSE]
        })))
    }
    ys <- do.call(c, lapply(seq_along(comps), function(c)
      drop(rot[[c]]$R %*% c(Y[comps[[c]], tr, drop = FALSE]))))
    Xs <- do.call(cbind, rotate_terms(X))
    lambda <- do.call(c, lapply(rot, `[[`, "lambda"))
    V_of <- function(dl) matrix(lambda + dl, ncol = 1L)
    dl <- if (!is.null(delta)) {
      if (mode == "multi_trait") delta else delta[[tr]]
    } else {
      .indri_gwas_reml_delta(matrix(ys, ncol = 1L), Xs, V_of)
    }
    V <- V_of(dl)
    null <- .indri_gwas_null(matrix(ys, ncol = 1L), Xs, V)
    # term groups: a, [d], gxe_1..gxe_m, each with one column per trait
    group <- function(g) (g - 1L) * p + seq_len(p)
    i_a <- group(1L)
    i_d <- if (has_d) group(2L) else integer()
    i_g <- unlist(lapply(seq_len(m_env), function(k) group(1L + has_d + k)))
    subsets <- list(a = i_a)
    if (has_d) subsets$ad <- c(i_a, i_d)
    if (m_env) subsets$ag <- c(i_a, i_g)
    if (has_d && m_env) subsets$adg <- c(i_a, i_d, i_g)
    store <- lapply(subsets, function(s)
      list(q = rep(NA_real_, M), ok = rep(FALSE, M),
           beta = matrix(NA_real_, M, length(s))))
    for (ch in .indri_gwas_chunks(M, chunk_size)) {
      g <- genotypes[, ch, drop = FALSE]
      terms <- c(rotate_terms(g),
                 if (has_d) rotate_terms(dom[, ch, drop = FALSE]),
                 unlist(lapply(seq_len(m_env), function(k)
                   rotate_terms(g * env[, k])), recursive = FALSE))
      fits <- .indri_gwas_nested(null, Xs, terms, V, subsets)
      for (s in names(subsets)) {
        store[[s]]$q[ch] <- fits[[s]]$q[1L, ]
        store[[s]]$ok[ch] <- fits[[s]]$ok[1L, ]
        store[[s]]$beta[ch, ] <- fits[[s]]$beta[1L, , ]
      }
    }
    d_ok <- if (has_d) estimable & store$ad$ok else rep(FALSE, M)
    g_ok <- if (!m_env) rep(FALSE, M) else
      ifelse(d_ok, if (has_d) store$adg$ok else FALSE, store$ag$ok)
    fit_of <- ifelse(d_ok, ifelse(g_ok, "adg", "ad"), ifelse(g_ok, "ag", "a"))
    q_of <- function(which) vapply(seq_len(M), function(j)
      store[[which[j]]]$q[j], numeric(1))
    q_fit <- q_of(fit_of)
    q_no_g <- q_of(ifelse(d_ok, "ad", "a"))
    q_no_d <- q_of(ifelse(g_ok, "ag", "a"))
    n <- length(ys)
    n_terms <- 1 + d_ok + m_env * g_ok
    coef <- function(j, cols) store[[fit_of[j]]]$beta[j, cols]
    pos <- function(j, what) {
      if (what == "a") return(seq_len(p))
      if (what == "d") return(p + seq_len(p))
      offset <- p * (1L + d_ok[j])
      offset + seq_len(p * m_env)
    }
    list(
      p = stats::pchisq(n * log(null$q0 / q_fit), p * n_terms,
                        lower.tail = FALSE),
      p_g = ifelse(g_ok, stats::pchisq(n * log(q_no_g / q_fit), p * m_env,
                                       lower.tail = FALSE), NA_real_),
      p_d = ifelse(d_ok, stats::pchisq(n * log(q_no_d / q_fit), p,
                                       lower.tail = FALSE), NA_real_),
      a = t(vapply(seq_len(M), function(j) coef(j, pos(j, "a")), numeric(p))),
      d = t(vapply(seq_len(M), function(j)
        if (d_ok[j]) coef(j, pos(j, "d")) else rep(NA_real_, p), numeric(p))),
      g = vapply(seq_len(M), function(j)
        if (g_ok[j]) coef(j, pos(j, "g")) else rep(NA_real_, p * m_env),
        numeric(p * m_env)),
      d_ok = d_ok, g_ok = g_ok, delta = dl)
  })

  P <- length(traits)
  na <- monomorphic
  shape <- function(x) {
    x <- matrix(x, M)
    x[na, ] <- NA_real_
    dimnames(x) <- list(snps, if (ncol(x) == P) traits else NULL)
    x
  }
  if (mode == "multi_trait") {
    r <- results[[1L]]
    p <- r$p; p_g <- r$p_g; p_d <- r$p_d
    p[na] <- p_g[na] <- p_d[na] <- NA_real_
    out <- list(pvalues = stats::setNames(p, snps), betas = shape(r$a),
                delta = unname(r$delta))
    gxe <- if (m_env) array(NA_real_, c(M, P, m_env),
                            dimnames = list(snps, traits, colnames(env)))
    if (m_env) for (k in seq_len(m_env))
      gxe[, , k] <- shape(t(matrix(r$g, P * m_env)[(k - 1L) * P + seq_len(P), ,
                                                    drop = FALSE]))
    if (m_env) {
      out$pvalues_gxe <- stats::setNames(p_g, snps)
      out$betas_gxe <- gxe
      out$gxe_estimable <- stats::setNames(r$g_ok, snps)
    }
    if (has_d) {
      out$pvalues_dominance <- stats::setNames(p_d, snps)
      out$betas_dominance <- shape(r$d)
      out$dominance_estimable <- stats::setNames(r$d_ok, snps)
    }
  } else {
    col <- function(field) shape(vapply(results, function(r) r[[field]],
                                        numeric(M)))
    out <- list(pvalues = col("p"), betas = col("a"),
                delta = stats::setNames(vapply(results, `[[`, numeric(1),
                                               "delta"), traits))
    if (m_env) {
      out$pvalues_gxe <- col("p_g")
      gxe <- array(NA_real_, c(M, P, m_env),
                   dimnames = list(snps, traits, colnames(env)))
      for (k in seq_len(P)) gxe[, k, ] <- t(matrix(results[[k]]$g, m_env))
      gxe[na, , ] <- NA_real_
      out$betas_gxe <- gxe
      out$gxe_estimable <- stats::setNames(
        Reduce(`&`, lapply(results, `[[`, "g_ok")), snps)
    }
    if (has_d) {
      out$pvalues_dominance <- col("p_d")
      out$betas_dominance <- col("d")
      out$dominance_estimable <- stats::setNames(
        Reduce(`&`, lapply(results, `[[`, "d_ok")), snps)
    }
  }
  out$n_components <- length(comps)
  out
}

# The generalized rotation of one connected component: Vg and Ve as in
# .indri_gwas_env(), for the individuals of the component and traits `tr`
# (trait-major), and R = U' L^-1 with the eigenvalues lambda.
.indri_env_rotation <- function(A, env, Gs, Es, tr) {
  n <- nrow(A)
  p <- length(tr)
  Vg <- Ve <- matrix(0, n * p, n * p)
  for (k in seq_len(p)) for (l in seq_len(p)) {
    sk <- paste(tr[k], env, sep = ":")
    sl <- paste(tr[l], env, sep = ":")
    rows <- (k - 1L) * n + seq_len(n)
    cols <- (l - 1L) * n + seq_len(n)
    Vg[rows, cols] <- A * Gs[sk, sl, drop = FALSE]
    Ve[rows, cols] <- diag(Es[cbind(sk, sl)], n)
  }
  L <- t(chol(Ve))
  Li <- forwardsolve(L, diag(n * p))      # L^-1
  W <- Li %*% Vg %*% t(Li)
  e <- eigen((W + t(W)) / 2, symmetric = TRUE)
  list(R = crossprod(e$vectors, Li), lambda = pmax(e$values, 0))
}

# Connected components of the relationship graph (nonzero entries of A).
.indri_components <- function(A) {
  N <- nrow(A)
  linked <- A != 0
  comp <- integer(N)
  current <- 0L
  for (start in seq_len(N)) {
    if (comp[start]) next
    current <- current + 1L
    queue <- start
    comp[start] <- current
    while (length(queue)) {
      i <- queue[1L]
      queue <- queue[-1L]
      new <- which(linked[i, ] & comp == 0L)
      comp[new] <- current
      queue <- c(queue, new)
    }
  }
  split(seq_len(N), comp)
}

# Nested weighted regressions with several SNP terms: for each rotated
# column c and each subset s of the terms, the reduction in the weighted
# residual sum of squares when the terms of s join Xs (Frisch-Waugh), their
# coefficients, and `ok`, FALSE where the terms are (numerically) collinear.
.indri_gwas_nested <- function(null, Xs, terms, V, subsets) {
  C <- ncol(null$r0)
  m <- ncol(terms[[1L]])
  K <- length(terms)
  res <- lapply(subsets, function(s)
    list(q = matrix(NA_real_, C, m), ok = matrix(FALSE, C, m),
         beta = array(NA_real_, c(C, m, length(s)))))
  for (c in seq_len(C)) {
    w <- null$W[, c]
    wr <- w * null$r0[, c]
    Tt <- lapply(terms, function(T)
      T - Xs %*% (null$XtWX_inv[[c]] %*% crossprod(Xs, T * w)))
    S <- array(0, c(K, K, m))
    r <- matrix(0, K, m)
    for (k in seq_len(K)) {
      r[k, ] <- colSums(Tt[[k]] * wr)
      for (l in k:K) {
        v <- colSums(Tt[[k]] * Tt[[l]] * w)
        S[k, l, ] <- v
        S[l, k, ] <- v
      }
    }
    for (si in seq_along(subsets)) {
      s <- subsets[[si]]
      for (j in seq_len(m)) {
        Ss <- matrix(S[s, s, j], length(s))
        ds <- diag(Ss)
        if (any(!(ds > 0))) next
        scale <- 1 / sqrt(ds)
        ev <- eigen(Ss * outer(scale, scale), symmetric = TRUE,
                    only.values = TRUE)$values
        if (min(ev) < 1e-10 * max(ev)) next
        b <- solve(Ss, r[s, j])
        res[[si]]$q[c, j] <- null$q0[c] - sum(b * r[s, j])
        res[[si]]$ok[c, j] <- TRUE
        res[[si]]$beta[c, j, ] <- b
      }
    }
  }
  names(res) <- names(subsets)
  res
}

# REML for the variance ratio, sigma^2 profiled: V_of(delta) gives the
# per-entry variances of the rotated response for that ratio.
.indri_gwas_reml_delta <- function(Ys, Xs, V_of) {
  n <- length(Ys) - ncol(Ys) * ncol(Xs)
  criterion <- function(log_delta) {
    V <- V_of(exp(log_delta))
    W <- 1 / V
    q <- 0
    logdet_A <- 0
    for (c in seq_len(ncol(Ys))) {
      A <- crossprod(Xs, Xs * W[, c])
      R <- chol(A)
      b <- backsolve(R, backsolve(R, crossprod(Xs, Ys[, c] * W[, c]),
                                  transpose = TRUE))
      r <- Ys[, c] - Xs %*% b
      q <- q + sum(W[, c] * r^2)
      logdet_A <- logdet_A + 2 * sum(log(diag(R)))
    }
    0.5 * (sum(log(V)) + logdet_A + n * log(q / n))
  }
  exp(stats::optimize(criterion, c(-12, 12), tol = 1e-9)$minimum)
}

.indri_gwas_rotation <- function(A) {
  e <- eigen((A + t(A)) / 2, symmetric = TRUE)
  if (min(e$values) < -1e-8 * max(abs(e$values)))
    stop("`A` must be positive semi-definite.", call. = FALSE)
  list(Q = e$vectors, lambda = pmax(e$values, 0))
}

.indri_gwas_rotate <- function(Q, x) if (is.null(Q)) unname(x) else crossprod(Q, x)

.indri_gwas_chunks <- function(M, chunk_size) {
  split(seq_len(M), ceiling(seq_len(M) / chunk_size))
}

# ---------------------------------------------------------------------------
# Input handling.

.indri_gwas_matrix <- function(x, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x))
    stop("`", name, "` must be numeric.", call. = FALSE)
  storage.mode(x) <- "double"
  x
}

# Put the rows of `x` in the order of the samples of Y: by name when both are
# named, otherwise by position.
.indri_gwas_rows <- function(x, ids, N, name) {
  if (!is.null(ids) && !is.null(rownames(x))) {
    missing <- setdiff(ids, rownames(x))
    if (length(missing))
      stop("`", name, "` has no rows for ", length(missing),
           " sample(s) of `Y`: ",
           paste(utils::head(missing, 5L), collapse = ", "),
           if (length(missing) > 5L) ", ...", call. = FALSE)
    return(x[ids, , drop = FALSE])
  }
  if (nrow(x) != N)
    stop("`", name, "` has ", nrow(x), " rows but `Y` has ", N,
         "; name the rows to match samples by name.", call. = FALSE)
  x
}

.indri_gwas_trait_matrix <- function(x, traits, name) {
  x <- .indri_gwas_matrix(x, name)
  P <- length(traits)
  if (!is.null(rownames(x)) && !is.null(colnames(x))) {
    if (!setequal(rownames(x), traits) || !setequal(colnames(x), traits))
      stop("`", name, "` does not describe the traits of `Y` (",
           paste(utils::head(traits, 5L), collapse = ", "),
           if (P > 5L) ", ...", ").", call. = FALSE)
    x <- x[traits, traits, drop = FALSE]
  } else if (nrow(x) != P || ncol(x) != P) {
    stop("`", name, "` must be ", P, " x ", P, ", one row and column per ",
         "trait of `Y`.", call. = FALSE)
  }
  if (any(!is.finite(x)) || max(abs(x - t(x))) > 1e-8 * max(abs(x), 1))
    stop("`", name, "` must be a finite symmetric matrix.", call. = FALSE)
  x <- (x + t(x)) / 2
  dimnames(x) <- list(traits, traits)
  x
}

.indri_gwas_chol <- function(x, name) {
  L <- tryCatch(t(chol(x)), error = function(e) NULL)
  if (is.null(L))
    stop("`", name, "` must be positive definite.", call. = FALSE)
  L
}

.indri_gwas_delta <- function(delta, mode, traits) {
  if (!is.numeric(delta) || any(!is.finite(delta)) || any(delta <= 0))
    stop("`delta` must be positive and finite.", call. = FALSE)
  if (mode == "multi_trait") {
    if (length(delta) != 1L)
      stop("`delta` must be a single number for the multi-trait test.",
           call. = FALSE)
    return(unname(delta))
  }
  if (length(delta) != length(traits))
    stop("`delta` must have one value per trait (", length(traits), ").",
         call. = FALSE)
  if (!is.null(names(delta))) {
    if (!setequal(names(delta), traits))
      stop("The names of `delta` must be the traits of `Y`.", call. = FALSE)
    return(delta[traits])
  }
  stats::setNames(delta, traits)
}

.indri_gwas_impute <- function(genotypes) {
  if (!anyNA(genotypes)) return(genotypes)
  for (j in which(colSums(is.na(genotypes)) > 0)) {
    g <- genotypes[, j]
    g[is.na(g)] <- if (all(is.na(g))) 0 else mean(g, na.rm = TRUE)
    genotypes[, j] <- g
  }
  genotypes
}

# ---------------------------------------------------------------------------
# Methods.

#' @export
print.indri_gwas <- function(x, ...) {
  what <- if (x$mode == "multi_trait") "multi-trait" else "single-trait"
  cat("indri association scan:", what, toupper(x$model),
      if (isTRUE(x$dominance)) "with additive and dominance effects",
      if (!is.null(x$gxe_terms))
        paste0("with interactions with ", paste(x$gxe_terms, collapse = ", ")),
      "\n")
  M <- NROW(x$pvalues)
  cat(sprintf("  %d SNPs, %d traits, %d individuals\n", M, length(x$traits),
              x$n_individuals))
  if (!is.null(x$delta)) {
    if (length(x$delta) == 1L) cat(sprintf("  null variance ratio delta: %.3g\n",
                                           x$delta))
    else cat(sprintf("  null variance ratio delta: %.3g to %.3g over traits\n",
                     min(x$delta), max(x$delta)))
  }
  if (!is.null(x$n_tests))
    cat(sprintf("  p-values adjusted for %.3g tests\n", x$n_tests))
  if (!is.null(x$n_permutations))
    cat(sprintf("  empirical p-values from %d permutations\n", x$n_permutations))
  if (!is.null(x$fdr))
    cat(.indri_gwas_fdr_line(x$fdr, x$fdr_threshold))
  if (all(is.na(x$pvalues))) {
    cat("  no SNP could be tested\n")
  } else {
    cat(sprintf("  smallest p-value: %.3g\n", min(x$pvalues, na.rm = TRUE)))
  }
  invisible(x)
}

#' @export
as.data.frame.indri_gwas <- function(x, row.names = NULL, optional = FALSE,
                                     ...) {
  info <- x$snp_info
  if (!"snp" %in% names(info))
    info <- cbind(snp = rownames(x$a), info, stringsAsFactors = FALSE)
  if (x$mode == "multi_trait") {
    a <- as.data.frame(x$a)
    names(a) <- paste0("a_", colnames(x$a))
    out <- cbind(info, pvalue = unname(x$pvalues), a)
    if (!is.null(x$pvalues_empirical)) {
      out$pvalue_empirical <- unname(x$pvalues_empirical)
      out$qvalue_empirical <- unname(x$qvalues_empirical)
    }
    if (!is.null(x$gxe_terms)) {
      out$pvalue_gxe <- unname(x$pvalues_gxe)
      for (e in x$gxe_terms) for (k in colnames(x$a))
        out[[paste0("gxe_", e, "_", k)]] <- unname(x$gxe[, k, e])
      out$gxe_estimable <- unname(x$gxe_estimable)
    }
    if (isTRUE(x$dominance)) {
      dom <- as.data.frame(x$d)
      names(dom) <- paste0("d_", colnames(x$d))
      out <- cbind(out, pvalue_dominance = unname(x$pvalues_dominance), dom,
                   dominance_estimable = unname(x$dominance_estimable))
    }
  } else {
    out <- do.call(rbind, lapply(colnames(x$pvalues), function(k) {
      part <- cbind(info, trait = k, pvalue = unname(x$pvalues[, k]),
                    a = unname(x$a[, k]), stringsAsFactors = FALSE)
      if (!is.null(x$pvalues_adjust))
        part$pvalue_adjust <- unname(x$pvalues_adjust[, k])
      if (!is.null(x$pvalues_empirical)) {
        part$pvalue_empirical <- unname(x$pvalues_empirical[, k])
        part$qvalue_empirical <- unname(x$qvalues_empirical)
      }
      if (!is.null(x$gxe_terms)) {
        part$pvalue_gxe <- unname(x$pvalues_gxe[, k])
        for (e in x$gxe_terms) part[[paste0("gxe_", e)]] <- unname(x$gxe[, k, e])
        part$gxe_estimable <- unname(x$gxe_estimable)
      }
      if (isTRUE(x$dominance)) {
        part$pvalue_dominance <- unname(x$pvalues_dominance[, k])
        part$d <- unname(x$d[, k])
        part$dominance_estimable <- unname(x$dominance_estimable)
      }
      part
    }))
  }
  rownames(out) <- row.names
  out
}

.indri_gwas_fdr_line <- function(fdr, threshold) {
  if (is.na(threshold))
    sprintf("  no SNP passes FDR %.3g (permutation estimate)\n", fdr)
  else sprintf("  FDR %.3g (permutation estimate): p <= %.3g\n", fdr, threshold)
}
