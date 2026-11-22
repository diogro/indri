#' Association scan with GEMMA
#'
#' Runs the linear mixed model association tests of GEMMA (Zhou and Stephens
#' 2012, 2014) from R and returns them as an `indri_gwas` scan, so that the
#' `plot()`, `summary()` and [indri_locus_G()] methods work on its results as
#' on those of [indri_gwas()]. GEMMA fits the genetic and residual covariances
#' itself, by REML under the null model, for each relationship matrix; no
#' [fit_G()] is needed.
#'
#' @param Y Numeric matrix, individuals by traits, with no missing values.
#' @param genotypes Numeric matrix, individuals by SNPs: dosages (0 to 2) of
#'   the allele the effects refer to. GEMMA replaces missing values by the
#'   SNP's mean.
#' @param A Relationship matrix of the individuals (in the row order of `Y`),
#'   or, with `loco = TRUE`, a list of them named by chromosome. `NULL` (the
#'   default) computes GEMMA's centred genomic relationship matrix from
#'   `genotypes`, \eqn{Z Z' / M} with each SNP centred on its mean (`-gk 1`);
#'   with `loco = TRUE`, one per chromosome from the SNPs on the other
#'   chromosomes.
#' @param X Fixed-effects design with full column rank, including the
#'   intercept (GEMMA adds none when covariates are given); defaults to an
#'   intercept.
#' @param snp_info Data frame with one row per column of `genotypes` and
#'   columns `snp`, `chrom` and `pos`; required for `loco = TRUE`.
#' @param mode `"multi_trait"` for GEMMA's multivariate LMM (one joint test of
#'   all traits per SNP), `"single_trait"` for a univariate LMM per trait.
#' @param test GEMMA's test: `"wald"` (`-lmm 1`), `"lrt"` (`-lmm 2`) or
#'   `"score"` (`-lmm 3`).
#' @param loco If `TRUE`, test the SNPs of each chromosome with a relationship
#'   matrix that leaves that chromosome out.
#' @param maf,miss GEMMA's SNP filters: minimum minor allele frequency and
#'   maximum missingness. SNPs that fail them get `NA`.
#' @param gemma Path to the GEMMA executable; by default
#'   `getOption("indri.gemma")`, then `gemma` on the PATH.
#' @param cores Number of GEMMA runs at a time (forked; one on Windows).
#' @param keep_dir Directory in which to keep GEMMA's input and output files,
#'   one subdirectory per run; by default they are written to a temporary
#'   directory and removed.
#'
#' @return An object of class `indri_gwas`, as from [indri_gwas()]: `pvalues`
#'   (a vector over SNPs, or SNPs by traits), `a` (the effects, SNPs by traits,
#'   in the units of `Y` per copy of the allele `genotypes` counts),
#'   `snp_info`, `mode`, `model` (`"lmm"`), `traits`, `n_individuals` and
#'   `call`; and from GEMMA, `engine = "gemma"`, `test`, `se` (standard errors
#'   of `a`) and `null`, GEMMA's REML estimates of the genetic (`Vg`) and
#'   residual (`Ve`) covariances under the null model, one per run.
#'
#' @details
#' Each run (one per chromosome with `loco = TRUE`, one per trait in
#' single-trait mode) writes its own BIMBAM files to its own directory, so runs
#' in parallel cannot read each other's results. Effects are per copy of the
#' allele that `genotypes` counts: the BIMBAM files list it as GEMMA's effect
#' allele.
#'
#' @references Zhou X, Stephens M (2012) Genome-wide efficient mixed-model
#'   analysis for association studies. Nature Genetics 44: 821-824.
#'   Zhou X, Stephens M (2014) Efficient multivariate linear mixed model
#'   algorithms for genome-wide association studies. Nature Methods 11: 407-409.
#' @seealso [indri_gwas()] for the same tests with the background covariances
#'   from [fit_G()]; [indri_gemma_kinship()] for the relationship matrices.
#' @export
indri_gemma <- function(Y, genotypes, A = NULL, X = NULL, snp_info = NULL,
                        mode = c("multi_trait", "single_trait"),
                        test = c("wald", "lrt", "score"), loco = FALSE,
                        maf = 0.01, miss = 0.05, gemma = NULL, cores = 1L,
                        keep_dir = NULL) {
  call <- match.call()
  mode <- match.arg(mode)
  test <- match.arg(test)
  cores <- .indri_count(cores, "cores")
  binary <- .indri_gemma_binary(gemma)

  Y <- .indri_gwas_matrix(Y, "Y")
  if (anyNA(Y))
    stop("`Y` has missing values; indri_gemma() needs complete phenotypes ",
         "(GEMMA's multivariate model would drop those individuals).", call. = FALSE)
  N <- nrow(Y)
  P <- ncol(Y)
  traits <- colnames(Y) %||% paste0("trait", seq_len(P))
  colnames(Y) <- traits
  ids <- rownames(Y) %||% as.character(seq_len(N))
  if (mode == "multi_trait" && P < 2L)
    stop("`mode = \"multi_trait\"` needs at least two traits.", call. = FALSE)

  genotypes <- .indri_gwas_matrix(genotypes, "genotypes")
  if (nrow(genotypes) != N)
    stop("`genotypes` must have one row per row of `Y`.", call. = FALSE)
  M <- ncol(genotypes)
  snps <- colnames(genotypes) %||% paste0("snp", seq_len(M))
  colnames(genotypes) <- snps
  if (!is.null(rownames(genotypes)) && !is.null(rownames(Y)) &&
      !identical(rownames(genotypes), rownames(Y)))
    stop("The rows of `genotypes` are not in the order of the rows of `Y`.",
         call. = FALSE)

  if (is.null(X)) X <- matrix(1, N, 1L, dimnames = list(NULL, "intercept"))
  X <- .indri_gwas_matrix(X, "X")
  if (nrow(X) != N) stop("`X` must have one row per row of `Y`.", call. = FALSE)
  if (anyNA(X)) stop("`X` has missing values.", call. = FALSE)
  if (qr(X)$rank < ncol(X))
    stop("`X` does not have full column rank.", call. = FALSE)
  if (qr(cbind(X, 1))$rank > ncol(X))
    stop("`X` cannot represent a constant; GEMMA adds no intercept when ",
         "covariates are given, so include one.", call. = FALSE)

  if (!is.null(snp_info)) {
    snp_info <- as.data.frame(snp_info, stringsAsFactors = FALSE)
    if (nrow(snp_info) != M)
      stop("`snp_info` must have one row per column of `genotypes`.", call. = FALSE)
    if (is.null(snp_info$snp)) snp_info$snp <- snps
  }
  if (loco && (is.null(snp_info) || is.null(snp_info$chrom)))
    stop("`loco = TRUE` needs `snp_info` with a `chrom` column.", call. = FALSE)

  # relationship matrices, one per run group (chromosome, or all SNPs)
  groups <- if (loco) as.character(snp_info$chrom) else rep("all", M)
  group_names <- unique(groups)
  kin <- if (is.null(A)) {
    if (loco) indri_gemma_kinship(genotypes, groups)
    else list(all = indri_gemma_kinship(genotypes))
  } else if (loco) {
    if (!is.list(A) || !all(group_names %in% names(A)))
      stop("With `loco = TRUE`, `A` must be a list of relationship matrices ",
           "named by chromosome (", paste(utils::head(group_names, 3L), collapse = ", "),
           ", ...).", call. = FALSE)
    A[group_names]
  } else {
    list(all = A)
  }
  kin <- lapply(kin, function(K) {
    K <- as.matrix(K)
    if (!identical(dim(K), c(N, N)))
      stop("Relationship matrices must be N x N, in the row order of `Y`.",
           call. = FALSE)
    if (!is.null(rownames(K)) && !is.null(rownames(Y)) &&
        !identical(rownames(K), rownames(Y)))
      stop("The relationship matrix is not in the row order of `Y`.", call. = FALSE)
    K
  })

  runs <- expand.grid(group = group_names,
                      trait = if (mode == "multi_trait") "all" else traits,
                      stringsAsFactors = FALSE)
  base <- keep_dir %||% tempfile("indri-gemma-")
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  if (is.null(keep_dir)) on.exit(unlink(base, recursive = TRUE), add = TRUE)

  one <- function(r) {
    g <- runs$group[r]; tr <- runs$trait[r]
    cols <- which(groups == g)
    dir <- file.path(base, paste0("run_", r))
    .indri_gemma_run(binary, dir, Y[, if (tr == "all") traits else tr, drop = FALSE],
                     genotypes[, cols, drop = FALSE], kin[[g]], X,
                     snp_info[cols, , drop = FALSE], test, maf, miss)
  }
  res <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(runs)), one, mc.cores = cores)
  } else {
    lapply(seq_len(nrow(runs)), one)
  }
  failed <- vapply(res, inherits, logical(1), "try-error")
  if (any(failed))
    stop("A GEMMA run failed: ", as.character(res[[which(failed)[1]]]), call. = FALSE)

  pcol <- c(wald = "p_wald", lrt = "p_lrt", score = "p_score")[[test]]
  if (mode == "multi_trait") {
    pvalues <- stats::setNames(rep(NA_real_, M), snps)
    a <- se <- matrix(NA_real_, M, P, dimnames = list(snps, traits))
    for (r in seq_along(res)) {
      d <- res[[r]]$assoc
      i <- match(d$rs, snps)
      pvalues[i] <- d[[pcol]]
      a[i, ] <- as.matrix(d[, paste0("beta_", seq_len(P))])
      se[i, ] <- sqrt(as.matrix(d[, paste0("Vbeta_", seq_len(P), "_", seq_len(P))]))
    }
  } else {
    pvalues <- a <- se <- matrix(NA_real_, M, P, dimnames = list(snps, traits))
    for (r in seq_along(res)) {
      d <- res[[r]]$assoc
      i <- match(d$rs, snps); k <- runs$trait[r]
      pvalues[i, k] <- d[[pcol]]
      if (!is.null(d$beta)) a[i, k] <- d$beta
      if (!is.null(d$se)) se[i, k] <- d$se
    }
  }
  null <- lapply(res, `[[`, "null")
  names(null) <- if (mode == "multi_trait") runs$group
                 else paste(runs$group, runs$trait, sep = ":")

  structure(list(pvalues = pvalues, a = a, se = se, delta = NA_real_,
                 pvalues_adjust = NULL, n_tests = NULL,
                 snp_info = snp_info %||% data.frame(snp = snps),
                 dominance = FALSE, gxe_terms = NULL, background = "kronecker",
                 mode = mode, model = "lmm", n_individuals = N, traits = traits,
                 engine = "gemma", test = test, loco = loco, null = null,
                 call = call),
            class = "indri_gwas")
}

#' Genomic relationship matrices as GEMMA computes them
#'
#' The centred relationship matrix of GEMMA's `-gk 1`: \eqn{Z Z' / M}, where
#' `Z` holds the `M` SNP dosages each centred on its mean (missing dosages
#' set to the mean). With `groups`, one matrix per group (chromosome) from the
#' SNPs of all other groups, the leave-one-chromosome-out (LOCO) matrices.
#'
#' @param genotypes Numeric matrix, individuals by SNPs (dosages).
#' @param groups Optional vector, one value per SNP (typically the
#'   chromosome).
#' @return A relationship matrix (individuals by individuals), or with
#'   `groups` a list of them named by group.
#' @export
indri_gemma_kinship <- function(genotypes, groups = NULL) {
  genotypes <- .indri_gwas_matrix(genotypes, "genotypes")
  Z <- sweep(genotypes, 2L, colMeans(genotypes, na.rm = TRUE))
  Z[is.na(Z)] <- 0
  keep <- colSums(Z^2) > 0          # monomorphic SNPs carry no information
  Z <- Z[, keep, drop = FALSE]
  ids <- rownames(genotypes)
  if (is.null(groups)) {
    K <- tcrossprod(Z) / ncol(Z)
    dimnames(K) <- list(ids, ids)
    return(K)
  }
  if (length(groups) != ncol(genotypes))
    stop("`groups` must have one value per column of `genotypes`.", call. = FALSE)
  groups <- as.character(groups)[keep]
  names_g <- unique(as.character(groups))
  parts <- lapply(names_g, function(g) tcrossprod(Z[, groups == g, drop = FALSE]))
  total <- Reduce(`+`, parts)
  m <- vapply(names_g, function(g) sum(groups == g), numeric(1))
  out <- lapply(seq_along(names_g), function(j) {
    K <- (total - parts[[j]]) / (sum(m) - m[j])
    dimnames(K) <- list(ids, ids)
    K
  })
  stats::setNames(out, names_g)
}

# The GEMMA executable, or an error naming where indri looked.
.indri_gemma_binary <- function(path = NULL) {
  path <- path %||% getOption("indri.gemma") %||% unname(Sys.which("gemma"))
  if (!nzchar(path) || !file.exists(path))
    stop("indri_gemma() needs the GEMMA program (https://github.com/",
         "genetics-statistics/GEMMA). Put it on the PATH, or give its location ",
         "with `gemma = <path>` or `options(indri.gemma = <path>)`.", call. = FALSE)
  normalizePath(path)
}

# One GEMMA run in its own directory: BIMBAM files, the LMM, and its results.
.indri_gemma_run <- function(binary, dir, Y, genotypes, K, X, snp_info, test,
                             maf, miss) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  f <- function(x) file.path(dir, x)
  num <- function(x) formatC(x, digits = 10, format = "g")
  # BIMBAM mean genotypes: SNP, effect allele, other allele, dosages
  g <- t(genotypes)
  g_txt <- matrix(ifelse(is.na(g), "NA", num(g)), nrow(g))
  writeLines(paste(rownames(g), "E", "O", apply(g_txt, 1L, paste, collapse = " ")),
             f("geno.txt"))
  utils::write.table(Y, f("pheno.txt"), row.names = FALSE, col.names = FALSE)
  utils::write.table(X, f("covar.txt"), row.names = FALSE, col.names = FALSE)
  utils::write.table(K, f("kinship.txt"), row.names = FALSE, col.names = FALSE)
  if (!is.null(snp_info) && all(c("chrom", "pos") %in% names(snp_info))) {
    utils::write.table(data.frame(rownames(g), snp_info$pos, snp_info$chrom),
                       f("annot.txt"), sep = ",", row.names = FALSE,
                       col.names = FALSE, quote = FALSE)
    annot <- c("-a", f("annot.txt"))
  } else annot <- character(0)
  lmm <- c(wald = "1", lrt = "2", score = "3")[[test]]
  args <- c("-g", f("geno.txt"), "-p", f("pheno.txt"), "-c", f("covar.txt"),
            "-k", f("kinship.txt"), annot, "-n", seq_len(ncol(Y)),
            "-lmm", lmm, "-maf", maf, "-miss", miss,
            "-outdir", f("output"), "-o", "result")
  log <- suppressWarnings(system2(binary, args, stdout = TRUE, stderr = TRUE))
  assoc_file <- f("output/result.assoc.txt")
  if (!file.exists(assoc_file))
    stop("GEMMA wrote no results in ", dir, ":\n",
         paste(utils::tail(log, 10L), collapse = "\n"), call. = FALSE)
  assoc <- utils::read.table(assoc_file, header = TRUE, stringsAsFactors = FALSE)
  list(assoc = assoc,
       null = .indri_gemma_null(f("output/result.log.txt"), colnames(Y)))
}

# GEMMA's REML Vg and Ve under the null model, from its log.
.indri_gemma_null <- function(logfile, traits) {
  if (!file.exists(logfile)) return(NULL)
  l <- readLines(logfile)
  P <- length(traits)
  if (P == 1L) {
    grab <- function(tag) {
      i <- grep(tag, l, fixed = TRUE)
      if (!length(i)) return(NA_real_)
      as.numeric(sub(".*=", "", l[i[1]]))
    }
    vg <- grab("vg estimate in the null model")
    ve <- grab("ve estimate in the null model")
    return(list(Vg = matrix(vg, 1, 1, dimnames = list(traits, traits)),
                Ve = matrix(ve, 1, 1, dimnames = list(traits, traits))))
  }
  grab <- function(tag) {
    i <- grep(tag, l, fixed = TRUE)
    if (!length(i)) return(NULL)
    v <- as.numeric(unlist(strsplit(trimws(sub("^##", "", l[i[1] + seq_len(P)])), "\\s+")))
    m <- matrix(0, P, P)
    m[upper.tri(m, diag = TRUE)] <- v[seq_len(P * (P + 1) / 2)]
    m <- t(m)
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    dimnames(m) <- list(traits, traits)
    m
  }
  list(Vg = grab("REMLE estimate for Vg in the null model"),
       Ve = grab("REMLE estimate for Ve in the null model"))
}
