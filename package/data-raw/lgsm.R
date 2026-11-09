# Build the LG/SM F3 data sets shipped in data/ from the mouseGrowthQTLs
# repository (a sibling checkout of indri). Run from the indri repository
# root:
#   Rscript package/data-raw/lgsm.R [path/to/mouseGrowthQTLs]
#
# Writes package/data/lgsm_f3.rda, lgsm_genotypes.rda, lgsm_markers.rda,
# lgsm_traits.rda, lgsm_pedigree.rda and lgsm_A.rda. Values are copied from the source files unchanged; only
# names are tidied (SexAN -> SEX) and rows put in ID order. The growth traits
# are the weekly gains of the original analyses (read_mouse_data.R), computed
# from the raw weekly weights, which are shipped too.

args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[1] else "../mouseGrowthQTLs"
data_dir <- file.path(src, "data")
stopifnot(dir.exists(data_dir))

## Phenotypes ---------------------------------------------------------------

raw <- utils::read.csv(file.path(data_dir, "growth traits",
  "F3Phenotypes_further corrected family data_corrected litter sizes.csv"),
  check.names = FALSE, stringsAsFactors = FALSE)
areas <- utils::read.csv(file.path(data_dir, "area traits", "areasF2F3.csv"))

weight <- paste0("WEEK", 1:10)
# every weekly gain the ten weights give (the original analyses used the
# first seven, growth12 to growth78)
growth <- paste0("growth", 1:9, 2:10)         # growth12, ..., growth910
necropsy <- c("FATPAD", "HEART", "KIDNEY", "SPLEEN", "LIVER")
area <- paste0("area", 1:7)
design <- c("ID", "FAMILY", "Sire", "Dam", "NURSE", "SexAN", "LSB", "LSW",
            "COHORT", "xfostpair")

lgsm_f3 <- raw[, c(design, weight, necropsy)]
for (k in 1:9)
  lgsm_f3[[growth[k]]] <- lgsm_f3[[weight[k + 1]]] - lgsm_f3[[weight[k]]]
lgsm_f3 <- lgsm_f3[, c(design, growth, weight, necropsy)]
names(lgsm_f3)[names(lgsm_f3) == "SexAN"] <- "SEX"
lgsm_f3 <- merge(lgsm_f3, areas[, c("ID", area)], by = "ID", all.x = TRUE)
lgsm_f3 <- lgsm_f3[order(lgsm_f3$ID), ]
rownames(lgsm_f3) <- NULL
stopifnot(!anyDuplicated(lgsm_f3$ID), nrow(lgsm_f3) == nrow(raw))

lgsm_traits <- list(growth = growth, weight = weight, area = area,
                    necropsy = necropsy)

## Markers -----------------------------------------------------------------

positions <- utils::read.csv(file.path(data_dir, "markers",
                                       "marker_positions.csv"),
                             check.names = FALSE, stringsAsFactors = FALSE)
lgsm_markers <- data.frame(
  snp = positions$Marker,
  chrom = as.integer(positions$chr),
  cM = as.numeric(positions[["Map Pos (cM)"]]),
  # one marker (rs3684333, chromosome 3) has "?" for its physical position
  bp = as.numeric(ifelse(positions[["Phys Pos (bp)"]] == "?", NA,
                         gsub(",", "", positions[["Phys Pos (bp)"]]))),
  stringsAsFactors = FALSE)

chroms <- lapply(1:19, function(chr)
  utils::read.csv(file.path(data_dir, "markers", sprintf("chrom%d.csv", chr))))
ids <- sort(Reduce(intersect, lapply(chroms, `[[`, "ID")))
coding <- function(prefix) {
  m <- do.call(cbind, lapply(chroms, function(x) {
    cols <- grep(paste0("^", prefix, "[0-9]+$"), names(x), value = TRUE)
    cols <- cols[order(as.integer(sub(prefix, "", cols)))]
    as.matrix(x[match(ids, x$ID), cols])
  }))
  storage.mode(m) <- "integer"
  dimnames(m) <- list(as.character(ids), lgsm_markers$snp)
  m
}
# The source files hold, per locus k, columns Ak, Dk and Ik in map order; the
# loci per chromosome match the rows of marker_positions.csv.
stopifnot(sapply(chroms, function(x) sum(grepl("^A[0-9]+$", names(x)))) ==
            as.vector(table(lgsm_markers$chrom)))
lgsm_genotypes <- list(A = coding("A"), D = coding("D"), I = coding("I"))
stopifnot(all(lgsm_genotypes$D == (lgsm_genotypes$A == 0L)))

## Missing sires, from the markers -----------------------------------------

# Four F3 mice have no recorded sire: 1133, 1139 and 1140 (family 199, dam
# 135) and 1109 (dam 207). Their sires were found by Mendelian exclusion:
# among the genotyped sires, the one with no locus where no allele of the
# dam plus one of the sire gives the offspring's genotype. On 25 mice with a
# known sire the true sire had 0 exclusions in every case and the next best
# candidate at least 12; here the inferred sire has 0 and the next best 12
# (family 199) and 23 (1109). Sire 165 also sired families 46 and 200, and
# sire 229 family 71.
dosage <- lgsm_genotypes$A + 1L
exclusions <- function(off, dam, sire) {
  g <- dosage[as.character(c(off, dam, sire)), ]
  from <- list(0L, 0:1, 1L)
  sum(!vapply(seq_len(ncol(g)), function(l)
    any(outer(from[[g[2, l] + 1L]], from[[g[3, l] + 1L]], `+`) == g[1, l]),
    logical(1)))
}
inferred <- c(`1133` = 165L, `1139` = 165L, `1140` = 165L, `1109` = 229L)
rows <- match(as.integer(names(inferred)), lgsm_f3$ID)
stopifnot(all(is.na(lgsm_f3$Sire[rows])),
          sum(is.na(lgsm_f3$Sire)) == length(inferred))
for (i in seq_along(rows))
  stopifnot(exclusions(lgsm_f3$ID[rows[i]], lgsm_f3$Dam[rows[i]],
                       inferred[[i]]) == 0L)
lgsm_f3$sire_inferred <- FALSE
lgsm_f3$Sire[rows] <- unname(inferred)
lgsm_f3$sire_inferred[rows] <- TRUE

## Pedigree and relationship matrix ----------------------------------------

# The cross: P are the inbred LG/J and SM/J strains, the F1 are all
# heterozygous (so identical), and the F2 parents of the F3 are unrelated
# and not inbred. The F2 are therefore the founders: F3 mice are related
# only through shared F2 parents.
parents <- sort(unique(c(lgsm_f3$Sire, lgsm_f3$Dam)))
stopifnot(!length(intersect(lgsm_f3$Sire, lgsm_f3$Dam)),
          !any(lgsm_f3$ID %in% parents), !anyNA(lgsm_f3$Sire))
lgsm_pedigree <- rbind(
  data.frame(ID = parents, Sire = NA_integer_, Dam = NA_integer_,
             generation = "F2", stringsAsFactors = FALSE),
  data.frame(ID = lgsm_f3$ID, Sire = lgsm_f3$Sire, Dam = lgsm_f3$Dam,
             generation = "F3", stringsAsFactors = FALSE))
lgsm_pedigree$sire_inferred <- c(rep(FALSE, length(parents)),
                                 lgsm_f3$sire_inferred)

sire <- as.character(lgsm_f3$Sire)
dam <- as.character(lgsm_f3$Dam)
shared <- function(a, b) outer(a, b, "==")
# a quarter of the parent pairs the two mice share (founders unrelated)
lgsm_A <- (shared(sire, sire) + shared(dam, dam) + shared(sire, dam) +
             shared(dam, sire)) / 4
diag(lgsm_A) <- 1
dimnames(lgsm_A) <- list(as.character(lgsm_f3$ID), as.character(lgsm_f3$ID))

dir.create("package/data", showWarnings = FALSE)
save(lgsm_pedigree, file = "package/data/lgsm_pedigree.rda", compress = "xz")
save(lgsm_A, file = "package/data/lgsm_A.rda", compress = "xz")
save(lgsm_f3, file = "package/data/lgsm_f3.rda", compress = "xz")
save(lgsm_traits, file = "package/data/lgsm_traits.rda", compress = "xz")
save(lgsm_markers, file = "package/data/lgsm_markers.rda", compress = "xz")
save(lgsm_genotypes, file = "package/data/lgsm_genotypes.rda", compress = "xz")
message(sprintf("lgsm_f3: %d mice; lgsm_genotypes: %d mice x %d loci",
                nrow(lgsm_f3), nrow(lgsm_genotypes$A), ncol(lgsm_genotypes$A)))
