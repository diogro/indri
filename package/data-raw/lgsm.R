# Build the LG/SM F3 data sets shipped in data/ from the mouseGrowthQTLs
# repository (a sibling checkout of indri). Run from the indri repository
# root:
#   Rscript package/data-raw/lgsm.R [path/to/mouseGrowthQTLs]
#
# Writes package/data/lgsm_f3.rda, lgsm_genotypes.rda, lgsm_markers.rda and
# lgsm_traits.rda. Values are copied from the source files unchanged; only
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
growth <- paste0("growth", 1:7, 2:8)          # growth12, ..., growth78
necropsy <- c("FATPAD", "HEART", "KIDNEY", "SPLEEN", "LIVER")
area <- paste0("area", 1:7)
design <- c("ID", "FAMILY", "Sire", "Dam", "NURSE", "SexAN", "LSB", "LSW",
            "COHORT", "xfostpair")

lgsm_f3 <- raw[, c(design, weight, necropsy)]
for (k in 1:7)
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

dir.create("package/data", showWarnings = FALSE)
save(lgsm_f3, file = "package/data/lgsm_f3.rda", compress = "xz")
save(lgsm_traits, file = "package/data/lgsm_traits.rda", compress = "xz")
save(lgsm_markers, file = "package/data/lgsm_markers.rda", compress = "xz")
save(lgsm_genotypes, file = "package/data/lgsm_genotypes.rda", compress = "xz")
message(sprintf("lgsm_f3: %d mice; lgsm_genotypes: %d mice x %d loci",
                nrow(lgsm_f3), nrow(lgsm_genotypes$A), ncol(lgsm_genotypes$A)))
