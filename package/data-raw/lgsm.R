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

## Pedigree and relationship matrix ----------------------------------------

# The cross: P are the inbred LG/J and SM/J strains, the F1 are all
# heterozygous (so identical), and the F2 parents of the F3 are unrelated
# and not inbred. The F2 are therefore the founders: F3 mice are related
# only through shared F2 parents. Four F3 mice have no recorded sire (three
# in family 199, and 1109); every dam mated once, so offspring of a dam are
# taken to share its unknown mate.
parents <- sort(unique(c(lgsm_f3$Sire, lgsm_f3$Dam)))
stopifnot(!length(intersect(lgsm_f3$Sire, lgsm_f3$Dam)),
          !any(lgsm_f3$ID %in% parents),
          all(table(unique(lgsm_f3[c("Sire", "Dam")])$Dam) == 1L))
lgsm_pedigree <- rbind(
  data.frame(ID = parents, Sire = NA_integer_, Dam = NA_integer_,
             generation = "F2", stringsAsFactors = FALSE),
  data.frame(ID = lgsm_f3$ID, Sire = lgsm_f3$Sire, Dam = lgsm_f3$Dam,
             generation = "F3", stringsAsFactors = FALSE))

sire <- ifelse(is.na(lgsm_f3$Sire), paste0("mate_of_", lgsm_f3$Dam),
               as.character(lgsm_f3$Sire))
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
