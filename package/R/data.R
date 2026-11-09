#' LG/SM F3 intercross: phenotypes
#'
#' Phenotypes of 1611 mice from the F3 generation of an intercross of the
#' LG/J and SM/J inbred strains, with their pedigree and rearing. The traits
#' fall in four classes, listed in [lgsm_traits]: growth (weekly weight
#' gains), weight (the raw weekly body weights), area (skull areas) and
#' necropsy (organ and fat pad weights).
#' Genotypes are in [lgsm_genotypes], and the markers in [lgsm_markers].
#'
#' @format A data frame with one row per mouse:
#' \describe{
#'   \item{ID}{mouse identifier; matches the row names of the
#'     [lgsm_genotypes] matrices (as character).}
#'   \item{FAMILY}{full-sib family.}
#'   \item{Sire, Dam}{parents (F2 mice, genotyped in [lgsm_genotypes]).}
#'   \item{NURSE}{the dam that reared the mouse (pups were cross-fostered,
#'     see `xfostpair`).}
#'   \item{SEX}{`"F"` or `"M"`.}
#'   \item{LSB, LSW}{litter size at birth and at weaning.}
#'   \item{COHORT}{rearing cohort.}
#'   \item{xfostpair}{cross-fostering indicator, as in the source.}
#'   \item{growth12, ..., growth78}{weekly weight gains, `growth12 = WEEK2 -
#'     WEEK1` to `growth78 = WEEK8 - WEEK7`: the growth traits of the
#'     original QTL analyses (growth).}
#'   \item{WEEK1, ..., WEEK10}{raw body weight at weeks 1 to 10 (weight).}
#'   \item{FATPAD, HEART, KIDNEY, SPLEEN, LIVER}{necropsy weights; a few
#'     missing.}
#'   \item{area1, ..., area7}{skull areas, `NA` for the 95 mice without
#'     them.}
#' }
#' Values are as in the source files, apart from the gains computed from
#' them: no transformation, no correction for the rearing covariates.
#'
#' @source Provenance to be documented. Built by `data-raw/lgsm.R` from the
#'   mouseGrowthQTLs analysis repository.
#' @seealso [lgsm_traits], [lgsm_genotypes], [lgsm_markers]
"lgsm_f3"

#' LG/SM F3 intercross: trait classes
#'
#' The columns of [lgsm_f3] in each trait class.
#'
#' @format A named list of character vectors: `growth` (`growth12` to
#'   `growth78`, weekly gains), `weight` (`WEEK1` to `WEEK10`, raw weekly
#'   weights), `area` (`area1` to `area7`) and `necropsy` (`FATPAD`,
#'   `HEART`, `KIDNEY`, `SPLEEN`, `LIVER`).
#' @source Built by `data-raw/lgsm.R`.
"lgsm_traits"

#' LG/SM F3 intercross: genotypes
#'
#' Genotypes of the LG/SM intercross mice at 353 markers, in the additive,
#' dominance and imprinting design of the original QTL analyses.
#'
#' @format A list of three integer matrices, mice by markers, with mouse IDs
#'   as row names and marker names (see [lgsm_markers]) as column names:
#' \describe{
#'   \item{A}{additive: -1 and 1 for the two homozygotes, 0 for
#'     heterozygotes. `A + 1` is a 0/1/2 dosage, as [indri_gwas()] expects.}
#'   \item{D}{dominance: 1 for heterozygotes, 0 for homozygotes (equal to
#'     `A == 0`).}
#'   \item{I}{imprinting: -1 or 1 for heterozygotes by the parent of origin
#'     of the alleles, 0 for homozygotes.}
#' }
#' The 1934 mice are the 1552 genotyped F3 mice of [lgsm_f3] and 382 of their
#' F2 parents. Which strain's allele is coded `+1` is as in the source.
#'
#' @source Provenance to be documented. Built by `data-raw/lgsm.R` from the
#'   mouseGrowthQTLs analysis repository.
#' @examples
#' dosage <- lgsm_genotypes$A + 1
#' f3 <- rownames(dosage) %in% lgsm_f3$ID
#' dim(dosage[f3, ])
"lgsm_genotypes"

#' LG/SM F3 intercross: marker map
#'
#' The 353 markers of [lgsm_genotypes], in column order.
#'
#' @format A data frame: `snp` (marker name), `chrom` (1 to 19), `cM`
#'   (genetic map position) and `bp` (physical position; `NA` for
#'   rs3684333, unknown in the source). Pass it as `snp_info` to
#'   [indri_gwas()] (with `pos = bp` or `pos = cM`) to place the markers on
#'   a Manhattan plot.
#' @source Provenance to be documented. Built by `data-raw/lgsm.R`.
"lgsm_markers"
