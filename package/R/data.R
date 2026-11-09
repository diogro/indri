#' LG/SM F3 intercross: phenotypes
#'
#' Phenotypes of 1611 mice from the F3 generation of an intercross of the
#' LG/J and SM/J inbred strains, with their pedigree and rearing. The traits
#' fall in four classes, listed in [lgsm_traits]: growth (weekly weight
#' gains), weight (the raw weekly body weights), area (mandible areas) and
#' necropsy (organ and fat pad weights).
#' Genotypes are in [lgsm_genotypes], the markers in [lgsm_markers], and the
#' relationship matrix in [lgsm_A].
#'
#' The cross: the parental (P) strains are fully inbred, the F1 are all
#' heterozygous, and the F2 parents of these F3 mice are unrelated. The F2
#' are the founders of the pedigree ([lgsm_pedigree]).
#'
#' @format A data frame with one row per mouse, in this column order:
#' \describe{
#'   \item{ID}{mouse identifier; matches the row names of the [lgsm_genotypes] matrices and of [lgsm_A] (as character).}
#'   \item{FAMILY}{full-sib family.}
#'   \item{Sire}{father (an F2 mouse, genotyped in [lgsm_genotypes]); four missing in the source were inferred from the markers, see `sire_inferred`.}
#'   \item{Dam}{mother (an F2 mouse).}
#'   \item{NURSE}{the dam that reared the mouse (pups were cross-fostered, see `xfostpair`).}
#'   \item{SEX}{`"F"` or `"M"`.}
#'   \item{LSB}{litter size at birth.}
#'   \item{LSW}{litter size at weaning.}
#'   \item{COHORT}{rearing cohort.}
#'   \item{xfostpair}{cross-fostering indicator, as in the source.}
#'   \item{growth12}{weight gain from week 1 to week 2, `WEEK2 - WEEK1` (growth).}
#'   \item{growth23}{weight gain from week 2 to week 3, `WEEK3 - WEEK2` (growth).}
#'   \item{growth34}{weight gain from week 3 to week 4, `WEEK4 - WEEK3` (growth).}
#'   \item{growth45}{weight gain from week 4 to week 5, `WEEK5 - WEEK4` (growth).}
#'   \item{growth56}{weight gain from week 5 to week 6, `WEEK6 - WEEK5` (growth).}
#'   \item{growth67}{weight gain from week 6 to week 7, `WEEK7 - WEEK6` (growth).}
#'   \item{growth78}{weight gain from week 7 to week 8, `WEEK8 - WEEK7` (growth).}
#'   \item{growth89}{weight gain from week 8 to week 9, `WEEK9 - WEEK8` (growth).}
#'   \item{growth910}{weight gain from week 9 to week 10, `WEEK10 - WEEK9` (growth).}
#'   \item{WEEK1}{body weight at week 1 (weight).}
#'   \item{WEEK2}{body weight at week 2 (weight).}
#'   \item{WEEK3}{body weight at week 3 (weight).}
#'   \item{WEEK4}{body weight at week 4 (weight).}
#'   \item{WEEK5}{body weight at week 5 (weight).}
#'   \item{WEEK6}{body weight at week 6 (weight).}
#'   \item{WEEK7}{body weight at week 7 (weight).}
#'   \item{WEEK8}{body weight at week 8 (weight).}
#'   \item{WEEK9}{body weight at week 9 (weight).}
#'   \item{WEEK10}{body weight at week 10 (weight).}
#'   \item{FATPAD}{fat pad weight at necropsy (necropsy); a few missing.}
#'   \item{HEART}{heart weight at necropsy (necropsy); a few missing.}
#'   \item{KIDNEY}{kidney weight at necropsy (necropsy); a few missing.}
#'   \item{SPLEEN}{spleen weight at necropsy (necropsy); a few missing.}
#'   \item{LIVER}{liver weight at necropsy (necropsy); a few missing.}
#'   \item{area1}{mandible area 1 (area); `NA` for the 95 mice without areas.}
#'   \item{area2}{mandible area 2 (area); `NA` for the 95 mice without areas.}
#'   \item{area3}{mandible area 3 (area); `NA` for the 95 mice without areas.}
#'   \item{area4}{mandible area 4 (area); `NA` for the 95 mice without areas.}
#'   \item{area5}{mandible area 5 (area); `NA` for the 95 mice without areas.}
#'   \item{area6}{mandible area 6 (area); `NA` for the 95 mice without areas.}
#'   \item{area7}{mandible area 7 (area); `NA` for the 95 mice without areas.}
#'   \item{sire_inferred}{`TRUE` for the four mice whose sire is missing in the source and was found by Mendelian exclusion: 1133, 1139 and 1140 (family 199, dam 135) are sired by 165, and 1109 (dam 207) by 229. Each inferred sire is consistent with the offspring and dam at every marker; the next best candidate fails at 12 or more (on mice with a recorded sire, the recorded one always had none and the next best 12 or more).}
#' }
#' Values are as in the source files, apart from the gains computed from
#' them: no transformation, no correction for the rearing covariates.
#'
#' @source Melo, D., Marroig, G. and Wolf, J. B. (2019). Genomic perspective
#'   on multivariate variation, pleiotropy, and evolution. Journal of
#'   Heredity 110, 479-493. \doi{10.1093/jhered/esz011}.
#'   Built by `data-raw/lgsm.R` from that study's analysis files.
#' @seealso [lgsm_traits], [lgsm_genotypes], [lgsm_markers]
"lgsm_f3"

#' LG/SM F3 intercross: trait classes
#'
#' The columns of [lgsm_f3] in each trait class.
#'
#' @format A named list of character vectors: `growth` (`growth12` to
#'   `growth910`, the nine weekly gains; the original QTL analyses used the
#'   first seven), `weight` (`WEEK1` to `WEEK10`, raw weekly
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
#' @source Melo, D., Marroig, G. and Wolf, J. B. (2019). Genomic perspective
#'   on multivariate variation, pleiotropy, and evolution. Journal of
#'   Heredity 110, 479-493. \doi{10.1093/jhered/esz011}.
#'   Built by `data-raw/lgsm.R` from that study's analysis files.
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
#' @source Melo, D., Marroig, G. and Wolf, J. B. (2019). Genomic perspective
#'   on multivariate variation, pleiotropy, and evolution. Journal of
#'   Heredity 110, 479-493. \doi{10.1093/jhered/esz011}.
#'   Built by `data-raw/lgsm.R` from that study's analysis files.
"lgsm_markers"

#' LG/SM F3 intercross: pedigree
#'
#' The F3 mice of [lgsm_f3] and their F2 parents. The P strains are fully
#' inbred and the F1 all heterozygous, so the F2 are unrelated founders and
#' have no recorded parents.
#'
#' @format A data frame: `ID`, `Sire`, `Dam` (`NA` for the 391 F2
#'   founders), `generation` (`"F2"` or `"F3"`) and `sire_inferred` (the four
#'   sires found from the markers; see [lgsm_f3]).
#' @source Built by `data-raw/lgsm.R` from [lgsm_f3].
#' @seealso [lgsm_A]
"lgsm_pedigree"

#' LG/SM F3 intercross: relationship matrix
#'
#' The additive relationship matrix of the 1611 F3 mice of [lgsm_f3], from
#' [lgsm_pedigree] with the F2 parents as unrelated, non-inbred founders:
#' 1 on the diagonal, 1/2 between full sibs, 1/4 between half sibs (sires
#' with more than one mate) and 0 otherwise. The four sires missing in the
#' source come from the markers (`lgsm_f3$sire_inferred`).
#'
#' @format A 1611 x 1611 numeric matrix with the F3 IDs (as character) as
#'   row and column names, in the order of `lgsm_f3$ID`. It is positive
#'   definite (smallest eigenvalue 1/2). Take the rows and columns of the
#'   mice in an analysis, as [fit_G()] and [indri_gwas()] match by name.
#' @source Built by `data-raw/lgsm.R` from [lgsm_pedigree].
#' @examples
#' table(lgsm_A[upper.tri(lgsm_A)])
"lgsm_A"
