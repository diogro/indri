# LG/SM F3 intercross: phenotypes

Phenotypes of 1611 mice from the F3 generation of an intercross of the
LG/J and SM/J inbred strains, with their pedigree and rearing. The
traits fall in four classes, listed in
[lgsm_traits](https://diogro.github.io/indri/reference/lgsm_traits.md):
growth (weekly weight gains), weight (the raw weekly body weights), area
(mandible areas) and necropsy (organ and fat pad weights). Genotypes are
in
[lgsm_genotypes](https://diogro.github.io/indri/reference/lgsm_genotypes.md),
the markers in
[lgsm_markers](https://diogro.github.io/indri/reference/lgsm_markers.md),
and the relationship matrix in
[lgsm_A](https://diogro.github.io/indri/reference/lgsm_A.md).

## Usage

``` r
lgsm_f3
```

## Format

A data frame with one row per mouse, in this column order:

- ID:

  mouse identifier; matches the row names of the
  [lgsm_genotypes](https://diogro.github.io/indri/reference/lgsm_genotypes.md)
  matrices and of
  [lgsm_A](https://diogro.github.io/indri/reference/lgsm_A.md) (as
  character).

- FAMILY:

  full-sib family.

- Sire:

  father (an F2 mouse, genotyped in
  [lgsm_genotypes](https://diogro.github.io/indri/reference/lgsm_genotypes.md));
  four missing in the source were inferred from the markers, see
  `sire_inferred`.

- Dam:

  mother (an F2 mouse).

- NURSE:

  the dam that reared the mouse (pups were cross-fostered, see
  `xfostpair`).

- SEX:

  `"F"` or `"M"`.

- LSB:

  litter size at birth.

- LSW:

  litter size at weaning.

- COHORT:

  rearing cohort.

- xfostpair:

  cross-fostering indicator, as in the source.

- growth12:

  weight gain from week 1 to week 2, `WEEK2 - WEEK1` (growth).

- growth23:

  weight gain from week 2 to week 3, `WEEK3 - WEEK2` (growth).

- growth34:

  weight gain from week 3 to week 4, `WEEK4 - WEEK3` (growth).

- growth45:

  weight gain from week 4 to week 5, `WEEK5 - WEEK4` (growth).

- growth56:

  weight gain from week 5 to week 6, `WEEK6 - WEEK5` (growth).

- growth67:

  weight gain from week 6 to week 7, `WEEK7 - WEEK6` (growth).

- growth78:

  weight gain from week 7 to week 8, `WEEK8 - WEEK7` (growth).

- growth89:

  weight gain from week 8 to week 9, `WEEK9 - WEEK8` (growth).

- growth910:

  weight gain from week 9 to week 10, `WEEK10 - WEEK9` (growth).

- WEEK1:

  body weight at week 1 (weight).

- WEEK2:

  body weight at week 2 (weight).

- WEEK3:

  body weight at week 3 (weight).

- WEEK4:

  body weight at week 4 (weight).

- WEEK5:

  body weight at week 5 (weight).

- WEEK6:

  body weight at week 6 (weight).

- WEEK7:

  body weight at week 7 (weight).

- WEEK8:

  body weight at week 8 (weight).

- WEEK9:

  body weight at week 9 (weight).

- WEEK10:

  body weight at week 10 (weight).

- FATPAD:

  fat pad weight at necropsy (necropsy); a few missing.

- HEART:

  heart weight at necropsy (necropsy); a few missing.

- KIDNEY:

  kidney weight at necropsy (necropsy); a few missing.

- SPLEEN:

  spleen weight at necropsy (necropsy); a few missing.

- LIVER:

  liver weight at necropsy (necropsy); a few missing.

- area1:

  mandible area 1 (area); `NA` for the 95 mice without areas.

- area2:

  mandible area 2 (area); `NA` for the 95 mice without areas.

- area3:

  mandible area 3 (area); `NA` for the 95 mice without areas.

- area4:

  mandible area 4 (area); `NA` for the 95 mice without areas.

- area5:

  mandible area 5 (area); `NA` for the 95 mice without areas.

- area6:

  mandible area 6 (area); `NA` for the 95 mice without areas.

- area7:

  mandible area 7 (area); `NA` for the 95 mice without areas.

- sire_inferred:

  `TRUE` for the four mice whose sire is missing in the source and was
  found by Mendelian exclusion: 1133, 1139 and 1140 (family 199,
  dam 135) are sired by 165, and 1109 (dam 207) by 229. Each inferred
  sire is consistent with the offspring and dam at every marker; the
  next best candidate fails at 12 or more (on mice with a recorded sire,
  the recorded one always had none and the next best 12 or more).

Values are as in the source files, apart from the gains computed from
them: no transformation, no correction for the rearing covariates.

## Source

Melo, D., Marroig, G. and Wolf, J. B. (2019). Genomic perspective on
multivariate variation, pleiotropy, and evolution. Journal of Heredity
110, 479-493.
[doi:10.1093/jhered/esz011](https://doi.org/10.1093/jhered/esz011) .
Built by `data-raw/lgsm.R` from that study's analysis files.

## Details

The cross: the parental (P) strains are fully inbred, the F1 are all
heterozygous, and the F2 parents of these F3 mice are unrelated. The F2
are the founders of the pedigree
([lgsm_pedigree](https://diogro.github.io/indri/reference/lgsm_pedigree.md)).

## See also

[lgsm_traits](https://diogro.github.io/indri/reference/lgsm_traits.md),
[lgsm_genotypes](https://diogro.github.io/indri/reference/lgsm_genotypes.md),
[lgsm_markers](https://diogro.github.io/indri/reference/lgsm_markers.md)
