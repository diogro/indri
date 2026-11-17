# LG/SM F3 intercross: pedigree

The F3 mice of
[lgsm_f3](https://diogro.github.io/indri/reference/lgsm_f3.md) and their
F2 parents. The P strains are fully inbred and the F1 all heterozygous,
so the F2 are unrelated founders and have no recorded parents.

## Usage

``` r
lgsm_pedigree
```

## Format

A data frame: `ID`, `Sire`, `Dam` (`NA` for the 391 F2 founders),
`generation` (`"F2"` or `"F3"`) and `sire_inferred` (the four sires
found from the markers; see
[lgsm_f3](https://diogro.github.io/indri/reference/lgsm_f3.md)).

## Source

Built by `data-raw/lgsm.R` from
[lgsm_f3](https://diogro.github.io/indri/reference/lgsm_f3.md).

## See also

[lgsm_A](https://diogro.github.io/indri/reference/lgsm_A.md)
