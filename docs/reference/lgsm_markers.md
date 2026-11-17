# LG/SM F3 intercross: marker map

The 353 markers of
[lgsm_genotypes](https://diogro.github.io/indri/reference/lgsm_genotypes.md),
in column order.

## Usage

``` r
lgsm_markers
```

## Format

A data frame: `snp` (marker name), `chrom` (1 to 19), `cM` (genetic map
position) and `bp` (physical position; `NA` for rs3684333, unknown in
the source). Pass it as `snp_info` to
[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)
(with `pos = bp` or `pos = cM`) to place the markers on a Manhattan
plot.

## Source

Melo, D., Marroig, G. and Wolf, J. B. (2019). Genomic perspective on
multivariate variation, pleiotropy, and evolution. Journal of Heredity
110, 479-493.
[doi:10.1093/jhered/esz011](https://doi.org/10.1093/jhered/esz011) .
Built by `data-raw/lgsm.R` from that study's analysis files.
