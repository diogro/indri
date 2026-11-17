# Plot an association scan

A Manhattan plot (`-log10(p)` along the genome) or a QQ plot (observed
against expected `-log10(p)` under the null), drawn with ggplot2 and
returned as a `ggplot` object to print, save or extend.

## Usage

``` r
# S3 method for class 'indri_gwas'
plot(
  x,
  type = c("manhattan", "qq"),
  which = c("joint", "dominance", "gxe", "empirical"),
  trait = NULL,
  threshold = NULL,
  highlight = NULL,
  ...
)
```

## Arguments

- x:

  An `indri_gwas` scan.

- type:

  `"manhattan"` or `"qq"`.

- which:

  The p-values to draw: `"joint"` (`pvalues`, the default),
  `"dominance"` (`pvalues_dominance`), `"gxe"` (`pvalues_gxe`) or
  `"empirical"` (`pvalues_empirical`).

- trait:

  For single-trait scans, the trait to draw; `NULL` (the default) draws
  each SNP's smallest p-value over traits, adjusted when the scan was
  (`pvalues_adjust`). The QQ plot of a single-trait scan uses every SNP
  and trait when `trait` is `NULL`.

- threshold:

  Significance line on the Manhattan plot, as a p-value; defaults to
  Bonferroni over the SNPs drawn (`0.05 / M`). `NA` for none.

- highlight:

  Names of SNPs to mark and label.

- ...:

  Unused.

## Value

A `ggplot` object.

## Details

The Manhattan plot places SNPs by `chrom` and `pos` from `snp_info` when
both are there, chromosomes end to end; otherwise in column order. SNPs
with `NA` p-values (monomorphic, or without a dominance estimate) are
left out.

## See also

[`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)
