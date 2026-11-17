# Summarise an association scan

The SNPs with the smallest p-values and the genomic inflation factor
`lambda_gc` (the median test statistic over its expectation under the
null; near 1 for a calibrated scan of mostly null SNPs).

## Usage

``` r
# S3 method for class 'indri_gwas'
summary(object, n = 10L, ...)
```

## Arguments

- object:

  An `indri_gwas` scan.

- n:

  Number of top SNPs (or SNP-trait pairs) to list.

- ...:

  Unused.

## Value

An object of class `summary.indri_gwas` with `top` (a data frame),
`lambda_gc` (one value, or one per trait for single-trait scans),
`n_tested` and the scan's description.
