# Package index

## Estimating G

Split the traits into overlapping subsets, fit each, reassemble.

- [`fit_G()`](https://diogro.github.io/indri/reference/fit_G.md) :
  Estimate the additive genetic covariance matrix
- [`fit_G_env()`](https://diogro.github.io/indri/reference/fit_G_env.md)
  : Environment-specific genetic covariance (experimental)
- [`indri_engine()`](https://diogro.github.io/indri/reference/indri_engine.md)
  [`indri_engines()`](https://diogro.github.io/indri/reference/indri_engine.md)
  : Block engines
- [`indri_locus_G()`](https://diogro.github.io/indri/reference/indri_locus_G.md)
  : Genetic covariance matrices from locus effects

## Genome scans

- [`indri_gwas()`](https://diogro.github.io/indri/reference/indri_gwas.md)
  : Genome-wide association scan with the assembled covariances
- [`summary(`*`<indri_gwas>`*`)`](https://diogro.github.io/indri/reference/summary.indri_gwas.md)
  : Summarise an association scan
- [`plot(`*`<indri_gwas>`*`)`](https://diogro.github.io/indri/reference/plot.indri_gwas.md)
  : Plot an association scan
- [`indri_effective_tests()`](https://diogro.github.io/indri/reference/indri_effective_tests.md)
  : Effective number of independent traits

## Building blocks

The steps fit_G() runs, for custom pipelines.

- [`indri_subsets()`](https://diogro.github.io/indri/reference/indri_subsets.md)
  : Trait subsets for the block estimator
- [`indri_a_forms()`](https://diogro.github.io/indri/reference/indri_a_forms.md)
  : Relationship matrices in the forms the engines want
- [`indri_inflate()`](https://diogro.github.io/indri/reference/indri_inflate.md)
  : Place a subset covariance matrix into the full trait space
- [`indri_assemble()`](https://diogro.github.io/indri/reference/indri_assemble.md)
  : Average the subset estimates into one matrix
- [`indri_assemble_draws()`](https://diogro.github.io/indri/reference/indri_assemble_draws.md)
  : Assemble block draws into one array
- [`indri_project_psd()`](https://diogro.github.io/indri/reference/indri_project_psd.md)
  : Fit the assembled blocks to the closest positive semi-definite
  matrix

## Data

The LG/SM F3 mouse intercross.

- [`lgsm_A`](https://diogro.github.io/indri/reference/lgsm_A.md) : LG/SM
  F3 intercross: relationship matrix
- [`lgsm_f3`](https://diogro.github.io/indri/reference/lgsm_f3.md) :
  LG/SM F3 intercross: phenotypes
- [`lgsm_genotypes`](https://diogro.github.io/indri/reference/lgsm_genotypes.md)
  : LG/SM F3 intercross: genotypes
- [`lgsm_markers`](https://diogro.github.io/indri/reference/lgsm_markers.md)
  : LG/SM F3 intercross: marker map
- [`lgsm_pedigree`](https://diogro.github.io/indri/reference/lgsm_pedigree.md)
  : LG/SM F3 intercross: pedigree
- [`lgsm_traits`](https://diogro.github.io/indri/reference/lgsm_traits.md)
  : LG/SM F3 intercross: trait classes

## Package

- [`indri`](https://diogro.github.io/indri/reference/indri-package.md)
  [`indri-package`](https://diogro.github.io/indri/reference/indri-package.md)
  : indri: scalable additive genetic covariance estimation
