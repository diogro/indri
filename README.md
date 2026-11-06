# indri

Scalable estimation of the additive genetic covariance matrix **G** for
high-dimensional traits, by combining trait subsampling (LiMMBo) with small
per-subset fits from a pluggable engine — Stan via stanimal, REML via sommer, or
REML via ASReml-R.

Status: repo scaffold. The design is settled enough to code from, but no
implementation has landed yet — see [NOTES.md](NOTES.md) for the estimator, the
engines, the architecture, the Stan delta and the milestones.

## Scope

In scope: `G` (traits × traits additive genetic covariance) from a phenotype
matrix `Y`, a design `X` and a known relationship matrix `A`, using one of
`stanimal` (Bayesian), `sommer` (REML) or `asreml` (REML, sparse `A`) as the
block engine.

Out of scope for now: `E`/`P` reporting, association testing, breeding values,
missing-data imputation, model comparison.
