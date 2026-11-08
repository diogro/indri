# limmbo reference fixtures

Association results from limmbo (LIMIX 1.x) on limmbo's own test data
(`limmbo/io/test/data/`), restricted to the first 250 samples, used by
`test-gwas-limmbo.R`.

Regenerate from a limmbo checkout with its pixi environment:

```sh
cd limmbo
PYTHONPATH=. pixi run python /path/to/make_limmbo_fixtures.py <outdir>
```

then convert `K.csv`, the relationship matrix, to `A.rds`, the name indri
uses for it (`saveRDS(A, "A.rds", compress = "xz")`).
`log_delta.csv` holds the null variance ratio LIMIX fitted for each test;
the tests pass it to `indri_gwas(delta =)`.
