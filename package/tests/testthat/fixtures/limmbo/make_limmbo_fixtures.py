"""Reference association results from limmbo for indri's GWAS tests.

Run from the limmbo checkout:  pixi run python make_limmbo_fixtures.py OUTDIR
Uses limmbo's own test data, the first N_SAMPLES samples (to keep the
relatedness fixture small), covariates (limmbo's reader appends the intercept column),
and the shipped Cg/Cn. Writes every input and output as CSV.
"""
import os
import sys

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from limmbo.io.reader import ReadData
from limmbo.io.input import InputData
from limmbo.core.gwas import GWAS
from limmbo.utils.utils import effectiveTests

N_SAMPLES = 250
out = sys.argv[1]
os.makedirs(out, exist_ok=True)


def data(name):
    return resource_filename('limmbo', 'io/test/data/' + name)


rd = ReadData(verbose=False)
rd.getPhenotypes(file_pheno=data('pheno.csv'))
rd.getCovariates(file_covariates=data('covs.csv'))
rd.getRelatedness(file_relatedness=data('relatedness.csv'), delim=',')
rd.getGenotypes(file_genotypes=data('genotypes.csv'))
rd.getVarianceComponents(file_Cg=data('Cg.csv'), file_Cn=data('Cn.csv'))

# getCovariates() appends a column of ones itself, so the covariates already
# carry the intercept indri requires.
covs = rd.covariates

ind = InputData(verbose=False)
ind.addPhenotypes(phenotypes=rd.phenotypes, phenotype_ID=rd.phenotype_ID,
                  pheno_samples=np.asarray(rd.phenotypes.index))
ind.addCovariates(covariates=covs)
ind.addRelatedness(relatedness=rd.relatedness,
                   relatedness_samples=rd.relatedness_samples)
ind.addGenotypes(genotypes=rd.genotypes, geno_samples=np.asarray(rd.genotypes.index),
                 genotypes_info=rd.genotypes_info)
ind.addVarianceComponents(Cg=rd.Cg, Cn=rd.Cn)
samples = np.asarray(rd.phenotypes.index)[:N_SAMPLES]
ind.commonSamples(samplelist=samples)

ids = list(ind.samples)
traits = list(ind.phenotype_ID)
snps = list(ind.genotypes_info.index)


def frame(x, rows, cols):
    return pd.DataFrame(np.asarray(x), index=rows, columns=cols)


frame(ind.phenotypes, ids, traits).to_csv(os.path.join(out, 'Y.csv'))
covnames = ['cov%d' % i for i in range(1, covs.shape[1])] + ['intercept']
frame(ind.covariates, ids, covnames).to_csv(os.path.join(out, 'X.csv'))
frame(ind.genotypes, ids, snps).to_csv(os.path.join(out, 'genotypes.csv'))
frame(ind.relatedness, ids, ids).to_csv(os.path.join(out, 'K.csv'))
frame(ind.Cg, traits, traits).to_csv(os.path.join(out, 'Cg.csv'))
frame(ind.Cn, traits, traits).to_csv(os.path.join(out, 'Cn.csv'))
info = ind.genotypes_info.copy()
info.index.name = 'snp'
info.to_csv(os.path.join(out, 'snp_info.csv'))

gwas = GWAS(datainput=ind, verbose=False)
ldelta = {}
for setup in ('lmm', 'lm'):
    res = gwas.runAssociationAnalysis(setup=setup, mode='multitrait')
    tag = 'multi_trait_%s' % setup
    # The null variance ratio LIMIX fitted: V = Cg (x) K + delta Cn (x) I.
    ldelta[tag] = [float(np.asarray(res['lm'].getLdelta0()).ravel()[0])]
    pd.DataFrame({'snp': snps, 'pvalue': np.asarray(res['pvalues']).ravel()}
                 ).to_csv(os.path.join(out, tag + '_pvalues.csv'), index=False)
    frame(np.asarray(res['betas']).T, snps, traits).to_csv(
        os.path.join(out, tag + '_betas.csv'))
    for adjust in (None, 'bonferroni', 'effective'):
        gwas.adjustBy = None
        res = gwas.runAssociationAnalysis(setup=setup, mode='singletrait',
                                          adjustSingleTrait=adjust)
        tag = 'single_trait_%s' % setup
        if adjust is None:
            # Per trait, V = K + delta I (K = I for the LM).
            ldelta[tag] = list(np.asarray(res['lm'].ldelta_0)[:, 0])
            frame(np.asarray(res['pvalues']).T, snps, traits).to_csv(
                os.path.join(out, tag + '_pvalues.csv'))
            frame(np.asarray(res['betas']).T, snps, traits).to_csv(
                os.path.join(out, tag + '_betas.csv'))
        else:
            frame(np.asarray(res['pvalues_adjust']).T, snps, traits).to_csv(
                os.path.join(out, '%s_padjust_%s.csv' % (tag, adjust)))

# The covariance limmbo's multi-trait LM uses for Cn, and the effective
# number of tests, both from the phenotypes as passed.
frame(np.cov(np.asarray(ind.phenotypes).T), traits, traits).to_csv(
    os.path.join(out, 'Cn_lm.csv'))
pd.DataFrame({'effective_tests': [effectiveTests(np.asarray(ind.phenotypes))]}
             ).to_csv(os.path.join(out, 'effective_tests.csv'), index=False)
pd.DataFrame({'test': [t for t in ldelta for _ in ldelta[t]],
              'trait': [traits[i] if len(ldelta[t]) > 1 else ''
                        for t in ldelta for i in range(len(ldelta[t]))],
              'log_delta': [v for t in ldelta for v in ldelta[t]]}
             ).to_csv(os.path.join(out, 'log_delta.csv'), index=False)
print('wrote', sorted(os.listdir(out)))
