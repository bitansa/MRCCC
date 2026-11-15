# MRCCC

Bayesian Mendelian randomization for causal cell-cell communication.

MRCCC implements MR-CCC, a Bayesian Mendelian randomization framework for
inferring causal ligand-receptor signalling between cell types from
population-scale single-cell expression and genotype data. For each
ligand-receptor-pathway triplet, cis-eQTL genotypes instrument ligand
expression in the sender cell type and receptor expression in the receiver
cell type, and a working model with a ligand-by-receptor interaction is fitted
to downstream pathway activity by a blocked Gibbs sampler written in
RcppArmadillo. A spike-and-slab prior on the causal block yields a posterior
inclusion probability (PIP) that serves as the communication score. The
retained draws support credible intervals, posterior sign probabilities, the
receptor level at which the ligand effect reverses sign, Monte Carlo error and
convergence diagnostics, and Bayesian false discovery rate control across
many triplets.

## Installation

Once the package is available from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MRCCC")
```

The development version can be installed from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("bitansa/MRCCC")
```

The sampler is compiled from C++ through RcppArmadillo, so a working C++
toolchain is required (Rtools on Windows, Xcode command line tools on macOS,
`build-essential` and `gfortran` on Debian-based Linux).

## Quick start

MRCCC works on donor-level (pseudo-bulk) expression of a sender and a
receiver cell type, stored as `SummarizedExperiment` objects with donors in
columns, and on the genotypes of the same donors, stored as a
`RangedSummarizedExperiment` of allele dosages (0, 1 or 2 copies of the
alternative allele) with SNP positions. `simulateMRCCC()` produces data in
this layout.

```r
library(MRCCC)

sim <- simulateMRCCC(n = 300, seed = 1)

# Select cis-eQTL instruments, match donors, read covariates from colData
# and score the pathway, all in one step.
dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
                    ligand = "LIG1", receptor = "REC1",
                    pathway = sim$pathwayGenes,
                    covariates = c("sex", "age", "PC1", "PC2"))

fit <- mrccc(dat, seed = 1)
fit
pip(fit)
effectSizes(fit)
convergence(fit)
credibleIntervals(fit)
signProbabilities(fit)
signReversal(fit)
plotEffectCurve(fit)
bayesFDR(c(0.95, 0.80, 0.55, 0.20))
```

`selectInstruments()` implements the instrument selection of the paper:
SNPs within 200 kb of the gene's transcription start site are screened by a
covariate-adjusted regression of the gene's expression on dosage in the same
donors, ranked by p-value and LD-clumped (r^2 <= 0.8) to at most ten.

## Chain length and convergence

`mrccc()` runs `nChains` chains (default 4) of `nIter` iterations (default
20,000, of which the first 2,000 are discarded as burn-in) from dispersed
starting values and pools the retained draws. Per-chain Monte Carlo standard errors combine
in quadrature, so `C` chains of `nIter` iterations carry about the same
standard error as one chain of `C * nIter`. Several chains also give the
Gelman-Rubin statistic, reported by `convergence()`. `mrccc()` warns when
any R-hat exceeds 1.01.

`instrumentStrength()` reports the first-stage partial F statistic of each
instrument set, and a message is issued when either falls below 10. The
one-sample Bayesian model does not require F > 10, but posterior inclusion
probabilities are less decisive under weak instruments and longer chains may
be needed.

The posterior inclusion probability is the mean of an autocorrelated binary
chain, so its Monte Carlo standard error, reported by `pipMCSE()`, matters
when the PIP is near the discovery threshold. `mrccc()` warns when the PIP
lies within two standard errors of it. The real-data analysis in the paper
uses four chains of 100,000 iterations per triplet. Its simulation studies
use one chain per fit, of 20,000 iterations when the instruments are strong
and of 400,000 iterations when they are weak or otherwise misspecified.

## Reproducing the paper

The scripts that reproduce the simulation studies and the real-data analysis
of the paper are maintained separately at
<https://github.com/bitansa/MR-CCC> and archived on Zenodo under DOI
[10.5281/zenodo.23167575](https://doi.org/10.5281/zenodo.23167575). The
processed data used by those scripts are archived on Zenodo under DOI
[10.5281/zenodo.19675075](https://doi.org/10.5281/zenodo.19675075).

## Citation

Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for Causal
Cell-Cell Communication. arXiv:2604.23917.
<https://arxiv.org/abs/2604.23917>

```r
citation("MRCCC")
```

## License

MIT. See the `LICENSE` file.
