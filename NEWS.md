# MRCCC 0.99.7

* `onek1kExample()` provides a small real-data example from the OneK1K
  cohort (NK cells to monocytes, interleukin signalling), derived from the
  public data of the paper by `inst/scripts/make-onek1k-example.R`. The
  vignette fits and interprets the IL-18 triplet on it.
* `mrcccScreen()` fits many triplets, in parallel through BiocParallel, and
  `screenTable()` summarises them with Bayesian false discovery rate
  control. `bayesFDR()` and `bfdrSelect()` accept a list of fits.
* `mrccc()` gains a `BPPARAM` argument. Each chain runs under its own seed,
  derived from `seed`, so the result does not depend on the back end.
* `MRCCCData` extends `S4Vectors::Annotated`, so the instrument-selection
  tables are read with `metadata()`.
* New accessors `firstStageEffects()` (posterior effect of each instrument
  on expression), `posteriorMeans()`, `inputScale()` and `fitSettings()`.
* `credibleIntervals()` gains an `included` argument. With
  `included = TRUE` the intervals use only the draws with `gamma = 1`, as
  reported for the discoveries in the paper.
* `buildTriplet()` accepts the name of a numeric `colData` column of the
  receiver as a precomputed pathway score. The covariate design is built on
  the retained donors, so a factor level found only among dropped donors
  leaves no empty column.
* Input checks are clearer. Missing covariate values and repeated donor
  names are reported for every component. `selectInstruments()` stops when
  the gene's expression is the same in every donor, and its message for an
  empty result gives the number of SNPs in the window and after the carrier
  filter. `initScale` must be positive when `nChains > 1`, and `ridge` must
  be positive.
* The borderline-PIP warning suggests a chain length that pools about
  400,000 draws, as in the paper, for the number of chains used.
* The fitting code of `mrccc()` is split into small internal helpers.
* In the first iteration of a chain with dispersed starting values, the
  update of the outcome variance now uses the drawn starting values of
  `beta_X` and `beta_XZ` rather than zero. This affected only that one
  iteration, which is always within the burn-in. Matrix inversions that the
  sampler computed but never used are removed.
* The default `minCandidates` of `selectInstruments()` is now 1, as in the
  analysis of the paper: a gene keeps its instruments when at least one SNP
  in its window can be screened.
* The vignette has an installation section, a comparison with related
  software, examples of building the inputs from a `SingleCellExperiment`
  and from a VCF file, examples of `MRCCCData()`, the accessors and
  `mrcccHyperparameters()`, a screen of several simulated triplets with the
  code for the full real-data screen, and references.
* The help pages describe every element of the stored settings, summaries
  and diagnostics, explain the choice of R-hat threshold, and point to the
  vignette. The National Institutes of Health are listed as funder with the
  grant numbers.
* New tests check that the sampler recovers simulated signals and nulls and
  that the credible intervals cover the generating values. They also check
  the diagnostics against hand-computed values, and the convergence
  warnings.

# MRCCC 0.99.6

* The interface is now built on S4 classes. `MRCCCData` holds one
  ligand-receptor-pathway triplet with a validity method that checks donor
  alignment, allele dosages in [0, 2] and missing values. `mrccc()` fits the
  model and returns an `MRCCCFit` object, read with the accessors `pip()`,
  `pipMCSE()`, `effectSizes()`, `convergence()`, `instrumentStrength()` and
  `posteriorDraws()`.
* `buildTriplet()` assembles a triplet from sender and receiver
  `SummarizedExperiment` objects and a genotype `RangedSummarizedExperiment`:
  donors are matched by column name, covariates are read from `colData`,
  and the pathway is scored from the receiver's expression of its genes.
* `selectInstruments()` selects cis-eQTL instruments as in the paper: a
  window of 200 kb around the transcription start site, a carrier-frequency
  filter, a covariate-adjusted screen in the same donors, and greedy LD
  clumping at r^2 <= 0.8 to at most ten SNPs.
* `simulateMRCCC()` returns simulated data in this layout.
* Exported functions use camelCase: `credibleIntervals()`,
  `signProbabilities()`, `signReversal()`, `plotEffectCurve()`,
  `bayesFDR()`, `bfdrSelect()` and `mrcccHyperparameters()`. `mr_ccc()`,
  `se_row()` and `simulate_mrccc()` are replaced by the functions above,
  and the compiled sampler is no longer exported.
* The vignette describes the input data, including the units of the
  genotypes and the definition of the instruments, and its example fit
  converges at the default chain length.
* The help page of `mrccc()` and the vignette state the chain settings used
  in the paper, for analyses that compare with them.

# MRCCC 0.99.5

* `mr_ccc_gibbs()` gains the argument `save_first_stage` (default `FALSE`).
  When `TRUE`, the retained draws of the first-stage instrument effects are
  returned as `Pi_X_draws` and `Pi_Z_draws`, so that credible intervals for
  the instrument effects can be computed inside the model. Storing them
  uses no random numbers. All other outputs are unchanged.

# MRCCC 0.99.4

* The CITATION entry includes the DOI of the preprint.

# MRCCC 0.99.3

* `V` is now a required argument of `mr_ccc()`.
* `sign_reversal()` uses a scale-free cutoff for near-zero interaction
  draws: a draw is dropped when its standardised interaction coefficient is
  at most `1e-8` in absolute value.
* Additional input-validation tests.

# MRCCC 0.99.2

* `mr_ccc()` now centres every column of the instrument matrices `G`, `H`
  and the covariate matrix `V`, as it already did for `X`, `Z` and `Y`.
  Neither first stage of the model has an intercept, so uncentred designs
  (for example genotype dosages coded 0/1/2) forced each first-stage fit
  through the origin and shifted the projections `X*` and `Z*` and the fixed
  prior scale. The documentation of `mr_ccc_gibbs()`, which does not centre
  its inputs, now states that all six inputs must be centred.
* `mr_ccc()` stops with an explicit message when `cbind(G, V)` or
  `cbind(H, V)` is rank deficient after centring.
* `mr_ccc()` warns when any Gelman-Rubin R-hat exceeds 1.01, and issues a
  message when only one chain is run or when the Monte Carlo standard error
  of the PIP cannot be estimated.
* `mr_ccc()` and `simulate_mrccc()` restore the caller's random number state
  after running under a supplied `seed`. The fitted object records the
  g-prior scale used in `settings$g`.
* `mr_ccc_gibbs()` validates its arguments, reports numerically singular
  matrices with an informative error, and can be interrupted from R.

# MRCCC 0.99.1

* The priors on the second-stage coefficients `(beta_X, beta_XZ)` and
  `beta_Z` are now scaled by a fixed diagonal computed once from the
  least-squares first stage, instead of by the current `X*`, `Z*` at every
  sweep. Because `X*` and `Z*` are functions of the first-stage parameters,
  the prior used in 0.99.0 entered their full conditionals, which the sampler did not
  account for. With the scale fixed the first-stage updates are exact as
  written and every update remains conjugate. The scale is returned in
  `settings$prior_scale`.

# MRCCC 0.99.0

* Initial Bioconductor submission.
* `mr_ccc()` fits the MR-CCC model to one ligand-receptor-pathway triplet
  with one or more Gibbs chains, returning the posterior inclusion
  probability with its Monte Carlo standard error, effect estimates on the
  raw and standardised scales, pooled posterior draws, the joint
  log-likelihood of every retained draw, and the Gelman-Rubin statistic.
* `mr_ccc()` also reports first-stage instrument strength: the partial F of
  each instrument block after the covariates, with a message when either
  block is weak by the conventional two-stage least squares threshold.
* `credible_intervals()`, `sign_probabilities()`, `sign_reversal()` and
  `plot_effect_curve()` summarise the posterior draws.
* `bayes_fdr()` and `bfdr_threshold_set()` provide Bayesian false discovery
  rate control across many triplets.
* `simulate_mrccc()` generates data from the working model with an
  unmeasured confounder.
* `mr_ccc_gibbs()` exposes the compiled RcppArmadillo sampler.
