#' @include AllClasses.R AllGenerics.R MRCCCData.R
NULL

#' Fit the MR-CCC model to one ligand-receptor-pathway triplet
#'
#' Fits the Bayesian Mendelian randomization model for causal cell-cell
#' communication to one triplet stored in an [MRCCCData-class] object. The
#' cis-eQTL genotypes instrument ligand expression in the sender cell type
#' and receptor expression in the receiver cell type. The model asks whether
#' the genetically predicted ligand level, modulated by the genetically
#' predicted receptor level, shifts pathway activity in the receiver. Several
#' Gibbs chains are run and pooled, and the posterior inclusion probability
#' (PIP) is returned with its Monte Carlo standard error and convergence
#' diagnostics.
#'
#' @param x An [MRCCCData-class] object, usually from [buildTriplet()].
#' @param nIter Total number of Gibbs iterations per chain.
#' @param burnIn Number of initial iterations discarded per chain. Must be
#'   smaller than `nIter`.
#' @param thin Thinning interval: every `thin`-th post-burn-in iteration is
#'   retained.
#' @param nChains Number of independent chains, pooled for all posterior
#'   summaries. With two or more chains the Gelman-Rubin statistic is
#'   reported and the chains start from dispersed values. Pooling `C` chains
#'   of `nIter` iterations gives roughly the Monte Carlo precision of one
#'   chain of `C * (nIter - burnIn)` retained draws. Each chain repeats the
#'   burn-in, which is the only extra cost.
#' @param initScale Standard deviation of the normal distribution from which
#'   \eqn{\beta_X} and \eqn{\beta_{XZ}} are initialised, on the centred
#'   input scale, when `nChains > 1`. The default 1 gives standard normal
#'   starting values. With a single chain the sampler starts at zero.
#' @param pipThreshold Threshold in \eqn{(0, 1)} that defines a discovery
#'   (`PIP > pipThreshold`). The default of 0.5 corresponds to the median
#'   probability model. Used for the borderline warning and by
#'   [summary()][summary,MRCCCFit-method].
#' @param hyper Prior hyperparameters, as returned by
#'   [mrcccHyperparameters()].
#' @param seed Optional integer. If supplied, the chains are run under this
#'   seed, which makes the whole multi-chain fit reproducible, and the
#'   caller's random number state is restored afterwards.
#' @param verbose Logical. If `TRUE`, each chain prints a progress line every
#'   1000 iterations.
#' @param ... Not used.
#'
#' @details
#' **Model.** For donor \eqn{i}, with \eqn{X^* = E[X \mid G, V]} and
#' \eqn{Z^* = E[Z \mid H, V]} the genetically predicted ligand and receptor
#' levels, the outcome model is
#' \deqn{Y = \mu + \gamma(\beta_X X^* + \beta_{XZ} X^* Z^*) + \beta_Z Z^*
#'   + \alpha_Y^\top V + \tilde\varepsilon_Y,}
#' where \eqn{X} is ligand expression, \eqn{Z} receptor expression, \eqn{Y}
#' pathway activity, \eqn{G} and \eqn{H} the instrument dosages, \eqn{V} the
#' covariates, and \eqn{\gamma \in \{0, 1\}} a spike-and-slab inclusion
#' indicator for the causal block \eqn{(\beta_X, \beta_{XZ})}. Under
#' \eqn{\gamma = 0} the prior variance of that block is multiplied by `nu1`
#' (\eqn{10^{-4}} by default), so the coefficients are concentrated tightly
#' at zero. The PIP \eqn{P(\gamma = 1 \mid \mathrm{data})} is the
#' communication score. The receptor-modulated ligand effect is
#' \eqn{\beta_X + \beta_{XZ} z}, which changes sign at
#' \eqn{\tau = -\beta_X / \beta_{XZ}}. See [signReversal()].
#'
#' **Centring and scale.** Neither first stage carries an intercept, so the
#' expression, pathway, genotype and covariate data are all centred
#' internally. Effects are reported on the centred input scale and on a
#' standardised scale (per standard deviation of \eqn{X} and \eqn{Z}, in
#' standard deviations of \eqn{Y}). See [effectSizes()].
#'
#' **Priors.** The first-stage blocks and the covariate block carry Zellner
#' g-priors on their designs. The second-stage blocks are scaled by a fixed
#' diagonal computed once from the least-squares first stage, so every
#' Gibbs update is an exact conjugate step. Unless overridden in `hyper`,
#' every scale \eqn{g} is \eqn{\min(n, 100)}. See [mrcccHyperparameters()].
#'
#' **Instrument strength.** The first-stage partial F statistic of each
#' instrument set, after the covariates, is reported by
#' [instrumentStrength()]. A message is issued when either F is below 10.
#' MR-CCC does not require F > 10, but PIPs are less decisive with weak
#' instruments and longer chains may be needed.
#'
#' **Chains, Monte Carlo error and convergence.** Odd-numbered chains start
#' with \eqn{\gamma = 1} and even-numbered chains with \eqn{\gamma = 0}.
#' The Monte Carlo standard error (MCSE) of the PIP is estimated by
#' non-overlapping batch means within each chain and combined across chains.
#' It is `NA` when the indicator never changes value in any chain, so that
#' the PIP is exactly 0 or 1. A warning is issued when
#' \eqn{|\mathrm{PIP} - t| < 2\,\mathrm{MCSE}} for threshold \eqn{t},
#' because the discovery status is then not resolved at the current chain
#' length, and when any Gelman-Rubin \eqn{\hat R} exceeds
#' 1.01. \eqn{\hat R} is reported for \eqn{\beta_X}, \eqn{\beta_{XZ}},
#' \eqn{\gamma} and the joint log-likelihood. It requires at least two
#' chains. The default chain length is a starting point. Reported analyses
#' should use chains long enough that neither warning is issued.
#'
#' **Settings used in the paper.** The real-data analysis ran four chains of
#' 100,000 iterations per triplet with a burn-in of 2,000. The simulation
#' studies ran one chain per fit, of 20,000 iterations when the instruments
#' are strong (scenarios S1 to S3) and of 400,000 iterations under weak or
#' otherwise misspecified instruments (scenarios S4 to S7), with a burn-in of
#' 2,000. A comparison with MR-CCC should use chains of these lengths, or
#' longer, for the corresponding setting.
#'
#' @return An [MRCCCFit-class] object. See [MRCCCFit-accessors].
#'
#' @seealso [buildTriplet()], [MRCCCFit-accessors], [credibleIntervals()],
#'   [signProbabilities()], [signReversal()], [plotEffectCurve()],
#'   [bayesFDR()]
#'
#' @references
#' Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for
#' Causal Cell-Cell Communication. arXiv:2604.23917.
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' fit
#' pip(fit)
#' effectSizes(fit)
#'
#' @rdname mrccc
#' @export
setMethod("mrccc", "MRCCCData",
    function(x, nIter = 20000L, burnIn = 2000L, thin = 1L, nChains = 4L,
             initScale = 1, pipThreshold = 0.5,
             hyper = mrcccHyperparameters(), seed = NULL, verbose = FALSE,
             ...) {
    .fitMRCCC(
        X = x@ligandExpression, Z = x@receptorExpression,
        Y = x@pathwayActivity, G = x@ligandGenotypes,
        H = x@receptorGenotypes, V = x@covariates,
        nIter = nIter, burnIn = burnIn, thin = thin, nChains = nChains,
        initScale = initScale, pipThreshold = pipThreshold, hyper = hyper,
        seed = seed, verbose = verbose, labels = x@labels)
})

# The fitting engine. Arguments are numeric vectors (X, Z, Y) and matrices
# (G, H, V) already aligned by donor, as guaranteed by the MRCCCData validity
# method.
.fitMRCCC <- function(X, Z, Y, G, H, V, nIter, burnIn, thin, nChains,
                      initScale, pipThreshold, hyper, seed, verbose,
                      labels = .makeLabels(NULL)) {
    X <- as_column(X, "ligandExpression")
    Z <- as_column(Z, "receptorExpression")
    Y <- as_column(Y, "pathwayActivity")
    G <- as_design(G, "ligandGenotypes")
    H <- as_design(H, "receptorGenotypes")
    V <- as_design(V, "covariates")
    n <- nrow(X)

    check_design(G, "ligandGenotypes")
    check_design(H, "receptorGenotypes")
    check_design(V, "covariates")
    p_total <- ncol(G) + ncol(H) + ncol(V)
    if (n <= p_total + 4L) {
        stop("Too few donors for the model: ", n, " donors against ",
             p_total, " instrument and covariate columns plus four ",
             "second-stage coefficients. Supply more donors or fewer ",
             "instruments.", call. = FALSE)
    }

    check_scalar(nIter, "nIter")
    check_scalar(burnIn, "burnIn")
    check_scalar(thin, "thin")
    check_scalar(nChains, "nChains")
    check_scalar(initScale, "initScale")
    check_scalar(pipThreshold, "pipThreshold")
    if (nIter != round(nIter) || nIter < 1) {
        stop("'nIter' must be a positive integer.", call. = FALSE)
    }
    if (burnIn != round(burnIn) || burnIn < 0) {
        stop("'burnIn' must be a non-negative integer.", call. = FALSE)
    }
    if (nIter <= burnIn) {
        stop("'nIter' (", nIter, ") must be greater than 'burnIn' (",
             burnIn, ").", call. = FALSE)
    }
    if (thin != round(thin) || thin < 1) {
        stop("'thin' must be an integer greater than or equal to 1.",
             call. = FALSE)
    }
    if (nChains != round(nChains) || nChains < 1) {
        stop("'nChains' must be an integer greater than or equal to 1.",
             call. = FALSE)
    }
    if (initScale < 0) {
        stop("'initScale' must be non-negative.", call. = FALSE)
    }
    if (pipThreshold <= 0 || pipThreshold >= 1) {
        stop("'pipThreshold' must lie strictly between 0 and 1.",
             call. = FALSE)
    }
    if (!methods::is(hyper, "mrcccHyperparameters")) {
        stop("'hyper' must be created by mrcccHyperparameters().",
             call. = FALSE)
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(seed)) {
        check_scalar(seed, "seed")
        if (seed != round(seed)) {
            stop("'seed' must be a whole number.", call. = FALSE)
        }
    }
    nKeepExpected <- (nIter - burnIn) %/% thin
    if (nKeepExpected < 2L) {
        stop("The chosen 'nIter' (", nIter, "), 'burnIn' (", burnIn,
             ") and 'thin' (", thin, ") retain only ", nKeepExpected,
             " draw(s) per chain. At least two are required.",
             call. = FALSE)
    }
    nIter <- as.integer(nIter)
    burnIn <- as.integer(burnIn)
    thin <- as.integer(thin)
    nChains <- as.integer(nChains)

    # ---- Centring and scale ----
    sdX <- stats::sd(X[, 1L])
    sdZ <- stats::sd(Z[, 1L])
    sdY <- stats::sd(Y[, 1L])
    if (any(!is.finite(c(sdX, sdZ, sdY))) || any(c(sdX, sdZ, sdY) <= 0)) {
        stop("Ligand expression, receptor expression and pathway activity ",
             "must each vary across donors.", call. = FALSE)
    }
    meanX <- mean(X[, 1L])
    meanZ <- mean(Z[, 1L])
    meanY <- mean(Y[, 1L])
    Xc <- X - meanX
    Zc <- Z - meanZ
    Yc <- Y - meanY
    zRangeStd <- range(Zc[, 1L] / sdZ)
    Gc <- center_columns(G)
    Hc <- center_columns(H)
    Vc <- center_columns(V)
    check_rank(cbind(Gc, Vc), "ligand instruments and covariates")
    check_rank(cbind(Hc, Vc), "receptor instruments and covariates")

    # ---- Instrument strength ----
    instruments <- data.frame(
        exposure = c("ligand", "receptor"),
        nInstruments = c(ncol(G), ncol(H)),
        F = c(first_stage_F(Xc[, 1L], G, V), first_stage_F(Zc[, 1L], H, V)),
        stringsAsFactors = FALSE
    )
    instruments$weak <- is.finite(instruments$F) & instruments$F < 10
    if (any(instruments$weak)) {
        wk <- instruments[instruments$weak, , drop = FALSE]
        weakText <- sprintf("%s (F = %s)", wk$exposure,
                            vapply(wk$F, fmt, character(1)))
        message("First-stage F below 10 for the ",
                gsub(", ", " and ", toString(weakText), fixed = TRUE),
                ": the instruments are weak by the two-stage least-squares ",
                "heuristic. MR-CCC does not require F > 10, but posterior ",
                "inclusion probabilities are less decisive in this regime ",
                "and longer chains may be needed. See instrumentStrength().")
    }

    gVal <- if (is.null(hyper$g)) min(n, 100) else hyper$g
    initRaw <- if (nChains > 1L) initScale else 0

    # ---- Run the chains ----
    if (!is.null(seed)) withr::local_seed(as.integer(seed))
    fits <- vector("list", nChains)
    for (cc in seq_len(nChains)) {
        fits[[cc]] <- mr_ccc_gibbs(
            Xc, Zc, Yc, Gc, Hc, Vc,
            n_iter = nIter, burn_in = burnIn, thin = thin,
            a_sigma = hyper$aSigma, b_sigma = hyper$bSigma,
            a_rho = hyper$aRho, b_rho = hyper$bRho, nu1 = hyper$nu1,
            gG = gVal, gV = gVal, gH = gVal, gZ = gVal, gBeta = gVal,
            ridge = hyper$ridge,
            init_gamma = if (cc %% 2L == 0L) 0L else 1L,
            init_scale = initRaw, verbose = verbose)
    }
    nKeep <- as.integer(fits[[1L]]$n_keep)

    chainOf <- function(name) lapply(fits, function(f) as.numeric(f[[name]]))
    bXch <- chainOf("Beta_X_draws")
    bXZch <- chainOf("Beta_XZ_draws")
    gch <- chainOf("gamma_draws")
    llch <- chainOf("loglik_draws")
    draws <- data.frame(
        chain = rep(seq_len(nChains), each = nKeep),
        beta_X = unlist(bXch, use.names = FALSE),
        beta_XZ = unlist(bXZch, use.names = FALSE),
        beta_Z = unlist(chainOf("Beta_Z_draws"), use.names = FALSE),
        gamma = unlist(gch, use.names = FALSE),
        mu = unlist(chainOf("mu_draws"), use.names = FALSE),
        loglik = unlist(llch, use.names = FALSE)
    )

    pipVal <- mean(draws$gamma)
    mc <- mcse_pooled(lapply(gch, mcse_batch))
    pipMcse <- unname(mc[["mcse"]])
    pipEss <- unname(mc[["ess"]])

    diagPars <- c("beta_X", "beta_XZ", "gamma", "loglik")
    diagnostics <- data.frame(
        parameter = diagPars,
        rhat = if (nChains >= 2L) {
            vapply(list(bXch, bXZch, gch, llch), gelman_rubin, numeric(1))
        } else {
            rep(NA_real_, length(diagPars))
        },
        stringsAsFactors = FALSE
    )

    rh <- diagnostics$rhat
    bad <- !is.na(rh) & (is.infinite(rh) | rh > 1.01)
    if (any(bad)) {
        warning("Gelman-Rubin R-hat exceeds 1.01 for ",
                toString(diagnostics$parameter[bad]),
                ": the chains have not mixed. Increase 'nIter' and 'burnIn'.",
                call. = FALSE)
    }
    if (nChains == 1L) {
        message("Only one chain was run, so R-hat is unavailable. Use ",
                "nChains >= 2 to assess convergence.")
    }
    # When every chain holds the indicator at the same value throughout, the
    # PIP is 0 or 1 with no observed Monte Carlo variation and the MCSE is
    # reported as NA without comment.
    constantGamma <- length(unique(draws$gamma)) == 1L
    if (!is.finite(pipMcse) && !constantGamma) {
        message("The Monte Carlo uncertainty of the PIP could not be ",
                "estimated (too few retained draws per chain, or a constant ",
                "indicator chain), so the borderline check was skipped.")
    }
    if (is.finite(pipMcse) && abs(pipVal - pipThreshold) < 2 * pipMcse) {
        suggestion <- if (nIter < 100000L) {
            paste("Increase 'nIter' to 100000, or to 400000 if the PIP",
                  "remains borderline")
        } else if (nIter < 400000L) {
            "Increase 'nIter' to 400000"
        } else {
            paste("Increase 'nChains', which pools additional draws at",
                  "proportional cost")
        }
        warning("The posterior inclusion probability (", fmt(pipVal),
                ") lies within 2 * MCSE (MCSE = ", fmt(pipMcse),
                ") of the discovery threshold (", pipThreshold,
                "): the discovery status of this triplet is not resolved at ",
                "this chain length. ", suggestion, ".", call. = FALSE)
    }

    methods::new("MRCCCFit",
        pip = pipVal, pipMCSE = pipMcse, pipESS = pipEss,
        draws = draws,
        scale = list(sdX = sdX, sdZ = sdZ, sdY = sdY, meanX = meanX,
                     meanZ = meanZ, meanY = meanY, zRangeStd = zRangeStd),
        diagnostics = diagnostics,
        instruments = instruments,
        settings = list(nIter = nIter, burnIn = burnIn, thin = thin,
                        nChains = nChains, nKeep = nKeep,
                        pipThreshold = pipThreshold, hyper = hyper,
                        g = gVal, initScale = initScale,
                        priorScale = fits[[1L]]$prior_scale),
        labels = labels)
}
