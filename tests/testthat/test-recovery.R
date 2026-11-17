# Checks of the results themselves, beyond their structure.

test_that("a strong simulated signal gives a high PIP and the right sign", {
    sim <- simulateMRCCC(n = 300, betaX = 0.3, betaXZ = 0.3, seed = 11)
    dat <- makeTriplet(n = 300, seed = 11, betaX = 0.3, betaXZ = 0.3)
    fit <- suppressWarnings(suppressMessages(
        mrccc(dat, nIter = 4000, burnIn = 1000, nChains = 2, seed = 1)))
    expect_gt(pip(fit), 0.9)
    raw <- effectSizes(fit, standardized = FALSE)
    expect_gt(raw$mean[raw$parameter == "beta_XZ"], 0)
    expect_lt(abs(raw$mean[raw$parameter == "beta_XZ"] - 0.3), 0.2)
    expect_lt(abs(raw$mean[raw$parameter == "beta_X"] - 0.3), 0.2)
    # The 95% credible intervals cover the generating values.
    ci <- credibleIntervals(fit, standardized = FALSE)
    truth <- c(beta_X = sim$truth$betaX, beta_XZ = sim$truth$betaXZ)
    ci <- ci[match(names(truth), ci$parameter), ]
    expect_true(all(ci$lower <= truth & truth <= ci$upper))
})

test_that("a null simulated triplet gives a low PIP", {
    dat <- makeTriplet(n = 300, seed = 12, betaX = 0, betaXZ = 0)
    fit <- suppressWarnings(suppressMessages(
        mrccc(dat, nIter = 4000, burnIn = 1000, nChains = 2, seed = 1)))
    expect_lt(pip(fit), 0.5)
})

test_that("gelman_rubin matches hand-computed values", {
    gr <- MRCCC:::gelman_rubin
    # Means 2 and 3, within-chain variances 1, n = 3:
    # W = 1, B = 3 * var(c(2, 3)) = 1.5, R = sqrt((2/3) * 1 + 1.5 / 3).
    expect_equal(gr(list(c(1, 2, 3), c(2, 3, 4))), sqrt(7 / 6))
    # Identical chains: B = 0, so R = sqrt((n - 1) / n).
    expect_equal(gr(list(c(1, 2, 3), c(1, 2, 3))), sqrt(2 / 3))
    expect_identical(gr(list(c(1, 1, 1), c(2, 2, 2))), Inf)
    expect_true(is.na(gr(list(c(1, 1, 1), c(1, 1, 1)))))
    expect_true(is.na(gr(list(c(1, 2, 3)))))
})

test_that("mcse_batch and mcse_pooled match hand-computed values", {
    x <- c(rep(0, 30), rep(1, 30))
    # Two batches of 30 with means 0 and 1: sigma2 = 30 * 0.5 = 15.
    mc <- MRCCC:::mcse_batch(x, n_batch = 2L)
    expect_equal(mc[["mcse"]], sqrt(15 / 60))
    expect_equal(mc[["ess"]], 60 * stats::var(x) / 15)
    # A chain whose batch means are all equal has no estimable error.
    expect_true(is.na(MRCCC:::mcse_batch(rep(c(0, 1), 30))[["mcse"]]))
    pooled <- MRCCC:::mcse_pooled(list(c(mcse = 0.1, ess = 10),
                                       c(mcse = 0.2, ess = 20)))
    expect_equal(pooled[["mcse"]], sqrt(0.1^2 + 0.2^2) / 2)
    expect_equal(pooled[["ess"]], 30)
})

test_that("convergence warnings and messages are issued when due", {
    msgs <- MRCCC:::.convergenceMessages
    base <- list(
        diagnostics = data.frame(parameter = c("beta_X", "beta_XZ", "gamma",
                                               "loglik"),
                                 rhat = c(1, 1, 1, 1)),
        nChains = 4L, pip = 0.9, mcse = 0.01, constantGamma = FALSE)
    set <- list(pipThreshold = 0.5, nIter = 20000L)
    expect_silent(msgs(base, set))

    notMixed <- base
    notMixed$diagnostics$rhat[1L] <- 1.05
    expect_warning(msgs(notMixed, set), "R-hat exceeds 1.01 for beta_X")

    borderline <- base
    borderline$pip <- 0.52
    borderline$mcse <- 0.02
    expect_warning(msgs(borderline, set), "within 2 \\* MCSE")
    expect_warning(msgs(borderline, set), "100000")

    oneChain <- base
    oneChain$nChains <- 1L
    oneChain$diagnostics$rhat <- NA_real_
    expect_message(msgs(oneChain, set), "Only one chain was run")

    noMcse <- base
    noMcse$mcse <- NA_real_
    expect_message(msgs(noMcse, set), "could not be estimated")
    noMcse$constantGamma <- TRUE
    expect_silent(msgs(noMcse, set))
})

test_that("first-stage effects and posterior means are reported", {
    dat <- makeTriplet(n = 150, seed = 1)
    fit <- suppressWarnings(suppressMessages(
        mrccc(dat, nIter = 1500, burnIn = 300, nChains = 2, seed = 1)))
    fs <- firstStageEffects(fit)
    expect_equal(fs$snp, c(colnames(ligandGenotypes(dat)),
                           colnames(receptorGenotypes(dat))))
    expect_true(all(fs$lower <= fs$mean & fs$mean <= fs$upper))
    pm <- posteriorMeans(fit)
    expect_equal(pm$covariateEffects$covariate, colnames(covariates(dat)))
    expect_named(pm$residualVariances, c("ligand", "receptor", "pathway"))
    expect_true(all(pm$residualVariances > 0))
    expect_true(pm$inclusionProbabilityPrior > 0 &&
                    pm$inclusionProbabilityPrior < 1)
})

test_that("chains give the same result on a parallel back end", {
    skip_on_os("windows")
    dat <- makeTriplet(n = 150, seed = 1)
    serial <- suppressWarnings(suppressMessages(
        mrccc(dat, nIter = 600, burnIn = 100, nChains = 2, seed = 3)))
    parallel <- suppressWarnings(suppressMessages(
        mrccc(dat, nIter = 600, burnIn = 100, nChains = 2, seed = 3,
              BPPARAM = BiocParallel::MulticoreParam(2))))
    expect_identical(posteriorDraws(parallel), posteriorDraws(serial))
    expect_error(mrccc(dat, BPPARAM = "serial"), "BiocParallelParam")
})

test_that("mrcccScreen and screenTable summarise several triplets", {
    triplets <- list(a = makeTriplet(n = 150, seed = 1),
                     b = makeTriplet(n = 150, seed = 2, betaX = 0,
                                     betaXZ = 0))
    fits <- suppressWarnings(suppressMessages(
        mrcccScreen(triplets, nIter = 1500, burnIn = 300, nChains = 2,
                    seed = 1)))
    expect_named(fits, c("a", "b"))
    expect_true(all(vapply(fits, methods::is, logical(1),
                           class2 = "MRCCCFit")))
    again <- suppressWarnings(suppressMessages(
        mrcccScreen(triplets, nIter = 1500, burnIn = 300, nChains = 2,
                    seed = 1)))
    expect_identical(posteriorDraws(again$b), posteriorDraws(fits$b))
    tab <- screenTable(fits, alpha = 0.2)
    expect_equal(nrow(tab), 2L)
    expect_true(all(diff(tab$pip) <= 0))
    expect_equal(tab$bFDR, sort(tab$bFDR))
    expect_equal(sort(tab$pip), sort(vapply(fits, pip, numeric(1)),
                                     decreasing = FALSE), ignore_attr = TRUE)
    expect_error(mrcccScreen(list(1)), "list of MRCCCData")
    expect_error(screenTable(list(1)), "list of MRCCCFit")
})
