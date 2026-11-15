dat <- makeTriplet(n = 150, seed = 2, betaX = 0.3, betaXZ = -0.3)
fit <- suppressWarnings(mrccc(dat, nIter = 1500, burnIn = 300, seed = 1))

test_that("bayesFDR reproduces the nested-set definition in input order", {
    expect_equal(bayesFDR(c(0.9, 0.8, 0.5)), c(0.1, 0.15, 0.2666667),
                 tolerance = 1e-6)
    expect_equal(bayesFDR(c(0.5, 0.9, 0.8)), c(0.2666667, 0.1, 0.15),
                 tolerance = 1e-6)
    expect_length(bayesFDR(numeric(0)), 0L)
    expect_error(bayesFDR(c(0.5, NA)), "NA")
    expect_error(bayesFDR(c(0.5, 1.2)), "\\[0, 1\\]")
})

test_that("bfdrSelect selects the largest qualifying nested set", {
    p <- c(0.9, 0.8, 0.5, 0.2)
    expect_equal(bfdrSelect(p, alpha = 0.2), c(TRUE, TRUE, FALSE, FALSE))
    expect_equal(bfdrSelect(p, alpha = 0.3), c(TRUE, TRUE, TRUE, FALSE))
    expect_equal(bfdrSelect(p, alpha = 0.05), rep(FALSE, 4))
    expect_equal(bfdrSelect(c(0.2, 0.9, 0.5, 0.8), alpha = 0.2),
                 c(FALSE, TRUE, FALSE, TRUE))
    expect_error(bfdrSelect(p, alpha = 1), "'alpha' must lie")
})

test_that("bayesFDR accepts a list of fits", {
    expect_equal(bayesFDR(list(fit, fit)), bayesFDR(rep(pip(fit), 2)))
})

test_that("credibleIntervals are ordered and respect the scale", {
    ci <- credibleIntervals(fit)
    expect_equal(ci$parameter, c("beta_X", "beta_XZ", "beta_Z"))
    expect_true(all(ci$lower <= ci$median & ci$median <= ci$upper))
    ciRaw <- credibleIntervals(fit, standardized = FALSE)
    s <- fit@scale
    sf <- c(s$sdX / s$sdY, s$sdX * s$sdZ / s$sdY, s$sdZ / s$sdY)
    expect_equal(ci$median, ciRaw$median * sf)
    ci90 <- credibleIntervals(fit, level = 0.9)
    expect_true(all(ci90$lower >= ci$lower & ci90$upper <= ci$upper))
    expect_error(credibleIntervals(fit, level = 1), "'level' must lie")
})

test_that("signProbabilities returns valid probabilities", {
    sp <- signProbabilities(fit)
    expect_equal(nrow(sp), 1L)
    expect_equal(sp$pip, pip(fit))
    expect_true(sp$pPositive + sp$pNegative <= 1 + 1e-12)
    if (sp$nIncluded > 0) {
        expect_true(sp$pPositiveIncluded + sp$pNegativeIncluded <= 1 + 1e-12)
    }
})

test_that("signReversal uses paired draws on the standardised scale", {
    sr <- signReversal(fit)
    expect_equal(nrow(sr), 1L)
    expect_true(sr$nDraws > 0)
    expect_true(sr$tauLower <= sr$tauMedian && sr$tauMedian <= sr$tauUpper)
    expect_true(sr$pInRange >= 0 && sr$pInRange <= 1)
    d <- posteriorDraws(fit)
    s <- fit@scale
    ok <- abs(d$beta_XZ * s$sdX * s$sdZ / s$sdY) > 1e-8
    tau <- (-d$beta_X[ok] / d$beta_XZ[ok]) / s$sdZ
    expect_equal(sr$tauMedian, stats::median(tau))
    expect_equal(c(sr$zMin, sr$zMax), s$zRangeStd)
})

test_that("summaries reject objects that are not fits", {
    expect_error(credibleIntervals(list()), "unable to find an inherited")
    expect_error(signReversal("a"), "unable to find an inherited")
})

test_that("mrcccHyperparameters validates its inputs", {
    h <- mrcccHyperparameters()
    expect_s3_class(h, "mrcccHyperparameters")
    expect_null(h$g)
    expect_equal(mrcccHyperparameters(g = 50)$g, 50)
    expect_error(mrcccHyperparameters(nu1 = 1), "'nu1' must lie")
    expect_error(mrcccHyperparameters(aSigma = 0), "must be positive")
    expect_error(mrcccHyperparameters(aRho = 0),
                 "'aRho' and 'bRho' must be positive")
    expect_error(mrcccHyperparameters(g = -1), "'g' must be positive")
    expect_error(mrcccHyperparameters(ridge = -1e-8), "non-negative")
})
