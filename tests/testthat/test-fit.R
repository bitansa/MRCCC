dat <- makeTriplet(n = 150, seed = 1)
fit <- suppressWarnings(mrccc(dat, nIter = 1500, burnIn = 300, seed = 1))

test_that("mrccc returns a well-formed MRCCCFit", {
    expect_s4_class(fit, "MRCCCFit")
    expect_true(methods::validObject(fit))
    expect_true(pip(fit) >= 0 && pip(fit) <= 1)
    expect_named(pipMCSE(fit), c("mcse", "ess"))
    expect_equal(fitSettings(fit)$nChains, 4L)
    expect_identical(tripletLabels(fit), tripletLabels(dat))
})

test_that("draws are pooled across chains", {
    d <- posteriorDraws(fit)
    expect_equal(nrow(d), 4L * 1200L)
    expect_equal(sort(unique(d$chain)), 1:4)
    expect_true(all(d$gamma %in% c(0, 1)))
    expect_equal(mean(d$gamma), pip(fit))
    expect_true(all(is.finite(as.matrix(d))))
    expect_gt(stats::var(d$loglik), 0)
})

test_that("effect sizes are reported on both scales", {
    raw <- effectSizes(fit, standardized = FALSE)
    std <- effectSizes(fit)
    expect_equal(std$parameter, c("beta_X", "beta_XZ", "beta_Z"))
    s <- inputScale(fit)
    sf <- c(s$sdX / s$sdY, s$sdX * s$sdZ / s$sdY, s$sdZ / s$sdY)
    expect_equal(std$mean, raw$mean * sf)
    dStd <- posteriorDraws(fit, standardized = TRUE)
    expect_equal(mean(dStd$beta_XZ), std$mean[2L])
})

test_that("strong simulated instruments are reported as strong", {
    fi <- instrumentStrength(fit)
    expect_equal(fi$exposure, c("ligand", "receptor"))
    expect_equal(fi$nInstruments,
                 c(ncol(ligandGenotypes(dat)), ncol(receptorGenotypes(dat))))
    expect_true(all(fi$F > 10))
    expect_false(any(fi$weak))
})

test_that("weak instruments trigger a message", {
    inp <- makeInputs(n = 150, seed = 5)
    withr::local_seed(6)
    noise <- matrix(stats::rbinom(150 * 3, 2, 0.4), 150, 3,
                    dimnames = list(names(inp$x), paste0("n", 1:3)))
    weak <- MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                      pathwayActivity = inp$y, ligandGenotypes = noise,
                      receptorGenotypes = inp$H, covariates = inp$V)
    expect_message(
        fw <- suppressWarnings(mrccc(weak, nIter = 400, burnIn = 100,
                                     nChains = 1, seed = 1)),
        "First-stage F below 10 for the ligand")
    expect_true(instrumentStrength(fw)$weak[1L])
})

test_that("a single chain gives NA R-hat and a message", {
    expect_message(
        f1 <- suppressWarnings(mrccc(dat, nIter = 600, burnIn = 100,
                                     nChains = 1, seed = 1)),
        "Only one chain was run")
    expect_true(all(is.na(convergence(f1)$rhat)))
})

test_that("a seed makes the fit reproducible and restores the RNG state", {
    f2 <- suppressWarnings(mrccc(dat, nIter = 1500, burnIn = 300, seed = 1))
    expect_identical(pip(f2), pip(fit))
    expect_identical(posteriorDraws(f2), posteriorDraws(fit))

    withr::local_seed(99)
    before <- stats::runif(1)
    withr::local_seed(99)
    invisible(suppressWarnings(mrccc(dat, nIter = 300, burnIn = 100,
                                     seed = 1)))
    expect_identical(stats::runif(1), before)
})

test_that("shifting covariates leaves the fit unchanged", {
    shifted <- MRCCCData(
        ligandExpression = ligandExpression(dat),
        receptorExpression = receptorExpression(dat),
        pathwayActivity = pathwayActivity(dat),
        ligandGenotypes = ligandGenotypes(dat),
        receptorGenotypes = receptorGenotypes(dat),
        covariates = covariates(dat) + 5)
    fs <- suppressWarnings(mrccc(shifted, nIter = 1500, burnIn = 300,
                                 seed = 1))
    expect_equal(pip(fs), pip(fit), tolerance = 1e-8)
    expect_equal(posteriorDraws(fs)$beta_XZ, posteriorDraws(fit)$beta_XZ,
                 tolerance = 1e-6)
})

test_that("thinning reduces the number of retained draws", {
    ft <- suppressWarnings(suppressMessages(
        mrccc(dat, nIter = 700, burnIn = 100, thin = 3, nChains = 1,
              seed = 1)))
    expect_equal(fitSettings(ft)$nKeep, 200L)
    expect_equal(nrow(posteriorDraws(ft)), 200L)
})

test_that("show, summary and plot methods work", {
    expect_output(show(fit), "MRCCCFit: LIG1 -> REC1")
    sm <- summary(fit)
    expect_s4_class(sm, "MRCCCSummary")
    expect_output(show(sm), "Posterior inclusion probability")
    tab <- as.data.frame(sm)
    expect_s3_class(tab, "data.frame")
    expect_equal(tab$parameter, c("beta_X", "beta_XZ", "beta_Z"))
    expect_s3_class(plot(fit), "ggplot")
    expect_s3_class(plotEffectCurve(fit, standardized = FALSE), "ggplot")
})

test_that("the compiled sampler honours nKeep, verbose and its checks", {
    gibbs <- MRCCC:::mr_ccc_gibbs
    ctr <- function(v) matrix(v - mean(v), ncol = 1)
    ctrCols <- function(M) scale(M, center = TRUE, scale = FALSE)
    X <- ctr(ligandExpression(dat))
    Z <- ctr(receptorExpression(dat))
    Y <- ctr(pathwayActivity(dat))
    Gc <- ctrCols(ligandGenotypes(dat))
    Hc <- ctrCols(receptorGenotypes(dat))
    Vc <- ctrCols(covariates(dat))
    out <- gibbs(X, Z, Y, Gc, Hc, Vc, n_iter = 300, burn_in = 100, thin = 2)
    expect_equal(out$n_keep, 100L)
    expect_length(out$gamma_draws, 100L)
    expect_named(out$prior_scale, c("d_X", "d_XZ", "d_Z"))
    expect_output(
        gibbs(X, Z, Y, Gc, Hc, Vc, n_iter = 1000, burn_in = 100, thin = 10,
              verbose = TRUE),
        "Iteration 1000 / 1000")
    run <- function(save) {
        withr::with_seed(11, gibbs(X, Z, Y, Gc, Hc, Vc, n_iter = 400,
                                   burn_in = 100, thin = 3,
                                   save_first_stage = save))
    }
    off <- run(FALSE)
    on <- run(TRUE)
    expect_identical(on[names(off)], off)
    expect_equal(dim(on$Pi_X_draws), c(off$n_keep, ncol(Gc)))
    expect_error(gibbs(X, Z, Y, Gc, Hc, Vc, n_iter = 300, burn_in = 100,
                       thin = 0), "thin")
    expect_error(gibbs(X, Z, Y, Gc[-1, , drop = FALSE], Hc, Vc,
                       n_iter = 300, burn_in = 100), "same number of rows")
})
