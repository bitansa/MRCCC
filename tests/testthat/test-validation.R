inp <- makeInputs(n = 60, seed = 3)
dat <- dataFrom(inp)

smallFit <- function(..., nIter = 50, burnIn = 10) {
    mrccc(dat, nIter = nIter, burnIn = burnIn, ...)
}

test_that("MCMC settings are validated", {
    expect_error(smallFit(nIter = 100, burnIn = 100),
                 "'nIter' \\(100\\) must be greater than 'burnIn'")
    expect_error(smallFit(nIter = 50.5), "'nIter' must be a positive integer")
    expect_error(smallFit(burnIn = -1), "'burnIn' must be a non-negative")
    expect_error(smallFit(thin = 0), "'thin' must be an integer")
    expect_error(smallFit(nChains = 0), "'nChains' must be an integer")
    expect_error(smallFit(initScale = -1), "'initScale' must be non-negative")
    expect_error(smallFit(initScale = 0, nChains = 2),
                 "'initScale' must be positive when nChains > 1")
    expect_error(smallFit(nIter = 20, burnIn = 10, thin = 40), "retain only")
    expect_error(smallFit(seed = 1.5), "'seed' must be a whole number")
    expect_error(smallFit(verbose = NA), "'verbose' must be TRUE or FALSE")
    expect_error(smallFit(hyper = list(aSigma = 3)),
                 "'hyper' must be created by mrcccHyperparameters")
    for (t in c(0, 1, 1.5)) {
        expect_error(smallFit(pipThreshold = t), "'pipThreshold' must lie")
    }
})

test_that("constant and duplicated instrument columns are rejected", {
    Gconst <- inp$G
    Gconst[, 2] <- 1
    d1 <- MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                    pathwayActivity = inp$y, ligandGenotypes = Gconst,
                    receptorGenotypes = inp$H, covariates = inp$V)
    expect_error(mrccc(d1, nIter = 50, burnIn = 10),
                 "'ligandGenotypes' has constant column")

    Gdup <- inp$G
    Gdup[, 3] <- Gdup[, 1]
    d2 <- MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                    pathwayActivity = inp$y, ligandGenotypes = Gdup,
                    receptorGenotypes = inp$H, covariates = inp$V)
    expect_error(mrccc(d2, nIter = 50, burnIn = 10), "duplicated column")

    Vconst <- cbind(inp$V, c3 = 1)
    d3 <- MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                    pathwayActivity = inp$y, ligandGenotypes = inp$G,
                    receptorGenotypes = inp$H, covariates = Vconst)
    expect_error(mrccc(d3, nIter = 50, burnIn = 10),
                 "already includes an intercept")
})

test_that("rank-deficient designs are rejected after centring", {
    Gbad <- cbind(inp$G, flip = 2 - inp$G[, 1])
    d <- MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                   pathwayActivity = inp$y, ligandGenotypes = Gbad,
                   receptorGenotypes = inp$H, covariates = inp$V)
    expect_error(mrccc(d, nIter = 50, burnIn = 10), "rank deficient")
})

test_that("too few donors for the number of columns is rejected", {
    idx <- seq_len(12)
    d <- MRCCCData(ligandExpression = inp$x[idx],
                   receptorExpression = inp$z[idx],
                   pathwayActivity = inp$y[idx],
                   ligandGenotypes = inp$G[idx, ],
                   receptorGenotypes = inp$H[idx, ],
                   covariates = inp$V[idx, ])
    expect_error(mrccc(d, nIter = 50, burnIn = 10),
                 "Too few donors for the model")
})

test_that("mrccc requires an MRCCCData object", {
    expect_error(mrccc(list()), "unable to find an inherited method")
})
