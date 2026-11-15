# Shared test data. testthat sources helper files before the tests.

testCovariates <- c("sex", "age", "PC1", "PC2")

makeTriplet <- function(n = 150, seed = 1, ...) {
    sim <- simulateMRCCC(n = n, seed = seed, ...)
    suppressMessages(
        buildTriplet(sim$sender, sim$receiver, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = sim$pathwayGenes,
                     covariates = testCovariates))
}

# Small aligned inputs for MRCCCData(), with dosage genotypes.
makeInputs <- function(n = 60, seed = 3) {
    withr::local_seed(seed)
    ids <- sprintf("D%03d", seq_len(n))
    G <- matrix(stats::rbinom(n * 3, 2, 0.35), n, 3,
                dimnames = list(ids, paste0("g", 1:3)))
    H <- matrix(stats::rbinom(n * 3, 2, 0.35), n, 3,
                dimnames = list(ids, paste0("h", 1:3)))
    V <- matrix(stats::rnorm(n * 2), n, 2,
                dimnames = list(ids, c("c1", "c2")))
    x <- drop(scale(G) %*% rep(0.6, 3)) + stats::rnorm(n)
    z <- drop(scale(H) %*% rep(0.6, 3)) + stats::rnorm(n)
    y <- 0.3 * x + 0.3 * x * z + stats::rnorm(n)
    list(x = stats::setNames(x, ids), z = stats::setNames(z, ids),
         y = stats::setNames(y, ids), G = G, H = H, V = V)
}

dataFrom <- function(inp) {
    MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
              pathwayActivity = inp$y, ligandGenotypes = inp$G,
              receptorGenotypes = inp$H, covariates = inp$V)
}
