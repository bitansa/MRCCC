sim <- simulateMRCCC(n = 150, seed = 1)

test_that("simulateMRCCC returns data in the layout of real data", {
    expect_s4_class(sim$sender, "RangedSummarizedExperiment")
    expect_s4_class(sim$receiver, "RangedSummarizedExperiment")
    expect_s4_class(sim$genotypes, "RangedSummarizedExperiment")
    expect_identical(colnames(sim$sender), colnames(sim$genotypes))
    expect_true(all(c("LIG1", "REC1") %in%
                        c(rownames(sim$sender), rownames(sim$receiver))))
    d <- SummarizedExperiment::assay(sim$genotypes, "dosage")
    expect_true(all(d %in% c(0, 1, 2)))
    expect_true(all(testCovariates %in%
                        colnames(SummarizedExperiment::colData(sim$receiver))))
    expect_error(simulateMRCCC(n = 40, seed = 1.5), "whole number")
})

test_that("selectInstruments follows the window, screen and clumping rules", {
    sel <- selectInstruments(sim$sender, "LIG1", sim$genotypes,
                             covariates = testCovariates)
    expect_s4_class(sel, "DataFrame")
    expect_true(nrow(sel) >= 1L && nrow(sel) <= 10L)
    expect_true(all(diff(sel$pvalue) >= 0))

    # Every selected SNP lies within 200 kb of the LIG1 transcription start
    # site (1 Mb on chr1).
    rr <- SummarizedExperiment::rowRanges(sim$genotypes)[sel$snp]
    pos <- methods::slot(IRanges::ranges(rr), "start")
    expect_true(all(abs(pos - 1e6) <= 2e5))

    # No pair of selected SNPs exceeds the LD bound.
    D <- SummarizedExperiment::assay(sim$genotypes)[sel$snp, , drop = FALSE]
    if (nrow(D) > 1L) {
        r2 <- stats::cor(t(D))^2
        expect_true(max(r2[upper.tri(r2)]) <= 0.8 + 1e-12)
    }

    # The screen reproduces lm() with the covariates.
    V <- MRCCC:::.resolveCovariates(testCovariates, sim$sender)
    y <- SummarizedExperiment::assay(sim$sender)["LIG1", ]
    g <- D[1L, ]
    p <- summary(stats::lm(y ~ g + V))$coefficients["g", 4]
    expect_equal(sel$pvalue[1L], p, tolerance = 1e-8)

    meta <- S4Vectors::metadata(sel)
    expect_equal(meta$gene, "LIG1")
    expect_equal(meta$nDonors, 150L)
})

test_that("selectInstruments honours maxSNPs and r2Max", {
    one <- selectInstruments(sim$sender, "LIG1", sim$genotypes,
                             covariates = testCovariates, maxSNPs = 1)
    expect_equal(nrow(one), 1L)
    clumped <- selectInstruments(sim$sender, "LIG1", sim$genotypes,
                                 covariates = testCovariates)
    unclumped <- selectInstruments(sim$sender, "LIG1", sim$genotypes,
                                   covariates = testCovariates, r2Max = 1)
    expect_true(nrow(unclumped) >= nrow(clumped))
    expect_identical(unclumped$snp[1L], clumped$snp[1L])
})

test_that("selectInstruments explains an empty result", {
    # The same SNPs with chromosomes named "1" and "2" instead of "chr1" and
    # "chr2".
    rr <- SummarizedExperiment::rowRanges(sim$genotypes)
    rr2 <- GenomicRanges::GRanges(
        seqnames = ifelse(startsWith(names(rr), "rsL"), "1", "2"),
        ranges = IRanges::ranges(rr))
    names(rr2) <- names(rr)
    geno <- SummarizedExperiment::SummarizedExperiment(
        assays = SummarizedExperiment::assays(sim$genotypes),
        rowRanges = rr2)
    expect_message(
        sel <- selectInstruments(sim$sender, "LIG1", geno,
                                 covariates = testCovariates),
        "name chromosomes")
    expect_equal(nrow(sel), 0L)
    expect_error(selectInstruments(sim$sender, "NOGENE", sim$genotypes),
                 "not a row name")
    flat <- sim$sender
    SummarizedExperiment::assay(flat)["LIG1", ] <- 1
    expect_error(selectInstruments(flat, "LIG1", sim$genotypes),
                 "same in every donor")
})

test_that("factor levels of dropped donors leave no empty covariate", {
    recv <- sim$receiver
    batch <- rep(c("b1", "b2"), length.out = ncol(recv))
    batch[1:5] <- "b3"
    recv$batch <- batch
    dat <- suppressMessages(
        buildTriplet(sim$sender[, -(1:5)], recv, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = sim$pathwayGenes,
                     covariates = c(testCovariates, "batch")))
    expect_false("batchb3" %in% colnames(covariates(dat)))
    expect_true(all(apply(covariates(dat), 2L, stats::var) > 0))
})

test_that("buildTriplet assembles a valid, donor-aligned triplet", {
    dat <- makeTriplet()
    expect_s4_class(dat, "MRCCCData")
    expect_true(methods::validObject(dat))
    expect_length(donors(dat), 150L)
    expect_identical(rownames(ligandGenotypes(dat)), donors(dat))
    expect_true("sexmale" %in% colnames(covariates(dat)))
    expect_equal(unname(tripletLabels(dat)[c("ligand", "receptor")]),
                 c("LIG1", "REC1"))
    pwExpr <- SummarizedExperiment::assay(sim$receiver)[sim$pathwayGenes, ]
    expect_equal(unname(pathwayActivity(dat)), unname(colMeans(pwExpr)))
    expect_output(show(dat), "MRCCCData")
})

test_that("buildTriplet matches donors by name, not position", {
    dat <- makeTriplet()
    shuffled <- sim$receiver[, rev(colnames(sim$receiver))]
    dat2 <- suppressMessages(
        buildTriplet(sim$sender, shuffled, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = sim$pathwayGenes,
                     covariates = testCovariates))
    expect_identical(receptorExpression(dat2), receptorExpression(dat))
    expect_identical(pathwayActivity(dat2), pathwayActivity(dat))

    fewer <- sim$receiver[, -(1:5)]
    expect_message(
        dat3 <- buildTriplet(sim$sender, fewer, sim$genotypes,
                             ligand = "LIG1", receptor = "REC1",
                             pathway = sim$pathwayGenes,
                             covariates = testCovariates),
        "donors shared")
    expect_length(donors(dat3), 145L)
})

test_that("buildTriplet accepts explicit instruments and precomputed scores", {
    sel <- selectInstruments(sim$receiver, "REC1", sim$genotypes,
                             covariates = testCovariates)
    score <- stats::setNames(stats::rnorm(150), colnames(sim$receiver))
    dat <- suppressMessages(
        buildTriplet(sim$sender, sim$receiver, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1", pathway = score,
                     covariates = testCovariates,
                     ligandSNPs = sim$truth$ligandCausal,
                     receptorSNPs = sel))
    expect_identical(colnames(ligandGenotypes(dat)), sim$truth$ligandCausal)
    expect_identical(colnames(receptorGenotypes(dat)), sel$snp)
    expect_equal(unname(pathwayActivity(dat)), unname(score))
})

test_that("the pc1 pathway score is oriented with mean expression", {
    dat <- suppressMessages(
        buildTriplet(sim$sender, sim$receiver, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = sim$pathwayGenes,
                     covariates = testCovariates, pathwayMethod = "pc1"))
    pwExpr <- SummarizedExperiment::assay(sim$receiver)[sim$pathwayGenes, ]
    expect_gt(stats::cor(pathwayActivity(dat), colMeans(pwExpr)), 0)
})

test_that("buildTriplet reports unusable inputs", {
    expect_error(
        buildTriplet(sim$sender, sim$receiver, sim$genotypes,
                     ligand = "NOPE", receptor = "REC1",
                     pathway = sim$pathwayGenes, covariates = "sex"),
        "is not a row name of 'sender'")
    expect_error(
        buildTriplet(sim$sender, sim$receiver, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = sim$pathwayGenes, covariates = "height"),
        "not a column of colData")
    expect_error(
        buildTriplet(sim$sender, sim$receiver, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = sim$pathwayGenes),
        "'covariates' is required")
    expect_error(
        buildTriplet(sim$sender, sim$receiver, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = c("X1", "X2"), covariates = "sex"),
        "None of the pathway genes")
})

test_that("MRCCCData aligns donors by name", {
    inp <- makeInputs()
    dat <- dataFrom(inp)
    rev_ <- rev(names(inp$z))
    dat2 <- MRCCCData(ligandExpression = inp$x,
                      receptorExpression = inp$z[rev_],
                      pathwayActivity = inp$y,
                      ligandGenotypes = inp$G[rev_, ],
                      receptorGenotypes = inp$H,
                      covariates = inp$V[rev_, ])
    expect_identical(dat2, dat)
})

test_that("MRCCCData rejects misaligned or invalid inputs", {
    inp <- makeInputs()
    expect_error(
        MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z[-1],
                  pathwayActivity = inp$y, ligandGenotypes = inp$G,
                  receptorGenotypes = inp$H, covariates = inp$V),
        "lacks 1 of the 60 donors")
    expect_error(
        MRCCCData(ligandExpression = unname(inp$x),
                  receptorExpression = inp$z, pathwayActivity = inp$y,
                  ligandGenotypes = inp$G, receptorGenotypes = inp$H,
                  covariates = inp$V),
        "Name the donors of every component")
    Gbad <- inp$G
    Gbad[1, 1] <- 3
    expect_error(
        MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                  pathwayActivity = inp$y, ligandGenotypes = Gbad,
                  receptorGenotypes = inp$H, covariates = inp$V),
        "allele dosages")
    ybad <- inp$y
    ybad[2] <- NA
    expect_error(
        MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                  pathwayActivity = ybad, ligandGenotypes = inp$G,
                  receptorGenotypes = inp$H, covariates = inp$V),
        "contains NA")
    Vna <- inp$V
    Vna[3, 1] <- NA
    expect_error(
        MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                  pathwayActivity = inp$y, ligandGenotypes = inp$G,
                  receptorGenotypes = inp$H, covariates = Vna),
        "'covariates' contains missing values")
    Vdup <- inp$V
    rownames(Vdup)[2] <- rownames(Vdup)[1]
    expect_error(
        MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                  pathwayActivity = inp$y, ligandGenotypes = inp$G,
                  receptorGenotypes = inp$H, covariates = Vdup),
        "Donor names are not unique in covariates")
    expect_error(
        MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                  pathwayActivity = inp$y, ligandGenotypes = inp$G,
                  receptorGenotypes = inp$H, covariates = inp$V,
                  labels = c(gene = "A")),
        "'labels' must be named")
    idx <- 1:8
    expect_error(
        MRCCCData(ligandExpression = inp$x[idx],
                  receptorExpression = inp$z[idx],
                  pathwayActivity = inp$y[idx],
                  ligandGenotypes = inp$G[idx, ],
                  receptorGenotypes = inp$H[idx, ],
                  covariates = inp$V[idx, ]),
        "at least 10 donors")
})

test_that("unnamed inputs are taken in order", {
    inp <- makeInputs()
    strip <- function(v) unname(v)
    dat <- MRCCCData(ligandExpression = strip(inp$x),
                     receptorExpression = strip(inp$z),
                     pathwayActivity = strip(inp$y),
                     ligandGenotypes = unname(inp$G),
                     receptorGenotypes = unname(inp$H),
                     covariates = unname(inp$V))
    expect_identical(donors(dat), paste0("donor", 1:60))
})

test_that("data frame covariates are expanded to indicators", {
    inp <- makeInputs()
    df <- data.frame(sex = rep(c("female", "male"), 30),
                     age = seq(20, 79), row.names = names(inp$x))
    dat <- MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                     pathwayActivity = inp$y, ligandGenotypes = inp$G,
                     receptorGenotypes = inp$H, covariates = df)
    expect_identical(colnames(covariates(dat)), c("sexmale", "age"))
    df$sex <- "female"
    expect_error(
        MRCCCData(ligandExpression = inp$x, receptorExpression = inp$z,
                  pathwayActivity = inp$y, ligandGenotypes = inp$G,
                  receptorGenotypes = inp$H, covariates = df),
        "single value")
})

test_that("buildTriplet reads a precomputed pathway score from colData", {
    rec <- sim$receiver
    score <- stats::rnorm(ncol(rec))
    cd <- SummarizedExperiment::colData(rec)
    cd$pathwayScore <- score
    SummarizedExperiment::colData(rec) <- cd
    dat <- suppressMessages(
        buildTriplet(sim$sender, rec, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = "pathwayScore",
                     covariates = testCovariates))
    expect_equal(unname(pathwayActivity(dat)), score)
    expect_match(S4Vectors::metadata(dat)$pathwayMethod, "colData column")
    cd$pathwayText <- "a"
    SummarizedExperiment::colData(rec) <- cd
    expect_error(
        buildTriplet(sim$sender, rec, sim$genotypes,
                     ligand = "LIG1", receptor = "REC1",
                     pathway = "pathwayText", covariates = testCovariates),
        "must be numeric")
})

test_that("the instrument tables are available through metadata()", {
    dat <- makeTriplet()
    md <- S4Vectors::metadata(dat)
    expect_s4_class(md$ligandInstruments, "DataFrame")
    expect_identical(as.character(md$ligandInstruments$snp),
                     colnames(ligandGenotypes(dat)))
    expect_identical(md$pathwayMethod, "mean")
})

test_that("the OneK1K example has the documented layout", {
    ex <- onek1kExample()
    expect_named(ex, c("sender", "receiver", "genotypes", "triplets"))
    expect_s4_class(ex$genotypes, "RangedSummarizedExperiment")
    expect_identical(colnames(ex$sender), colnames(ex$receiver))
    cd <- SummarizedExperiment::colData(ex$receiver)
    expect_true(all(c("sex", "ageBin", "PC1", "PC2", "PC3",
                      "interleukinPC1", "interleukinMean") %in%
                        colnames(cd)))
    expect_true(all(ex$triplets$ligand %in% rownames(ex$sender)))
    expect_true(all(ex$triplets$receptor %in% rownames(ex$receiver)))
})
