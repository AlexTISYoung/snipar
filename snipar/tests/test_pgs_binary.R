#!/usr/bin/env Rscript
# Run with: Rscript snipar/tests/test_pgs_binary.R
# Numerical checks and complete CLI exercises; no Monte Carlo calibration run.
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(numDeriv))

script_argument <- grep("^--file=", commandArgs(), value = TRUE)[1L]
test_path <- normalizePath(sub("^--file=", "", script_argument))
repository_root <- dirname(dirname(dirname(test_path)))
cli_path <- file.path(repository_root, "snipar", "scripts", "pgs.R")
source(cli_path)

assert_close <- function(actual, expected, tolerance, label) {
    if (!identical(dim(actual), dim(expected)) ||
        length(actual) != length(expected) ||
        any(!is.finite(actual)) || any(!is.finite(expected)) ||
        max(abs(actual - expected)) > tolerance) {
        stop(label, " differs by more than ", tolerance)
    }
}

expect_error <- function(expression, pattern) {
    message <- tryCatch({ force(expression); NULL },
                        error = function(error) conditionMessage(error))
    if (is.null(message) || !grepl(pattern, message, ignore.case = TRUE)) {
        stop("Expected an error matching '", pattern, "'; received: ", message)
    }
}

fit_test_model <- function(formula, data, nAGQ = 1L) {
    glmer(formula, data = data, family = binomial(), nAGQ = nAGQ,
          control = glmerControl(optimizer = "bobyqa", tolPwrss = 1e-9,
              optCtrl = list(maxfun = 200000L),
              check.conv.singular = "ignore"))
}

check_parser_guards <- function() {
    arguments <- c("output", "--pgs", "scores.tsv", "--phenofile", "phenotypes.tsv")
    options <- parse_pgs_binary_options(arguments)
    stopifnot(is.null(options$phen), is.null(options$covar), is.null(options$covar_cols))
    expect_error(parse_pgs_binary_options(c(arguments, "--covar_cols", "age")),
                 "requires --covar")
    expect_error(parse_pgs_binary_options(c(arguments, "--unknownoption")),
                 "Unknown option")
    for (value in c("0", "-1", "1.5", "NaN", "garbage")) {
        expect_error(parse_pgs_binary_options(c(arguments, "--maxfun", value)),
                     "maxfun.*positive integer")
    }
    cat("PASS: CLI optional arguments, covariate dependency and invalid-option guards\n")
}

check_table_width <- function() {
    path <- tempfile()
    on.exit(unlink(path))
    writeLines(c("FID IID case_status", "11 01 81 1", "12 02 82 0"), path)
    expect_error(read_pgs_binary_table(path, "NA"), "field count")
    writeLines(c("FID IID case_status", "11 01 1", "12 02"), path)
    expect_error(read_pgs_binary_table(path, "NA"), "field count")
    writeLines(c("FID IID case_status", "", "11\t01\t1", "12 02 0"), path)
    table <- read_pgs_binary_table(path, "NA")
    stopifnot(identical(table$FID, c("11", "12")), identical(table$IID, c("01", "02")))
}

generate_example <- function(n_families = 300L) {
    n <- 2L * n_families
    proband <- rep(rnorm(n_families), each = 2L) * sqrt(.5) +
        rnorm(n) * sqrt(.5)
    innovation <- rep(rnorm(n_families), each = 2L) * sqrt(.5) +
        rnorm(n) * sqrt(.5)
    parental <- .65 * proband + innovation
    intercept <- rep(rnorm(n_families, sd = .8), each = 2L)
    y <- rbinom(n, 1L, plogis(-.5 + .6 * proband + .4 * parental + intercept))
    data.frame(y = y, proband = proband, parental = parental,
               FID = factor(rep(seq_len(n_families), each = 2L)))
}

check_numerics <- function() {
    set.seed(20260908L)
    example <- generate_example()
    fits <- list(
        population = fit_test_model(y ~ proband + (1 | FID), example),
        direct = fit_test_model(y ~ proband + parental + (1 | FID), example)
    )
    influences <- list()
    for (name in names(fits)) {
        fit <- fits[[name]]
        parameters <- c(getME(fit, "theta"), fixef(fit))
        stopifnot(parameters[[1L]] > 1e-4)
        X <- getME(fit, "X")
        y <- getME(fit, "y")
        groups <- droplevels(getME(fit, "flist")[[1L]])
        group_index <- as.integer(groups)
        family_loglik <- function(par) laplace_family_loglik(par, X, y, group_index)
        loglik <- function(par) sum(family_loglik(par))
        assert_close(loglik(parameters), as.numeric(logLik(fit)), 1e-5,
                     paste(name, "Laplace log likelihood"))
        information <- numDeriv::hessian(function(par) -loglik(par), parameters)
        lme4_information <- fit@optinfo$derivs$Hessian / 2
        stopifnot(max(abs(information - lme4_information)) /
                      max(abs(lme4_information)) < .005)
        score <- numDeriv::grad(loglik, parameters)
        stopifnot(max(abs(score) / sqrt(diag(information))) < .005)
        scores <- numDeriv::jacobian(family_loglik, parameters)
        scores <- sweep(scores, 2L, colMeans(scores), FUN = "-")
        expected <- scores %*% solve(information)[, -1L, drop = FALSE]
        actual <- logistic_family_influence(fit)
        stopifnot(identical(colnames(actual$influences), colnames(X)),
                  identical(rownames(actual$influences), levels(groups)))
        assert_close(unname(actual$influences), unname(expected), 2e-6,
                     paste(name, "all fixed-effect family influences"))
        assert_close(unname(actual$influence),
                     unname(actual$influences[, "proband"]), 1e-12,
                     paste(name, "selected proband influence"))
        influences[[name]] <- actual$influences
    }
    result <- logistic_direct_population_ratio(fits$population, fits$direct)
    population <- influences$population
    direct <- influences$direct[rownames(population), , drop = FALSE]
    correction <- nrow(population) / (nrow(population) - 1L)
    assert_close(result$population_vcov, crossprod(population) * correction,
                 1e-12, "Population sandwich covariance")
    assert_close(result$direct_vcov, crossprod(direct) * correction,
                 1e-12, "Direct sandwich covariance")
    assert_close(result$cross_vcov, crossprod(population, direct) * correction,
                 1e-12, "Cross-model sandwich covariance")
    stacked <- cbind(population, direct)
    stopifnot(min(eigen(crossprod(stacked) * correction,
                         symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
    joint <- result$joint_vcov[c("population", "direct"), c("population", "direct")]
    ratio <- result$ratio
    gradient <- c(-ratio$direct / ratio$population^2, 1 / ratio$population)
    expected_se <- sqrt(as.numeric(t(gradient) %*% joint %*% gradient))
    assert_close(ratio$se, expected_se, 1e-12, "Ratio delta-method SE")
    assert_close(ratio$population_se^2, joint["population", "population"],
                 1e-12, "Population ratio marginal variance")
    assert_close(ratio$direct_se^2, joint["direct", "direct"],
                 1e-12, "Direct ratio marginal variance")
    identical_result <- logistic_direct_population_ratio(fits$population, fits$population)
    assert_close(identical_result$ratio$ratio, 1, 1e-12, "Identical-model ratio")
    assert_close(identical_result$ratio$se, 0, 1e-10, "Identical-model ratio SE")

    scaled <- example
    scaled$proband <- 2 * scaled$proband
    scaled_fit <- fit_test_model(y ~ proband + parental + (1 | FID), scaled)
    expect_error(logistic_direct_population_ratio(fits$population, scaled_fit),
                 "same|scale|proband|score")
    subset_fit <- fit_test_model(y ~ proband + parental + (1 | FID), example[-1L, ])
    expect_error(logistic_direct_population_ratio(fits$population, subset_fit),
                 "same|observation|famil")
    non_laplace <- fit_test_model(y ~ proband + (1 | FID), example, nAGQ = 0L)
    expect_error(logistic_family_influence(non_laplace), "nAGQ|Laplace")
    cat("PASS: Laplace likelihood, Hessian, full influence matrices, ratio covariance and guards\n")
}

check_zero_variance <- function() {
    set.seed(20260909L)
    data <- data.frame(y = rep(c(0L, 1L), 100L),
                       FID = factor(rep(seq_len(100L), each = 2L)))
    data$proband <- .5 * (2 * data$y - 1) + rnorm(nrow(data))
    fit <- fit_test_model(y ~ proband + (1 | FID), data)
    stopifnot(getME(fit, "theta") < 1e-5)
    X <- getME(fit, "X")
    y <- getME(fit, "y")
    group <- droplevels(getME(fit, "flist")[[1L]])
    probability <- plogis(as.vector(X %*% fixef(fit)))
    parameters <- c(0, fixef(fit))
    expected_loglik <- sum(dbinom(y, 1L, probability, log = TRUE))
    actual_loglik <- sum(laplace_family_loglik(parameters, X, y, as.integer(group)))
    assert_close(actual_loglik, expected_loglik, 1e-10, "Zero-variance GLM likelihood")
    scores <- rowsum(X * (y - probability), as.integer(group), reorder = TRUE)
    scores <- sweep(scores, 2L, colMeans(scores), FUN = "-")
    information <- crossprod(X, X * (probability * (1 - probability)))
    expected <- scores %*% solve(information)
    actual <- logistic_family_influence(fit)
    assert_close(unname(actual$influences), unname(expected), 1e-5,
                 "Zero-variance full GLM cluster influences")
    stopifnot(actual$diagnostics$boundary)
    cat("PASS: zero random-intercept variance agrees with the full GLM sandwich\n")
}

check_pirls_precision <- function() {
    # Mixed singleton/sibling families with a strong PGS effect expose a stale
    # Laplace log-determinant at tolPwrss=1e-9 in lme4 1.1-32.
    set.seed(20261128L)
    sizes <- sample(c(1L, 2L), 600L, replace = TRUE, prob = c(.33, .67))
    groups <- rep(seq_along(sizes), sizes)
    proband <- rnorm(length(sizes), sd = .3)[groups] +
        rnorm(length(groups), sd = .3)
    phenotype <- rbinom(length(groups), 1L,
                        plogis(4.5 * proband + rnorm(length(sizes), sd = 1.5)[groups]))
    data <- data.frame(phenotype = phenotype, proband = proband, FID = factor(groups))
    fit <- fit_pgs_binary_model(data, "proband", 200000L, 1L)
    influence <- logistic_family_influence(fit)
    stopifnot(influence$diagnostics$pirls_tolerance <= 1e-12,
              abs(influence$diagnostics$loglik_difference) <= max(1e-5, nrow(data) * 1e-7),
              influence$diagnostics$max_scaled_score < .005,
              all(is.finite(influence$influences)))
    # The independent likelihood guard must still reject an inconsistent fit.
    inconsistent <- fit
    inconsistent@devcomp$cmp[["dev"]] <- inconsistent@devcomp$cmp[["dev"]] + 1
    expect_error(logistic_family_influence(inconsistent), "Laplace likelihood disagrees")
    cat("PASS: tight PIRLS fit agrees with family likelihood and retains discrepancy guard\n")
}

write_input <- function(data, path) {
    write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

run_cli <- function(arguments, expected_status = 0L) {
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
        shQuote(c(cli_path, arguments)), stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    if ((expected_status == 0L && status != 0L) ||
        (expected_status != 0L && status == 0L)) {
        stop("Unexpected CLI status ", status, ":\n", paste(output, collapse = "\n"))
    }
    output
}

read_output <- function(path) {
    read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
               stringsAsFactors = FALSE, comment.char = "", quote = "")
}

read_covariance <- function(path) {
    as.matrix(read.table(path, header = TRUE, row.names = 1L, sep = "\t",
                         check.names = FALSE, comment.char = "", quote = ""))
}

check_cli_outputs <- function(prefix, expected_n, expected_cases, expected_families,
                              expected_direct_terms) {
    estimates <- read_output(paste0(prefix, ".effects.txt"))
    stopifnot(identical(names(estimates), c("statistic", "est", "SE")),
              !anyDuplicated(estimates$statistic),
              sum(estimates$statistic == "direct_population_ratio") == 1L,
              !any(file.exists(paste0(prefix, c(".1.effects.txt", ".2.effects.txt",
                                                ".direct_population_ratio.txt")))))
    effects <- lapply(1:2, function(generation) {
        rows <- startsWith(estimates$statistic, paste0(generation, "."))
        output <- estimates[rows, ]
        output$term <- substring(output$statistic, 3L)
        output
    })
    covariance <- lapply(1:2, function(generation) {
        read_covariance(paste0(prefix, ".", generation, ".vcov.txt"))
    })
    expected_statistics <- c(paste0("1.", rownames(covariance[[1L]])),
                             paste0("2.", rownames(covariance[[2L]])),
                             "direct_population_ratio")
    stopifnot(identical(estimates$statistic, expected_statistics))
    for (generation in 1:2) {
        stopifnot(identical(effects[[generation]]$term, rownames(covariance[[generation]])),
                  identical(rownames(covariance[[generation]]), colnames(covariance[[generation]])),
                  all(is.finite(effects[[generation]]$est)))
        assert_close(effects[[generation]]$SE^2, unname(diag(covariance[[generation]])),
                     1e-10, paste("CLI generation", generation, "SE/covariance agreement"))
        stopifnot(min(eigen(covariance[[generation]], symmetric = TRUE,
                             only.values = TRUE)$values) > -1e-9)
        model_covariance <- read_covariance(paste0(prefix, ".", generation, ".model_vcov.txt"))
        stopifnot(identical(dimnames(model_covariance), dimnames(covariance[[generation]])),
                  all(is.finite(model_covariance)))
    }
    stopifnot(all(c("proband", "age", "PC1") %in% effects[[1L]]$term),
              all(expected_direct_terms %in% effects[[2L]]$term))
    ratio <- estimates[estimates$statistic == "direct_population_ratio", ]
    stopifnot(all(is.finite(unlist(ratio[c("est", "SE")]))))
    population <- effects[[1L]]$est[match("proband", effects[[1L]]$term)]
    direct <- effects[[2L]]$est[match("proband", effects[[2L]]$term)]
    assert_close(ratio$est, direct / population, 1e-10, "CLI ratio")
    cross <- read_covariance(paste0(prefix, ".1_2.cross_vcov.txt"))
    stopifnot(identical(rownames(cross), effects[[1L]]$term),
              identical(colnames(cross), effects[[2L]]$term))
    joint <- read_covariance(paste0(prefix, ".direct_population.vcov.txt"))
    joint <- joint[c("population", "direct"), c("population", "direct")]
    gradient <- c(-direct / population^2, 1 / population)
    assert_close(ratio$SE, sqrt(as.numeric(t(gradient) %*% joint %*% gradient)),
                 1e-10, "CLI delta-method SE reconstruction")
    assert_close(joint["population", "population"],
                 covariance[[1L]]["proband", "proband"], 1e-10,
                 "CLI joint population variance")
    assert_close(joint["direct", "direct"],
                 covariance[[2L]]["proband", "proband"], 1e-10,
                 "CLI joint direct variance")
    assert_close(joint["population", "direct"], cross["proband", "proband"],
                 1e-10, "CLI cross-model proband covariance")
    combined <- rbind(cbind(covariance[[1L]], cross),
                      cbind(t(cross), covariance[[2L]]))
    stopifnot(min(eigen(combined, symmetric = TRUE, only.values = TRUE)$values) > -1e-9)
    info <- read_output(paste0(prefix, ".model_info.tsv"))
    stopifnot(nrow(info) == 2L,
              all(c("n", "cases", "controls", "n_families", "ratio_status",
                    "ratio_method", "ratio_sampling_correlation") %in% names(info)),
              identical(info$generation, 1:2),
              all(info$n == expected_n), all(info$cases == expected_cases),
              all(info$controls == expected_n - expected_cases),
              all(info$n_families == expected_families),
              all(info$convergence == "OK"), all(info$nAGQ == 1L),
              all(info$optimizer == "bobyqa"),
              all(info$pirls_tolerance <= 1e-12),
              all(is.finite(info$max_scaled_score)), all(info$max_scaled_score <= .05),
              length(unique(info$ratio_status)) == 1L,
              all(info$ratio_status %in% c("ok", "weak_denominator")),
              all(info$ratio_method == "Laplace_family_cluster_sandwich_delta"))
    assert_close(info$sandwich_se, sqrt(diag(joint)), 1e-10,
                 "CLI proband diagnostic SEs")
    expected_correlation <- joint["population", "direct"] / sqrt(prod(diag(joint)))
    assert_close(info$ratio_sampling_correlation, rep(expected_correlation, 2L),
                 1e-10, "CLI ratio sampling-correlation diagnostic")
    expected_status <- if (abs(population) < 1.96 * sqrt(joint[1L, 1L])) "weak_denominator" else "ok"
    stopifnot(all(info$ratio_status == expected_status))
    invisible(ratio)
}

check_cli <- function() {
    work <- tempfile("snipar-binary-cli-")
    dir.create(work)
    on.exit(unlink(work, recursive = TRUE), add = TRUE)
    set.seed(731204L)
    n_families <- 500L
    n <- 2L * n_families
    FID <- rep(sprintf("%05d", seq_len(n_families)), each = 2L)
    IID <- rep(c("001", "002"), n_families)
    paternal <- rep(rnorm(n_families), each = 2L)
    maternal <- rep(rnorm(n_families), each = 2L)
    proband <- .5 * (paternal + maternal) + rnorm(n, sd = sqrt(.5))
    sibling <- .5 * (paternal + maternal) + rnorm(n, sd = sqrt(.5))
    age <- rnorm(n)
    PC1 <- rnorm(n)
    random_intercept <- rep(rnorm(n_families, sd = 1), each = 2L)
    outcome <- rbinom(n, 1L, plogis(-.4 + .7 * proband + .25 * paternal +
                                    .15 * maternal + .25 * age - .2 * PC1 +
                                    random_intercept))
    pgs <- data.frame(FID = FID, IID = IID, proband = proband,
                      paternal = paternal, maternal = maternal)
    phenotypes <- data.frame(FID = FID, IID = IID, unused = 1 - outcome, case = outcome)
    covariates <- data.frame(FID = FID, IID = IID, age = age, PC1 = PC1)
    pgs$paternal[c(3L, 20L)] <- NA_real_
    phenotypes$case[25L] <- NA_integer_
    covariates$PC1[32L] <- NA_real_
    eligible <- complete.cases(pgs, phenotypes$case, covariates)
    expected_n <- sum(eligible)
    expected_cases <- sum(outcome[eligible])
    expected_families <- length(unique(FID[eligible]))
    pgs_path <- file.path(work, "scores.tsv")
    phen_path <- file.path(work, "phenotypes.tsv")
    covar_path <- file.path(work, "covariates.tsv")
    write_input(pgs[sample(n), ], pgs_path)
    write_input(phenotypes[sample(n), ], phen_path)
    write_input(covariates[sample(n), ], covar_path)
    prefix <- file.path(work, "separate")
    arguments <- c(prefix, "--pgs", pgs_path, "--phenofile", phen_path,
                   "--phen", "case", "--covar", covar_path,
                   "--covar_cols", "PC1,age")
    run_cli(arguments)
    check_cli_outputs(prefix, expected_n, expected_cases, expected_families,
                      c("proband", "paternal", "maternal", "PC1", "age"))

    # With parental-only scores, verify the parental-sum option, sibling term,
    # covariate-column defaults, phenotype indexing and common score scaling.
    pgs_sum <- data.frame(FID = FID, IID = IID, proband = proband,
                          parental = pgs$paternal + pgs$maternal, sibling = sibling)
    write_input(pgs_sum[sample(n), ], pgs_path)
    sum_prefix <- file.path(work, "parental_sum")
    run_cli(c(sum_prefix, "--pgs", pgs_path, "--phenofile", phen_path,
              "--phen_index", "2", "--covar", covar_path,
              "--parsum", "--fit_sib", "--scale_pgs"))
    check_cli_outputs(sum_prefix, expected_n, expected_cases, expected_families,
                      c("proband", "parental", "sibling", "PC1", "age"))

    duplicate_path <- file.path(work, "duplicates.tsv")
    write_input(rbind(pgs, pgs[1L, ]), duplicate_path)
    duplicate_prefix <- file.path(work, "duplicate_failure")
    message <- run_cli(c(duplicate_prefix, "--pgs", duplicate_path,
                         "--phenofile", phen_path, "--phen", "case",
                         "--covar", covar_path), expected_status = 1L)
    stopifnot(any(grepl("duplicate|unique", message, ignore.case = TRUE)),
              length(Sys.glob(paste0(duplicate_prefix, ".*"))) == 0L)

    # Invalid inputs must fail before any fits or output files are created.
    write_input(pgs, pgs_path)
    invalid_input <- function(label, phenofile = phen_path, covar = covar_path,
                              extra = character(), pattern) {
        failed_prefix <- file.path(work, label)
        output <- run_cli(c(failed_prefix, "--pgs", pgs_path, "--phenofile", phenofile,
                            "--phen", "case", "--covar", covar, extra),
                          expected_status = 1L)
        stopifnot(any(grepl(pattern, output, ignore.case = TRUE)),
                  length(Sys.glob(paste0(failed_prefix, ".*"))) == 0L)
    }
    duplicated_phenotype <- file.path(work, "duplicate_phenotype.tsv")
    write_input(rbind(phenotypes, phenotypes[1L, ]), duplicated_phenotype)
    invalid_input("duplicate_phenotype_failure", phenofile = duplicated_phenotype,
                  pattern = "duplicate|unique")
    nonbinary_path <- file.path(work, "nonbinary.tsv")
    nonbinary <- phenotypes
    nonbinary$case[100L] <- 2L
    write_input(nonbinary, nonbinary_path)
    invalid_input("nonbinary_failure", phenofile = nonbinary_path,
                  pattern = "0.*1|Bernoulli|binary")
    constant_path <- file.path(work, "constant_covariate.tsv")
    constant <- covariates
    constant$age <- 1
    write_input(constant, constant_path)
    invalid_input("constant_covariate_failure", covar = constant_path,
                  pattern = "constant|rank|variance|varia|singular")
    invalid_input("trailing_comma_failure", extra = c("--covar_cols", "age,"),
                  pattern = "covar_cols")
    cat("PASS: CLI covariates, string IDs, matched samples, all covariance files, parental sums and sibling effects\n")
}

check_parser_guards()
check_table_width()
check_numerics()
check_zero_variance()
check_pirls_precision()
check_cli()
cat("All binary PGS tests passed.\n")
