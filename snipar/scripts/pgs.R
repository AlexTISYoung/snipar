#!/usr/bin/env Rscript
# Binary-outcome family PGS regression and direct/population effect ratios.
# Run with --help. This script is self-contained and can also be sourced in R.

# Cross-model covariance for two random-intercept logistic mixed models.
# Each cluster contributes its Laplace-approximated marginal log likelihood.
# Differentiating in BOTH the fixed effects and random-intercept SD includes
# uncertainty from estimating the variance component.  The observed joint
# Hessian supplies the bread; stacked family scores supply the sandwich meat.

laplace_family_loglik <- function(parameters, X, y, family_index) {
    theta <- parameters[[1L]]
    beta <- parameters[-1L]
    eta_fixed <- as.numeric(X %*% beta)
    family_index <- as.integer(family_index)
    n_families <- max(family_index)
    family_sum <- function(values) {
        as.numeric(rowsum(values, family_index, reorder = TRUE))
    }
    u <- numeric(n_families)
    converged <- FALSE
    for (iteration in seq_len(100L)) {
        probability <- plogis(eta_fixed + theta * u[family_index])
        curvature <- 1 + theta^2 * family_sum(probability * (1 - probability))
        step <- (theta * family_sum(y - probability) - u) / curvature
        # The conditional objective is strictly concave; capping Newton steps
        # keeps large starting displacements stable for separated families.
        step <- pmax(-2, pmin(2, step))
        u <- u + step
        if (max(abs(step)) < 1e-11) {
            converged <- TRUE
            break
        }
    }
    if (!converged) stop("Conditional random-effect modes failed to converge")
    eta <- eta_fixed + theta * u[family_index]
    probability <- plogis(eta)
    curvature <- 1 + theta^2 * family_sum(probability * (1 - probability))
    # log(1 + exp(eta)), evaluated without overflow.
    log_normalizer <- pmax(eta, 0) + log1p(exp(-abs(eta)))
    family_sum(y * eta - log_normalizer) - 0.5 * u^2 - 0.5 * log(curvature)
}


logistic_family_influence <- function(fit, term = "proband") {
    if (!inherits(fit, "glmerMod") || family(fit)$family != "binomial" ||
        family(fit)$link != "logit") {
        stop("Expected a binomial-logit glmer model")
    }
    if (lme4::getME(fit, "devcomp")$dims[["nAGQ"]] != 1L) {
        stop("Laplace score covariance requires glmer nAGQ=1")
    }
    theta <- lme4::getME(fit, "theta")
    Zt <- lme4::getME(fit, "Zt")
    if (length(theta) != 1L || any(Matrix::colSums(Zt != 0) != 1) ||
        any(Zt@x != 1)) {
        stop("Exactly one unweighted random intercept per observation is required")
    }
    if (!is.null(model.offset(model.frame(fit))) ||
        any(stats::weights(fit, type = "prior") != 1)) {
        stop("Offset and weighted logistic models are not supported")
    }
    X <- lme4::getME(fit, "X")
    y <- lme4::getME(fit, "y")
    if (any(!y %in% c(0, 1))) stop("Individual Bernoulli outcomes are required")
    group <- droplevels(lme4::getME(fit, "flist")[[1L]])
    family_index <- as.integer(group)
    n_families <- nlevels(group)
    if (n_families < 3L) stop("At least three independent families are required")
    if (!(term %in% colnames(X))) stop("Coefficient not found: ", term)
    parameters <- c(theta, lme4::fixef(fit))
    loglik <- laplace_family_loglik(parameters, X, y, family_index)
    fitted_loglik <- as.numeric(logLik(fit))
    discrepancy <- sum(loglik) - fitted_loglik
    if (abs(discrepancy) > max(1e-5, length(y) * 1e-7)) {
        stop("Family Laplace likelihood disagrees with lme4 by ", discrepancy)
    }
    boundary <- as.numeric(theta) <= 1e-4
    keep <- if (boundary) seq.int(2L, length(parameters)) else seq_along(parameters)
    scores <- matrix(NA_real_, nrow = n_families, ncol = length(keep))
    for (column in seq_along(keep)) {
        index <- keep[[column]]
        step <- .Machine$double.eps^(1 / 3) * (1 + abs(parameters[[index]]))
        upper <- lower <- parameters
        upper[[index]] <- upper[[index]] + step
        lower[[index]] <- lower[[index]] - step
        scores[, column] <- (
            laplace_family_loglik(upper, X, y, family_index) -
                laplace_family_loglik(lower, X, y, family_index)
        ) / (2 * step)
    }
    hessian <- fit@optinfo$derivs$Hessian
    if (is.null(hessian) || !all(dim(hessian) == length(parameters))) {
        stop("The joint lme4 Hessian for theta and fixed effects is unavailable")
    }
    # lme4 stores the Hessian of -2 log L, not -log L.
    information <- hessian[keep, keep, drop = FALSE] / 2
    inverse_information <- chol2inv(chol(information))
    coefficient_indices <- match(seq_len(ncol(X)) + 1L, keep)
    scaled_score <- colSums(scores) / sqrt(diag(information))
    if (max(abs(scaled_score)) > 0.05) {
        stop("Laplace score is not sufficiently close to zero at fitted parameters: ",
             max(abs(scaled_score)))
    }
    scores <- sweep(scores, 2L, colMeans(scores), FUN = "-")
    influences <- scores %*% inverse_information[, coefficient_indices, drop = FALSE]
    dimnames(influences) <- list(levels(group), colnames(X))
    influence <- influences[, term]
    names(influence) <- levels(group)
    list(
        influence = influence,
        influences = influences,
        estimate = as.numeric(lme4::fixef(fit)[[term]]),
        model_se = sqrt(as.matrix(vcov(fit))[term, term]),
        n = length(y),
        response = y,
        row_names = rownames(model.frame(fit)),
        groups = as.character(group),
        diagnostics = data.frame(
            n_families = n_families,
            random_intercept_sd = as.numeric(theta),
            boundary = boundary,
            laplace_loglik = sum(loglik),
            lme4_loglik = fitted_loglik,
            loglik_difference = discrepancy,
            max_scaled_score = max(abs(scaled_score)),
            hessian_condition = kappa(information),
            model_based_se = sqrt(as.matrix(vcov(fit))[term, term]),
            stringsAsFactors = FALSE
        )
    )
}


logistic_direct_population_ratio <- function(fit_1, fit_2, term = "proband") {
    population <- logistic_family_influence(fit_1, term)
    direct <- logistic_family_influence(fit_2, term)
    if (!identical(population$row_names, direct$row_names) ||
        !identical(population$response, direct$response) ||
        !identical(population$groups, direct$groups) ||
        !identical(as.numeric(lme4::getME(fit_1, "X")[, term]),
                   as.numeric(lme4::getME(fit_2, "X")[, term]))) {
        stop("Both logistic models must use the same observations, families and proband scores")
    }
    family_ids <- rownames(population$influences)
    if (!setequal(family_ids, rownames(direct$influences))) {
        stop("The two models do not contain identical families")
    }
    influence_1 <- population$influences
    influence_2 <- direct$influences[family_ids, , drop = FALSE]
    n_families <- nrow(influence_1)
    correction <- n_families / (n_families - 1)
    covariance_1 <- crossprod(influence_1) * correction
    covariance_2 <- crossprod(influence_2) * correction
    cross_covariance <- crossprod(influence_1, influence_2) * correction
    influence <- cbind(population = influence_1[, term], direct = influence_2[, term])
    joint <- crossprod(influence) * correction
    correlation <- joint[1L, 2L] / sqrt(prod(diag(joint)))
    ratio <- se <- NA_real_
    status <- "ok"
    if (any(!is.finite(joint)) || any(diag(joint) <= 0)) {
        status <- "invalid_covariance"
    }
    if (!is.finite(population$estimate) || population$estimate == 0) {
        status <- "undefined_zero_denominator"
    } else {
        ratio <- direct$estimate / population$estimate
        if (status == "ok") {
            gradient <- c(-ratio / population$estimate, 1 / population$estimate)
            # Evaluate via influences to avoid cancellation for nearly equal fits.
            variance <- sum(as.numeric(influence %*% gradient)^2) * correction
            if (!is.finite(variance)) stop("Non-finite delta-method variance")
            se <- sqrt(variance)
            if (abs(population$estimate) < 1.96 * sqrt(joint[1L, 1L])) {
                status <- "weak_denominator"
            }
        }
    }
    if (status != "ok") warning("Direct/population PGS ratio: ", status)
    diagnostics <- rbind(
        cbind(generation = 1L, population$diagnostics),
        cbind(generation = 2L, direct$diagnostics)
    )
    diagnostics$sandwich_se <- sqrt(diag(joint))
    list(
        ratio = data.frame(
            ratio = ratio, se = se,
            direct = direct$estimate, population = population$estimate,
            direct_se = sqrt(joint[2L, 2L]), population_se = sqrt(joint[1L, 1L]),
            covariance = joint[1L, 2L], sampling_correlation = correlation,
            n = direct$n, n_families = n_families, status = status,
            method = "Laplace_family_cluster_sandwich_delta",
            stringsAsFactors = FALSE
        ),
        joint_vcov = joint,
        population_vcov = covariance_1,
        direct_vcov = covariance_2,
        cross_vcov = cross_covariance,
        diagnostics = diagnostics
    )
}


pgs_binary_usage <- function() {
    cat(paste0(
        "Usage: Rscript pgs.R OUT --pgs FILE --phenofile FILE [options]\n\n",
        "Fit 1- and 2-generation logistic models with a family random intercept.\n",
        "Both models use the same complete cases; the ratio and SE are automatic.\n\n",
        "Inputs (whitespace-delimited, with headers):\n",
        "  --pgs FILE          FID IID proband paternal maternal PGS columns\n",
        "  --phenofile FILE    FID IID and binary 0/1 phenotype column(s)\n",
        "  --phen NAME         Phenotype column name (default: first after FID/IID)\n",
        "  --phen_index INT    1-based phenotype column index, excluding FID/IID\n",
        "  --covar FILE        Optional FID IID and numeric covariate columns\n",
        "  --covar_cols LIST   Comma-separated covariates (default: all in FILE)\n",
        "  --parsum            Use parental PGS, or sum paternal and maternal PGS\n",
        "  --fit_sib           Include sibling PGS in the 2-generation model\n",
        "  --scale_pgs         Divide all PGS by complete-case proband PGS SD\n",
        "  --missing_char STR  Additional missing-value token (default: NA)\n",
        "  --maxfun INT        Optimizer evaluation limit (default: 200000)\n",
        "  --help              Show this help\n\n",
        "OUT is the output prefix (--out OUT is also accepted). R package lme4\n",
        "is required. Coefficients are log odds; ratio is the 2-gen/1-gen\n",
        "proband coefficient ratio. Sampling covariance files and the ratio SE\n",
        "use a family-cluster sandwich of Laplace scores (nAGQ=1).\n"
    ))
}


parse_pgs_binary_options <- function(arguments) {
    options <- list(out = NULL, pgs = NULL, phenofile = NULL, phen = NULL,
                    covar = NULL, covar_cols = NULL,
                    phen_index = 1L, maxfun = 200000L, missing_char = "NA",
                    parsum = FALSE, fit_sib = FALSE, scale_pgs = FALSE)
    flags <- c("parsum", "fit_sib", "scale_pgs")
    valued <- c("out", "pgs", "phenofile", "phen", "phen_index", "covar",
                "covar_cols", "missing_char", "maxfun")
    seen <- character()
    index <- 1L
    while (index <= length(arguments)) {
        argument <- arguments[[index]]
        if (argument %in% c("--help", "-h")) {
            pgs_binary_usage()
            return(NULL)
        }
        if (!startsWith(argument, "--")) {
            key <- "out"
            value <- argument
        } else {
            key <- substring(argument, 3L)
            if (key %in% flags) {
                value <- TRUE
            } else if (key %in% valued) {
                index <- index + 1L
                if (index > length(arguments) || startsWith(arguments[[index]], "--")) {
                    stop("Missing value for --", key)
                }
                value <- arguments[[index]]
            } else {
                stop("Unknown option: ", argument)
            }
        }
        if (key %in% seen) stop("Option supplied more than once: ", key)
        seen <- c(seen, key)
        options[[key]] <- value
        index <- index + 1L
    }
    for (key in c("out", "pgs", "phenofile")) {
        if (is.null(options[[key]]) || !nzchar(options[[key]])) {
            stop("Missing required argument: ", key, "; use --help for usage")
        }
    }
    if (all(c("phen", "phen_index") %in% seen)) {
        stop("Specify only one of --phen and --phen_index")
    }
    for (key in c("phen_index", "maxfun")) {
        value <- suppressWarnings(as.numeric(options[[key]]))
        if (length(value) != 1L || !is.finite(value) || value < 1 ||
            value != floor(value) || value > .Machine$integer.max) {
            stop("--", key, " must be a positive integer")
        }
        options[[key]] <- as.integer(value)
    }
    if (!is.null(options$covar_cols) && is.null(options$covar)) {
        stop("--covar_cols requires --covar")
    }
    options
}


read_pgs_binary_table <- function(path, missing_char) {
    data <- read.table(path, header = TRUE, colClasses = "character",
                       stringsAsFactors = FALSE, check.names = FALSE,
                       comment.char = "", quote = "", row.names = NULL,
                       na.strings = unique(c("NA", "NaN", "nan", "N/A", missing_char)))
    if (anyDuplicated(names(data))) stop(path, " contains duplicate column names")
    if (!all(c("FID", "IID") %in% names(data))) {
        stop(path, " must contain FID and IID columns")
    }
    if (nrow(data) == 0L) stop(path, " contains no observations")
    if (anyNA(data$FID) || anyNA(data$IID) ||
        any(!nzchar(data$FID)) || any(!nzchar(data$IID))) {
        stop(path, " contains missing FID or IID values")
    }
    if (anyDuplicated(data[c("FID", "IID")])) {
        stop(path, " contains duplicate FID/IID pairs")
    }
    data
}


pgs_binary_key <- function(data) paste(data$FID, data$IID, sep = "\r")


pgs_binary_numeric <- function(values, name) {
    numeric <- suppressWarnings(as.numeric(values))
    if (any(!is.na(values) & is.na(numeric))) {
        stop("Column ", name, " must contain numeric values or declared missing values")
    }
    numeric
}


prepare_pgs_binary_data <- function(options) {
    scores <- read_pgs_binary_table(options$pgs, options$missing_char)
    parents <- if (options$parsum) "parental" else c("paternal", "maternal")
    if (options$parsum && !("parental" %in% names(scores))) {
        if (!all(c("paternal", "maternal") %in% names(scores))) {
            stop("--parsum requires parental, or both paternal and maternal PGS columns")
        }
        scores$parental <- pgs_binary_numeric(scores$paternal, "paternal") +
            pgs_binary_numeric(scores$maternal, "maternal")
    }
    score_names <- c("proband", if (options$fit_sib) "sibling", parents)
    absent <- setdiff(score_names, names(scores))
    if (length(absent)) stop("PGS file is missing columns: ", paste(absent, collapse = ", "))
    analysis <- scores[, c("FID", "IID", score_names), drop = FALSE]
    phenotypes <- read_pgs_binary_table(options$phenofile, options$missing_char)
    phenotype_names <- setdiff(names(phenotypes), c("FID", "IID"))
    phenotype <- options$phen
    if (is.null(phenotype)) {
        if (options$phen_index > length(phenotype_names)) {
            stop("--phen_index exceeds the number of phenotype columns")
        }
        phenotype <- phenotype_names[[options$phen_index]]
    }
    if (!(phenotype %in% phenotype_names)) stop("Phenotype column not found: ", phenotype)
    phenotype_matches <- match(pgs_binary_key(analysis), pgs_binary_key(phenotypes))
    analysis$phenotype <- phenotypes[[phenotype]][phenotype_matches]
    covariate_names <- character()
    if (!is.null(options$covar)) {
        covariates <- read_pgs_binary_table(options$covar, options$missing_char)
        covariate_names <- setdiff(names(covariates), c("FID", "IID"))
        if (!is.null(options$covar_cols)) {
            covariate_names <- trimws(strsplit(options$covar_cols, ",", fixed = TRUE)[[1L]])
            if (!length(covariate_names) || any(!nzchar(covariate_names)) ||
                endsWith(trimws(options$covar_cols), ",") ||
                anyDuplicated(covariate_names) ||
                !all(covariate_names %in% setdiff(names(covariates), c("FID", "IID")))) {
                stop("--covar_cols must list distinct covariate columns present in --covar")
            }
        }
        reserved <- c("FID", "IID", "phenotype", "proband", "paternal", "maternal",
                      "parental", "sibling", "(Intercept)")
        if (any(covariate_names %in% reserved) || any(grepl("`", covariate_names, fixed = TRUE))) {
            stop("Covariate names must not reuse model/ID column names or contain backticks")
        }
        covariate_matches <- match(pgs_binary_key(analysis), pgs_binary_key(covariates))
        for (name in covariate_names) analysis[[name]] <- covariates[[name]][covariate_matches]
    }
    numeric_names <- c("phenotype", score_names, covariate_names)
    for (name in numeric_names) analysis[[name]] <- pgs_binary_numeric(analysis[[name]], name)
    valid <- Reduce(`&`, lapply(analysis[numeric_names], is.finite))
    input_n <- nrow(analysis)
    analysis <- analysis[valid, , drop = FALSE]
    if (!nrow(analysis) || !identical(sort(unique(analysis$phenotype)), c(0, 1))) {
        stop("The complete-case phenotype must contain both 0 and 1 and no other values")
    }
    score_sd <- sd(analysis$proband)
    if (!is.finite(score_sd) || score_sd <= 0) stop("Proband PGS must have positive finite SD")
    if (options$scale_pgs) {
        for (name in score_names) analysis[[name]] <- analysis[[name]] / score_sd
    }
    analysis$FID <- factor(analysis$FID)
    families <- table(analysis$FID)
    if (length(families) < 3L || !any(families > 1L)) {
        stop("At least three families and at least one multi-person family are required")
    }
    list(data = analysis, phenotype = phenotype, covariates = covariate_names,
         score_names = score_names, parents = parents, input_n = input_n,
         excluded_n = input_n - nrow(analysis), proband_sd = score_sd)
}


fit_pgs_binary_model <- function(analysis, terms, maxfun, generation) {
    quoted <- paste0("`", terms, "`")
    fixed_formula <- reformulate(quoted, response = "phenotype")
    design <- model.matrix(fixed_formula, analysis)
    if (qr(design)$rank < ncol(design)) {
        stop(generation, "-generation design is rank deficient; check constant or collinear PGS/covariates")
    }
    formula <- as.formula(paste("phenotype ~", paste(quoted, collapse = " + "), "+ (1 | FID)"))
    message("Fitting ", generation, "-generation logistic mixed model")
    fit <- lme4::glmer(
        formula, data = analysis, family = binomial(link = "logit"), nAGQ = 1L,
        control = lme4::glmerControl(
            optimizer = "bobyqa", optCtrl = list(maxfun = maxfun),
            calc.derivs = TRUE, tolPwrss = 1e-9, check.conv.singular = "ignore"
        )
    )
    optimizer_code <- fit@optinfo$conv$opt
    messages <- fit@optinfo$conv$lme4$messages
    if ((length(optimizer_code) && any(optimizer_code != 0)) || length(messages)) {
        stop(generation, "-generation model failed convergence checks: ",
             paste(c(optimizer_code, messages), collapse = "; "))
    }
    fit
}


write_pgs_binary_matrix <- function(matrix, path) {
    output <- data.frame(term = rownames(matrix), matrix, check.names = FALSE)
    write.table(output, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}


pgs_binary_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
    options <- parse_pgs_binary_options(arguments)
    if (is.null(options)) return(invisible(NULL))
    if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("R package lme4 is required; install it with install.packages('lme4')")
    }
    prepared <- prepare_pgs_binary_data(options)
    analysis <- prepared$data
    common <- c("proband", prepared$covariates)
    full <- c(common, if (options$fit_sib) "sibling", prepared$parents)
    fit_1 <- fit_pgs_binary_model(analysis, common, options$maxfun, 1L)
    fit_2 <- fit_pgs_binary_model(analysis, full, options$maxfun, 2L)
    result <- logistic_direct_population_ratio(fit_1, fit_2)
    # Complete all statistical calculations before writing any result files.
    output_directory <- dirname(options$out)
    if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)
    fits <- list(fit_1, fit_2)
    covariances <- list(result$population_vcov, result$direct_vcov)
    for (generation in seq_len(2L)) {
        prefix <- paste0(options$out, ".", generation)
        fit <- fits[[generation]]
        covariance <- covariances[[generation]]
        estimates <- data.frame(
            term = names(lme4::fixef(fit)), estimate = as.numeric(lme4::fixef(fit)),
            se = sqrt(diag(covariance)), stringsAsFactors = FALSE
        )
        write.table(estimates, paste0(prefix, ".effects.txt"), sep = "\t",
                    row.names = FALSE, quote = FALSE, na = "NA")
        write_pgs_binary_matrix(covariance, paste0(prefix, ".vcov.txt"))
        write_pgs_binary_matrix(as.matrix(vcov(fit)), paste0(prefix, ".model_vcov.txt"))
    }
    write_pgs_binary_matrix(result$cross_vcov, paste0(options$out, ".1_2.cross_vcov.txt"))
    write_pgs_binary_matrix(result$joint_vcov, paste0(options$out, ".direct_population.vcov.txt"))
    write.table(result$ratio, paste0(options$out, ".direct_population_ratio.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    diagnostics <- result$diagnostics
    diagnostics$phenotype <- prepared$phenotype
    diagnostics$n <- nrow(analysis)
    diagnostics$cases <- sum(analysis$phenotype == 1)
    diagnostics$controls <- sum(analysis$phenotype == 0)
    diagnostics$input_n <- prepared$input_n
    diagnostics$excluded_n <- prepared$excluded_n
    diagnostics$proband_sd_before_scaling <- prepared$proband_sd
    diagnostics$pgs_scaled <- options$scale_pgs
    diagnostics$coefficient_scale <- if (options$scale_pgs) "log_odds_per_proband_PGS_SD" else "log_odds_per_PGS_unit"
    diagnostics$sampling_covariance <- "family_cluster_sandwich"
    diagnostics$convergence <- "OK"
    diagnostics$optimizer <- "bobyqa"
    diagnostics$nAGQ <- 1L
    write.table(diagnostics, paste0(options$out, ".model_info.tsv"), sep = "\t",
                row.names = FALSE, quote = FALSE, na = "NA")
    cat("Complete-case sample:", nrow(analysis), "individuals in", nlevels(analysis$FID),
        "families; excluded:", prepared$excluded_n, "\n")
    cat("Direct/population PGS ratio:", result$ratio$ratio, "(SE", result$ratio$se,
        ");", result$ratio$status, "\n")
    cat("Results written with prefix:", options$out, "\n")
    invisible(result)
}


if (sys.nframe() == 0L) {
    tryCatch(pgs_binary_main(), error = function(error) {
        message("Error: ", conditionMessage(error))
        quit(status = 1L)
    })
}
