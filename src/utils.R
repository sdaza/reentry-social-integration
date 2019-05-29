# texreg function for brms extraction
library(texreg)

# general brms extract method
extract.brms = function(model, include.r2 = TRUE, include.loo = FALSE, ...) {
    s = summary(model)
    # fixed coefficients
    coefficient.names = names(s$fixed[,1])
    coefficients = s$fixed[,1]
    ci.low = s$fixed[, "l-95% CI"]
    ci.upper = s$fixed[, "u-95% CI"]

    # random
    if ('random' %in% names(s)) {
        r = s$random[[1]]
        random.names = stringr::str_replace_all(rownames(r), stringr::fixed('\\_'), "\\_")
        random.estimates = r[,1]
        random.lower = r[,'l-95% CI']
        random.upper = r[, 'u-95% CI']

        coefficient.names = c(coefficient.names, random.names)
        coefficients = c(coefficients, random.estimates)
        ci.low = c(ci.low, random.lower)
        ci.upper = c(ci.upper, random.upper)
    }

    gof = numeric()
    gof.names = character()
    gof.decimal = logical()
    gof = c(gof, s$nobs)
    gof.names = c(gof.names, "Num.\ obs.")
    gof.decimal = c(gof.decimal, FALSE)

    if ('ngrps' %in% names(s)) {
        for (i in seq_along(s$ngrps)) {
            gof = c(gof, s$ngrps[[i]])
            gof.names = c(gof.names, paste('Num.\ obs. ',
            stringr::str_replace_all(names(s$ngrps[i]), stringr::fixed('_'), '\\_')))
            gof.decimal = c(gof.decimal, FALSE)
        }
    }

    if (include.loo == TRUE) {
        loo = loo::loo(model, reloo=TRUE)
        gof = c(gof, loo$estimates[3,1])
        gof.names = c(gof.names, "LOO Information Criterion")
        gof.decimal = c(gof.decimal, FALSE)
    }

    if (include.r2 == TRUE) {
        r2 = brms::bayes_R2(model)
        gof = c(gof, round(r2[1, 1], 2))
        gof.names = c(gof.names, "Bayes $R^2$")
        gof.decimal = c(gof.decimal, TRUE)
    }

    tr = texreg::createTexreg(coef.names = coefficient.names,
                              coef = coefficients,
                              ci.low = ci.low,
                              ci.up = ci.upper,
                              gof.names = gof.names,
                              gof = gof,
                              gof.decimal = gof.decimal
                              )
    return(tr)
}

setMethod("extract",
          signature = className("brmsfit_multiple", "brms"),
          definition = extract.brms)


# extract function that selects coefficients for multivariate models
extract.brms.select_coeff = function(model, include.r2 = TRUE, include.loo = FALSE,
                                     coeff_pattern = NULL,
                                     iteration = NULL) {
    s = summary(model)

    # fixed coefficients
    select = grepl(coeff_pattern, names(s$fixed[, 1]))
    coefficient.names = names(s$fixed[,1])[select]
    coefficients = s$fixed[, 1][select]
    ci.low = s$fixed[, "l-95% CI"][select]
    ci.upper = s$fixed[, "u-95% CI"][select]

    # rename coefficients
    coefficient.names = gsub(coeff_pattern, '', coefficient.names)

    # random
    if ('random' %in% names(s)) {
        r = s$random[[1]]
        random.names = stringr::str_replace_all(rownames(r), stringr::fixed('\\_'), "\\_")
        random.estimates = r[,1]
        random.lower = r[,'l-95% CI']
        random.upper = r[, 'u-95% CI']

        coefficient.names = c(coefficient.names, random.names)
        coefficients = c(coefficients, random.estimates)
        ci.low = c(ci.low, random.lower)
        ci.upper = c(ci.upper, random.upper)
    }

    gof = numeric()
    gof.names = character()
    gof.decimal = logical()

    if (iteration == 1) {

        gof = c(gof, s$nobs)
        gof.names = c(gof.names, "Num.\ obs.")
        gof.decimal = c(gof.decimal, FALSE)

        if ('ngrps' %in% names(s)) {
            for (i in seq_along(s$ngrps)) {
                gof = c(gof, s$ngrps[[i]])
                gof.names = c(gof.names, paste('Num.\ obs. ',
                stringr::str_replace_all(names(s$ngrps[i]), stringr::fixed('_'), '\\_')))
                gof.decimal = c(gof.decimal, FALSE)
            }
        }

        if (include.loo == TRUE) {
            loo = loo::loo(model, reloo=TRUE)
            gof = c(gof, loo$estimates[3,1])
            gof.names = c(gof.names, "LOO Information Criterion")
            gof.decimal = c(gof.decimal, FALSE)
        }

        if (include.r2 == TRUE) {
            r2 = brms::bayes_R2(model)
            gof = c(gof, round(r2[1, 1], 2))
            gof.names = c(gof.names, "Bayes $R^2$")
            gof.decimal = c(gof.decimal, TRUE)
        }

    }

    tr = texreg::createTexreg(coef.names = coefficient.names,
                              coef = coefficients,
                              ci.low = ci.low,
                              ci.up = ci.upper,
                              gof.names = gof.names,
                              gof = gof,
                              gof.decimal = gof.decimal
                              )
    return(tr)
}

# auxiliary function to create texreg objects
create_texreg_multivariate = function(model, dependent_variables_regex,
                                      rhat_min = 0.90,
                                      rhat_max = 1.05,
                                      include.r2 = FALSE,
                                      include.loo = FALSE) {

    check_convergence_mi(model, low = rhat_min, high = rhat_max)
    texreg_objs = list()
    for (i in seq_along(dependent_variables_regex)) {
        print(paste0('::::: create table number ', i))
        texreg_objs[[i]] = extract.brms.select_coeff(model,
                           coeff_pattern = dependent_variables_regex[i],
                           include.r2 = include.r2,
                           include.loo = include.loo,
                           iteration = i)
    }
    return(texreg_objs)
}

# check convergence multiple imputation
check_convergence_mi = function(model, low=.90, high=1.05) {
    total = sum(model$rhats > high | model$rhats < low)
    if (total > 0) { stop('Convergence problems, please check model!') }
    else { print('Checking model convergence: no problems, go ahead and ignore warnings!') }
}

# end
