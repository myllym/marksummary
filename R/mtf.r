# Mark test functions.
#
# If you would like to add support for a new mark test function, e.g. 'xyz',
# do the following:
# - Create a function named create_mtf_func_xyz. Follow examples.
# - Add the calculation of the scaling coefficient into mark_distr_stats.
#   Follow examples.
# - If the xyz mark test function can only produce nonnegative values, add
#   'K_xyz' into the appropriate places in all_besags_l.
# - Add 'xyz' into the vector called available in check_mtf.
# - If suitable, add 'xyz' into the default argument of the parameter
#   mtf_name in summ_func_random_labelling.

#' Check that the given mark test function names are valid.
check_mtf <- function(mtf_name) {
    available <- c('1', 'm', 'mm', 'gamma', 'gammaAbs', 'mor', 'morAbs')

    matched <- mtf_name %in% available
    n_match <- sum(matched)
    no_match_idx <- which(!matched)
    if (length(no_match_idx) > 0L) {
        stop('The following mark test function names were not recognized: ',
             paste(mtf_name[no_match_idx], collapse = ', '),
             '.')
    }
    if (n_match < 1L) {
        stop('No mark test functions were given.')
    }
}



create_mtf_func_1 <- function(mark_stat_l) {
    return(function(m1, m2) { rep.int(1, length(m1)) })
}

create_mtf_func_m <- function(mark_stat_l) {
    return(function(m1, m2) { m1 })
}

create_mtf_func_mm <- function(mark_stat_l) {
    return(function(m1, m2) { m1 * m2 })
}

create_mtf_func_gamma <- function(mark_stat_l) {
    return(function(m1, m2) {
               mark_diff <- m1 - m2
               0.5 * (mark_diff * mark_diff)
           })
}

create_mtf_func_gammaAbs <- function(mark_stat_l) {
    return(function(m1, m2) { abs(m1 - m2) })
}

create_mtf_func_mor <- function(mark_stat_l) {
    mean_mark <- mark_stat_l[['mean']]
    return(function(m1, m2) { (m1 - mean_mark) * (m2 - mean_mark) })
}

create_mtf_func_morAbs <- function(mark_stat_l) {
    mean_mark <- mark_stat_l[['mean']]
    return(function(m1, m2) { abs((m1 - mean_mark) * (m2 - mean_mark)) })
}

#' Given the names of mark test functions and the required mark statistics,
#' produce mark test functions as closures.
create_mark_test_funcs <- function(mtf_name, mark_stat_l) {
    mtf_func_constructor <- paste('create_mtf_func_', mtf_name, sep = '')
    mtf_func_l <- lapply(mtf_func_constructor, function(constructor) {
                         get(constructor)(mark_stat_l)
                  })
}

#' Adds only those mark statistics on the result list which are needed by
#' the K_f functions to be calculated from the given f.
mark_distr_stats <- function(marks, mtf_name) {
    add_scaling_coeff <- function(res_l, mtf, coeff) {
        res_l[[paste(mtf, '_coeff', sep = '')]] <- coeff
        res_l
    }

    res <- list()

    n_mark <- length(marks)

    if (is.element('1', mtf_name)) {
        res <- add_scaling_coeff(res, '1', 1)
    }

    if (is.element('m', mtf_name)) {
        if (length(res[['mean']]) < 1L) {
            res[['mean']] <- mean(marks)
        }
        res <- add_scaling_coeff(res, 'm', 1 / res[['mean']])
    }

    if (is.element('mm', mtf_name)) {
        if (length(res[['mean']]) < 1L) {
            res[['mean']] <- mean(marks)
        }
        mean_mark <- res[['mean']]
        res <- add_scaling_coeff(res, 'mm', 1 / (mean_mark * mean_mark))
    }

    if (is.element('gamma', mtf_name)) {
        if (length(res[['var']]) < 1L) {
            res[['var']] <- var(marks)
        }
        res <- add_scaling_coeff(res, 'gamma', 1 / res[['var']])
    }

    if (is.element('mor', mtf_name)) {
        if (length(res[['mean']]) < 1L) {
            res[['mean']] <- mean(marks)
        }
        if (length(res[['var']]) < 1L) {
            res[['var']] <- var(marks)
        }
        res <- add_scaling_coeff(res, 'mor', 1 / res[['var']])
    }

    if (is.element('morAbs', mtf_name)) {
        if (length(res[['mean']]) < 1L) {
            res[['mean']] <- mean(marks)
        }
        if (length(res[['mor_abs_scaling']]) < 1L) {
            # == 1 / (n_mark * n_mark) * sum(marks - mean(marks))
            res[['mor_abs_scaling']] <- n_mark * mean(marks - res[['mean']])
        }
        res <- add_scaling_coeff(res, 'morAbs',
                                 1 / res[['mor_abs_scaling']])
    }

    if (is.element('gammaAbs', mtf_name)) {
        if (length(res[['gamma_abs_scaling']]) < 1L) {
            res[['gamma_abs_scaling']] <-
                sum(colSums(abs(outer(marks, marks, '-')))) /
                (n_mark * n_mark)
        }
        res <- add_scaling_coeff(res, 'gammaAbs',
                                 1 / res[['gamma_abs_scaling']])
    }

    if (length(res) < 1L) {
        stop('None of the given mark test function names were recognized.')
    }
    res
}



#' Create an individual weight vector for each mark test function.
#'
#' Create an individual weight vector for each mark test function. Bake in
#' 1 / lambda ^ 2 and 1 / c_f so that all weighing can be handled by one
#' multiplication per pair of points per mark test function.
#' @return A matrix of weights. Each column is for a different mark test
#'   function. Each row is for a different _ordered_ point pair where
#'   (x1, x2) is different from (x2, x1).
specialise_weight <- function(weight_vec, mark_stat, mtf_name,
                              one_per_lambda2) {
    mtf_coeff_vec <- sapply(mark_stat[paste(mtf_name, '_coeff', sep = '')],
                            identity)
    outer(weight_vec, mtf_coeff_vec * one_per_lambda2, '*')
}

create_mtfs_and_weights <- function(marks, mtf_name, edge_corr,
                                    one_per_lambda2) {
    check_mtf(mtf_name)

    mark_distr_stat_l <- mark_distr_stats(marks, mtf_name)
    mtf_func_l <- create_mark_test_funcs(mtf_name, mark_distr_stat_l)

    weight_m <- specialise_weight(edge_corr, mark_distr_stat_l, mtf_name,
                                  one_per_lambda2)

    list(mtf_func_l = mtf_func_l, weight_m = weight_m)
}
