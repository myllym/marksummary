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
#'
#' @param mtf_name A string giving the name of the mark test function.
#'                 Possible names are '1', 'm', 'mm', 'gamma', 'gammaAbs', 'morAbs'.
check_mtf <- function(mtf_name) {
    available <- c('1', 'm', 'mm', 'gamma', 'gammaAbs', 'morAbs')

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

create_mtf_func_morAbs <- function(mark_stat_l) {
    mean_mark <- mark_stat_l[['mean']]
    return(function(m1, m2) { abs((m1 - mean_mark) * (m2 - mean_mark)) })
}

#' Create mark test functions
#'
#' Given the names of mark test functions and the required mark statistics,
#' produce mark test functions as closures.
#'
#' @param mtf_name The name of mark test function.
#' @param mark_stat_l List of required mark statistics.
create_mark_test_funcs <- function(mtf_name, mark_stat_l) {
    mtf_func_constructor <- paste('create_mtf_func_', mtf_name, sep = '')
    mtf_func_l <- lapply(mtf_func_constructor, function(constructor) {
                         get(constructor)(mark_stat_l)
                  })
}

#' Statistics of the mark distribution
#'
#' Adds only those mark statistics on the result list which are needed by
#' the K_f functions to be calculated from the given f.
#'
#' @param marks A vector of marks of a point pattern.
#' @param mtf_name The name of mark test function.
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

    if (is.element('morAbs', mtf_name)) {
        if (length(res[['mean']]) < 1L) {
            res[['mean']] <- mean(marks)
        }
        if (length(res[['mor_abs_scaling']]) < 1L) {
            diff_abs_sum <- sum(abs(marks - res[['mean']]))
            res[['mor_abs_scaling']] <- 1 / (n_mark * n_mark) *
                                        (diff_abs_sum * diff_abs_sum)
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

    res
}



#' Individualize weight vectors for each mark test function.
#'
#' Create an individual weight vector for each mark test function. The
#' weight of the point pair (x_1, x_2) is
#' w_{1, 2} = e_{1, 2} / (lambda^2 * c_f),
#' where e_{1, 2} is the edge correction term, lambda^2 is the estimate of
#' the intensity squared and c_f is the scaling term for the mark test
#' function. This way each pair of points has exactly one coefficient
#' tailored for them.
#'
#' @param weight_vec A vector of weights.
#' @param mark_stat Mark statistics.
#' @param mtf_name A vector containing the names of the mark test functions.
#' @param one_per_lambda2 An estimate for the inverse of the estimate of
#'                        the intensity squared
#' @return A matrix of weights. Each column is for a different mark test
#'   function. Each row is for a different _ordered_ point pair where
#'   (x1, x2) is different from (x2, x1).
specialise_weight <- function(weight_vec, mark_stat, mtf_name,
                              one_per_lambda2) {
    mtf_coeff_vec <- sapply(mark_stat[paste(mtf_name, '_coeff', sep = '')],
                            identity)
    outer(weight_vec, mtf_coeff_vec * one_per_lambda2, '*')
}

#' Create mark test functions and weights
#'
#' Create a list of mark test functions and a matrix of weights.
#'
#' @param marks A vector of marks on which the mark test functions and their
#'   scaling coefficients will be based.
#' @param mtf_name A vector containing the names of the mark test functions.
#' @param edge_corr A vector of edge correction coefficients for point
#'   pairs.
#' @param one_per_lambda2 The value of 1 / lambda^2 where lambda^2 is the
#'   estimate of the intensity squared. A scalar.
#' @return A list containing a list of mark tests functions (mtf_func_l) and
#'   a weight matrix (weight_m). Each mark test function takes two numeric
#'   vectors containing the marks of the points x_1 and the points x_2,
#'   respectively, from point pairs (x_1, x_2). The functions are in the
#'   same order as mtf_name. Each column of the weight matrix is reserved
#'   for use with the corresponding mark test function. The columns are in
#'   the same order as mtf_name. Each row corresponds to one pair of points.
create_mtfs_and_weights <- function(marks, mtf_name, edge_corr,
                                    one_per_lambda2) {
    check_mtf(mtf_name)

    mark_distr_stat_l <- mark_distr_stats(marks, mtf_name)
    mtf_func_l <- create_mark_test_funcs(mtf_name, mark_distr_stat_l)

    weight_m <- specialise_weight(edge_corr, mark_distr_stat_l, mtf_name,
                                  one_per_lambda2)

    list(mtf_func_l = mtf_func_l, weight_m = weight_m)
}



#' Is L-transform valid?
#'
#' Ask for which mark test functions Besag's L-transformation can be done
#' given a mark distribution.
#'
#' Given names of mark test functions and information on the (one-point)
#' mark distribution, returns the names of those mark test functions which
#' have a valid Besag's L-transform for the specified mark distribution.
#'
#' @inheritParams summ_func_random_labelling
#' @param is_any_mark_neg A logical scalar value describing whether any of
#'   the marks in the mark distribution have a negative value.
#' @param is_any_mark_pos A logical scalar value describing whether any of
#'   the marks in the mark distribution have a positive value.
#' @return A vector of the names of those mark test function names which
#'   have a valid Besag's L-transform. The names will be in the same order
#'   as they appear in mtf_name.
#' @export
besags_L_valid <- function(mtf_name, is_any_mark_neg, is_any_mark_pos) {
    check_mtf(mtf_name)

    if (length(is_any_mark_neg) != 1L || !is.logical(is_any_mark_neg)) {
        stop('is_any_mark_neg must be a logical scalar vector.')
    }
    if (length(is_any_mark_pos) != 1L || !is.logical(is_any_mark_pos)) {
        stop('is_any_mark_pos must be a logical scalar vector.')
    }

    if (is_any_mark_neg) {
        if (is_any_mark_pos) {
            # Positive and negative marks.
            transformable <- c('1', 'gamma', 'morAbs', 'gammaAbs')
        } else {
            # Only nonpositive marks.
            transformable <- c('1', 'mm', 'gamma', 'morAbs', 'gammaAbs')
        }
    } else {
        if (is_any_mark_pos) {
            # Only nonnegative marks.
            transformable <- c('1', 'm', 'mm', 'gamma', 'morAbs',
                               'gammaAbs')
        } else {
            # Only zero marks. A bit silly.
            transformable <- c('1', 'm', 'mm', 'gamma', 'morAbs',
                               'gammaAbs')
        }
    }

    intersect(mtf_name, transformable)
}
