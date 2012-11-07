#' Returns a weight matrix for all point pairs using translational
#' edge-correction.
#'
#' Calculates a matrix of weights for all point pairs using the
#' translational edge-correction. Notice that the inverse has already been
#' taken in the returned matrix. Each matrix element w_{ij}^{-1} has value
#' 1.0 / ((win.x - abs(xi - xj)) * (win.y - abs(yi - yj))).
#'
#' @references
#' [1] J. Illian, A. Penttinen, H. Stoyan, and D. Stoyan, Statistical
#'     Analysis and Modelling of Spatial Point Patterns, 1st ed. p. 484.
#'     John Wiley & Sons, Ltd, 2008.
#'
#' @param window A window object from \code{\link[spatstat]{owin}} where the
#'   points reside. The window has to have type "rectangle".
#' @param x A vector containing the x coordinates of the points.
#' @param y A vector containing the y coordinates of the points.
#' @return A matrix containing the weights. Notice that the inverse has
#'   already been taken. The order of the rows and columns is the same as in
#'   the given coordinate vectors.
translational_correction <- function(pattern) {
    weigh_one_dim <- function(win_len, pos) {
        win_len - abs(outer(pos, pos, '-'))
    }

    window <- pattern[['window']]
    x <- pattern[['x']]
    y <- pattern[['y']]
    with(window, 1.0 / (weigh_one_dim(diff(xrange), x) *
                        weigh_one_dim(diff(yrange), y)))
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

#' Checks that the given pattern is valid.
#'
#' @importFrom spatstat is.ppp
#' @importFrom spatstat is.empty.ppp
check_pattern <- function(pattern) {
    if (!spatstat::is.ppp(pattern)) {
        stop('The pattern must be a ppp object.')
    }
    if (spatstat::is.empty.ppp(pattern)) {
        stop('The pattern must not be empty.')
    }
    marks <- pattern[['marks']]
    if (length(marks) < 1L) {
        stop('The pattern must have marks.')
    }
    if (pattern[['markformat']] != 'vector') {
        stop('The markformat must be vector.')
    }
    if (!is.numeric(marks)) {
        stop('The marks must be of type numeric or integer.')
    }
    if (pattern[['window']][['type']] != 'rectangle') {
        stop('The window must have type \"rectangle\".')
    }
    if (any(duplicated(matrix(c(pattern[['x']], pattern[['y']]), ncol = 2L),
                       MARGIN = 1L))) {
        stop('The pattern must be simple i.e. without duplicate points ',
             'location-wise.')
    }
}

#' Estimate 1 / \lambda ^ 2.
#'
#' Notice that this estimator is biased.
#' TODO: Consider unbiased estimator.
#'
#' @importFrom spatstat area.owin
one_per_lambda_to_the_second <- function(pattern) {
    area <- spatstat::area.owin(pattern[['window']])
    n <- pattern[['n']]
    one_per_lambda2 <- area * area / (n * n)
    #one_per_lambda2 <- area * area / (n * (n - 1L))
}




#' Create an individual weight vector for each mark test function.
#'
#' Create an individual weight vector for each mark test function. Bake in
#' 1 / \lambda ^ 2 and 1 / c_f so that all weighing can be handled by one
#' multiplication per pair of points per mark test function.
#' @return A matrix of weights. Each column is for different mark test
#'   functions. The columns are named. Each row is for a different _ordered_
#'   point pair where (x1, x2) is different from (x2, x1).
specialise_weight <- function(weight_vec, mark_stat, mtf_name,
                              one_per_lambda2) {
    mtf_coeff_vec <- sapply(mark_stat[paste(mtf_name, '_coeff', sep = '')],
                            identity)
    outer(weight_vec, mtf_coeff_vec * one_per_lambda2, '*')
}



create_mtf_func_1 <- function(mark_stat) {
    return(function(m1, m2) { rep.int(1, length(m1)) })
}

create_mtf_func_m <- function(mark_stat) {
    return(function(m1, m2) { m1 })
}

create_mtf_func_mm <- function(mark_stat) {
    return(function(m1, m2) { m1 * m2 })
}

create_mtf_func_gamma <- function(mark_stat) {
    return(function(m1, m2) {
               mark_diff <- m1 - m2
               0.5 * (mark_diff * mark_diff)
           })
}

create_mtf_func_gammaAbs <- function(mark_stat) {
    return(function(m1, m2) { abs(m1 - m2) })
}

create_mtf_func_mor <- function(mark_stat) {
    mean_mark <- mark_stat[['mean']]
    return(function(m1, m2) { (m1 - mean_mark) * (m2 - mean_mark) })
}

create_mtf_func_morAbs <- function(mark_stat) {
    mean_mark <- mark_stat[['mean']]
    return(function(m1, m2) { abs((m1 - mean_mark) * (m2 - mean_mark)) })
}



#' Given the names of mark test functions and the required mark statistics,
#' produce mark test functions as closures.
prepare_mark_test_funcs <- function(mtf_name, mark_stat) {
    mtf_func_constructor <- paste('create_mtf_func_', mtf_name, sep = '')
    mtf_func_l <- lapply(mtf_func_constructor, function(constructor) {
                         get(constructor)(mark_stat)
                  })
}



#' Summary functions for one pattern.
#'
#' Uses the function \code{\link{summ_func_random_labelling}}. Eats the
#' number of permutations the user has given and replaces it with zero.
#' @export
summ_func <- function(..., n_perm = 0L) {
    summ_func_random_labelling(..., n_perm = 0L)
}

#' Calculates summary functions for a pattern and its simulations under the
#' hypothesis of random labelling.
#'
#' @return An array containing the summary functions for the original
#'   pattern and the randomly labelled patterns. Dimensions: orig_and_perm,
#'   summ_func, r.
#' @importFrom spatstat pairdist.ppp
#' @importFrom abind abind
#' @export
summ_func_random_labelling <-
    function(pattern, edge_corr_func = translational_correction,
             mtf_name = c('1', 'm', 'mm', 'gamma', 'gammaAbs', 'mor',
                          'morAbs'),
             n_perm = 999L, r_max = NULL, r_vec = NULL, do_besags_L = TRUE,
             method = 'permute', do_order_pairs = TRUE, ...) {

    # Check input.
    #
    # Handle r_max and r_vec.
    # Unmarked stuff.
    # - Interesting point pairs (needs r_max).
    # - Which r-vector bin each point pair belongs to?
    # Mark distribution stats.
    # Unmarked stuff.
    # - Weights (need mark distribution stats).
    # Prepare MTFs (need mark distribution stats).
    # For each marked pattern:
    #     Choose and order marks (need unmarked stuff).
    # Throw into K_f.
    #
    # Rearrange, name results, etc. The boring stuff.

   check_pattern(pattern)

    if (length(mtf_name) < 1L) {
        stop('At least one test mark function name needs to be given.')
    }

    if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 0L) {
        stop('Number of permutations has to be a scalar, finite, ',
             'non-negative number.')
    }

    radius_l <- consider_radius(pattern, r_max = r_max, r_vec = r_vec,
                                do_order_pairs = do_order_pairs)
    r_vec <- radius_l[['r_vec']]
    n_bin <- radius_l[['n_bin']]
    nearby_arr_idx <- radius_l[['nearby_arr_idx']]
    bin_idx <- radius_l[['bin_idx']]

    # Function that takes pattern, nearby_arr_idx, mtf_name, returns:
    # - mark_stats?
    # - mtf_func_vec
    # - weight?
    orig_marks <- pattern[['marks']]
    orig_mark_distr_stat_l <- mark_distr_stats(orig_marks, mtf_name)

    # TODO: Handle unmarked point pattern in a function and return only what
    #       is needed.
    # what marks_within_radius needs: nearby_arr_idx, bin_order_idx
    # what K_f_native needs: nearby_arr_idx, bin_order_idx, n_bin, weights


    # TODO: only for interesting pairs
    edge_corr_m <- edge_corr_func(pattern)
    edge_corr <- edge_corr_m[nearby_arr_idx]

    orig_mark_l <- marks_within_radius(orig_marks, nearby_arr_idx)
    orig_mark1 <- orig_mark_l[['mark1']]
    orig_mark2 <- orig_mark_l[['mark2']]

    one_per_lambda2 <- one_per_lambda_to_the_second(pattern)
    orig_weight_m <- specialise_weight(edge_corr, orig_mark_distr_stat_l,
                                       mtf_name, one_per_lambda2)

    orig_mtf_func_l <- prepare_mark_test_funcs(mtf_name,
                                               orig_mark_distr_stat_l)

    orig_summ_func_m <- K_f(orig_mark1, orig_mark2, bin_idx, orig_weight_m,
                            orig_mtf_func_l, mtf_name, n_bin)

    if (method == 'permute') {
        sim_summ_func_a <-
            replicate(n_perm, {
                      perm_marks <- sample(orig_marks)

                      perm_mark_l <- marks_within_radius(perm_marks,
                                                         nearby_arr_idx)
                      perm_mark1_vec <- perm_mark_l[['mark1']]
                      perm_mark2_vec <- perm_mark_l[['mark2']]

                      perm_summ_func_m <- K_f(perm_mark1_vec,
                                              perm_mark2_vec, bin_idx,
                                              orig_weight_m,
                                              orig_mtf_func_l, mtf_name,
                                              n_bin)
            })
    } else {
        # TODO: Implement sampling from empirical mark distribution and
        #       calculate new mark statistics, mark test functions, combined
        #       weight matrix etc. for each set of marks.
        stop('Only the method \"permute\" has been implemented thusfar.')
    }



    # FIXME: This part needs work.



    # Combine summary functions from the original pattern and simulations.
    if (length(sim_summ_func_a) < 1L) {
        all_summ_func_a <-
            array(orig_summ_func_m, c(dim(orig_summ_func_m), 1L),
                  dimnames = list(r = NULL,
                                  summ_func =
                                  dimnames(orig_summ_func_m)[[2]],
                                  orig_and_perm = 'original'))
    } else {
        all_summ_func_a <- abind::abind(orig_summ_func_m, sim_summ_func_a,
                                        along = 3L)
    }

    # After this line, the dimension order should be: orig_and_perm,
    # summ_func, r.
    all_summ_func_a <- aperm(all_summ_func_a, c(3L, 2L, 1L))

    # Keep the names.
    dimnames(all_summ_func_a) <-
        list(orig_and_perm = c('original', rep.int('', n_perm)),
             summ_func = dimnames(all_summ_func_a)[[2]],
             r = NULL)

    if (do_besags_L) {
        l_summ_a <- all_besags_l(all_summ_func_a)
        all_summ_func_a <- abind::abind(all_summ_func_a, l_summ_a,
                                        along = 2L)
    }

    # Keep the names.
    # FIXME: Remove this repetition. The K_f names have to exist for besags_L and after all is done, all names have to exist.
    dimnames(all_summ_func_a) <-
        list(orig_and_perm = c('original', rep.int('', n_perm)),
             summ_func = dimnames(all_summ_func_a)[[2]],
             r = NULL)

    list(r = r_vec, a = all_summ_func_a)
}

#' Make sure that only the non-negative functions will get Besag's
#' L-transform.
#'
#' FIXME: Consider negative marks.
#'
#' @param all_k_a An array containing all K_f for all patterns. Dimensions:
#'   orig_and_perm, summ_func, r.
all_besags_l <- function(all_k_a) {
    available_names <- dimnames(all_k_a)[['summ_func']]
    transformable_names <- c(
                             'K',
                             'K_m',
                             'K_mm',
                             'K_gamma',
                             'K_morAbs',
                             'K_gammaAbs'
                            )
    matched <- intersect(available_names, transformable_names)
    all_l_a <- besags_l(all_k_a[, matched, , drop = FALSE])

    # Create L-names.
    # First row 'K', second row mtf.
    split_name_m <- matrix(sapply(strsplit(matched, '_', fixed = TRUE),
                                  identity),
                           ncol = length(matched))
    split_name_m[1, ] <- 'L'
    l_names <- apply(split_name_m, 2, paste, collapse = '_')
    dimnames(all_l_a)[2] <- list(summ_func = l_names)

    all_l_a
}

#' Besag's L-transform.
besags_l <- function(k) {
    sqrt(1 / pi * k)
}



K_f <- function(...) {
    kappa_f_m <- rho_f(...)
    k_f_m <- apply(kappa_f_m, 2, cumsum)
    colnames(k_f_m) <- paste('K_', colnames(kappa_f_m), sep = '')
    k_f_m
}

#' Sum and order the vector containing weighted mark values for each pair.
#'
#' @useDynLib marksummary
sum_vector <- function(sValue, sIdxVec, sNBin) {
    require(Rcpp)
    .Call("sum_vector", sValue, sIdxVec, sNBin, PACKAGE = 'marksummary')
}

rho_f <- function(mark1, mark2, bin_idx, weight_m, mtf_func_l, mtf_name,
                  n_bin) {
    val_m <- sapply(seq_along(mtf_func_l), function(mtf_idx) {
                    mtf_func <- mtf_func_l[[mtf_idx]]

                    # Calculate result of mark test function, weigh and sum
                    # bins together.
                    val <- mtf_func(mark1, mark2) *
                           weight_m[, mtf_idx, drop = TRUE]

                    #vec <- rep.int(0, n_bin)
                    #for (bin_idx in seq_along(sorted_bin_idx)) {
                    #    bin <- sorted_bin_idx[bin_idx]
                    #    vec[bin] <- vec[bin] + val[bin_idx]
                    #}

                    #sum_vec <- sapply(split(val, sorted_bin_idx), sum)
                    #vec <- rep.int(0, n_bin)
                    #vec[sorted_uniq_bin_idx] <- sum_vec

                    vec <- sum_vector(val, bin_idx, n_bin)

                    vec
             })
    colnames(val_m) <- mtf_name
    val_m
}
