#' Checks that the given pattern is valid.
#'
#' @importFrom spatstat is.ppp
#' @importFrom spatstat is.empty.ppp
check_pattern <- function(pattern) {
    if (length(pattern) < 1L) {
        stop('Pattern was not given.')
    }
    if (!is.ppp(pattern)) {
        stop('The pattern must be a ppp object.')
    }
    if (is.empty.ppp(pattern)) {
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

#' Check that the given number of permutations makes sense.
check_n_perm <- function(n_perm) {
    if (length(n_perm) != 1L || !is.finite(n_perm) || n_perm < 0L) {
        stop('Number of permutations has to be a scalar, finite, ',
             'non-negative number.')
    }
}

#' Estimate 1 / lambda ^ 2.
#'
#' Notice that this estimator is biased.
#'
#' @importFrom spatstat area.owin
one_per_lambda_squared <- function(pattern, use_biased = FALSE) {
    area <- area.owin(pattern[['window']])
    n_point <- pattern[['n']]
    if (use_biased) {
        one_per_lambda2 <- area * area / (n_point * n_point)
    } else {
        one_per_lambda2 <- area * area / (n_point * (n_point - 1L))
    }
    one_per_lambda2
}

#' Summary functions for one pattern.
#'
#' Uses the function \code{\link{summ_func_random_labelling}}. Eats the
#' number of permutations the user has given and replaces it with zero.
#'
#' @seealso summ_func_random_labelling
#' @return A matrix with dimensions: summ_func, r.
#' @export
summ_func <- function(..., n_perm = 0L) {
    res <- summ_func_random_labelling(..., n_perm = 0L)
    res[['a']] <- res[['a']]['original', , , drop = TRUE]
    res
}

#' Calculates summary functions for a pattern and its simulations under the
#' hypothesis of random labelling.
#'
#' @param pattern A \code{\link[spatstat]{ppp}} object as the simple marked
#'   point pattern to be analysed. The marks need to be in the form of a
#'   numeric vector. The window has to have the type "rectangle".
#' @param edge_corr_func The name of the edge correction to be used. Options
#'   are 'translate' and 'none'.
#' @param mtf_name A vector of mark test function names. "1" stands for the
#'   unmarked K-function.
#' @param n_perm The number of permutations.
#' @param r_max A positive scalar value representing the maximum radius that
#'   should be considered. r_vec overrides r_max. By default, r_max is NULL
#'   and will get a sensible default.
#' @param r_vec A monotonically increasing vector of non-negative r-values
#'   to act as the endpoints of the bins for the K_f-functions. r_vec
#'   overrides r_max. The bins are exclusive on the left and inclusive on
#'   the right. If the first vector element has value zero, it will be
#'   regarded as the collapsed bin [0, 0], and the next bin will start from
#'   and exclude 0.
#' @param do_besags_L A boolean describing whether Besag's L-function should
#'   also be returned where available.
#' @param method The name of the method to create simulations under the null
#'   hypothesis.
#' @param ... Currently unused.
#' @return An array containing the summary function estimates for the
#'   original pattern and the randomly labelled patterns. Dimensions:
#'   orig_and_perm, summ_func, r. The estimates for the given marked point
#'   pattern are on the row named "original".
#' @importFrom spatstat pairdist.ppp
#' @importFrom abind abind
#' @export
summ_func_random_labelling <-
    function(pattern, edge_correction = 'translate',
             mtf_name = c('1', 'm', 'mm', 'gamma', 'gammaAbs', 'morAbs'),
             n_perm = 999L, r_max = NULL, r_vec = NULL, do_besags_L = TRUE,
             method = 'permute', ...) {
    check_pattern(pattern)
    check_n_perm(n_perm)

    # Handle things related to distances of points.
    radius_l <- consider_radius(pattern, r_max, r_vec)
    r_vec <- radius_l[['r_vec']]
    n_bin <- radius_l[['n_bin']]
    nearby_arr_idx <- radius_l[['nearby_arr_idx']]
    bin_idx <- radius_l[['bin_idx']]

    # Handle other things related to the unmarked point pattern.
    one_per_lambda2 <- one_per_lambda_squared(pattern)
    edge_corr <- do_edge_correction(pattern, edge_correction,
                                    nearby_arr_idx)

    orig_marks <- pattern[['marks']]

    # Create mark test functions and the corresponding weights for pairs.
    mtf_l <- create_mtfs_and_weights(orig_marks, mtf_name, edge_corr,
                                     one_per_lambda2)
    orig_mtf_func_l <- mtf_l[['mtf_func_l']]
    orig_weight_m <- mtf_l[['weight_m']]

    # Pick relevant marks.
    orig_mark_l <- marks_within_radius(orig_marks, nearby_arr_idx)
    orig_mark1 <- orig_mark_l[['mark1']]
    orig_mark2 <- orig_mark_l[['mark2']]

    # Calculate summary function estimates.
    orig_summ_func_m <- K_f(orig_mark1, orig_mark2, bin_idx, orig_weight_m,
                            orig_mtf_func_l, mtf_name, n_bin)

    # FIXME: K_1 gets calculated again and again. It is enough to do it for
    # the original marks. Then remove '1' from mtf_name and add it back
    # later for the simulations.

    # Simulate from the null hypothesis and estimate summary functions.
    if (method == 'permute') {
        # Result dimensions: permutation, summ_func, r.
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



    # FIXME: This part needs cleaning up.



    # Combine summary functions from the original pattern and simulations.
    if (length(sim_summ_func_a) < 1L) {
        all_summ_func_a <-
            array(orig_summ_func_m, c(dim(orig_summ_func_m), 1L),
                  dimnames = list(r = NULL,
                                  summ_func =
                                  dimnames(orig_summ_func_m)[[2]],
                                  orig_and_perm = 'original'))
    } else {
        all_summ_func_a <- abind(orig_summ_func_m, sim_summ_func_a,
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
        is_any_mark_neg <- any(orig_marks < 0)
        is_any_mark_pos <- any(orig_marks > 0)
        l_summ_a <- all_besags_l(all_summ_func_a, is_any_mark_neg,
                                 is_any_mark_pos)
        all_summ_func_a <- abind(all_summ_func_a, l_summ_a, along = 2L)
    }

    # Keep the names.
    # FIXME: Remove the silly repetition. The K_f names have to exist before
    #        Besag's L. Before returning all names have to be in order.
    dimnames(all_summ_func_a) <-
        list(orig_and_perm = c('original', rep.int('', n_perm)),
             summ_func = dimnames(all_summ_func_a)[[2]],
             r = NULL)

    list(r = r_vec, a = all_summ_func_a)
}



#' Besag's L-transform.
besags_l <- function(k) {
    sqrt(1 / pi * k)
}

#' Make sure that only the non-negative functions will get Besag's
#' L-transform.
#'
#' @param all_k_a An array containing all K_f for all patterns. Dimensions:
#'   orig_and_perm, summ_func, r.
all_besags_l <- function(all_k_a, is_any_mark_neg, is_any_mark_pos) {
    available <- dimnames(all_k_a)[['summ_func']]
    if (is_any_mark_neg) {
        if (is_any_mark_pos) {
            # Positive and negative marks.
            transformable <- c('K_1', 'K_gamma', 'K_morAbs', 'K_gammaAbs')
        } else {
            # Only nonpositive marks.
            transformable <- c('K_1', 'K_mm', 'K_gamma', 'K_morAbs',
                               'K_gammaAbs')
        }
    } else {
        if (is_any_mark_pos) {
            # Only nonnegative marks.
            transformable <- c('K_1', 'K_m', 'K_mm', 'K_gamma', 'K_morAbs',
                               'K_gammaAbs')
        } else {
            # Only zero marks. A bit silly.
            transformable <- c('K_1', 'K_m', 'K_mm', 'K_gamma', 'K_mor',
                               'K_morAbs', 'K_gammaAbs')
        }
    }
    matched <- intersect(available, transformable)
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



#' Calculates the cumulated K_f-estimates.
K_f <- function(...) {
    rho_f_m <- rho_f(...)
    # as.matrix is needed if rho_f_m has only one r-value.
    n_r <- nrow(rho_f_m)
    K_f_m <- as.matrix(apply(rho_f_m, 2, cumsum), nrow = n_r)
    if (n_r == 1L) {
        # apply flips the matrix over if there is only one r-value.
        K_f_m <- t(K_f_m)
    }
    colnames(K_f_m) <- paste('K_', colnames(rho_f_m), sep = '')
    K_f_m
}

#' Weigh and sum into bins.
#'
#' Weigh the value of the f-function for each pair, sum together results
#' of the same bin and place the sum in the correct bin.
#'
#' @param f_value A vector of values from the mark test function
#'   corresponding to each pair of points.
#' @param weight A vector of weights for each pair of points.
#' @param idx_vec A vector of bin indices to guide each pair of points into
#'   the correct bin.
#' @param n_bin The amount of bins that the returned vector should have.
#' @return A rho_f vector with length of n_bin.
#' @useDynLib marksummary
weigh_and_bin <- function(f_value, weight, idx_vec, n_bin) {
    .Call('weighAndBin', f_value, weight, idx_vec, n_bin,
          PACKAGE = 'marksummary')
}

#' Calculates the uncumulated rho_f-estimates.
rho_f <- function(mark1, mark2, bin_idx, weight_m, mtf_func_l, mtf_name,
                  n_bin) {
    val_m <- as.matrix(vapply(seq_along(mtf_func_l),
                  function(mtf_idx) {
                      mtf_func <- mtf_func_l[[mtf_idx]]
                      f_value <- mtf_func(mark1, mark2)
                      one_rho_f <-
                          weigh_and_bin(f_value,
                                        weight_m[, mtf_idx, drop = TRUE],
                                        bin_idx, n_bin)
                      one_rho_f
                  }, numeric(n_bin)))
    # vapply acts differently depending on whether the result is a vector or
    # a scalar.
    if (n_bin == 1L) {
        val_m <- t(val_m)
    }
    colnames(val_m) <- mtf_name
    val_m
}
