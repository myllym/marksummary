#' Summary functions for one pattern.
#'
#' Uses the function \code{\link{summ_func_random_labelling}}. Eats the
#' number of permutations the user has given and replaces it with zero.
#'
#' @seealso summ_func_random_labelling
#' @inheritParams summ_func_random_labelling
#' @param n_perm Ignored.
#' @return A matrix with dimensions: summ_func, r.
#' @importFrom abind adrop
#' @export
summ_func <- function(..., n_perm = 0L) {
    res <- summ_func_random_labelling(..., n_perm = 0L)
    # There are no simulations.
    res[['a']] <- adrop(res[['a']], drop = 1L)
    res
}

#' Summary functions for random labelling
#'
#' Calculates summary functions for a pattern and its simulations under the
#' hypothesis of random labelling.
#'
#' @param pattern A \code{\link[spatstat]{ppp}} object as the simple marked
#'   point pattern to be analysed. The marks need to be in the form of a
#'   numeric vector. The window has to have the type "rectangle".
#' @param edge_correction The name of the edge correction to be used. Options
#'   are 'translate' and 'none'.
#' @param mtf_name A vector of mark test function names. "1" stands for the
#'   unmarked K-function. Accepted values are '1', 'm', 'mm', 'gamma',
#'   'gammaAbs' and 'morAbs'.
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
#' @param use_biased_lambda2 A logical scalar on whether to use the biased
#'   or the unbiased (in the Poisson case) estimate of the intensity
#'   squared.
#' @param method The name of the method to create simulations under the null
#'   hypothesis. 'permute' results in permutations of the marks. Using
#'   'sample' will sample the marks from the empirical mark distribution
#'   with replacement. 'permute' is the default.
#' @param ... Currently unused.
#' @return A list containing the radius vector (r) and an array (a)
#'   containing the summary function estimates for the original pattern and
#'   the randomly labelled patterns with dimensions orig_and_perm,
#'   summ_func, r. The estimates for the given marked point pattern are on
#'   the row named "original".
#' @importFrom spatstat pairdist.ppp
#' @importFrom abind abind
#' @export
#' @examples
#' require(spatstat)
#' pp <- spruces
#' res <- summ_func_random_labelling(pp)
#' plot(res, mtf_name = "m", L=TRUE)
summ_func_random_labelling <-
    function(pattern, edge_correction = 'translate',
             mtf_name = c('1', 'm', 'mm', 'gamma', 'gammaAbs', 'morAbs'),
             n_perm = 999L, r_max = NULL, r_vec = NULL, do_besags_L = TRUE,
             method = 'permute', use_biased_lambda2 = FALSE, ...) {
    check_pattern(pattern)
    check_n_perm(n_perm)

    # Handle things related to distances of points.
    radius_l <- consider_radius(pattern, r_max, r_vec)
    r_vec <- radius_l[['r_vec']]
    n_bin <- radius_l[['n_bin']]
    nearby_arr_idx <- radius_l[['nearby_arr_idx']]
    bin_idx <- radius_l[['bin_idx']]

    # Handle other things related to the unmarked point pattern.
    one_per_lambda2 <-
        one_per_lambda_squared(pattern,
                               use_biased_lambda2 = use_biased_lambda2)
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
    # Dims: r, mtf
    orig_summ_func_m <- K_f(orig_mark1, orig_mark2, bin_idx, orig_weight_m,
                            orig_mtf_func_l, mtf_name, n_bin)

    # Change dims into: orig_and_perm, summ_func, r
    orig_summ_func_a <- t(orig_summ_func_m)
    dim(orig_summ_func_a) <- c(1L, dim(orig_summ_func_a))
    dimnames(orig_summ_func_a) <- c(list(orig_and_perm = c('original')),
                                    rev(dimnames(orig_summ_func_m)))

    # FIXME: K_1 gets calculated again and again. It is enough to do it for
    # the original marks if it is asked for by the user at all. Then remove
    # '1' from mtf_name, orig_weight_m, orig_mtf_func_l and add K_1 and
    # possibly L_1 back later for the simulations.

    if (n_perm > 0L) {
        # Simulate from the null hypothesis and estimate summary functions.
        if (method %in% 'permute') {
            # Result dimensions: r, summ_func, permutation.
            sim_summ_func_a <-
                replicate(n_perm, {
                          perm_marks <- sample(orig_marks, replace = FALSE)

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
        } else if (method %in% 'sample') {
            # Result dimensions: r, summ_func, permutation.
            sim_summ_func_a <-
                replicate(n_perm, {
                          perm_marks <- sample(orig_marks, replace = TRUE)

                          perm_mtf_l <-
                              create_mtfs_and_weights(perm_marks,
                                                      mtf_name,
                                                      edge_corr,
                                                      one_per_lambda2)
                          perm_mtf_func_l <- perm_mtf_l[['mtf_func_l']]
                          perm_weight_m <- perm_mtf_l[['weight_m']]

                          perm_mark_l <- marks_within_radius(perm_marks,
                                                             nearby_arr_idx)
                          perm_mark1_vec <- perm_mark_l[['mark1']]
                          perm_mark2_vec <- perm_mark_l[['mark2']]

                          perm_summ_func_m <- K_f(perm_mark1_vec,
                                                  perm_mark2_vec, bin_idx,
                                                  perm_weight_m,
                                                  perm_mtf_func_l, mtf_name,
                                                  n_bin)
                })
        } else {
            stop('method was not recognized.')
        }

        sim_summ_func_a <-
            array(sim_summ_func_a,
                  dim=c(n_bin, dim(orig_summ_func_a)[2], n_perm),
                  dimnames=list(r=NULL,
                                summ_func=dimnames(orig_summ_func_a)[[2]],
                                orig_and_perm=NULL))
        sim_summ_func_a <- aperm(sim_summ_func_a, c(3L, 2L, 1L))

        all_summ_func_a <- abind(orig_summ_func_a, sim_summ_func_a,
                                 along = 1L)
        names(dimnames(all_summ_func_a)) <- c('orig_and_perm', 'summ_func',
                                              'r')
    } else {
        all_summ_func_a <- orig_summ_func_a
    }

    if (do_besags_L) {
        is_any_mark_neg <- any(orig_marks < 0)
        is_any_mark_pos <- any(orig_marks > 0)
        L_summ_a <- all_besags_L(all_summ_func_a,
                                 is_any_mark_neg = is_any_mark_neg,
                                 is_any_mark_pos = is_any_mark_pos)

        all_summ_func_a <- abind(all_summ_func_a, L_summ_a, along = 2L)
        names(dimnames(all_summ_func_a)) <- c('orig_and_perm', 'summ_func',
                                              'r')
    }

    res <- list(r = r_vec, a = all_summ_func_a, call = match.call())
    class(res) <- 'all_summ_func_rl'
    res
}

plot.all_summ_func_rl <- function(x, mtf_name = NULL, L = TRUE, nrow = NULL,
                                  ncol = NULL, ...) {
    # FIXME
    # Make a data frame according to the arguments, i.e. choose which
    # functions to plot.
    # Plot using ggplot.
}

#' Checks if a window object is a circle
#'
#' @param digit The accuracy to check whether the points at the polygonial boundary
#'               are on the arch of a circle.
is.disc <- function(x, digits=6, ...) {
    if(x[['type']] != 'polygonal') return(FALSE)
    if(diff(x[['xrange']]) != diff(x[['yrange']])) return(FALSE)
    r <- diff(x[['xrange']])/2
    x0 <- mean(x[['xrange']])
    y0 <- mean(x[['yrange']])
    if ( sum( round((x[['bdry']][[1]][['x']]-x0)^2 + (x[['bdry']][[1]][['y']]-y0)^2, digits=digits) == round(r^2, digits=digits) ) != length(x[['bdry']][[1]][['x']]) ) return(FALSE)
    return(TRUE)
}

#' Finds the radius and centre of a circle defined by a polygonial boundary
#'
#' @param digit The accuracy to check whether the points at the polygonial boundary
#'               are on the arch of a circle.
disc_param <- function(x, ...) {
    if(!is.disc(x)) stop("x is not a circ.")
    r <- diff(x[['xrange']])/2
    x0 <- mean(x[['xrange']])
    y0 <- mean(x[['yrange']])
    list(r = r, centre = c(x0, y0))
}


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
one_per_lambda_squared <- function(pattern, use_biased_lambda2) {
    if (length(use_biased_lambda2) != 1L ||
        !is.logical(use_biased_lambda2) ||
        !is.finite(use_biased_lambda2)) {
        stop('use_biased_lambda2 must be either TRUE or FALSE.')
    }

    area <- area.owin(pattern[['window']])
    n_point <- pattern[['n']]
    if (use_biased_lambda2) {
        one_per_lambda2 <- area * area / (n_point * n_point)
    } else {
        one_per_lambda2 <- area * area / (n_point * (n_point - 1L))
    }
    one_per_lambda2
}



#' Besag's L-transform.
besags_L <- function(k) {
    sqrt(1 / pi * k)
}

#' Channel to Besag's L.
#'
#' Pick only those K_f that have a valid Besag's L-transform and return the
#' transformed summary functions.
#'
#' @param all_k_a An array containing all K_f for all patterns. Dimensions:
#'   orig_and_perm, summ_func, r.
#' @return An array of dimensions: orig_and_perm, summ_func, r. The size of
#'   the second dimension depends on how many K_f functions had a valid
#'   L-transform.
all_besags_L <- function(all_k_a, is_any_mark_neg, is_any_mark_pos,
                         sep_char = '_', L_char = 'L') {
    K_f_name <- dimnames(all_k_a)[['summ_func']]
    # First row of sapply result is all 'K', second row mtf.
    mtf_name <- sapply(strsplit(K_f_name, sep_char, fixed=TRUE),
                       identity)[2, , drop = TRUE]
    valid_mtf_name <- besags_L_valid(mtf_name,
                                     is_any_mark_neg = is_any_mark_neg,
                                     is_any_mark_pos = is_any_mark_pos)
    match.idx <- match(valid_mtf_name, mtf_name)

    if (length(match.idx) < 1L) {
        all_l_a <- NULL
    } else {
        all_l_a <- besags_L(all_k_a[, match.idx, , drop = FALSE])
        L_f_name <- paste(L_char, mtf_name, sep = sep_char)
        dimnames(all_l_a)[2] <- list(summ_func = L_f_name)
    }

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
    names(dimnames(K_f_m)) <- c('r', 'summ_func')
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

#' Calculates the uncumulated rho_f-estimates
#'
#' @param mark1 The marks of the first points of pairs in a vector.
#' @param mark2 The marks of the second points of pairs in a vector.
#' @param bin_idx A vector of indices which tell which radius bin each point
#'   pair belongs to.
#' @param weight_m A weight matrix for point pairs.
#' @param mtf_func_l A list of mark test functions all of which take two
#'   real marks and return a real value.
#' @param mtf_name The names of the mark test functions. Matches mtf_func_l
#'   in order.
#' @param n_bin The number of bins in the radius vector.
#' @return A matrix of rho_f values. Dimensions r, mtf.
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
