#' Return the length of the shorter window side.
min_window_side_length <- function(pattern) {
    with(pattern[['window']], min(diff(xrange), diff(yrange)))
}

#' Round non-negative numbers and convert them to integers.
round_nonneg_into_int <- function(x) {
    if (any(x < 0)) {
        stop('x should be non-negative.')
    }
    as.integer(x + 0.5)
}

#' Check that the given maximum radius makes sense.
check_r_max <- function(r_max, min_window_side) {
    if (length(r_max) != 1L || !is.finite(r_max) || r_max <= 0 ||
        r_max > min_window_side) {
        stop('If given, r_max should be a positive scalar between zero and',
             ' the minimum window side length.')
    }
}

#' Check that the given radius vector makes sense.
check_r_vec <- function(r_vec, min_window_side) {
    n_r <- length(r_vec)
    if (n_r < 1L) {
        stop('r_vec must not be NULL.')
    }
    if (any(!is.finite(r_vec))) {
        stop('r_vec must have only finite values.')
    }
    if (max(r_vec) > min_window_side) {
        stop('r_vec must not have a greater maximum than is the length of ',
             'the shortest window edge.')
    }
    if (is.unsorted(r_vec, strictly = TRUE)) {
        stop('r_vec must increase monotonically.')
    }
    n_nonpositive <- sum(r_vec <= 0)
    if (n_r == 1L) {
        if (n_nonpositive > 0L) {
            stop('A scalar r_vec must not have value zero.')
        }
    } else {
        if (n_nonpositive > 1L) {
            stop('If r_vec is non-scalar, only one zero in the beginning ',
                 'is allowed.')
        }
    }
}

#' Create the maximum radius given the pattern.
create_r_max <- function(min_window_side, recommended_r_max_ratio = 0.25) {
    recommended_r_max_ratio * min_window_side
}

#' Create the radius vector given the maximum radius and pattern.
create_r_vec <- function(r_max, min_window_side, n_max_break = 2048L) {
    n_break_num <- r_max / min_window_side * n_max_break
    # Use round or floor. It should not matter much. If we chose the
    # recommended r_max and if there was no numerical error, n_break_num
    # should be a whole number as closely as numeric can represent.
    n_break <- round_nonneg_into_int(n_break_num)
    # + 1L for the "bin" for the r-interval [0, 0] which should be empty for
    # simple point patterns.
    r_vec <- seq(0, r_max, length.out = n_break + 1L)
}

#' Returns a function which gives C++ array indices.
#'
#' Big fat warning: binning_func returns indices for C++ array indexing.
#' Indexing starts from zero in C++. If you need to use the result of
#' binning_func() in R, add + 1L to the result.
#'
#' The bins are exclusive from left and inclusive from right.
#'
#' If the first r_vec value is zero, the index for that "bin" for the
#' r-interval [0, 0] is never returned.
#'
#' Also varying step size in r_vec is handled correctly.
create_binning_func <- function(r_vec) {
    n_r <- length(r_vec)
    if (n_r == 1L) {
        binning_func <- function(pair_dist) {
            rep.int(0L, length(pair_dist))
        }
    } else {
        is_first_zero <- r_vec[1] == 0

        steps <- diff(r_vec)
        if (n_r == 2L) {
            is_step_const <- TRUE
        } else {
            is_step_const <- isTRUE(all.equal(rep.int(0, n_r - 2L),
                                              diff(steps)))
        }

        if (is_step_const) {
            step <- steps[1]
            if (is_first_zero) {
                binning_func <- function(pair_dist) {
                    as.integer(ceiling(pair_dist / step))
                }
            } else {
                binning_func <- function(pair_dist) {
                    as.integer(ceiling(pair_dist / step)) - 1L
                }
            }
        } else {
            if (is_first_zero) {
                # Assume simple point pattern: Avoid zero index by not
                # having - 1L.
                binning_func <- function(pair_dist) {
                    cut(pair_dist, r_vec, labels = FALSE)
                }
            } else {
                # Add zero to catch point pairs closer than r_vec[1].
                binning_func <- function(pair_dist) {
                    cut(pair_dist, c(0, r_vec), labels = FALSE) - 1L
                }
            }
        }
    }
    binning_func
}

#' Pick indices of point pairs that are within r_max of each other.
#'
#' (x_i, x_i) is not considered a pair.
#'
#' @return A matrix of indices, row and column values separately.
pairs_within_r_max <- function(dist_m, r_max) {
    n_point <- nrow(dist_m)
    leq_r_max <- dist_m <= r_max
    false_diag <- matrix(!as.logical(diag(n_point)), ncol = n_point)
    nearby_arr_idx <- which(leq_r_max & false_diag, arr.ind = TRUE)
}

#' Handle things considering the radius vector.
#'
#' @return A list containing the radius vector (r_vec), the number of bins
#' (n_bin), the indexing matrix to pick the point pairs within r_max
#' (nearby_arr_idx) and the indices of the bins in which those point pairs
#' belong to distance-wise (bin_idx).
#' @importFrom spatstat pairdist.ppp
consider_radius <- function(pattern, r_max = NULL, r_vec = NULL,
                            do_order_pairs = TRUE) {
    min_window_side <- min_window_side_length(pattern)

    if (length(r_vec) > 0L) {
            check_r_vec(r_vec, min_window_side)
            r_max <- max(r_vec)
    } else {
        if (length(r_max) < 1L) {
            r_max <- create_r_max(min_window_side)
        } else {
            check_r_max(r_max)
        }
        r_vec <- create_r_vec(r_max, min_window_side)
    }
    n_bin <- length(r_vec)
    binning_func <- create_binning_func(r_vec)

    dist_m <- spatstat::pairdist.ppp(pattern)
    nearby_arr_idx <- pairs_within_r_max(dist_m, r_max)

    if (do_order_pairs) {
        nearby_dist <- dist_m[nearby_arr_idx]
        order_idx <- order(nearby_dist)
        nearby_arr_idx <- nearby_arr_idx[order_idx, , drop = FALSE]
    }

    bin_idx <- binning_func(dist_m[nearby_arr_idx])

    list(r_vec = r_vec, n_bin = n_bin, nearby_arr_idx = nearby_arr_idx,
         bin_idx = bin_idx)
}

#' Choose those marks that are within the radius.
#'
#' @param marks A vector containing all the marks of the point pattern in
#'   the order of the points.
#' @param nearby_arr_idx A logical vector that picks the points within maximum
#'   radius. A picked mark will be picked twice to accommodate for the
#'   point pair ordering.
#' @return A vector of only the needed marks, ordered.
marks_within_radius <- function(marks, nearby_arr_idx) {
    mark1_mark2_m <- matrix(marks[nearby_arr_idx], ncol = 2L)
    mark1 <- mark1_mark2_m[, 1, drop = TRUE]
    mark2 <- mark1_mark2_m[, 2, drop = TRUE]
    list(mark1 = mark1, mark2 = mark2)
}
