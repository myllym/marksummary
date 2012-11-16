#' Length of the shorter window side.
#'
#' @param A \code{\link[spatstat]{ppp}} object that has a rectangular
#'   window.
#' @return A scalar value giving the length of the shorter edge of the
#'   window.
min_window_side_length <- function(pattern) {
    with(pattern[['window']], min(diff(xrange), diff(yrange)))
}

#' Round non-negative numbers and convert them to integers.
#'
#' @param x A numeric type vector with only non-negative values.
#' @return An integer type vector with the given values rounded.
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
#' FIXME
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
#' If the first r_vec value is zero, the index for that "bin" matching the
#' r-interval [0, 0] should never be given by the returned function.
#'
#' Also varying step size in r_vec is handled correctly.
#'
#' @param r_vec The radius vector which defines the endpoints of the bins.
#' @return A function that turns the distance between a pair of points into
#'   a bin index.
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
                    as.integer(ceiling(1 / step * pair_dist))
                }
            } else {
                binning_func <- function(pair_dist) {
                    as.integer(ceiling(1 / step * pair_dist)) - 1L
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

#' Pick close pairs.
#'
#' Pick indices of point pairs that are within r_max of each other.
#'
#' If close enough to each other, (x_i, x_j) and (x_j, x_i) are both
#' included as long as i != j.
#'
#' @param dist_m A distance matrix, e.g. from
#'   \code{\link[spatstat]{pairdist}}, containing distances between all
#'   points.
#' @param r_max Maximum radius as a scalar.
#' @return A matrix of indices, row and column values separately.
pairs_within_r_max <- function(dist_m, r_max) {
    n_point <- nrow(dist_m)
    leq_r_max <- dist_m <= r_max
    false_diag <- matrix(!as.logical(diag(n_point)), ncol = n_point)
    nearby_arr_idx <- which(leq_r_max & false_diag, arr.ind = TRUE)
}

#' Handle things considering the radius vector.
#'
#' @param pattern A \code{\link[spatstat]{ppp}} object.
#' @param r_max The maximum radius to consider. r_vec overrides r_max.
#' @param r_vec The radius vector to consider. r_vec overrides r_max.
#' @return A list containing the radius vector (r_vec), the number of bins
#'   (n_bin), the indexing matrix to pick the point pairs within r_max
#'   (nearby_arr_idx) and the indices of the bins in which those point pairs
#'   belong to distance-wise (bin_idx).
#'
#'   Assume n is the amount of points in the point pattern and p <= n, then
#'   nearby_arr_idx is an integer matrix of size p x 2 and bin_idx is an
#'   integer vector of size p. nearby_arr_idx can be used to index an n x n
#'   matrix describing some value for each point pair OR a vector of length
#'   n describing some value of each point. Using nearby_arr_idx will pick p
#'   values in the former case, and 2 * p in the latter case.
#'
#'   For example, distances between all points can be stored simply in an
#'   n x n matrix. The mark of each point can be stored simply in a vector
#'   of length n. Use nearby_arr_idx and bin_idx accordingly.
#' @importFrom spatstat pairdist.ppp
consider_radius <- function(pattern, r_max, r_vec) {
    min_window_side <- min_window_side_length(pattern)

    if (length(r_vec) > 0L) {
            check_r_vec(r_vec, min_window_side)
            r_max <- max(r_vec)
    } else {
        if (length(r_max) < 1L) {
            r_max <- create_r_max(min_window_side)
        } else {
            check_r_max(r_max, min_window_side)
        }
        r_vec <- create_r_vec(r_max, min_window_side)
    }
    n_bin <- length(r_vec)
    binning_func <- create_binning_func(r_vec)

    dist_m <- pairdist.ppp(pattern)
    nearby_arr_idx <- pairs_within_r_max(dist_m, r_max)

    bin_idx <- binning_func(dist_m[nearby_arr_idx])

    list(r_vec = r_vec, n_bin = n_bin, nearby_arr_idx = nearby_arr_idx,
         bin_idx = bin_idx)
}

#' Choose the marks of those pairs where the distance between the point pair
#' is smaller than the maximum radius.
#'
#' Well, actually we select those marks that the index matrix points to.
#'
#' @param marks A vector containing all the marks of the point pattern in
#'   the order of the points.
#' @param nearby_arr_idx An index matrix that picks the points within
#'   maximum radius. A picked mark will be picked twice to accommodate for
#'   the point pair ordering.
#' @return A list of two mark vectors, mark1 and mark2.
marks_within_radius <- function(marks, nearby_arr_idx) {
    # Picking from a vector of length n with a matrix of size p x 2.
    mark1_mark2_m <- matrix(marks[nearby_arr_idx], ncol = 2L)
    mark1 <- mark1_mark2_m[, 1, drop = TRUE]
    mark2 <- mark1_mark2_m[, 2, drop = TRUE]
    list(mark1 = mark1, mark2 = mark2)
}
