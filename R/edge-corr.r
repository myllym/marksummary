#' Translational edge-correction
#'
#' Calculates a vector of edge correction weights for given point pairs
#' using the translational edge-correction. Notice that the inverse has
#' already been taken in the returned vector. Each matrix element
#' e_{1, 2}^{-1} has value
#' 1.0 / ((win_x - abs(x1 - x2)) * (win_y - abs(y1 - y2))).
#'
#' @references
#' [1] J. Illian, A. Penttinen, H. Stoyan, and D. Stoyan, Statistical
#'     Analysis and Modelling of Spatial Point Patterns, 1st ed. p. 188.
#'     John Wiley & Sons, Ltd, 2008.
#'
#' @param window A window object from \code{\link[spatstat]{owin}} where the
#'   points reside. The window has to have type "rectangle".
#' @param x1 A vector containing the x-coordinates of the first points of
#'   the pairs. The order of the elements has to correspond with y1.
#' @param y1 A vector containing the y-coordinates of the first points of
#'   the pairs. The order of the elements has to correspond with x1.
#' @param x2 A vector containing the x-coordinates of the second points of
#'   the pairs. The order of the elements has to correspond with y2.
#' @param y2 A vector containing the y-coordinates of the second points of
#'   the pairs. The order of the elements has to correspond with x2.
#' @return A vector containing the edge correction coefficients. Notice that
#'   the inverse has already been taken. The order of the pairs is the same
#'   as in ((x1, y1), (x2, y2)).
translational_correction <- function(window, x1, y1, x2, y2) {
    with(window, 1.0 / ((diff(xrange) - abs(x1 - x2)) *
                        (diff(yrange) - abs(y1 - y2))))
}

#' Calculate the edge correction for those point pairs that matter.
do_edge_correction <- function(pattern, corr_name, nearby_arr_idx) {
    if (corr_name == 'translate') {
        window <- pattern[['window']]
        x <- pattern[['x']]
        y <- pattern[['y']]

        first_idx <- nearby_arr_idx[, 1, drop = TRUE]
        second_idx <- nearby_arr_idx[, 2, drop = TRUE]
        x1 <- x[first_idx]
        y1 <- y[first_idx]
        x2 <- x[second_idx]
        y2 <- y[second_idx]

        edge_corr <- translational_correction(window, x1, y1, x2, y2)
    } else if (corr_name == 'none') {
        edge_corr <- rep.int(1, nrow(nearby_arr_idx))
    } else {
        stop('This edge correction has not been implemented: ', corr_name)
    }

    edge_corr
}
