#' Translational edge-correction
#'
#' Calculates a vector of edge correction weights for given point pairs
#' using the translational edge-correction. Notice that the inverse has
#' already been taken in the returned vector. Each matrix element
#' e_{1, 2}^{-1} has value
#' 1.0 / ((win_x - abs(x1 - x2)) * (win_y - abs(y1 - y2)))
#' for a rectangular window. For a circular window, it is
#' 1.0 / ( 2.0 * R^2 * acos(d/(2*R)) - d/2 * sqrt(4*R^2 - d^2) )
#' where d is the distance between (x1, y1) and (x2, y2) and R is the
#' radius of the disc.
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
    if(window[['type']] == 'rectangle')
        w <- with(window, 1.0 / ((diff(xrange) - abs(x1 - x2)) *
                                   (diff(yrange) - abs(y1 - y2))))
    if(is.disc(window)) {
        R <- disc_param(window)[['R']]
        d <- sqrt((x1-x2)^2 + (y1-y2)^2)
        w <- 1.0 / ( 2.0 * R^2 * acos(d/(2.0*R)) - d/2.0 * sqrt(4.0*R^2 - d^2) )
    }
    w
}

#' Calculate the edge correction for those point pairs that matter.
#'
#' @param pattern A \code{\link[spatstat]{ppp}} object from which the point
#'   pairs will be picked.
#' @param corr_name The name of the edge correction. Options are 'translate'
#'   and 'none'.
#' @param nearby_arr_idx An array index matrix with two columns. This
#'   decides for which point pairs an edge correction is calculated.
#' @return A vector of edge correction coefficients for each point pair
#'   picked by nearby_arr_idx.
#' @importFrom spatstat area.owin
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
        coeff <- 1 / spatstat::area.owin(pattern[['window']])
        edge_corr <- rep.int(coeff, nrow(nearby_arr_idx))
    } else {
        stop('This edge correction has not been implemented: ', corr_name)
    }

    edge_corr
}
