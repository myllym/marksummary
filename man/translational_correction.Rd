\name{translational_correction}
\alias{translational_correction}
\title{Translational edge-correction.}
\usage{
  translational_correction(window, x1, y1, x2, y2)
}
\arguments{
  \item{window}{A window object from
  \code{\link[spatstat]{owin}} where the points reside. The
  window has to have type "rectangle".}

  \item{x1}{A vector containing the x-coordinates of the
  first points of the pairs. The order of the elements has
  to correspond with y1.}

  \item{y1}{A vector containing the y-coordinates of the
  first points of the pairs. The order of the elements has
  to correspond with x1.}

  \item{x2}{A vector containing the x-coordinates of the
  second points of the pairs. The order of the elements has
  to correspond with y2.}

  \item{y2}{A vector containing the y-coordinates of the
  second points of the pairs. The order of the elements has
  to correspond with x2.}
}
\value{
  A vector containing the edge correction coefficients.
  Notice that the inverse has already been taken. The order
  of the pairs is the same as in ((x1, y1), (x2, y2)).
}
\description{
  Calculates a vector of edge correction weights for given
  point pairs using the translational edge-correction.
  Notice that the inverse has already been taken in the
  returned vector. Each matrix element e_{1, 2}^{-1} has
  value 1.0 / ((win_x - abs(x1 - x2)) * (win_y - abs(y1 -
  y2))).
}
\references{
  [1] J. Illian, A. Penttinen, H. Stoyan, and D. Stoyan,
  Statistical Analysis and Modelling of Spatial Point
  Patterns, 1st ed. p. 188.  John Wiley & Sons, Ltd, 2008.
}

