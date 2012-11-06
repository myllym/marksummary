\name{translational_correction}
\alias{translational_correction}
\title{Returns a weight matrix for all point pairs using translational
edge-correction.}
\usage{
  translational_correction(pattern)
}
\arguments{
  \item{window}{A window object from
  \code{\link[spatstat]{owin}} where the points reside. The
  window has to have type "rectangle".}

  \item{x}{A vector containing the x coordinates of the
  points.}

  \item{y}{A vector containing the y coordinates of the
  points.}
}
\value{
  A matrix containing the weights. Notice that the inverse
  has already been taken. The order of the rows and columns
  is the same as in the given coordinate vectors.
}
\description{
  Calculates a matrix of weights for all point pairs using
  the translational edge-correction. Notice that the
  inverse has already been taken in the returned matrix.
  Each matrix element w_{ij}^{-1} has value 1.0 / ((win.x -
  abs(xi - xj)) * (win.y - abs(yi - yj))).
}
\references{
  [1] J. Illian, A. Penttinen, H. Stoyan, and D. Stoyan,
  Statistical Analysis and Modelling of Spatial Point
  Patterns, 1st ed. p. 484.  John Wiley & Sons, Ltd, 2008.
}

