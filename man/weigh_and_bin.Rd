\name{weigh_and_bin}
\alias{weigh_and_bin}
\title{Weigh and sum into bins.}
\usage{
  weigh_and_bin(f_value, weight, idx_vec, n_bin)
}
\arguments{
  \item{f_value}{A vector of values from the mark test
  function corresponding to each pair of points.}

  \item{weight}{A vector of weights for each pair of
  points.}

  \item{idx_vec}{A vector of bin indices to guide each pair
  of points into the correct bin.}

  \item{n_bin}{The amount of bins that the returned vector
  should have.}
}
\value{
  A rho_f vector with length of n_bin.
}
\description{
  Weigh the value of the f-function for each pair, sum
  together results of the same bin and place the sum in the
  correct bin.
}

