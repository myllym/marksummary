\name{create_binning_func}
\alias{create_binning_func}
\title{Returns a function which gives C++ array indices.}
\usage{
  create_binning_func(r_vec)
}
\arguments{
  \item{r_vec}{The radius vector which defines the
  endpoints of the bins.}
}
\value{
  A function that turns the distance between a pair of
  points into a bin index.
}
\description{
  Big fat warning: binning_func returns indices for C++
  array indexing. Indexing starts from zero in C++. If you
  need to use the result of binning_func() in R, add + 1L
  to the result.
}
\details{
  The bins are exclusive from left and inclusive from
  right.

  If the first r_vec value is zero, the index for that
  "bin" for the r-interval [0, 0] is never returned.

  Also varying step size in r_vec is handled correctly.
}

