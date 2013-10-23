\name{summ_func}
\alias{summ_func}
\title{Summary functions for one pattern.}
\usage{
  summ_func(..., n_sim = 0L)
}
\arguments{
  \item{n_sim}{Ignored.}

  \item{...}{Currently unused.}
}
\value{
  A matrix with dimensions: summ_func, r.
}
\description{
  Uses the function
  \code{\link{summ_func_random_labelling}}. Eats the number
  of permutations the user has given and replaces it with
  zero.
}
\seealso{
  summ_func_random_labelling
}

