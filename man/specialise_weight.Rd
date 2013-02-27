\name{specialise_weight}
\alias{specialise_weight}
\title{Individualize weight vectors for each mark test function.}
\usage{
  specialise_weight(weight_vec, mark_stat, mtf_name,
    one_per_lambda2)
}
\value{
  A matrix of weights. Each column is for a different mark
  test function. Each row is for a different _ordered_
  point pair where (x1, x2) is different from (x2, x1).
}
\description{
  Create an individual weight vector for each mark test
  function. The weight of the point pair (x_1, x_2) is
  w_{1, 2} = e_{1, 2} / (lambda^2 * c_f), where e_{1, 2} is
  the edge correction term, lambda^2 is the estimate of the
  intensity squared and c_f is the scaling term for the
  mark test function. This way each pair of points has
  exactly one coefficient tailored for them.
}

