\name{specialise_weight}
\alias{specialise_weight}
\title{Create an individual weight vector for each mark test function.}
\usage{
  specialise_weight(weight_vec, mark_stat, mtf_name,
    one_per_lambda2)
}
\value{
  A matrix of weights. Each column is for different mark
  test functions. The columns are named. Each row is for a
  different _ordered_ point pair where (x1, x2) is
  different from (x2, x1).
}
\description{
  Create an individual weight vector for each mark test
  function. Bake in 1 / \lambda ^ 2 and 1 / c_f so that all
  weighing can be handled by one multiplication per pair of
  points per mark test function.
}

