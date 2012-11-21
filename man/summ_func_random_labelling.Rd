\name{summ_func_random_labelling}
\alias{summ_func_random_labelling}
\title{Calculates summary functions for a pattern and its simulations under the
hypothesis of random labelling.}
\usage{
  summ_func_random_labelling(pattern,
    edge_correction = "translate",
    mtf_name = c("1", "m", "mm", "gamma", "gammaAbs", "morAbs"),
    n_perm = 999L, r_max = NULL, r_vec = NULL,
    do_besags_L = TRUE, method = "permute", ...)
}
\arguments{
  \item{pattern}{A \code{\link[spatstat]{ppp}} object as
  the simple marked point pattern to be analysed. The marks
  need to be in the form of a numeric vector. The window
  has to have the type "rectangle".}

  \item{edge_corr_func}{The name of the edge correction to
  be used. Options are 'translate' and 'none'.}

  \item{mtf_name}{A vector of mark test function names. "1"
  stands for the unmarked K-function.}

  \item{n_perm}{The number of permutations.}

  \item{r_max}{A positive scalar value representing the
  maximum radius that should be considered. r_vec overrides
  r_max. By default, r_max is NULL and will get a sensible
  default.}

  \item{r_vec}{A monotonically increasing vector of
  non-negative r-values to act as the endpoints of the bins
  for the K_f-functions. r_vec overrides r_max. The bins
  are exclusive on the left and inclusive on the right. If
  the first vector element has value zero, it will be
  regarded as the collapsed bin [0, 0], and the next bin
  will start from and exclude 0.}

  \item{do_besags_L}{A boolean describing whether Besag's
  L-function should also be returned where available.}

  \item{method}{The name of the method to create
  simulations under the null hypothesis.}

  \item{...}{Currently unused.}
}
\value{
  An array containing the summary function estimates for
  the original pattern and the randomly labelled patterns.
  Dimensions: orig_and_perm, summ_func, r. The estimates
  for the given marked point pattern are on the row named
  "original".
}
\description{
  Calculates summary functions for a pattern and its
  simulations under the hypothesis of random labelling.
}

