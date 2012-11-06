\name{summ_func_random_labelling}
\alias{summ_func_random_labelling}
\title{Calculates summary functions for a pattern and its simulations under the
hypothesis of random labelling.}
\usage{
  summ_func_random_labelling(pattern,
    edge_corr_func = translational_correction,
    mtf_name = c("1", "m", "mm", "gamma", "gammaAbs", "mor", "morAbs"),
    n_perm = 999L, r_max = NULL, r_vec = NULL,
    do_besags_L = TRUE, method = "permute", ...)
}
\value{
  An array containing the summary functions for the
  original pattern and the randomly labelled patterns.
  Dimensions: orig_and_perm, summ_func, r.
}
\description{
  Calculates summary functions for a pattern and its
  simulations under the hypothesis of random labelling.
}

