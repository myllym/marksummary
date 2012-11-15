\name{consider_radius}
\alias{consider_radius}
\title{Handle things considering the radius vector.}
\usage{
  consider_radius(pattern, r_max, r_vec)
}
\value{
  A list containing the radius vector (r_vec), the number
  of bins (n_bin), the indexing matrix to pick the point
  pairs within r_max (nearby_arr_idx) and the indices of
  the bins in which those point pairs belong to
  distance-wise (bin_idx).
}
\description{
  Handle things considering the radius vector.
}

