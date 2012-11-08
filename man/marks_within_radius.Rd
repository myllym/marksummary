\name{marks_within_radius}
\alias{marks_within_radius}
\title{Choose the marks of those pairs where the distance between the point pair
is smaller than the maximum radius.}
\usage{
  marks_within_radius(marks, nearby_arr_idx)
}
\arguments{
  \item{marks}{A vector containing all the marks of the
  point pattern in the order of the points.}

  \item{nearby_arr_idx}{An index matrix that picks the
  points within maximum radius. A picked mark will be
  picked twice to accommodate for the point pair ordering.}
}
\value{
  A list of two mark vectors, mark1 and mark2.
}
\description{
  Well, actually we select those marks that the index
  matrix points to.
}

