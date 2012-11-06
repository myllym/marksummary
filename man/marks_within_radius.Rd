\name{marks_within_radius}
\alias{marks_within_radius}
\title{Choose those marks that are within the radius.}
\usage{
  marks_within_radius(marks, nearby_arr_idx, ordering_idx)
}
\arguments{
  \item{marks}{A vector containing all the marks of the
  point pattern in the order of the points.}

  \item{nearby_arr_idx}{A logical vector that picks the
  points within maximum radius. A picked mark will be
  picked twice to accommodate for the ordered point pairs.}

  \item{ordering_idx}{An integer vector to order the picked
  marks in the required order.}
}
\value{
  A vector of only the needed marks, ordered.
}
\description{
  Choose those marks that are within the radius.
}

