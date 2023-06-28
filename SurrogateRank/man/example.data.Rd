\name{example.data}
\alias{example.data}
\docType{data}
\title{
	Example data
}
\description{
	Example data use to illustrate the functions
}
\usage{data("example.data")}
\format{
 A list with 4 elements representing 25 observations from a treatment group (group 1) and 25 observations from a control group (group 0):
  \describe{
    \item{\code{y1}}{the primary outcome,Y, in group 1}
    \item{\code{y0}}{the primary outcome, Y, in group 0}
    \item{\code{s1}}{the surrogate marker, S, in group 1}
    \item{\code{s0}}{the surrogate marker, S, in group 0}
}
}
\examples{
data(example.data)
}

