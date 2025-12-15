\name{delta.calculate}
\alias{delta.calculate}
\title{
Calculates the rank-based test statistic for Y and S and the difference, delta
}
\description{
Calculates the rank-based test statistic for Y and the rank-based test statistic for S and the difference, delta, along with corresponding standard error estimates
}
\usage{
delta.calculate(full.data = NULL, yone = NULL, yzero = NULL, sone = NULL, szero = NULL)
}
\arguments{
  \item{full.data}{either full.data or yone, yzero, sone, szero must be supplied; if full data is supplied it must be in the following format: one observation per row, Y is in the first column, S is in the second column, treatment group (0 or 1) is in the third column.
}
  \item{yone}{
	primary outcome, Y, in group 1
}
  \item{yzero}{
	primary outcome, Y, in group 0
}
  \item{sone}{
	surrogate marker, S, in group 1
}
  \item{szero}{
	surrogate marker, S, in group 0
}
}
\value{
  \item{u.y}{rank-based test statistic for Y}
  \item{u.s}{rank-based test statistic for S}
  \item{delta}{difference, u.y-u.s}
  \item{sd.u.y}{standard error estimate of u.y}
  \item{sd.u.s}{standard error estimate of u.s}
  \item{sd.delta}{standard error estimate of delta}
}
\author{
	Layla Parast
}
\examples{
data(example.data)
delta.calculate(yone = example.data$y1, yzero = example.data$y0, sone = example.data$s1, 
szero = example.data$s0)
}

