\name{test.surrogate}
\alias{test.surrogate}
\title{
	Tests whether the surrogate is valid
}
\description{
	Calculates the rank-based test statistic for Y and the rank-based test statistic for S and the difference, delta, along with corresponding standard error estimates, then tests whether the surrogate is valid
}
\usage{
test.surrogate(full.data = NULL, yone = NULL, yzero = NULL, sone = NULL, 
szero = NULL, epsilon = NULL, power.want.s = 0.7, u.y.hyp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{full.data}{
	either full.data or {yone, yzero, sone, szero} must be supplied; if full data is supplied it must be in the following format: one observation per row, Y is in the first column, S is in the second column, treatment group (0 or 1) is in the third column.
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
  \item{epsilon}{
	threshold to use for delta, default calculates epsilon as a function of desired power for S
}
  \item{power.want.s}{
   desired power for S, default is 0.7
}
  \item{u.y.hyp}{
	hypothesized value of u.y used in the calculation of epsilon, default uses estimated valued of u.y
}
}
\value{
  \item{u.y}{rank-based test statistic for Y}
  \item{u.s}{rank-based test statistic for S}
  \item{delta}{difference, u.y-u.s}
  \item{sd.u.y}{standard error estimate of u.y}
  \item{sd.u.s}{standard error estimate of u.s}
  \item{sd.delta}{standard error estimate of delta}
  \item{ci.delta}{1-sided confidence interval for delta}
  \item{epsilon.used}{the epsilon value used for the test}
  \item{is.surrogate}{logical, TRUE if test indicates S is a good surrogate, FALSE otherwise}
	}
\author{
	Layla Parast}
\examples{
data(example.data)
test.surrogate(yone = example.data$y1, yzero = example.data$y0, sone = example.data$s1, 
szero = example.data$s0)
}
