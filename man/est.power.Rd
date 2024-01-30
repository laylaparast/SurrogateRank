\name{est.power}
\alias{est.power}
\title{
	Estimated power to detect a valid surrogate
}
\description{
	Calculates the estimated power to detect a valid surrogate given a total sample size and specified alternative
}
\usage{
est.power(n.total, rho = 0.8, u.y.alt, delta.alt, power.want.s = 0.7)
}
\arguments{
  \item{n.total}{
	total sample size in study
}
  \item{rho}{
	rank correlation between Y and S in group 0, default is 0.8
}
  \item{u.y.alt}{
	specified alternative for u.y
}
  \item{delta.alt}{
	specified alternative for u.s
}
  \item{power.want.s}{
	desired power for u.s, default is 0.7
}
}
\value{
estimated power
}
\author{
	Layla Parast
}
\examples{
est.power(n.total = 50, rho = 0.8, u.y.alt=0.9, delta.alt = 0.1)
}