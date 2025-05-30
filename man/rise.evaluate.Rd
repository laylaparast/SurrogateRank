% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rise.evaluate.R
\name{rise.evaluate}
\alias{rise.evaluate}
\title{Performs the evaluation stage of RISE: Two-Stage Rank-Based Identification of
High-Dimensional Surrogate Markers}
\usage{
rise.evaluate(
  yone,
  yzero,
  sone,
  szero,
  alpha = 0.05,
  power.want.s = NULL,
  epsilon = NULL,
  u.y.hyp = NULL,
  p.correction = "BH",
  n.cores = 1,
  alternative = "less",
  paired = FALSE,
  return.all.evaluate = TRUE,
  return.plot.evaluate = TRUE,
  evaluate.weights = TRUE,
  screening.weights = NULL,
  markers = NULL
)
}
\arguments{
\item{yone}{numeric vector of primary response values in the treated group.}

\item{yzero}{numeric vector of primary response values in the untreated group.}

\item{sone}{matrix or dataframe of surrogate candidates in the treated group
with dimension \code{n1 x p} where n1 is the number of treated samples
and p the number of candidates. Sample ordering must match exactly
\code{yone}.}

\item{szero}{matrix or dataframe of surrogate candidates in the untreated group
with dimension \code{n0 x p} where n0 is the number of untreated samples
and p the number of candidates. Sample ordering must match exactly
\code{yzero}.}

\item{alpha}{significance level for determining surrogate candidates. Default is
\code{0.05}.}

\item{power.want.s}{numeric in (0,1) - power desired for a test of treatment effect based
on the surrogate candidate. Either this or \code{epsilon} argument must
be specified.}

\item{epsilon}{numeric in (0,1) - non-inferiority margin for determining surrogate
validity. Either this or \code{power.want.s} argument must be specified.}

\item{u.y.hyp}{hypothesised value of the treatment effect on the primary response on
the probability scale. If not given, it will be estimated based on the
observations.}

\item{p.correction}{character. Method for p-value adjustment (see \code{p.adjust()}
function). Defaults to the Benjamini-Hochberg method (\code{"BH"}).}

\item{n.cores}{numeric giving the number of cores to commit to parallel computation
in order to improve computational time through the \code{pbmcapply()}
function. Defaults to \code{1}.}

\item{alternative}{character giving the alternative hypothesis type. One of
\code{c("less","two.sided")}, where "less" corresponds to a
non-inferiority test and "two.sided" corresponds to a two one-sided test
procedure. Default is "less".}

\item{paired}{logical flag giving if the data is independent or paired. If
\code{FALSE} (default), samples are assumed independent. If \code{TRUE},
samples are assumed to be from a paired design. The pairs are specified
by matching the rows of \code{yone} and \code{sone} to the rows of
\code{yzero} and \code{szero}.}

\item{return.all.evaluate}{logical flag. If \code{TRUE} (default), a dataframe will be returned
giving the evaluation of each individual marker passed to the
evaluation stage.}

\item{return.plot.evaluate}{logical flag. If \code{TRUE} (default), a ggplot2 object will be
returned allowing the user to visualise the association between the
composite surrogate on the individual-scale.}

\item{evaluate.weights}{logical flag. If \code{TRUE} (default), the composite surrogate is
constructed with weights as the absolute value of the inverse of the
delta values of each candidate, such that surrogates which are
predicted to be stronger receive more weight.}

\item{screening.weights}{dataframe with columns \code{marker} and \code{weight} giving the weight
in for the evaluation. Typically this is taken directly from the
screening stage as the output from the \code{rise.screen()} function.
Must be given if \code{evaluate.weights} is \code{TRUE}.}

\item{markers}{a vector of marker names (column names of szero and sone) to evaluate.
If not given, will default to evaluating all markers in the dataframes.}
}
\value{
A list with:
\itemize{
  \item \code{individual.metrics} If \code{return.all.evaluate = TRUE}, a dataframe of
  evaluation results for each significant marker.
  \item \code{gamma.s} A list with elements \code{gamma.s.one} and \code{gamma.s.zero}, giving
  the combined surrogate marker in the treated and untreated groups, respectively.
  \item \code{gamma.s.evaluate} A dataframe giving the evaluation of \code{gamma.s}.
  \item \code{gamma.s.plot} A ggplot2 plot showing \code{gamma.s} against the primary response
  on the rank-scale.
}
}

\description{
A set of high-dimensional surrogate candidates are evaluated jointly. Strength of surrogacy
is assessed through a rank-based measure of the similarity in treatment effects on a candidate
surrogate and the primary response.
}
\examples{
# Load high-dimensional example data
data("example.data.highdim")
yone <- example.data.highdim$y1
yzero <- example.data.highdim$y0
sone <- example.data.highdim$s1
szero <- example.data.highdim$s0
\donttest{
rise.evaluate.result <- rise.evaluate(yone, yzero, sone, szero, power.want.s = 0.8)}
}
\author{
Arthur Hughes
}
