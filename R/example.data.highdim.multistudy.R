#' High‑dimensional, multi-study surrogate candidate example dataset
#'
#' A simulated high‑dimensional, multi-study dataset for demonstrating the RISE-meta methodology
#' implemented in \pkg{SurrogateRank}. The data contains primary response and
#' 1000 surrogate candidates from 10 studues of varying sample sizes.
#' 10% of the surrogate candidates are "valid", where the strength of the surrogacy depends on the study.
#'
#' @format A list containing :
#' \describe{
#'   \item{y1}{primary response in treated}
#'   \item{y0}{primary response in untreated}
#'   \item{s1}{\code{P} surrogate candidates in treated}
#'   \item{s0}{\code{P} surrogate candidates in untreated}
#'   \item{study}{vector of \code{K} study names}
#'   \item{hyp}{vector of \code{P} hypotheses where \code{null false} if the surrogate is valid}
#' }
#'
#' @usage data("example.data.highdim.multistudy", package = "SurrogateRank")
#'
#' @examples
#' data("example.data.highdim.multistudy", package = "SurrogateRank")
#' head(example.data.highdim.multistudy)
#'
#' @source Simulated for package examples.
#' @docType data
#' @keywords datasets
"example.data.highdim.multistudy"
