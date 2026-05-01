#' High‑dimensional multi-study individual participant surrogate candidate example dataset
#'
#' A simulated high‑dimensional dataset for demonstrating the RISE-Meta methodology
#' implemented in \pkg{SurrogateRank}. The data contains primary response and
#' 100 surrogate candidates from 25 treated individuals and 25 untreated individuals across
#' 5 different studies,
#' where 10% of the surrogate candidates are "valid".
#'
#' @format A list containing :
#' \describe{
#'   \item{y1}{primary response in treated}
#'   \item{y0}{primary response in untreated}
#'   \item{s1}{1000 surrogate candidates in treated}
#'   \item{s0}{1000 surrogate candidates in untreated}
#'   \item{study1}{study names for treated}
#'   \item{study0}{study names for untreated}
#'   \item{hyp}{for each surrogate, \code{null false} if the surrogate is valid}
#' }
#'
#' @usage data("example.data.highdim.multistudy.ipd", package = "SurrogateRank")
#'
#' @examples
#' data("example.data.highdim.multistudy.ipd", package = "SurrogateRank")
#' head(example.data.highdim.multistudy.ipd)
#'
#' @source Simulated for package examples.
#' @docType data
#' @keywords datasets
"example.data.highdim.multistudy.ipd"
