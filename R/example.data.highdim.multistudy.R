#' High-dimensional, multi-study surrogate candidate example dataset
#'
#' A simulated high-dimensional, multi-study dataset for demonstrating the RISE-meta methodology
#' implemented in \pkg{SurrogateRank}, generated with the generate.example.data.highdim.multistudy() function.
#' The data contains treatment effect measures on the primary endpoint and on 500 surrogate candidates,
#' where the first 50 of these candidates are "valid" surrogates.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{uy}{Numeric vector of length \code{M} containing treatment effects on the primary endpoint across trials.}
#'   \item{us}{Numeric matrix of dimension \code{M} times \code{J} containing treatment effects on each of the \code{J} candidate markers.}
#'   \item{hyp}{Vector of length \code{J} containing the truth of surrogate validity. \code{null false} corresponds to valid surrogates, whereas \code{null true} corresponds to invalid surrogates.}
#'   \item{epsilon}{Value of \code{epsilon} used to define surrogate validity.}
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
