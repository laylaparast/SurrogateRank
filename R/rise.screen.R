#' Function to perform the screening stage of RISE : Two-Stage Rank-Based Identification of
#' High-Dimensional Surrogate Markers
#'
#' @description
#' A set of high-dimensional surrogate candidates are screened one-by-one to identify strong
#' candidates. Strength of surrogacy is assessed through a rank-based measure of the similarity in
#' treatment effects on a candidate surrogate and the primary response. P-values corresponding to
#' hypothesis testing on this measure are corrected for the high number of statistical tests
#' performed.
#'
#' @param yone numeric vector of primary response values in the treated group.
#' @param yzero numeric vector of primary response values in the untreated group.
#' @param sone matrix or dataframe of surrogate candidates in the treated group with dimension
#'   \code{n1 x p} where n1 is the number of treated samples and p the number of candidates. Sample
#'   ordering must match exactly yone.
#' @param szero matrix or dataframe of surrogate candidates in the untreated group with dimension
#'   \code{n0 x p} where n0 is the number of untreated samples and p the number of candidates. Sample
#'   ordering must match exactly yzero.
#' @param alpha significance level for determining surrogate candidates. Default is \code{0.05}.
#' @param power.want.s numeric in (0,1) - power desired for a test of treatment effect based on the
#'   surrogate candidate. Either this or \code{epsilon} argument must be specified.
#' @param epsilon numeric in (0,1) - non-inferiority margin for determining surrogate validity. Either
#'   this or \code{power.want.s} argument must be specified.
#' @param u.y.hyp hypothesised value of the treatment effect on the primary response on the probability
#'   scale. If not given, it will be estimated based on the observations.
#' @param p.correction character. Method for p-value adjustment (see \code{p.adjust()} function).
#'   Defaults to the Benjamini-Hochberg method (\code{"BH"}).
#' @param n.cores numeric giving the number of cores to commit to parallel computation in order to
#'   improve computational time through the \code{pbmcapply()} function. Defaults to \code{1}.
#' @param alternative character giving the alternative hypothesis type. One of
#'   \code{c("less","two.sided")}, where "less" corresponds to a non-inferiority test and "two.sided"
#'   corresponds to a two one-sided test procedure. Default is "less".
#' @param paired logical flag giving if the data is independent or paired. If \code{FALSE} (default),
#'   samples are assumed independent. If \code{TRUE}, samples are assumed to be from a paired design.
#'   The pairs are specified by matching the rows of \code{yone} and \code{sone} to the rows of
#'   \code{yzero} and \code{szero}.
#' @param return.all.screen logical flag. If \code{TRUE} (default), a dataframe will be returned giving
#'   the screening results for all candidates. Else, only the significant candidates will be returned.
#'
#' @return a list with elements \itemize{
#'   \item \code{screening.metrics} : dataframe of screening results (for each candidate marker - delta,
#'     CI, sd, epsilon, p-values).
#'   \item \code{significant.markers}: character vector of markers with \code{p_adjusted < alpha}
#'   \item \code{screening.weights}: dataframe giving marker names and the inverse absolute value of the
#'     associated deltas.
#' }
#'
#' @import dplyr pbmcapply
#' @export
#' @author Arthur Hughes
#'
#' @examples
#' # Load high-dimensional example data
#' data("example.data.highdim")
#' yone <- example.data.highdim$y1
#' yzero <- example.data.highdim$y0
#' sone <- example.data.highdim$s1
#' szero <- example.data.highdim$s0
#' rise.screen.result <- rise.screen(yone, yzero, sone, szero, power.want.s = 0.8)
rise.screen <- function(yone,
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
                        return.all.screen = TRUE) {
  # Data formatting
  ## Convert dataframes to numeric matrices
  if (is.data.frame(sone) | is.data.frame(szero)) {
    sone <- as.matrix(sone)
    szero <- as.matrix(szero)
  }

  # If no column names on surrogate candidates, set them as the column indices
  if (is.null(colnames(sone))) {
    colnames(sone) <- paste0("marker", 1:ncol(sone))
    colnames(szero) <- paste0("marker", 1:ncol(szero))
  }

  # Validity checks

  ## Check same number of samples in primary response and surrogates
  n0 <- length(yzero)
  n1 <- length(yone)

  if (nrow(szero) != n0) {
    stop("szero does not have the same number of samples as yzero.")
  }

  if (nrow(sone) != n1) {
    stop("sone does not have the same number of samples as yone.")
  }

  ## if in paired mode, yone/sone must have exactly the same number of samples as yzero/szero
  if (paired) {
    if (length(yone) != length(yzero)) {
      stop("Paired mode is requested but the number of samples in yone does not match that of yzero.")
    } else if (length(sone) != length(szero)) {
      stop("Paired mode is requested but the number of samples in sone does not match that of szero.")
    }
  }

  ## Check that either epsilon or power.want.s is specified
  if (is.null(epsilon) & is.null(power.want.s)) {
    stop("Must specify either epsilon or power.want.s.")
  }

  # Screen markers by applying surrogate test in parallel
  ## First define a function that we can then apply in parallel
  .test_marker <- function(idx) {
    args <- list(
      yone = yone,
      yzero = yzero,
      sone = sone[, idx],
      szero = szero[, idx],
      alternative = alternative,
      paired = paired,
      power.want.s = power.want.s,
      epsilon = epsilon
    )

    res <- do.call(test.surrogate.extension, args)

    c(
      delta = res$delta.estimate,
      ci_lower = res$ci.delta[1],
      ci_upper = res$ci.delta[2],
      sd = res$sd.delta,
      epsilon = res$epsilon.used,
      p_unadjusted = res$p.delta
    )
  }

  P <- ncol(sone)
  raw_list <- pbmclapply(1:P, .test_marker, mc.cores = n.cores)

  # Process the list results into a dataframe
  results <- do.call(rbind, raw_list) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(
      marker        = colnames(sone),
      p_adjusted    = p.adjust(p_unadjusted, method = p.correction)
    ) %>%
    dplyr::select(marker, epsilon, delta, sd, ci_lower, ci_upper, p_unadjusted, p_adjusted)

  # Retreive names of significant markers
  significant_markers <- results %>%
    filter(p_adjusted < alpha) %>%
    pull(marker)

  if (!return.all.screen) {
    results <- results %>%
      filter(marker %in% significant_markers)
  }

  # For the weights, if delta is 0, we set it to the next closest value
  min.nonzero.delta <- min(results %>%
    filter(delta != 0) %>%
    pull(delta) %>%
    abs())

  screening.weights <- results %>%
    filter(marker %in% significant_markers) %>%
    mutate(nonzero.delta = ifelse(delta != 0, delta, min.nonzero.delta)) %>%
    mutate(weight = 1 / (abs(nonzero.delta))) %>%
    dplyr::select(marker, weight)

  return(
    list(
      screening.metrics   = results,
      significant.markers = significant_markers,
      screening.weights   = screening.weights
    )
  )
}
