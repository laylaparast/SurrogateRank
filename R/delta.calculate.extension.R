#' Calculate Delta: Difference in Rank-based Statistics for Two Outcomes
#'
#' This function estimates the difference (`delta`) between two rank-based statistics 
#' (e.g., Wilcoxon statistics or paired ranks) for a primary outcome and a surrogate, 
#' under either an independent or paired design.
#'
#' @description
#' This function calculates the difference in treatment effects on a univariate marker 
#' and on a continuous primary response. This extends the \code{delta.calculate()} function 
#' from the \code{SurrogateRank} package to the case where samples may be paired instead of 
#' independent, and where a two sided test is desired.
#' 
#' @param yone               numeric vector of primary response values in the treated group.
#' @param yzero              numeric vector of primary response values in the untreated group.
#' @param sone               matrix or dataframe of surrogate candidates in the treated group 
#'                           with dimension \code{n1 x p} where n1 is the number of treated samples 
#'                           and p the number of candidates. Sample ordering must match exactly 
#'                           \code{yone}. 
#' @param szero              matrix or dataframe of surrogate candidates in the untreated group 
#'                           with dimension \code{n0 x p} where n0 is the number of untreated 
#'                           samples and p the number of candidates. Sample ordering must match 
#'                           exactly \code{yzero}.
#' @param paired             logical flag giving if the data is independent or paired. If 
#'                           \code{FALSE} (default), samples are assumed independent. If 
#'                           \code{TRUE}, samples are assumed to be from a paired design. The 
#'                           pairs are specified by matching the rows of \code{yone} and 
#'                           \code{sone} to the rows of \code{yzero} and \code{szero}.
#'
#' @return A list with the following elements: 
#' \itemize{
#'   \item{\code{u.y}: Rank-based test statistic for the primary outcome}
#'   \item{\code{u.s}: Rank-based test statistic for the surrogate}
#'   \item{\code{delta.estimate}: Estimated difference between outcome and surrogate statistics}
#'   \item{\code{sd.u.y}: Standard deviation of the outcome statistic}
#'   \item{\code{sd.u.s}: Standard deviation of the surrogate statistic}
#'   \item{\code{sd.delta}: Standard error of the delta estimate}
#' }
#' 
#' @import dplyr
#' @export
#' @author Arthur Hughes, Layla Parast
#'
#' @examples
#' # Load data
#' data("example.data")
#' yone = example.data$y1
#' yzero = example.data$y0
#' sone = example.data$s1
#' szero = example.data$s0
#' delta.calculate.extension.result = delta.calculate.extension(
#'   yone, yzero, sone, szero, paired = TRUE
#' )
delta.calculate.extension = function(yone, 
                                     yzero,
                                     sone, 
                                     szero,
                                     paired = FALSE){
  
  # Validity checks
  ## Check same number of samples in primary response and surrogates
  n0 = length(yzero)
  n1 = length(yone)
  
  if (length(szero) != n0){
    stop("szero does not have the same number of samples as yzero.")
  }
  
  if (length(sone) != n1){
    stop("sone does not have the same number of samples as yone.")
  }
  
  ## if in paired mode, yone/sone must have exactly the same number of samples as yzero/szero
  if (paired){
    if (length(yone) != length(yzero)){
      stop("Paired mode is requested but the number of samples in yone does not match that of yzero.")
    } else if (length(sone) != length(szero)){
      stop("Paired mode is requested but the number of samples in sone does not match that of szero.")
    }
  }
  
  # Independent case
  if (!paired) {
    # Use Wilcoxon rank-sum test for independent samples
    test.y = wilcox.test(yone, yzero, exact = FALSE)
    test.s = wilcox.test(sone, szero, exact = FALSE)
    
    # Normalize test statistics to estimate U statistics
    n1.f = length(yone)
    n0.f = length(yzero)
    u.y = (n1.f * n0.f)^(-1) * test.y$statistic
    u.s = (n1.f * n0.f)^(-1) * test.s$statistic
    delta.estimate = u.y - u.s
    
    # Estimate the variance using Hajek projection variance estimators
    m.count = n1.f
    n.count = n0.f
    
    # Variance components from treated group comparisons
    V10.Xi.Y = sapply(yone, var.wil, b = yzero)
    V10.Xi.S = sapply(sone, var.wil, b = szero)
    
    # Variance components from control group comparisons (flip ranks)
    V01.Yj.Y = sapply(yzero, var.wil, b = yone, flip = TRUE)
    V01.Yj.S = sapply(szero, var.wil, b = sone, flip = TRUE)
    
    # Covariance matrix components from treated group
    s10.11.YY = var(V10.Xi.Y)
    s10.12.YS = cov(V10.Xi.Y, V10.Xi.S)
    s10.22.SS = var(V10.Xi.S)
    s10.21.SY = cov(V10.Xi.Y, V10.Xi.S)  # same as s10.12
    
    # Covariance matrix components from control group
    s01.11.YY = var(V01.Yj.Y)
    s01.12.YS = cov(V01.Yj.Y, V01.Yj.S)
    s01.22.SS = var(V01.Yj.S)
    s01.21.SY = cov(V01.Yj.Y, V01.Yj.S)
    
    # Combine into covariance matrices
    S10 = matrix(c(s10.11.YY, s10.12.YS, s10.21.SY, s10.22.SS), nrow = 2)
    S01 = matrix(c(s01.11.YY, s01.12.YS, s01.21.SY, s01.22.SS), nrow = 2)
    
    # Total variance matrix
    S.mat = (1/m.count) * S10 + (1/n.count) * S01
    
    # Standard deviations for outcome and surrogate
    sd.y = sqrt(S.mat[1, 1])
    sd.s = sqrt(S.mat[2, 2])
    
    # Linear combination L = (1, -1) to compute variance of delta = u.y - u.s
    L = t(as.matrix(c(1, -1)))
    sd.est = sqrt(L %*% S.mat %*% t(L))
    
  } else if (paired){
    # For paired samples, use within-pair comparisons
    n = length(yone)
    
    # Proportion of pairs where yone > yzero
    u.y = mean(ifelse(yone > yzero, 1, 0))
    u.s = mean(ifelse(sone > szero, 1, 0))
    delta.estimate = u.y - u.s
    
    # Estimate variances
    sd.y = sd(ifelse(yone > yzero, 1, 0))
    sd.s = sd(ifelse(sone > szero, 1, 0))
    
    # Difference in indicators for each pair
    di = ifelse(yone > yzero, 1, 0) - ifelse(sone > szero, 1, 0)
    var_delta = var(di) / n
    sd.est = sqrt(var_delta)
  }
  
  # Return a list of results
  return(list(
    u.y = as.numeric(u.y),
    u.s = as.numeric(u.s),
    delta.estimate = as.numeric(delta.estimate),
    sd.u.y = sd.y,
    sd.u.s = sd.s,
    sd.delta = as.numeric(sd.est)
  ))
}
