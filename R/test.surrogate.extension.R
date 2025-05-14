#' Function to test for trial-level surrogacy of a single marker extended to the paired, two sided test setting
#' 
#' @description
#' This function tests for surrogacy of a univariate marker with respect to a continuous primary response.
#' This extends the \code{test.surrogate()} function from the \code{SurrogateRank} package to the case where
#' samples may be paired instead of independent, and where a two sided test is desired.
#' 
#' @param yone               numeric vector of primary response values in the treated group.
#' @param yzero              numeric vector of primary response values in the untreated group.
#' @param sone               matrix or dataframe of surrogate candidates in the treated group 
#'                           with dimension \code{n1 x p} where n1 is the number of treated samples and p the number of candidates.
#'                           Sample ordering must match exactly yone. 
#' @param szero              matrix or dataframe of surrogate candidates in the untreated group 
#'                           with dimension \code{n0 x p} where n0 is the number of untreated samples and p the number of candidates.
#'                           Sample ordering must match exactly yzero.
#' @param alpha              significance level for determining surrogate candidates. Default is \code{0.05}.
#' @param power.want.s       numeric in (0,1) - power desired for a test of treatment effect based on the surrogate candidate. 
#'                           Either this or \code{epsilon} argument must be specified.
#' @param epsilon            numeric in (0,1) - non-inferiority margin for determining surrogate validity.
#'                           Either this or \code{power.want.s} argument must be specified.
#' @param u.y.hyp            hypothesised value of the treatment effect on the primary response on the probability
#'                           scale. If not given, it will be estimated based on the observations.
#' @param alternative        character giving the alternative hypothesis type. One of \code{c("less","two.sided")}, where "less" 
#'                           corresponds to a non-inferiority test and "two.sided" corresponds to a two one-sided test procedure. Default is 
#'                           "less".
#' @param paired             logical flag giving if the data is independent or paired. 
#'                           If \code{FALSE} (default), samples are assumed independent. 
#'                           If \code{TRUE}, samples are assumed to be from a paired design. The pairs are specified by matching the rows of 
#'                           \code{yone} and \code{sone} to the rows of \code{yzero} and \code{szero}. 
#' 
#'
#'@return A list containing:
#' \itemize{
#'   \item{\code{u.y}}{Estimated rank-based treatment effect on the outcome.}
#'   \item{\code{u.s}}{Estimated rank-based treatment effect on the surrogate.}
#'   \item{\code{delta.estimate}}{Estimated difference in treatment effects: \code{u.y - u.s}.}
#'   \item{\code{sd.u.y}}{Standard deviation of \code{u.y}.}
#'   \item{\code{sd.u.s}}{Standard deviation of \code{u.s}.}
#'   \item{\code{sd.delta}}{Standard deviation of \code{delta.estimate}.}
#'   \item{\code{ci.delta}}{One-sided confidence interval upper bound for \code{delta.estimate}.}
#'   \item{\code{p.delta}}{p-value for validity of trial-level surrogacy.}
#'   \item{\code{epsilon.used}}{Non-inferiority threshold used in the test.}
#'   \item{\code{is.surrogate}}{\code{TRUE} if the surrogate passes the test, else \code{FALSE}.}
#' }
#' 
#' @import dplyr 
#' @author Arthur Hughes, Layla Parast
#' 
#' @examples
#' # Load data
#' data("example.data")
#' yone = example.data$y1
#' yzero = example.data$y0
#' sone = example.data$s1
#' szero = example.data$s0
#' 
#' test.surrogate.extension.result = test.surrogate.extension(yone, yzero, sone, szero, power.want.s = 0.8, paired = TRUE, alternative = "two.sided")

test.surrogate.extension = function(yone, 
                                    yzero,
                                    sone, 
                                    szero,
                                    alpha = 0.05, 
                                    power.want.s = NULL,
                                    epsilon = NULL,
                                    u.y.hyp = NULL,
                                    alternative ="less",
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
  
  ## Check that either epsilon or power.want.s is specified
  if (is.null(epsilon) & is.null(power.want.s)){
    stop("Must specify either epsilon or power.want.s.")
  }
  
  # Compute treatment effects and standard errors using delta.calculate()
  dd = delta.calculate.extension(
    yone = yone, 
    yzero = yzero, 
    sone = sone,
    szero = szero, 
    paired = paired)
  
  # If epsilon is not provided, estimate it based on power.want.s
  if (is.null(epsilon) & !paired) {
    n1 = length(yone)
    n0 = length(yzero)
    
    # Null SD under Wilcoxon/Mann-Whitney assumptions
    sd.null = sqrt((n1 + n0 + 1) / (12 * n1 * n0))
    
    # Adjust quantiles for two-sided test and power calculation
    z.alpha.2 = qnorm(1 - (alpha / 2))
    u.s.power = 0.5 - (qnorm(1 - power.want.s) - z.alpha.2) * sd.null
    
    # Estimate epsilon as the difference between observed u.y (or hypothesized) and u.s under power
    if (is.null(u.y.hyp)) {
      epsilon = dd$u.y - u.s.power
    } else {
      epsilon = u.y.hyp - u.s.power
    }
  } else if (is.null(epsilon) & paired){
    # If epsilon is not provided, estimate it based on power.want.s
    n = length(yone)
    
    # Null SD under paired-setting assumptions
    ## First estimate the probability of ties in the responses 
    pi = mean(yone == yzero)
    ## Now estimate the null
    sd.null = sqrt((1-pi)/(4*n))
    
    # Adjust quantiles for two-sided test and power calculation
    z.alpha.2 = qnorm(1 - (alpha / 2))
    u.s.power = 0.5 - (qnorm(1 - power.want.s) - z.alpha.2) * sd.null
    
    # Estimate epsilon as the difference between observed u.y (or hypothesized) and u.s under power
    if (is.null(u.y.hyp)) {
      epsilon = dd$u.y - u.s.power
    } else {
      epsilon = u.y.hyp - u.s.power
    }
  }
  
  if (alternative == "less"){
    # Compute upper bound of one-sided (1 - alpha) CI for delta = u.y - u.s
    z.alpha = qnorm(1 - alpha)
    ci.delta = c(-1, dd$delta.estimate + z.alpha * dd$sd.delta)
    
    # Decision rule for non-inferiority:
    # If upper bound of CI is less than epsilon, surrogate is acceptable
    is.surrogate = ci.delta[2] < epsilon
    
    # Compute a corresponding p-value 
    p = pnorm(dd$delta.estimate, 
              epsilon, 
              dd$sd.delta)
    
  } else if (alternative == "two.sided"){
    n = length(yone)
    # Compute (1 - 2*alpha) CI for delta = u.y - u.s
    z.alpha = qnorm(1 - alpha)
    ci.delta = c(dd$delta.estimate - z.alpha * dd$sd.delta, 
                 dd$delta.estimate + z.alpha * dd$sd.delta)
    
    # Calculate p-value corresponding to null: delta > epsilon
    p1 = pnorm(dd$delta.estimate, 
               epsilon, 
               dd$sd.delta)
    
    # Calculate p-value corresponding to null: delta < -epsilon
    p2 = pnorm(dd$delta.estimate, 
               -epsilon, 
               dd$sd.delta, 
               lower.tail = F)
    
    
    
    p = max(p1,p2)
    
    is.surrogate = p < alpha
  }
  
  return(list(
    u.y = as.numeric(dd$u.y),
    u.s = as.numeric(dd$u.s),
    delta.estimate = as.numeric(dd$delta.estimate),
    sd.u.y = dd$sd.u.y,
    sd.u.s = dd$sd.u.s,
    sd.delta = as.numeric(dd$sd.delta),
    ci.delta = ci.delta,
    p.delta = p,
    epsilon.used = epsilon,
    is.surrogate = is.surrogate
  ))
  
}