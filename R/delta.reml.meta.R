#' Function to compute and test random-effects meta-analysis summary of delta value estimates for a single marker
#'
#' @param delta numeric vector of delta values per study
#' @param sd.delta numeric vector of standard error of delta values per study
#' @param epsilon numeric non-inferiority margin for testing cross-study validity
#' @param alpha numeric significance level of test. Note : using the two-one-sided test (\code{alternative = "two.sided"})
#'   produces a (1-2\code{alpha})*100% confidence interval, so you may consider halving your 
#'   desired \code{alpha} if using this option.
#' @param alternative character giving the alternative hypothesis type for testing the summary effect.
#'   One of \code{c("less","two.sided")}, where "less" corresponds to a non-inferiority test and "two.sided"
#'   corresponds to a two one-sided test procedure. Default is "two.sided".
#' @param sample.sizes numeric vector of sample sizes. These should be 
#'   provided in case the estimated standard error for a given study is exactly 0. In order to compute 
#'   estimates of cochrane's Q and I-squared, a conservative bound for the variance is used. This bound
#'   is rooted in U-statistic theory and is taken by ignoring the covariance term when estmating the variance 
#'   for delta. 
#' @param tol numeric convergence tolerance for finding a root of the score equation
#' @param verbose logical flag indicating whether messages should be printed, defaults to \code{FALSE}
#'
#' @return a list with elements \itemize{
#'   \item \code{n.studies} : numeric, number of studies considered
#'   \item \code{tau2} : numeric, estimated tau-squared (between-study heterogeneity)
#'   \item \code{mu.delta} : numeric, estimated mean of distribution of delta
#'   \item \code{se.delta} : numeric, standard error of delta summary estimate with Hartung-Knapp adjustment
#'   \item \code{ci.delta.upper} : numeric, upper confidence interval for mean of delta. 
#'          Note : if using the non-inferiority test (i.e. \code{alternative = "less"}), 
#'          these bounds correspond to a (1-\code{alpha})*100% confidence interval,
#'          whereas the two-one-sided test (i.e. \code{alternative = "two.sided"})
#'          corresponds to a (1-2\code{alpha})*100% interval.
#'   \item \code{ci.delta.lower} : numeric, lower confidence interval for mean of delta
#'   \item \code{p.lower} : numeric, if \code{alternative} is \code{"two.sided"}, gives the p-value corresponding to
#'   testing the null hypothesis that \code{delta} is greater than \code{-epsilon}.
#'   Value is \code{NULL} if \code{alternative} is \code{"less"}.
#'   \item \code{p.upper} : numeric, if \code{alternative} is \code{"two.sided"}, gives the p-value corresponding to
#'   testing the null hypothesis that \code{delta} is less than \code{epsilon}.
#'   Value is \code{NULL} if \code{alternative} is \code{"less"}.
#'   \item \code{p} : numeric, consensus p-value for hypothesis test for either the two-one-sided test or
#'   the non-inferiorty test.
#'   \item \code{Q} : numeric, Cochran's Q-statistic for heterogeneity between studies
#'   \item \code{I2} : numeric, Higgins-Thompson I-squared statistic representing the total percentage of variation
#'   attributable to between-study heterogeneity
#'   \item \code{weights.tau} : numeric vector of raw study weights for the summary measure
#'   \item \code{weights.tau.relative} : numeric vector of relative study weights for the summary measure,
#'   such that each weight is a percentage adding to 100%
#'   \item \code{weights.tau.sum} : numeric, sum of \code{weights.tau}
#' }
#' @export
#' @author Arthur Hughes
#'
delta.reml.meta <- function(delta = NULL,
                            sd.delta = NULL,
                            epsilon = NULL,
                            alpha = 0.05,
                            alternative = "two.sided",
                            sample.sizes = NULL,
                            tol = 1e-10,
                            verbose = FALSE) {
  # Validity checks
  n.studies = length(delta)
  
  if (length(sd.delta) != n.studies) {
    stop("sd.delta not the same length as delta")
  }
  
  if (is.null(epsilon)){
    stop("epsilon (equivalence margin) must be supplied")
  }
  
  # Remove values with exactly 0 standard deviation
  if (any(sd.delta == 0)) {
    delta = delta[-which(sd.delta == 0)]
    sd.delta = sd.delta[-which(sd.delta == 0)]
    # if (verbose) {
    #   message(
    #     "Note : one or multiple studies have exactly 0 estimated standard error.
    #           Removing this study from the analysis."
    #   )
    # }
  }
  
  # Reset the value of n.studies
  n.studies = length(delta)
  
  vi <- sd.delta^2
  
  # -------------------------
  # REML estimation of tau^2
  # -------------------------
  # weighted mean function for a given tau2
  mu.delta.of <- function(tau2) {
    w <- 1 / (vi + tau2)
    sum(w * delta) / sum(w)
  }
  
  # restricted score function: root gives REML tau^2
  score.fn <- function(tau2) {
    v <- vi + tau2
    w <- 1 / v
    mu.delta <- sum(w * delta) / sum(w)
    sum((delta - mu.delta)^2 / (v^2)) - sum(1 / v)
  }
  
  # negative restricted log-likelihood (minimize this if no root found)
  neg.restricted.loglik <- function(tau2) {
    v <- vi + tau2
    wsum <- sum(1 / v)
    mu.delta  <- sum((1 / v) * delta) / wsum
    0.5 * (sum(log(v)) + log(wsum) + sum((delta - mu.delta)^2 / v))
  }
  
  # choose an upper bound for tau2 search
  s2 <- var(delta)
  upper <- max(abs(s2 - mean(vi)), s2, max(vi), 1e-6)
  upper <- abs(upper) * 10 + 1e-6
  if (upper <= 0){
    upper <- 1e-6
  }
  
  # try uniroot on [0, upper] if a sign change exists
  tau2.hat <- NA
  
  s0 <- score.fn(0)
  su <- score.fn(upper)
  
  if (is.finite(s0) && is.finite(su) && s0 * su < 0) {
    ur <- tryCatch(
      uniroot(
        score.fn,
        lower = 0,
        upper = upper,
        tol = tol
      ),
      error = function(e)
        e
    )
    if (!inherits(ur, "error")) {
      tau2.hat <- max(0, ur$root)
      converged <- TRUE
    } else {
      converged <- FALSE
      if (verbose){
        message("uniroot failed; falling back to optimize on restricted log-likelihood")
      }
    }
  } else {
    converged = FALSE
  }
  
  if (!converged) {
    opt <- optimize(
      neg.restricted.loglik,
      lower = 0,
      upper = upper,
      tol = tol
    )
    tau2.hat <- max(0, opt$minimum)
  }
  
  # -------------------------
  # pooled estimates & Hartung-Knapp standard error
  # -------------------------
  w.tau <- 1 / (vi + tau2.hat)
  mu.delta.hat <- sum(w.tau * delta) / sum(w.tau)
  var.conv <- 1 / sum(w.tau)
  se.conv <- sqrt(var.conv)
  
  # Compute relative weights for studies
  pct.w <- 100 * w.tau / sum(w.tau)
  
  # Hartung-Knapp scaling factor q
  q <- (1 / (n.studies - 1)) * sum(w.tau * (delta - mu.delta.hat)^2)
  # ensure q non-negative (numerical)
  q <- max(0, q)
  se.HK <- sqrt(q) * sqrt(var.conv)
  
  # HK CI (uses t_{n.studies-1})
  tcrit <- qt(1 - alpha, df = n.studies - 1)
  ci.HK <- c(mu.delta.hat - tcrit * se.HK, mu.delta.hat + tcrit * se.HK)
  
  # -------------------------
  # Two-one-sided test (TOST) using HK SE
  # -------------------------
  if (alternative == "two.sided") {
    T.L <- (mu.delta.hat + epsilon) / se.HK
    T.U <- (mu.delta.hat - epsilon) / se.HK
    p.lower <- 1 - pt(T.L, df = n.studies - 1)   # p-value for H0L: mu.delta <= -epsilon
    p.upper <- pt(T.U, df = n.studies - 1)       # p-value for H0U: mu.delta >= +epsilon
    p.final <- max(p.lower, p.upper)
  } else {
    T.U <- (mu.delta.hat - epsilon) / se.HK
    p.final <- pt(T.U, df = n.studies - 1)
    p.lower = NULL
    p.upper = NULL
  }
  
  
  # -------------------------
  # Cochran's Q and I^2
  # -------------------------
  # Fixed-effect weights use vi only
  
  # In the case where variance exactly 0, replace this with a conservative estimate based on asymptotic U-statistic theory
  # when ignoring the covariance term
  # The variance of U-statistics is derived in the RISE paper
  
  if (any(vi == 0)){
    variance.bound = 1/(2*sample.sizes)
    
    vi[vi == 0] = variance.bound[vi == 0]
  }
  
  w0 <- 1 / vi
  w0.sum <- sum(w0)
  delta.FE <- sum(w0 * delta) / w0.sum
  Q <- sum(w0 * (delta - delta.FE)^2)
  # handle degenerate Q
  if (is.finite(Q) && Q > (n.studies - 1)) {
    I2 <- max(0, (Q - (n.studies - 1)) / Q) * 100
  } else {
    I2 <- max(0, (Q - (n.studies - 1)) / max(Q, 1e-12)) * 100
    I2 <- max(0, I2)
  }
  
  # -------------------------
  # Prediction intervals
  # -------------------------
  # predictive standard error: HK variance of the pooled mean + between-study variance
  var_pred <- (se.HK^2) + tau2.hat
  se.pred <- sqrt(var_pred)
  
  # degrees of freedom for PI: M-2 
  df_pi <- n.studies - 2
  tcrit_pi <- qt(1 - alpha / 2, df = df_pi)
  
  pi.lower <- mu.delta.hat - tcrit_pi * se.pred
  pi.upper <- mu.delta.hat + tcrit_pi * se.pred
  
  # -------------------------
  # prepare output
  # -------------------------
  results <- list(
    n.studies = n.studies,
    tau2 = tau2.hat,
    mu.delta = mu.delta.hat,
    se.delta = se.HK,
    ci.delta.upper = ci.HK[2],
    ci.delta.lower = ci.HK[1],
    pi.lower = pi.lower,
    pi.upper = pi.upper,
    p.lower = p.lower,
    p.upper = p.upper,
    p = p.final,
    Q = Q,
    I2 = I2,
    weights.tau = w.tau,
    weights.tau.relative = pct.w,
    weights.tau.sum = sum(w.tau)
  )
  
  return(list(results = results))
}
