#' Function to perform meta-analysis of summary statistics and hypothesis testing for a single marker
#'
#' @param delta numeric vector of delta values per study
#' @param sd.delta numeric vector of standard error of delta values per study
#' @param epsilon numeric non-inferiority margin for testing cross-study validity
#' @param alpha numeric significance level of test. Note : using the two-one-sided test (\code{alternative = "two.sided"})
#'   produces a (1-2\code{alpha})*100% confidence interval.
#' @param alternative character giving the alternative hypothesis type for testing the summary effect.
#'   One of \code{c("less","two.sided")}, where "less" corresponds to a non-inferiority test and
#'   "two.sided" corresponds to a two one-sided test procedure. Default is "two.sided".
#' @param tol numeric convergence tolerance for finding a root of the score equation
#' @param verbose logical flag indicating whether messages should be printed, defaults to \code{FALSE}
#' @param test character giving the type of test to be performed. The default is \code{knha}, corresponding to
#' variance estimation using the more conservative Hartung-Knapp estimator and performes tests with the t-distribution,
#'  whereas setting this argument to \code{z} estimates the variance with the conventional estimator and uses a normal approximation for testing.
#' @param meta.analysis.method character giving the meta-analysis method to be used. The default is \code{RE}, corresponding to
#' random-effects meta-analysis, whereas setting this argument to \code{FE} uses fixed-effects meta-analysis.
#'
#' @return a list with elements \itemize{
#'   \item \code{n.studies} : numeric, number of studies considered
#'   \item \code{tau2} : numeric, estimated tau-squared (between-study heterogeneity)
#'   \item \code{mu.delta} : numeric, estimated mean of distribution of delta
#'   \item \code{se.delta} : numeric, standard error of delta summary estimate
#'   \item \code{ci.delta.upper} : numeric, upper confidence interval for mean of delta.
#'          Note : if using the non-inferiority test (i.e. \code{alternative = "less"}),
#'          these bounds correspond to a (1-\code{alpha})*100% confidence interval,
#'          whereas the two-one-sided test (i.e. \code{alternative = "two.sided"})
#'          corresponds to a (1-2\code{alpha})*100% interval.
#'   \item \code{ci.delta.lower} : numeric, lower confidence interval for mean of delta
#'   \item \code{p.lower.ci} : numeric, if \code{alternative} is \code{"two.sided"}, gives the p-value corresponding to
#'   testing the null hypothesis that \code{delta} is less than \code{-epsilon}.
#'   Value is \code{NULL} if \code{alternative} is \code{"less"}.
#'   \item \code{p.upper.ci} : numeric, if \code{alternative} is \code{"two.sided"}, gives the p-value corresponding to
#'   testing the null hypothesis that \code{delta} is less than \code{epsilon}.
#'   Value is \code{NULL} if \code{alternative} is \code{"less"}.
#'   \item \code{p} : numeric, consensus p-value for hypothesis test for either the two-one-sided test or
#'   the non-inferiorty test.
#'   \item \code{pi.lower}, \code{pi.upper} : numeric, prediction interval for the \code{delta} of a new study.
#'          It uses the same critical value as the confidence interval, so it is a (1-2\code{alpha})*100% interval
#'          (e.g. 90% for \code{alpha = 0.05}) when \code{alternative = "two.sided"} and a one-sided
#'          (1-\code{alpha})*100% interval with \code{pi.lower = -1} when \code{alternative = "less"}. Hence
#'          \code{p.pi < alpha} if and only if the prediction interval lies within the equivalence margin.
#'   \item \code{p.lower.pi} : numeric, if \code{alternative} is \code{"two.sided"}, gives the p-value corresponding to
#'   testing the null hypothesis that the predictive distribution for a new study's \code{delta} is less than \code{-epsilon}.
#'   Value is \code{NULL} if \code{alternative} is \code{"less"}.
#'   \item \code{p.upper.pi} : numeric, if \code{alternative} is \code{"two.sided"}, gives the p-value corresponding to
#'   testing the null hypothesis that the predictive distribution for a new study's \code{delta} is greater than \code{epsilon}.
#'   Value is \code{NULL} if \code{alternative} is \code{"less"}.
#'   \item \code{p.pi} : numeric, consensus p-value for the two-one-sided test (or non-inferiority test) applied to the
#'   predictive distribution rather than the pooled mean.
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
                            tol = 1e-10,
                            verbose = FALSE,
                            test = "knha",
                            meta.analysis.method = "RE") {
  # Validity checks
  n.studies = length(delta)
  
  if (length(sd.delta) != n.studies) {
    stop("sd.delta not the same length as delta")
  }
  
  if (is.null(epsilon)) {
    stop("epsilon (equivalence margin) must be supplied")
  }
  
  if (!meta.analysis.method %in% c("RE", "FE")) {
    stop("meta.analysis.method must be one of c(\"RE\", \"FE\")")
  }
  
  # Remove values with exactly 0 standard deviation
  if (any(sd.delta == 0)) {
    delta = delta[-which(sd.delta == 0)]
    sd.delta = sd.delta[-which(sd.delta == 0)]
  }
  
  # Reset the value of n.studies
  n.studies = length(delta)

  # At least two studies with non-degenerate standard errors are needed to pool
  if (n.studies < 2) {
    warning(sprintf(
      "Only %d study(ies) with a non-zero standard error; meta-analysis not performed and NA returned.",
      n.studies
    ))
    two.sided <- alternative == "two.sided"
    return(list(results = list(
      n.studies = n.studies, tau2 = NA_real_, mu.delta = NA_real_, se.delta = NA_real_,
      ci.delta.upper = NA_real_, ci.delta.lower = NA_real_,
      p.lower.ci = if (two.sided) NA_real_ else NULL,
      p.upper.ci = if (two.sided) NA_real_ else NULL,
      pi.lower = NA_real_, pi.upper = NA_real_,
      p.lower.pi = if (two.sided) NA_real_ else NULL,
      p.upper.pi = if (two.sided) NA_real_ else NULL,
      p.pi = NA_real_, p = NA_real_, Q = NA_real_, I2 = NA_real_,
      weights.tau = NA_real_, weights.tau.relative = NA_real_, weights.tau.sum = NA_real_
    )))
  }
  
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
  if (upper <= 0) {
    upper <- 1e-6
  }
  
  # try uniroot on [0, upper] if a sign change exists
  tau2.hat <- NA
  
  if (meta.analysis.method == "RE") {
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
        if (verbose) {
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
  } else {
    tau2.hat <- 0
  }
  
  # -------------------------
  # pooled estimates
  # -------------------------
  w.tau <- if (meta.analysis.method == "RE"){
    1 / (vi + tau2.hat)
  }
  else {
    1 / vi
  }
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
  
  # -------------------------
  # CI, SE, and tests (branch on test argument)
  # -------------------------
  if (meta.analysis.method == "FE") {
    se.final <- se.conv
    zcrit <- qnorm(1 - alpha)
    if(alternative == "two.sided") {
      ci.final <- c(max(-1, mu.delta.hat - zcrit * se.final), min(1, mu.delta.hat + zcrit * se.final))
    } else {
      ci.final <- c(-1, min(1, mu.delta.hat + zcrit * se.final))
    }
    
    if (alternative == "two.sided") {
      T.L <- (mu.delta.hat + epsilon) / se.final
      T.U <- (mu.delta.hat - epsilon) / se.final
      p.lower.ci <- 1 - pnorm(T.L)
      p.upper.ci <- pnorm(T.U)
      p.final <- max(p.lower.ci, p.upper.ci)
    } else {
      T.U <- (mu.delta.hat - epsilon) / se.final
      p.final <- pnorm(T.U)
      p.lower.ci <- NULL
      p.upper.ci <- NULL
    }
  } else if (meta.analysis.method == "RE" & test == "z") {
    se.final <- se.conv
    zcrit <- qnorm(1 - alpha)
    if(alternative == "two.sided") {
      ci.final <- c(max(-1, mu.delta.hat - zcrit * se.final), min(1, mu.delta.hat + zcrit * se.final))
    } else {
      ci.final <- c(-1, min(mu.delta.hat + zcrit * se.final, 1))
    }
    
    if (alternative == "two.sided") {
      T.L <- (mu.delta.hat + epsilon) / se.final
      T.U <- (mu.delta.hat - epsilon) / se.final
      p.lower.ci <- 1 - pnorm(T.L)
      p.upper.ci <- pnorm(T.U)
      p.final <- max(p.lower.ci, p.upper.ci)
    } else {
      T.U <- (mu.delta.hat - epsilon) / se.final
      p.final <- pnorm(T.U)
      p.lower.ci <- NULL
      p.upper.ci <- NULL
    }
  } else if (meta.analysis.method == "RE" & test == "knha"){
    # Default: Hartung-Knapp, t-distribution
    se.final <- se.HK
    tcrit <- qt(1 - alpha, df = n.studies - 1)
    
    if(alternative == "two.sided") {
      ci.final <- c(max(-1, mu.delta.hat - tcrit * se.final),
                    min(1, mu.delta.hat + tcrit * se.final))
    } else {
      ci.final <- c(-1,
                    min(1,mu.delta.hat + tcrit * se.final))
    }
    
    if (alternative == "two.sided") {
      T.L <- (mu.delta.hat + epsilon) / se.final
      T.U <- (mu.delta.hat - epsilon) / se.final
      p.lower.ci <- 1 - pt(T.L, df = n.studies - 1)
      p.upper.ci <- pt(T.U, df = n.studies - 1)
      p.final <- max(p.lower.ci, p.upper.ci)
    } else {
      T.U <- (mu.delta.hat - epsilon) / se.final
      p.final <- pt(T.U, df = n.studies - 1)
      p.lower.ci <- NULL
      p.upper.ci <- NULL
    }
  }
  
  # -------------------------
  # Cochran's Q and I^2
  # -------------------------
  # Fixed-effect weights use vi only
  
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
  # predictive standard error: variance of pooled mean + between-study variance
  var_pred <- (se.final^2) + tau2.hat
  se.pred <- sqrt(var_pred)
  
  # The prediction interval is built with the same critical value as the confidence interval
  # (one-sided level alpha), so that the TOST / non-inferiority test on the predictive
  # distribution rejects at level alpha exactly when the interval lies inside the margin:
  #   two.sided : (1 - 2*alpha)*100% two-sided interval, e.g. 90% for alpha = 0.05
  #   less      : (1 - alpha)*100% one-sided interval (-1, upper], e.g. 95% for alpha = 0.05
  if (test == "z" || meta.analysis.method == "FE") {
    crit_pi <- qnorm(1 - alpha)
    p.fn <- function(q, lower.tail = TRUE) pnorm(q, lower.tail = lower.tail)
  } else {
    df_pi <- n.studies - 1
    crit_pi <- qt(1 - alpha, df = df_pi)
    p.fn <- function(q, lower.tail = TRUE) pt(q, df = df_pi, lower.tail = lower.tail)
  }

  pi.upper <- mu.delta.hat + crit_pi * se.pred
  pi.lower <- if (alternative == "two.sided") mu.delta.hat - crit_pi * se.pred else -1

  # TOST / non-inferiority test on the predictive distribution
  T.U.pi <- (mu.delta.hat - epsilon) / se.pred
  if (alternative == "two.sided") {
    T.L.pi <- (mu.delta.hat + epsilon) / se.pred
    p.lower.pi <- p.fn(T.L.pi, lower.tail = FALSE)
    p.upper.pi <- p.fn(T.U.pi)
    p.pi <- max(p.lower.pi, p.upper.pi)
  } else {
    p.pi <- p.fn(T.U.pi)
    p.lower.pi <- NULL
    p.upper.pi <- NULL
  }

  # -------------------------
  # prepare output
  # -------------------------
  results <- list(
    n.studies = n.studies,
    tau2 = tau2.hat,
    mu.delta = mu.delta.hat,
    se.delta = se.final,
    ci.delta.upper = ci.final[2],
    ci.delta.lower = ci.final[1],
    p.lower.ci = p.lower.ci,
    p.upper.ci = p.upper.ci,
    pi.lower = pi.lower,
    pi.upper = pi.upper,
    p.lower.pi = p.lower.pi,
    p.upper.pi = p.upper.pi,  
    p.pi = p.pi,  
    p = p.final,
    Q = Q,
    I2 = I2,
    weights.tau = w.tau,
    weights.tau.relative = pct.w,
    weights.tau.sum = sum(w.tau)
  )
  
  return(list(results = results))
}
