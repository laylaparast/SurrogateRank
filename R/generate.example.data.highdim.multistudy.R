#' Generate high-dimensional multi-study surrogate marker trial-level effects
#'
#' Generates simulated trial-level treatment effects for multiple
#' surrogate markers across multiple studies, including both valid and invalid
#' surrogates. This function implements a hierarchical random-effects model:
#' true trial-level effects are drawn from marker-specific means with
#' between-trial heterogeneity, and observed trial-level effects include
#' additional within-study sampling error.
#'
#' @param epsilon Numeric in (0,1). Defines the region of validity for the
#'   surrogate marker means. Markers with mean discrepancy within
#'   \code{[-epsilon, epsilon]} are valid; others are invalid.
#' @param M Integer. Number of trials (studies) to simulate. Must be > 1.
#' @param sample_sizes Numeric vector of length \code{M}. Sample size for
#'   each trial. Used to compute within-study variances.
#' @param J Integer. Total number of markers to simulate (valid + invalid).
#' @param prop_valid Numeric, between 0 and 1. Proportion of markers that are valid.
#' @param u_tau_min Numeric >= 0. Lower bound of marker-specific between-trial
#'   heterogeneity variance (\eqn{\tau_j^2}).
#' @param u_tau_max Numeric >= u_tau_min. Upper bound of marker-specific
#'   between-trial heterogeneity variance (\eqn{\tau_j^2}).
#' @param u_nu_min Numeric > 0. Lower bound of marker-specific variance
#'   component (\eqn{\nu_j}) used to scale within-study sampling error.
#' @param u_nu_max Numeric >= u_nu_min. Upper bound of marker-specific variance
#'   component (\eqn{\nu_j}) used to scale within-study sampling error.
#' @param prop_invalid_under Numeric, between 0 and 1. Probability that an invalid
#'   marker underestimates the treatment effect on Y.
#' @param invalid_at_boundary default \code{FALSE}, meaning invalid surrogates are generated
#'   uniformly across the entire invalid region. If \code{TRUE}, generates invalid surrogates
#'   at the boundary values defined by epsilon. This is the worst-case scenario for invalid
#'   surrogates and thus is useful for checking the calibration of the method.
#' @param invalid_mean_discrete vector of discrete numeric values to sample true
#' means of valid surrogates at. These values must be greater or equal in absolute value than epsilon.
#' @param valid_mean_discrete vector of discrete numeric values to sample true
#' means of valid surrogates at. These values must be smaller in absolute value than epsilon.
#' @param seed numeric giving a seed for reproducibility
#'
#' @return A list with the following components:
#' \describe{
#'   \item{delta}{M x J matrix of observed trial-level discrepancies
#'        (\eqn{\hat{\delta}_{m,j}}) including sampling error.}
#'   \item{sd.delta}{M x J matrix of within-study standard deviations
#'        (\eqn{\sigma_{m,j}}).}
#'   \item{n}{Numeric vector of sample sizes for each trial.}
#'   \item{hyp}{Character vector of length J, "null true" for valid markers and
#'        "null false" for invalid markers.}
#'   \item{mu.true}{Numeric vector of true marker-level mean discrepancies
#'        (\eqn{\mu_{\delta,j}}).}
#'   \item{tau2.true}{Numeric vector of marker-specific between-trial
#'        heterogeneity variances (\eqn{\tau_j^2}).}
#' }
#'
#' @details
#' The function first draws marker-level parameters:
#' \eqn{\mu_{\delta,j}} from the validity or invalidity region, \eqn{\tau_j^2}
#' from a uniform distribution, and \eqn{\nu_j} from a uniform distribution.
#' Then, for each trial, true trial-level effects are drawn as
#' \eqn{\delta_{m,j}^{true} \sim N(\mu_{\delta,j}, \tau_j^2)}, and
#' observed effects include independent within-study sampling error
#' \eqn{\hat{\delta}_{m,j} \sim N(\delta_{m,j}^{true}, \nu_j / n_m)}.
#'
#' @examples
#' res <- generate.example.data.highdim.multistudy(
#'   epsilon = 0.2,
#'   M = 5,
#'   sample_sizes = c(25, 50, 100, 150, 250),
#'   J = 500,
#'   prop_valid = 0.1
#' )
#' dim(res$delta)       # 5 x 500
#' head(res$mu.true)
#'
#' @export
generate.example.data.highdim.multistudy <- function(epsilon = 0.2,
                                            M = 5,
                                            sample_sizes = c(25, 50, 100, 150, 250),
                                            J = 500,
                                            prop_valid = 0.1,
                                            u_tau_min = 0.01,
                                            # interpreted as variance lower bound
                                            u_tau_max = 0.1,
                                            # interpreted as variance upper bound
                                            u_nu_min = 0.01,
                                            u_nu_max = 0.1,
                                            prop_invalid_under = 0.5,
                                            invalid_at_boundary = FALSE,
                                            invalid_mean_discrete = NULL,
                                            valid_mean_discrete = NULL,
                                            seed = 12345) {
  # Set seed for reproducibility
  set.seed(seed)
  
  ## --- input checks
  if (prop_valid < 0 || prop_valid > 1) {
    stop("prop_valid must be between 0 and 1.")
  }
  
  if (epsilon <= 0 || epsilon >= 1) {
    stop("epsilon should be in (0,1).")
  }
  
  if (M <= 1) {
    stop("Number of trials M must be > 1.")
  }
  
  if (length(sample_sizes) != M) {
    stop("sample_sizes must be length M.")
  }
  
  if (J < 1) {
    stop("J must be >= 1.")
  }
  
  if (!is.null(valid_mean_discrete)) {
    if (any(abs(valid_mean_discrete) > epsilon)) {
      stop("All values in valid_mean_discrete must be smaller or equal in absolute value than epsilon.")
    }
  }
  
  if (!is.null(invalid_mean_discrete)) {
    if (any(abs(invalid_mean_discrete) < epsilon)) {
      stop("All values in invalid_mean_discrete must be greater or equal in absolute value than epsilon.")
    }
  }
  
  ## --- how many valid / invalid markers
  J_valid <- round(prop_valid * J)
  J_invalid <- J - J_valid
  
  ## --- draw marker-specific heterogeneity variances (tau^2) and convert to sd
  tau2_j <- runif(J, min = u_tau_min, max = u_tau_max)   # tau^2 for each marker
  tau_sd_j <- sqrt(tau2_j)                               # sd for rnorm
  
  ## --- draw marker-specific noise components nu_j
  nu_j <- runif(J, min = u_nu_min, max = u_nu_max)
  
  ## --- create within-study variances matrix: sigma_{m,j}^2 = nu_j / n_m
  inv_n <- 1 / sample_sizes
  sigma_m_j <- outer(inv_n, nu_j, FUN = "*")   # M x J matrix
  
  ## --- assign mu_delta for markers (first valid then invalid)
  mu_j <- numeric(J)
  if (J_valid > 0) {
    if (is.null(valid_mean_discrete)) {
      mu_j[1:J_valid] <- runif(J_valid, min = -epsilon, max = epsilon)
    } else {
      mu_j[1:J_valid] <- sample(valid_mean_discrete, J_valid, replace = TRUE)
    }
  }
  
  if (J_invalid > 0) {
    invalid_idx <- (J_valid + 1):J
    s_j <- rbinom(J_invalid, size = 1, prob = prop_invalid_under)  # 1 => underestimates (positive mu)
    
    if (any(s_j == 0)) {
      if (invalid_at_boundary) {
        mu_j[invalid_idx[s_j == 0]] = -epsilon
      } else if (!is.null(invalid_mean_discrete)) {
        mu_j[invalid_idx[s_j == 0]] <- sample(
          invalid_mean_discrete[invalid_mean_discrete < 0],
          size = sum(s_j == 0),
          replace = TRUE
        )
      } else {
        mu_j[invalid_idx[s_j == 0]] <- runif(sum(s_j == 0), min = -1, max = -epsilon)  # negative -> overestimates
      }
    }
    
    if (any(s_j == 1)) {
      if (invalid_at_boundary) {
        mu_j[invalid_idx[s_j == 1]] = epsilon
      } else if (!is.null(invalid_mean_discrete)) {
        mu_j[invalid_idx[s_j == 1]] <- sample(
          invalid_mean_discrete[invalid_mean_discrete > 0],
          size = sum(s_j == 1),
          replace = TRUE
        )
      } else {
        mu_j[invalid_idx[s_j == 1]] <- runif(sum(s_j == 1), min = epsilon, max = 1)    # positive -> underestimates
      }
    }
  }
  
  ## --- simulate true delta and observed delta
  delta_true <- matrix(NA_real_, nrow = M, ncol = J)
  delta_estimated <- matrix(NA_real_, nrow = M, ncol = J)
  
  for (j in seq_len(J)) {
    delta_true[, j] <- rnorm(n = M, mean = mu_j[j], sd = tau_sd_j[j])
    delta_estimated[, j] <- rnorm(n = M,
                                  mean = delta_true[, j],
                                  sd = sqrt(sigma_m_j[, j]))
  }
  
  ## --- assemble outputs
  colnames(delta_estimated) <- paste0("S", seq_len(J))
  rownames(delta_estimated) <- paste0("Trial ", seq_len(M))
  colnames(sigma_m_j) <- paste0("S", seq_len(J))
  rownames(sigma_m_j) <- paste0("Trial ", seq_len(M))
  
  hyp <- c(rep("null false", J_valid), rep("null true", J_invalid))
  sd_delta <- sqrt(sigma_m_j)
  
  res <- list(
    delta = delta_estimated,
    sd.delta = sd_delta,
    n = sample_sizes,
    hyp = hyp,
    mu.true = mu_j,
    tau2.true = tau2_j,
    invalid_at_boundary = invalid_at_boundary
  )
  
  return(res)
}
