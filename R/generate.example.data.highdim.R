#' Generate individual participant data for high-dimensional surrogate candidates and response
#'
#' Generates individual participant data for high-dimensional surrogate candidates using one
#' of two data generating processes, as described in Hughes A et al (2025) <doi:10.1002/sim.70241>.
#'
#' @param n1 positive numeric giving the sample size in the treated group
#' @param n0 positive numeric giving the sample size in the untreated group
#' @param p positive numeric giving the number of markers to generate
#' @param prop_valid numeric between 0 and 1 (inclusive) giving the proportion of surrogate
#' candidates to generate as valid.
#' @param valid_sigma non-negative numeric giving the standard deviation for valid candidates
#' @param corr non-negative numeric giving the correlation between the surrogate candidates
#' @param mode character taking values in c("simple", "complex"). If "simple", generates all variables with
#' (multivariate) normal distributions. Else, uses a more complex exponential distribution.
#' @param y0_mean numeric giving the mean of the primary endpoint in the untreated group
#' @param y0_sd non-negative numeric giving the standard deviation of the primary endpoint in the untreated group
#' @param y1_mean numeric giving the mean of the primary endpoint in the treated group
#' @param y1_sd non-negative numeric giving the standard deviation of the primary endpoint in the treated group
#' @param s0_mean numeric giving the mean of the surrogate candidates in the untreated group
#' @param s0_sd non-negative numeric giving the standard deviation of the surrogate candidates in the untreated group
#' @param s1_mean numeric giving the mean of the surrogate candidates in the treated group
#' @param s1_sd non-negative numeric giving the standard deviation of the surrogate candidates in the treated group
#' @param seed numeric giving a seed for reproducibility
#'
#' @import dplyr, MASS
#' @return A list with the following components:
#' \describe{
#'   \item{y1}{vector containing primary endpoint values in treated group}
#'   \item{y0}{vector containing primary endpoint values in untreated group}
#'   \item{s1}{n1 times p matrix containing surrogate candidate values in treated group}
#'   \item{s0}{n0 times p matrix containing surrogate candidate values in untreated group}
#'   \item{hyp}{character vector giving the truth behind the null hypothesis for each surrogate candidate}
#' }
#' 
#' @examples
#' res <- generate.example.data.highdim(n1 = 25, n0 = 25, p = 500, prop_valid = 1)
#' dim(res$s1)       # 25 x 500
#' @export
generate.example.data.highdim = function(n1,
                                         n0,
                                         p,
                                         prop_valid,
                                         valid_sigma = 1,
                                         corr = 0,
                                         mode = "simple",
                                         y0_mean = 0,
                                         y0_sd = 1,
                                         y1_mean = 3,
                                         y1_sd = 1,
                                         s0_mean = 0,
                                         s0_sd = 1,
                                         s1_mean = 0,
                                         s1_sd = 1,
                                         seed = 12345) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # If simple data generating process requested
  if (mode == "simple") {
    # Compute number of valid and invalid surrogates
    p_valid = prop_valid * p  %>%
      as.numeric()
    p_invalid = (1 - prop_valid) * p %>%
      as.numeric()
    
    # generate primary responses from multivariate normal distributions
    y1 = rnorm(n1, y1_mean, y1_sd) # treated response
    y0 = rnorm(n0, y0_mean, y0_sd) # untreated repsonse
    
    # generate candidate surrogates
    ## generate random invalid surrogate means from uniform distribution
    mm = runif(p_invalid, min = 0.5, max = 2.5)
    ## generate random invalid surrogate standard deviations from uniform distribution
    ss = runif(p_invalid, min = 0.5, max = 2)
    ## Construct variance-covariance matrix with correlation as off-diagonal elements
    Sigma_invalid = matrix(corr, nrow = p_invalid, ncol = p_invalid)
    ## Diagonal elements are standard deviations
    diag(Sigma_invalid) = ss
    
    
    if (prop_valid != 0) {
      # if valid surrogates required
      
      # ensure the covariance matrix is valid when all surrogates are requested valid
      if (p_valid == 1) {
        Sigma_valid = c(valid_sigma)
        s1.valid = y1 +  MASS::mvrnorm(n = n1,
                                 mu = 0,
                                 Sigma = Sigma_valid)
        s0.valid = y0 +  MASS::mvrnorm(n = n0,
                                 mu = 0,
                                 Sigma = Sigma_valid)
      } else if (p_valid > 1) {
        Sigma_valid = matrix(corr * valid_sigma, nrow = p_valid, ncol = p_valid)
        diag(Sigma_valid) = rep(valid_sigma, p_valid)
        
        # Generate valid surrogates in treated group by perturbing primary response
        s1.valid = matrix(y1,
                          nrow = n1,
                          ncol = p_valid,
                          byrow = TRUE) +
          MASS::mvrnorm(n = n1,
                  mu = rep(0, p_valid),
                  Sigma = Sigma_valid)
        
        # valid surrogates in untreated group by perturbing primary response
        s0.valid = matrix(y0,
                          nrow = n0,
                          ncol = p_valid,
                          byrow = TRUE) +
          MASS::mvrnorm(n = n0,
                  mu = rep(0, p_valid),
                  Sigma = Sigma_valid)
      }
      
      # invalid surrogates
      if (prop_valid != 1) {
        # if invalid surrogates required
        if (p_invalid > 1) {
          s1.invalid = MASS::mvrnorm(n = n1,
                               mu = mm,
                               Sigma = Sigma_invalid)
          s0.invalid = MASS::mvrnorm(n = n0,
                               mu = mm,
                               Sigma = Sigma_invalid)
        } else {
          mm = runif(1, min = 0.5, max = 2.5) # generate random invalid means, sds
          ss = runif(1, min = 0.5, max = 2)
          s1.invalid = rnorm(n = n1,
                             mean = mm,
                             sd = ss)
          s0.invalid = rnorm(n = n0,
                             mean = mm,
                             sd = ss)
        }
        # bind candidates together
        s1 = cbind(s1.valid, s1.invalid)
        s0 = cbind(s0.valid, s0.invalid)
      } else {
        s1 = s1.valid
        s0 = s0.valid
      }
      
    } else {
      # if no valid surrogates required
      # invalid surrogates
      s1 = MASS::mvrnorm(n = n1,
                   mu = mm,
                   Sigma = Sigma_invalid)
      s0 = MASS::mvrnorm(n = n0,
                   mu = mm,
                   Sigma = Sigma_invalid)
    }
    # store hypothesis truths
    hyp = c(rep("null false", p_valid), rep("null true", p_invalid))
    
    # return generated data as a list
    return(list(
      y1 = y1,
      y0 = y0,
      s1 = s1,
      s0 = s0,
      hyp = hyp
    ))
  } else if (mode == "complex") {
    # calculate the number of valid and invalid surrogates
    p_valid = prop_valid * p
    p_invalid = (1 - prop_valid) * p
    
    # generate primary responses from multivariate normal distributions
    y1 = rnorm(n1, y1_mean, y1_sd) # treated response
    y0 = rnorm(n0, y0_mean, y0_sd) # untreated repsonse
    
    # generate candidate surrogates
    lambda = runif(p_invalid, min = 0.5, max = 2.5)
    if (prop_valid != 0) {
      #if valid surrogates required
      # valid surrogate covariance matrix
      Sigma_valid = matrix(corr * valid_sigma, nrow = p_valid, ncol = p_valid)
      diag(Sigma_valid) = rep(valid_sigma, p_valid)
      s1.valid = matrix(y1^3,
                        nrow = n1,
                        ncol = p_valid,
                        byrow = TRUE) +
        MASS::mvrnorm(n = n1,
                mu = rep(0, p_valid),
                Sigma = Sigma_valid)
      
      s0.valid = matrix(y0^3,
                        nrow = n0,
                        ncol = p_valid,
                        byrow = TRUE) +
        MASS::mvrnorm(n = n0,
                mu = rep(0, p_valid),
                Sigma = Sigma_valid)
      
      if (prop_valid != 1) {
        #if invalid surrogates required
        s0.invalid = sapply(lambda, function(rate)
          rexp(n0, rate))
        s1.invalid = sapply(lambda, function(rate)
          rexp(n1, rate))
        s1 = cbind(s1.valid, s1.invalid)
        s0 = cbind(s0.valid, s0.invalid)
      } else {
        s1 = s1.valid
        s0 = s0.valid
      }
    } else {
      s0 = sapply(lambda, function(rate)
        rexp(n0, rate))
      s1 = sapply(lambda, function(rate)
        rexp(n1, rate))
    }
    hyp = c(rep("null false", p_valid), rep("null true", p_invalid))
    return(list(
      y1 = y1,
      y0 = y0,
      s1 = s1,
      s0 = s0,
      hyp = hyp
    ))
  }
}
