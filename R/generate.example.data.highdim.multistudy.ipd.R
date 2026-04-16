#' Generate multi-study individual participant data for high-dimensional surrogate candidates and response
#'
#' Generates individual participant data for high-dimensional surrogate candidates using one
#' of two data generating processes, as described in Hughes A et al (2025) <doi:10.1002/sim.70241>.
#'
#' @param M number of studies
#' @param n1 positive numeric giving the sample size in the treated groups
#' @param n0 positive numeric giving the sample size in the untreated groups
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
#'   \item{study1}{study names for treated samples}
#'   \item{study0}{study names for untreated samples}
#'   \item{hyp}{character vector giving the truth behind the null hypothesis for each surrogate candidate}
#' }
#' 
#' @examples
#' res <- generate.example.data.highdim.multistudy.ipd(
#' M = 5,
#' n1 = 25,
#' n0 = 25,
#' p = 500,
#' prop_valid = 1
#' )
#' dim(res$s1)       # (5 studies x 25 individuals = 125) x 500
#' @export
generate.example.data.highdim.multistudy.ipd = function(M,
                                                        n1,
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
  # Set global seed
  set.seed(seed)
  
  # Generate a different seed for each study
  seeds = sample(seq(1, 100000), size = M)
  
  data = list()
  for (m in 1:M) {
    # generate data for each study
    data[[m]] = generate.example.data.highdim(
      n1,
      n0,
      p,
      prop_valid,
      valid_sigma,
      corr,
      mode,
      y0_mean,
      y0_sd,
      y1_mean,
      y1_sd,
      s0_mean,
      s0_sd,
      s1_mean,
      s1_sd,
      seed = seeds[m]
    )
    
    data[[m]][["study1"]] = rep(paste0("study", m), n1)
    data[[m]][["study0"]] = rep(paste0("study", m), n0)
  }
  
  combined <- list(
    y1 = unlist(lapply(data, `[[`, "y1"), use.names = FALSE),
    y0 = unlist(lapply(data, `[[`, "y0"), use.names = FALSE),
    
    s1 = do.call(rbind, lapply(data, `[[`, "s1")),
    s0 = do.call(rbind, lapply(data, `[[`, "s0")),
    
    study1 = unlist(lapply(data, `[[`, "study1"), use.names = FALSE),
    study0 = unlist(lapply(data, `[[`, "study0"), use.names = FALSE),
    
    hyp = data[[1]]$hyp
  )
  
  colnames(combined$s0) = paste0("marker", seq(1:p))
  colnames(combined$s1) = paste0("marker", seq(1:p))
  
  return(combined)
}