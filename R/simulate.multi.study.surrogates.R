#' Title
#'
#' @param u_y_mean Mean of treatment effects on primary endpoint, bounded between 0 and 1
#' @param u_y_sd Standard deviation of treatment effects on primary endpoint
#' @param d_valid Average distance of treatment effect on valid surrogates and on primary outcome within each trial. 
#' Must satisfy \code{|d_valid| <= epsilon}.
#' @param d_invalid Average distance of treatment effect on invalid surrogates and on primary outcome within each trial.
#' Must satisfy \code{|d_invalid| > epsilon}.
#' @param u_s_sd_valid Standard deviation of distance of treatment effect on valid surrogates and on primary outcome.
#' @param u_s_sd_invalid Standard deviation of distance of treatment effect on invalid surrogates and on primary outcome.
#' @param epsilon Margin for definition of valid surrogates.
#' @param M Number of trials.
#' @param J Number of surrogate candidates.
#' @param prop_valid Proportion of surrogate candidates which are valid.
#'
#' @returns list containing elements \itemize{
#' \item \code{uy} : numeric vector of length \code{M} containing treatment effects on primary endpoint across trials
#' \item \code{us} : numeric matrix of dimension \code{M} times \code{J} containing treatment effects on each of \code{J} candidate markers. 
#' \item \code{hyp} : vector of length \code{J} containing the truth of surrogate validity. \code{null false} corresponds to valid surrogates
#' whereas \code{null true} corresponds to invalid surrogates. 
#' \item \code{epsilon} : value of epsilon used to define surrogate validity. 
#' }
#' @export
#'
#' @examples example.data.highdim.multistudy = simulate.multi.study.surrogates(u_y_mean = 0.8, u_y_sd = 0.1, d_valid = -0.1, d_invalid = -0.5, u_s_sd_valid = 0.1, u_s_sd_invalid = 0.1, epsilon = 0.2, M = 10, J = 500, prop_valid = 0.1)
simulate.multi.study.surrogates <- function(u_y_mean = 0.8,
                                            u_y_sd = 0.1,
                                            d_valid = -0.1,
                                            d_invalid = -0.5,
                                            u_s_sd_valid = 0.1,
                                            u_s_sd_invalid = 0.1,
                                            epsilon = 0.2,
                                            M = 10,
                                            J = 500,
                                            prop_valid = 0.1) {
  
  # Checks
  if (u_y_mean > 1 | u_y_mean < 0){
    stop(
      "u_y_mean must be between 0 and 1." 
    )
  }
  
  if (prop_valid > 0) {
    if (d_valid > 1 | d_valid < -1) {
      stop("d_valid  must be between -1 and 1.")
    }
  }
  
  if (d_invalid > 1 | d_invalid < -1){
    stop(
      "d_invalid must be between -1 and 1." 
    )
  }
  
  if (prop_valid > 1 | prop_valid < 0){
    stop(
      "prop_valid must be between 0 and 1." 
    )
  }
  
  if (epsilon > 1 | epsilon < 0){
    stop(
      "epsilon must be between 0 and 1." 
    )
  }
  
  if (M <= 1){
    stop(
      "The number of trials M must be greater than 1." 
    )
  }
  
  # Valid surrogate condition
  if (prop_valid > 0 && abs(d_valid) > epsilon) {
    stop(
      "You tried to generate some proportion of valid surrogates
         but |d_valid| is greater than epsilon.
         Either reduce d_valid or specify prop_valid = 0 to generate only invalid surrogates."
    )
  }
  
  # Invalid surrogate condition
  if (prop_valid < 1 && abs(d_invalid) <= epsilon) {
    stop(
      "You tried to generate some proportion of invalid surrogates
         but |d_invalid| is less than or equal to epsilon.
         Either increase d_invalid or specify prop_valid = 1 to generate only valid surrogates."
    )
  }
  
  # Number of valid / invalid surrogates
  J_valid <- round(prop_valid * J)
  J_invalid <- J - J_valid
  
  # Generate treatment effects for Y
  mu_nu <- qlogis(u_y_mean)
  nu_m <- rnorm(M, mean = mu_nu, sd = u_y_sd)
  u_y_m <- plogis(nu_m)
  
  # Helper function for bounding
  bound01 <- function(x)
    pmin(pmax(x, 0), 1)
  
  # Simulate valid surrogates
  us_valid <- replicate(J_valid, bound01(u_y_m + rnorm(M, d_valid, u_s_sd_valid)))
  
  # Simulate invalid surrogates
  us_invalid <- replicate(J_invalid, bound01(u_y_m + rnorm(M, d_invalid, u_s_sd_invalid)))
  
  # Combine
  us_all <- cbind(us_valid, us_invalid)
  
  colnames(us_all) <- paste0("S", 1:J)
  
  # Validity labels
  valid_label <- c(rep("null false", J_valid), rep("null true", J_invalid))
  
  # Output
  res = list(
    uy = u_y_m,
    us = us_all,
    hyp = valid_label,
    epsilon = epsilon
  )
  
  return(res)
}
