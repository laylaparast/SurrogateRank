#' Function to perform the evaluation stage of RISE : Two-Stage Rank-Based Identification of 
#' High-Dimensional Surrogate Markers
#' 
#' @description
#' A set of high-dimensional surrogate candidates are evaluated jointly. Strength of surrogacy 
#' is assessed through a rank-based measure of the similarity in treatment effects on a candidate 
#' surrogate and the primary response.
#' 
#' @param yone               numeric vector of primary response values in the treated group.
#' @param yzero              numeric vector of primary response values in the untreated group.
#' @param sone               matrix or dataframe of surrogate candidates in the treated group
#'                           with dimension \code{n1 x p} where n1 is the number of treated samples 
#'                           and p the number of candidates. Sample ordering must match exactly 
#'                           \code{yone}. 
#' @param szero              matrix or dataframe of surrogate candidates in the untreated group
#'                           with dimension \code{n0 x p} where n0 is the number of untreated samples 
#'                           and p the number of candidates. Sample ordering must match exactly 
#'                           \code{yzero}.
#' @param alpha              significance level for determining surrogate candidates. Default is 
#'                           \code{0.05}.
#' @param power.want.s       numeric in (0,1) - power desired for a test of treatment effect based 
#'                           on the surrogate candidate. Either this or \code{epsilon} argument must 
#'                           be specified.
#' @param epsilon            numeric in (0,1) - non-inferiority margin for determining surrogate 
#'                           validity. Either this or \code{power.want.s} argument must be specified.
#' @param u.y.hyp            hypothesised value of the treatment effect on the primary response on 
#'                           the probability scale. If not given, it will be estimated based on the 
#'                           observations.
#' @param p.correction       character. Method for p-value adjustment (see \code{p.adjust()} 
#'                           function). Defaults to the Benjamini-Hochberg method (\code{"BH"}).
#' @param n.cores            numeric giving the number of cores to commit to parallel computation 
#'                           in order to improve computational time through the \code{pbmcapply()} 
#'                           function. Defaults to \code{1}.
#' @param alternative        character giving the alternative hypothesis type. One of 
#'                           \code{c("less","two.sided")}, where "less" corresponds to a 
#'                           non-inferiority test and "two.sided" corresponds to a two one-sided test 
#'                           procedure. Default is "less".
#' @param paired             logical flag giving if the data is independent or paired. If 
#'                           \code{FALSE} (default), samples are assumed independent. If \code{TRUE}, 
#'                           samples are assumed to be from a paired design. The pairs are specified 
#'                           by matching the rows of \code{yone} and \code{sone} to the rows of 
#'                           \code{yzero} and \code{szero}.
#' @param return.all.evaluate logical flag. If \code{TRUE} (default), a dataframe will be returned 
#'                            giving the evaluation of each individual marker passed to the 
#'                            evaluation stage.
#' @param return.plot.evaluate logical flag. If \code{TRUE} (default), a ggplot2 object will be 
#'                              returned allowing the user to visualise the association between the 
#'                              composite surrogate on the individual-scale.
#' @param evaluate.weights    logical flag. If \code{TRUE} (default), the composite surrogate is 
#'                            constructed with weights as the absolute value of the inverse of the 
#'                            delta values of each candidate, such that surrogates which are 
#'                            predicted to be stronger receive more weight.
#' @param screening.weights   dataframe with columns \code{marker} and \code{weight} giving the weight 
#'                            in for the evaluation. Typically this is taken directly from the 
#'                            screening stage as the output from the \code{rise.screen()} function. 
#'                            Must be given if \code{evaluate.weights} is \code{TRUE}.
#'
#' @return a list with 
#' \itemize{
#'   \item \code{individual.metrics} if \code{return.all.evaluate}=\code{TRUE}, a dataframe of 
#'         evaluation results for each significant marker.
#'   \item \code{gamma.s} a list with elements \code{gamma.s.one} and  \code{gamma.s.zero}, giving 
#'         the combined surrogate marker in the treated and untreated groups, respectively.
#'   \item \code{gamma.s.evaluate} : a dataframe giving the evaluation of \code{gamma.s}
#'   \item \code{gamma.s.plot} : a ggplot2 plot showing \code{gamma.s} against the primary response 
#'         on the rank-scale.
#' }
#'
#' @import dplyr pbmcapply ggplot2
#' @export
#' @author Arthur Hughes
#'
#' @examples
#' # Load high-dimensional example data
#' data("example.data.highdim")
#' yone = example.data.highdim$y1
#' yzero = example.data.highdim$y0
#' sone = example.data.highdim$s1
#' szero = example.data.highdim$s0
#' rise.evaluate.result = rise.evaluate(yone, yzero, sone, szero, power.want.s = 0.8)
rise.evaluate = function(yone, 
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
                         paired  = FALSE, 
                         return.all.evaluate = TRUE,
                         return.plot.evaluate = TRUE,
                         evaluate.weights = TRUE,
                         screening.weights = NULL){
  
  # Data formatting
  ## Convert dataframes to numeric matrices
  if (is.data.frame(sone) | is.data.frame(szero)){
    sone <- as.matrix(sone)
    szero <- as.matrix(szero)
  }
  
  # If no column names on surrogate candidates, set them as the column indices
  if (is.null(colnames(sone))){
    colnames(sone) = paste0("marker",1:ncol(sone))
    colnames(szero) = paste0("marker",1:ncol(szero))
  } 
  
  # Validity checks
  
  ## Check same number of samples in primary response and surrogates
  n0 = length(yzero)
  n1 = length(yone)
  
  if (nrow(szero) != n0){
    stop("szero does not have the same number of samples as yzero.")
  }
  
  if (nrow(sone) != n1){
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
  
  ## Check that weights are provided
  if (evaluate.weights & is.null(screening.weights)){
    warning("Weighted combination requested but no weights provided. Using an unweighted sum instead.")
    evaluate.weights = FALSE
  }
  
  # Standardise and weight the candidates
  ## First combine the surrogate candidates together to standardise based on all samples
  s.combined = rbind(sone, szero)
  ## First standardise
  s.standardised <- scale(s.combined)
  ## Now weight
  if (evaluate.weights){
    # extract weights 
    w_vec <- setNames(screening.weights$weight, screening.weights$marker)
    
    # Now standardise by multiplying by weights
    s.standardised = sweep(
      s.standardised, 2, w_vec[colnames(s.standardised)], FUN = "*"
    )
  }
  
  # Form gamma.s by summing over the standardised, weighted predictors
  gamma.s = rowSums(s.standardised)
  
  # Now extract back the data by treatment group
  sone.standardised = s.standardised[c(1:nrow(sone)),]
  szero.standardised = s.standardised[-c(1:nrow(sone)),]
  gamma.s.one = gamma.s[c(1:nrow(sone))]
  gamma.s.zero = gamma.s[-c(1:nrow(sone))]
  
  # Test gamma.s as a surrogate 
  gamma.s.result = test.surrogate.extension(yone = yone,
                                            yzero = yzero,
                                            sone = gamma.s.one,
                                            szero = gamma.s.zero,
                                            alpha,power.want.s,
                                            epsilon,
                                            u.y.hyp,
                                            alternative,
                                            paired)
  
  gamma.s.evaluate = c(
    delta       = gamma.s.result$delta.estimate,
    ci_lower    = gamma.s.result$ci.delta[1],
    ci_upper    = gamma.s.result$ci.delta[2],
    sd          = gamma.s.result$sd.delta,
    epsilon     = gamma.s.result$epsilon.used,
    p_unadjusted= gamma.s.result$p.delta
  )
  
  ## Now make ggplot2 to visualise gamma.s on the individual scale
  if (return.plot.evaluate){
    # Make a treatment factor to colour the points
    treatment = as.factor(c(rep(1, length(yone)), 
                            rep(0, length(yzero)))
    )
    
    rank.df <- data.frame(
      treatment     = treatment,
      response.rank = rank(c(yone,yzero)),
      gamma.rank    = rank(gamma.s)
    )
    
    rho.val <- round(cor(rank.df$response.rank, rank.df$gamma.rank), 2)
    
    gamma.s.plot = rank.df %>%
      ggplot(aes(x = response.rank, y = gamma.rank, col = as.factor(treatment))) + 
      geom_point(size = 5) +
      xlab("Primary Response (rank)") +
      ylab(expression(gamma[S]~"(rank)")) +
      ggtitle("Ranks of primary response vs new surrogate in evaluation data") +
      guides(col = guide_legend(title = "Treatment")) +
      scale_color_manual(
        labels = levels(as.factor(treatment)),
        values = c("#1D8A99", "#c1121f")
      ) +
      theme_minimal() +
      geom_abline(slope = 1, col = "red", linetype = "dashed", linewidth = 1.4) +
      coord_fixed(ratio = 1) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = bquote(rho == .(round(cor(as.numeric(rank.df$response.rank), as.numeric(rank.df$gamma.rank)), 2))),
        vjust = 2, hjust = 3.3,
        color = "red",
        size = 15
      ) +
      theme_minimal(base_size = 22) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.background = element_rect(fill = 'white', color = 'white')
      )
  }
  
  # Now individually evaluate the candidate surrogates with the screen function
  if(return.all.evaluate){
    individual.metrics = rise.screen(
      yone, 
      yzero, 
      sone, 
      szero,
      alpha,
      power.want.s,
      epsilon, 
      u.y.hyp, 
      p.correction,
      n.cores, 
      alternative,
      paired, 
      return.all.screen = TRUE
    )[["screening.metrics"]]
  }
  
  return(list(
    individual.metrics        = list(sone.standardised = sone.standardised,
                                     szero.standardised = szero.standardised),
    gamma.s                   = list(gamma.s.one = gamma.s.one,
                                     gamma.s.zero = gamma.s.zero),
    gamma.s.evaluate          = gamma.s.evaluate,
    gamma.s.plot              = gamma.s.plot
  ))
}
