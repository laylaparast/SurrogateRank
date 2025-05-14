#' Function to perform RISE : Two-Stage Rank-Based Identification of High-Dimensional Surrogate Markers
#'
#' @description
#' RISE (Rank-Based Identification of High-Dimensional Surrogate Markers) is a two-stage method to identify 
#' and evaluate high-dimensional surrogate candidates of a continuous response.
#'
#' In the first stage (called screening), the high-dimensional candidates are screened one-by-one to identify
#' strong candidates. Strength of surrogacy is assessed through a rank-based measure of the similarity in 
#' treatment effects on a candidate surrogate and the primary response. P-values corresponding to hypothesis 
#' testing on this measure are corrected for the high number of statistical tests performed.
#'
#' In the second stage (called evaluation), candidates with an adjusted p-value below a given significance level are evaluated
#' by combining them into a single synthetic marker. The surrogacy of this marker is then assessed with the 
#' univariate test as described before.
#'
#' To avoid overfitting, the two stages are performed on separate data.
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
#' @param p.correction       character. Method for p-value adjustment (see \code{p.adjust()} function).
#'                           Defaults to the Benjamini-Hochberg method (\code{"BH"}).
#' @param n.cores            numeric giving the number of cores to commit to parallel computation in order to improve
#'                           computational time through the \code{pbmcapply()} function. Defaults to \code{1}.
#' @param alternative        character giving the alternative hypothesis type. One of \code{c("less","two.sided")}, where "less" 
#'                           corresponds to a non-inferiority test and "two.sided" corresponds to a two one-sided test procedure. Default is 
#'                           "less".
#' @param paired             logical flag giving if the data is independent or paired. 
#'                           If \code{FALSE} (default), samples are assumed independent. 
#'                           If \code{TRUE}, samples are assumed to be from a paired design. The pairs are specified by matching the rows of 
#'                           \code{yone} and \code{sone} to the rows of \code{yzero} and \code{szero}.
#' @param screen.proportion  numeric in (0,1) - proportion of data to be used for the screening stage. 
#'                           The default is \code{2/3}. If \code{1} is given, screening and evaluation will be performed on the same data.
#' @param return.all.screen  logical flag. If \code{TRUE} (default), a dataframe will be returned giving the screening results for 
#'                           all candidates. Else, only the significant candidates will be returned.
#' @param return.all.evaluate logical flag. If \code{TRUE} (default), a dataframe will be returned giving the evaluation of each
#'                           individual marker passed to the evaluation stage.
#' @param return.plot.evaluate logical flag. If \code{TRUE} (default), a ggplot2 object will be returned allowing the user to 
#'                           visualise the association between the composite surrogate on the individual-scale.
#' @param evaluate.weights   logical flag. If \code{TRUE} (default), the composite surrogate is constructed with weights as the 
#'                           absolute value of the inverse of the \delta values of each candidate, such that surrogates which are predicted
#'                           to be stronger receive more weight.
#'
#' @return a list with \itemize{
#'   \item \code{screening.results}: a list with \itemize{
#'     \item \code{screening.metrics} : dataframe of screening results (for each candidate marker - delta, CI, sd,epsilon, p-values).
#'     \item \code{significant_markers}: character vector of markers with \code{p_adjusted < alpha}.
#'   }
#'   \item \code{evaluate.results}: a list with \itemize{
#'     \item \code{individual.metrics} if \code{return.all.evaluate}=\code{TRUE}, a dataframe of evaluation results for each significant marker.
#'     \item \code{gamma.s} a list with elements \code{gamma.s.one} and  \code{gamma.s.zero}, giving the 
#'           combined surrogate marker in the treated and untreated groups, respectively.
#'     \item \code{gamma.s.evaluate} : a dataframe giving the evaluation of \code{gamma.s}
#'     \item \code{gamma.s.plot} : a ggplot2 plot showing \code{gamma.s} against the primary response on the rank-scale.
#'   }
#' }
#'
#' @import dplyr
#' @import pbmcapply
#' @import ggplot2
#' @author Arthur Hughes
#' 
#' 
#' 
#' @examples
#' # Load high-dimensional example data
#' data("example.data.highdim")
#' yone = example.data.highdim$y1
#' yzero = example.data.highdim$y0
#' sone = example.data.highdim$s1
#' szero = example.data.highdim$s0
#' 
#' rise.result = test.surrogate.rise(yone, yzero, sone, szero, power.want.s = 0.8)
#' 

test.surrogate.rise = function(yone, 
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
                               screen.proportion = 0.66, 
                               return.all.screen = TRUE,
                               return.all.evaluate = TRUE,
                               return.plot.evaluate = TRUE,
                               evaluate.weights = TRUE
){
  
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
  
  # STEP 0 : Sample splitting
  ## Assign participants to either screening or evaluation
  ### If in the independent sample setting, we must split the same proportion in both treatment groups
  if (!paired){
    split.indices.one = sample(seq_len(n1), 
                               size = screen.proportion * n1)
    split.indices.zero = sample(seq_len(n0), 
                                size = screen.proportion * n0)
    
    sone.screen = sone[split.indices.one,]
    szero.screen = szero[split.indices.zero,]
    yone.screen = yone[split.indices.one]
    yzero.screen = yzero[split.indices.zero]
    
    sone.evaluate = sone[-split.indices.one,]
    szero.evaluate = szero[-split.indices.zero,]
    yone.evaluate = yone[-split.indices.one]
    yzero.evaluate = yzero[-split.indices.zero]
    
  } else if (paired){
    split.indices = sample(seq_len(n1), 
                           size = screen.proportion * n1)
    
    sone.screen = sone[split.indices,]
    szero.screen = szero[split.indices,]
    yone.screen = yone[split.indices]
    yzero.screen = yzero[split.indices]
    
    sone.evaluate = sone[-split.indices,]
    szero.evaluate = szero[-split.indices,]
    yone.evaluate = yone[-split.indices]
    yzero.evaluate = yzero[-split.indices]
  }
  
  if(screen.proportion == 1){
    sone.screen = sone
    szero.screen = szero
    yone.screen = yone
    yzero.screen = yzero
    
    sone.evaluate = sone
    szero.evaluate = szero
    yone.evaluate = yone
    yzero.evaluate = yzero
  }
  
  # STEP 1 : SCREENING
  
  print(paste0("Screening step in progress."))
  
  screening.results = rise.screen(yone = yone.screen, 
                                  yzero = yzero.screen, 
                                  sone = sone.screen, 
                                  szero = szero.screen,
                                  alpha,
                                  power.want.s,
                                  epsilon, 
                                  u.y.hyp, 
                                  p.correction,
                                  n.cores, 
                                  alternative,
                                  paired, 
                                  return.all.screen
  )
  
  # STEP 2 : Evaluation
  ## Filter the data to only contain the strong candidates
  sone.evaluate = sone.evaluate %>% 
    as.data.frame() %>% 
    dplyr::select(any_of(screening.results$significant.markers))
  
  szero.evaluate = szero.evaluate %>% 
    as.data.frame() %>% 
    dplyr::select(any_of(screening.results$significant.markers))
  
  print(paste0("Evaluation step in progress."))
  
  evaluation.results = rise.evaluate(yone = yone.evaluate, 
                                     yzero = yzero.evaluate, 
                                     sone = sone.evaluate, 
                                     szero = szero.evaluate,
                                     alpha,
                                     power.want.s,
                                     epsilon, 
                                     u.y.hyp, 
                                     p.correction,
                                     n.cores, 
                                     alternative,
                                     paired, 
                                     return.all.evaluate,
                                     return.plot.evaluate,
                                     screening.weights = screening.results$screening.weights,
                                     evaluate.weights)
  
  return(list(screening.results = screening.results,
              evaluation.results = evaluation.results))
}


