#' Function to perform the evaluation stage of RISE-meta : Meta-Analysis of Rank-Based Identification of
#' High-Dimensional Surrogate Markers
#'
#' @param yone numeric vector of primary response values in the treated participants
#' @param yzero numeric vector of primary response values in the untreated participants
#' @param sone matrix or dataframe of surrogate candidates in the treated group with dimension
#'   \code{n1 x p} where n1 is the number of treated samples and p the number of candidates. Sample
#'   ordering must match exactly yone.
#' @param szero matrix or dataframe of surrogate candidates in the untreated group with dimension
#'   \code{n0 x p} where n0 is the number of treated samples and p the number of candidates. Sample
#'   ordering must match exactly yzero.
#' @param studyone character vector of length \code{n1} indicating the study corresponding to each treated sample.
#'   Ordering much match yone.
#' @param studyzero character vector of length \code{n0} indicating the study corresponding to each untreated sample.
#'   Ordering much match yzero.
#' @param alpha significance level for determining valid surrogates. Default is \code{0.05}.
#' @param power.want.s.study numeric in (0,1) - power desired for a test of treatment effect based on the
#'   surrogate candidate. If \code{return.all.evaluate = TRUE}, either this or \code{epsilon.study} argument must be specified.
#' @param epsilon.study numeric in (0,1) - non-inferiority margin for determining surrogate validity in the
#'   within-study screening phase. If \code{return.all.evaluate = TRUE}, either this or \code{power.want.s.study} argument must be specified.
#' @param epsilon.meta.mode character string specifying the mode to choose the value of the acceptable margin defined
#' by epsilon. By default, this is set to "user", where the value of epsilon is fixed by the user, defined by the
#' value of the argument \code{epsilon.meta}. The alternative is to set this as "mean.power", which corresponds to
#' taking the mean value of epsilon across studies such that the power to detect departures from the null within
#' each study is defined by the \code{power.want.s.study} argument.
#' @param epsilon.meta numeric in (0,1) - non-inferiority margin for determining surrogate validity
#'   in the meta-analysis stage. Must be specified.
#' @param u.y.hyp hypothesised value of the treatment effect on the primary response on the probability
#'   scale. If not given, it will be estimated based on the observations.
#' @param p.correction character. Method for p-value adjustment (see \code{p.adjust()} function).
#'   Defaults to the Benjamini-Hochberg method (\code{"BH"}).
#' @param n.cores numeric giving the number of cores to commit to parallel computation in order to
#'   improve computational time through the \code{pbmcapply()} function. Defaults to \code{1}.
#' @param alternative character giving the alternative hypothesis type. One of
#'   \code{c("less","two.sided")}, where "less" corresponds to a non-inferiority test and "two.sided"
#'   corresponds to a two one-sided test procedure. Default is "two.sided".
#' @param test character giving the type of test to be performed. The default is \code{knha}, corresponding to
#' variance estimation using the more conservative Hartung-Knapp estimator and performes tests with the t-distribution,
#'  whereas setting this argument to \code{z} estimates the variance with the conventional estimator and uses a normal approximation for testing.
#' @param paired.all logical flag giving if the data is independent or paired. If \code{FALSE} (default),
#'   samples are assumed independent. If \code{TRUE}, all samples are assumed to be from a paired design.
#'   The pairs are specified by matching the rows of \code{yone} and \code{sone} to the rows of
#'   \code{yzero} and \code{szero}, and all studies must be paired. If only some studies are paired and
#'   others have independent samples, one may specify the \code{paired.studies} argument instead.
#' @param paired.studies character vector specifying the names of the studies in \code{studyone} or \code{studyzero}
#' which are paired, in the case where some have paired designs and others do not. By default, this is \code{NULL},
#' indicating that study designs are all specified by the \code{paired.all} argument.
#' @param evaluate.weights    logical flag. If \code{TRUE} (default), the composite surrogate is
#'                            constructed with weights such that surrogates which are
#'                            predicted to be stronger receive more weight.
#' @param screening.weights   dataframe with columns \code{marker} and \code{weight} giving the weight
#'                            in for the evaluation. Typically this is taken directly from the
#'                            screening stage as the output from the \code{rise.screen.meta()} function.
#'                            Must be given if \code{evaluate.weights} is \code{TRUE}.
#' @param markers   a vector of marker names (column names of szero and sone) to evaluate.
#'                             If not given, will default to evaluating all markers in the dataframes.
#' @param return.all.evaluate logical flag. If \code{TRUE}, a dataframe will be returned
#'                            giving the meta-analysis evaluation of each individual marker passed to the
#'                            evaluation stage. Defaults to \code{FALSE} for computational time.
#' @param weight.mode character giving the type of weighting to return to be used in case \code{return.all.evaluate = TRUE}.
#'                            See \code{rise.screen.meta()} for detail.
#' @param return.forest.plot logical flag. If \code{TRUE} (default), a forest plot of the effect sizes for the
#' combined signature across studies, with its meta-analysis summary measure, will be included in the output.
#' @param return.fit.plot logical flag. If \code{TRUE} (default), a plot of the effects on the primary response
#' versus the effects on the combined surrogate signature for each study will be included in the output.
#' @param show.pooled.effect logical flag. If \code{TRUE} (default), the forest plot will show the pooled effect
#' estimate. Otherwise, it will just show the individual trial estimates.
#' @param meta.analysis.method character giving the meta-analysis method to be used. The default is \code{RE}, corresponding to
#' random-effects meta-analysis, whereas setting this argument to \code{FE} uses fixed-effects meta-analysis.
#'
#' @return a list with elements \itemize{
#'   \item \code{individual.metrics} : if \code{return.all.evaluate}=\code{TRUE}, a list containing
#'         dataframes \code{individual.metrics.study} (per-study results for individual markers) and
#'         \code{individual.metrics.meta} (meta-analysis results for individual markers).
#'   \item \code{evaluation.metrics.study} : study-level results for the combined marker, gamma.
#'   \item \code{evaluation.metrics.meta} : meta-analysis results for the combined marker, gamma.
#'   \item \code{gamma.s} : a list with elements \code{gamma.s.one} and \code{gamma.s.zero}, giving
#'         the values of the combined surrogate marker gamma in the treated and untreated groups, respectively.
#'   \item \code{gamma.s.plot} : if \code{return.forest.plot} and/or \code{return.fit.plot} are \code{TRUE},
#'   returns evaluation plots as a list
#' }
#' @import dplyr pbmcapply ggplot2 cowplot
#' @export
#' @author Arthur Hughes
#'
#' @examples
#' # Load high-dimensional example data
# data("example.data.highdim.multistudy")
# yone <- example.data.highdim.multistudy$y1
# yzero <- example.data.highdim.multistudy$y0
# sone <- example.data.highdim.multistudy$s1
# szero <- example.data.highdim.multistudy$s0
# studyone <- example.data.highdim.multistudy$study1
# studyzero <- example.data.highdim.multistudy$study0
# rise.meta.screen.result <- rise.screen.meta(yone, yzero, sone, szero, studyone, studyzero, epsilon.study = 0.2, epsilon.meta = 0.2, n.cores = 12)
# markers = rise.meta.screen.result[["significant.markers"]]
# screening.weights = rise.meta.screen.result[["screening.weights"]]
# rise.meta.evaluate.result <- rise.evaluate.meta(yone, yzero, sone, szero, studyone, studyzero, epsilon.meta = 0.2, markers = markers, screening.weights = screening.weights, return.all.evaluate = T, epsilon.study = 0.1, n.cores = 12)
rise.evaluate.meta = function(yone,
                              yzero,
                              sone,
                              szero,
                              studyone,
                              studyzero,
                              alpha = 0.05,
                              power.want.s.study = NULL,
                              epsilon.study = NULL,
                              epsilon.meta.mode = "user",
                              epsilon.meta = NULL,
                              u.y.hyp = NULL,
                              p.correction = "BH",
                              n.cores = 1,
                              alternative = "two.sided",
                              test = "knha",
                              paired.all = FALSE,
                              paired.studies = NULL,
                              evaluate.weights = TRUE,
                              screening.weights = NULL,
                              weight.mode = "diff.epsilon",
                              markers = NULL,
                              return.all.evaluate = FALSE,
                              return.forest.plot = TRUE,
                              return.fit.plot = TRUE,
                              show.pooled.effect = TRUE,
                              meta.analysis.method = "RE") {
  # DATA FORMATTING #
  ## Convert dataframes to numeric matrices
  if (is.data.frame(sone) | is.data.frame(szero)) {
    sone <- as.matrix(sone)
    szero <- as.matrix(szero)
  }
  
  # Ensure sone and szero are matrices with at least one column
  if (is.null(dim(sone))) {
    sone <- matrix(sone, ncol = 1)
  }
  if (is.null(dim(szero))) {
    szero <- matrix(szero, ncol = 1)
  }
  
  # If they have one column but dropped dimension (edge safety)
  if (ncol(sone) == 0) {
    sone <- matrix(sone, ncol = 1)
  }
  if (ncol(szero) == 0) {
    szero <- matrix(szero, ncol = 1)
  }
  
  if (!is.null(markers)) {
    sone = sone %>%
      as.data.frame() %>%
      dplyr::select(any_of(markers))
    
    szero = szero %>%
      as.data.frame() %>%
      dplyr::select(any_of(markers))
  }
  
  # If no column names on surrogate candidates, set them as the column indices
  if (is.null(colnames(sone))) {
    colnames(sone) <- paste0("marker", 1:ncol(sone))
    colnames(szero) <- paste0("marker", 1:ncol(szero))
  }
  
  # VALIDITY CHECKS #
  
  ## Check that either epsilon.study or power.want.s.study is specified
  if (is.null(epsilon.study) & is.null(power.want.s.study)) {
    stop("Must specify either epsilon.study or power.want.s.study.")
  }
  
  # Check that epsilon.meta is specified if the user mode is desired.
  if (epsilon.meta.mode == "user" & is.null(epsilon.meta)) {
    stop("You have requsted a user-defined epsilon.meta, but you have not provided it.")
  }
  
  # Check that the power is specified if the mean.power approach for deciding epsilon is desired.
  if (epsilon.meta.mode == "mean.power" &
      is.null(power.want.s.study)) {
    stop(
      "You have requested the epsilon to be computed based on the within study power,
         but you have not provided the argument power.want.s.study."
    )
  }
  
  # Print a message if both mean.power mode is selected but epsilon.meta also specified
  if (epsilon.meta.mode == "mean.power" &
      !is.null(epsilon.meta)) {
    message(
      "You have selected the mean power approach to compute epsilon meta,
            but you have also specified epsilon.meta. Note that this value will be
            overwritten and the mean epsilon across studies defined by the power will
            be used instead."
    )
  }
  
  # If both power.want.s.study and epsilon.study, prioritise power and put a message
  if (!is.null(epsilon.study) & !is.null(power.want.s.study)) {
    message(
      "You provided both epsilon.study and power.want.s.study.
            Using the power argument to compute within-study epsilon by default.
            If you wanted to fix the epsilon yourself, set power.want.s.study = NULL."
    )
  }
  
  ## Check same number of samples in primary response and surrogates
  n0 <- length(yzero) # Number of samples in treated group
  n1 <- length(yone) # Number of samples in untreated group
  
  if (nrow(szero) != n0) {
    stop("szero does not have the same number of samples as yzero.")
  }
  
  if (nrow(sone) != n1) {
    stop("sone does not have the same number of samples as yone.")
  }
  
  if (length(studyzero) != n0) {
    stop("studyzero does not have the same number of samples as yzero.")
  }
  
  if (length(studyone) != n1) {
    stop("studyone does not have the same number of samples as yone.")
  }
  
  ## Check that weights are provided
  if (evaluate.weights & is.null(screening.weights)) {
    warning(
      "Weighted combination requested but no weights provided. Using an unweighted sum instead."
    )
    evaluate.weights <- FALSE
  }
  
  ## Retrieve the study names
  all.study.names = unique(studyone)
  
  ## if in paired mode, specified studies must have the same number of treated and untreated samples
  
  if (paired.all) {
    paired.studies = unique(studyone) # if in all-paired mode, set this argument to all study names
  }
  
  # Check all specified paired studies have the same number of treated and untreated samples
  for (study in paired.studies) {
    if (length(yone[paired.studies == study]) != length(yzero[paired.studies == study])) {
      stop(
        paste0(
          "Paired mode is requested but the number of samples in yone does not match that of yzero in study '",
          study,
          "'"
        )
      )
    }
  }
  
  # Form gamma as the weighted standardised sum of the markers within each study
  
  # Initialise results storage
  rise.evaluate.results.allstudies =  vector("list", length(all.study.names)) # prepare storage for study results
  ix <- 1L # initialise index
  
  gamma.s = vector("list", length(all.study.names))
  
  for (study in all.study.names) {
    yzero.study = yzero[studyzero == study] # Extract study samples
    yone.study = yone[studyone == study]
    
    sone.study = sone[studyone == study, , drop = FALSE]
    szero.study = szero[studyzero == study, , drop = FALSE]
    
    ## First combine the surrogate candidates together to standardise based on all samples
    s.combined.study <- rbind(sone.study, szero.study)
    ## First standardise
    s.standardised.study <- scale(s.combined.study)
    ## Now weight
    if (evaluate.weights) {
      # extract weights
      w_vec <- setNames(screening.weights$weight, screening.weights$marker)
      
      # Now standardise by multiplying by weights
      s.standardised.study <- sweep(s.standardised.study, 2, w_vec[colnames(s.standardised.study)], FUN = "*")
    }
    
    # Form gamma.s by summing over the standardised, weighted predictors
    gamma.s.study <- rowSums(s.standardised.study)
    
    # Now extract back the data by treatment group
    sone.standardised.study <- s.standardised.study[c(1:nrow(sone.study)), ]
    szero.standardised.study <- s.standardised.study[-c(1:nrow(sone.study)), ]
    gamma.s.one.study <- gamma.s.study[c(1:nrow(sone.study))]
    gamma.s.zero.study <- gamma.s.study[-c(1:nrow(sone.study))]
    
    gamma.s[[ix]] = list("gamma.s.one" = gamma.s.one.study, "gamma.s.zero" = gamma.s.zero.study)
    
    # Check if this study is paired
    if (is.null(paired.studies)) {
      paired.study = FALSE
    } else {
      paired.study = ifelse(study %in% paired.studies, TRUE, FALSE)
    }
    
    epsilon_arg <- if (is.null(power.want.s.study)) {
      epsilon.meta
    } else {
      NULL
    }
    
    # Test gamma.s as a surrogate
    gamma.s.result.study <- test.surrogate.extension(
      yone = yone.study,
      yzero = yzero.study,
      sone = gamma.s.one.study,
      szero = gamma.s.zero.study,
      alpha,
      epsilon = epsilon_arg,
      power.want.s = power.want.s.study,
      u.y.hyp,
      alternative,
      paired = paired.study
    )
    
    # Extract relevant screening results
    rise.evaluate.results.allstudies[[ix]] <- data.frame(
      "study" = study,
      "epsilon" = gamma.s.result.study$epsilon.used,
      "marker" = "gamma",
      "n" = length(yone.study) + length(yzero.study),
      "u.y" = gamma.s.result.study$u.y,
      "u.s" = gamma.s.result.study$u.s,
      "delta" = gamma.s.result.study$delta.estimate,
      "ci.lower" = gamma.s.result.study$ci.delta[1],
      "ci.upper" = gamma.s.result.study$ci.delta[2],
      "sd" = gamma.s.result.study$sd.delta,
      "p.unadjusted" = gamma.s.result.study$p.delta
    )
    
    # Increase index
    ix <- ix + 1L
    
  }
  
  names(gamma.s) = all.study.names
  
  # Bind the per-study results into a dataframe
  evaluation.metrics.study <- bind_rows(rise.evaluate.results.allstudies) %>%
    mutate(p.adjusted = p.adjust(p.unadjusted, method = p.correction))
  
  # If some studies have exactly 0 standard error, report this
  if (any(evaluation.metrics.study$sd == 0)) {
    sd0.studies = evaluation.metrics.study %>%
      filter(sd == 0) %>%
      pull(study) %>%
      unique()
    
    message(
      paste0(
        "Note: studies '",
        paste(sd0.studies, collapse = ", "),
        "' have degenerate (exactly 0) estimation of standard error",
        " for the combined marker gamma. ",
        "This is likely due to small sample size and perfect pairwise",
        " concordance between gamma and the primary endpoint.",
        " These studies will be removed from the meta-analysis."
      )
    )
    
    evaluation.metrics.study = evaluation.metrics.study %>%
      filter(sd != 0)
  }
  
  # Calculate average epsilon for meta-analysis
  if (epsilon.meta.mode == "mean.power") {
    
    epsilon.df = evaluation.metrics.study %>%
      dplyr::select(study, epsilon) %>%
      distinct()
    
    epsilon.meta = mean(epsilon.df$epsilon)
    
    message(
      "Using value ",
      round(epsilon.meta, 3),
      " for epsilon.meta, based on the mean of
            epsilon values across studies requiring specified power."
    )
    
    # Now repeat the meta-analysis using this value of epsilon across all studies
    
    # Initialise results storage
    rise.evaluate.results.allstudies =  vector("list", length(all.study.names)) # prepare storage for study results
    ix <- 1L # initialise index
    
    gamma.s = vector("list", length(all.study.names))
    
    for (study in all.study.names) {
      yzero.study = yzero[studyzero == study] # Extract study samples
      yone.study = yone[studyone == study]
      
      sone.study = sone[studyone == study, , drop = FALSE]
      szero.study = szero[studyzero == study, , drop = FALSE]
      
      ## First combine the surrogate candidates together to standardise based on all samples
      s.combined.study <- rbind(sone.study, szero.study)
      ## First standardise
      s.standardised.study <- scale(s.combined.study)
      ## Now weight
      if (evaluate.weights) {
        # extract weights
        w_vec <- setNames(screening.weights$weight, screening.weights$marker)
        
        # Now standardise by multiplying by weights
        s.standardised.study <- sweep(s.standardised.study, 2, w_vec[colnames(s.standardised.study)], FUN = "*")
      }
      
      # Form gamma.s by summing over the standardised, weighted predictors
      gamma.s.study <- rowSums(s.standardised.study)
      
      # Now extract back the data by treatment group
      sone.standardised.study <- s.standardised.study[c(1:nrow(sone.study)), ]
      szero.standardised.study <- s.standardised.study[-c(1:nrow(sone.study)), ]
      gamma.s.one.study <- gamma.s.study[c(1:nrow(sone.study))]
      gamma.s.zero.study <- gamma.s.study[-c(1:nrow(sone.study))]
      
      gamma.s[[ix]] = list("gamma.s.one" = gamma.s.one.study, "gamma.s.zero" = gamma.s.zero.study)
      
      # Check if this study is paired
      if (is.null(paired.studies)) {
        paired.study = FALSE
      } else {
        paired.study = ifelse(study %in% paired.studies, TRUE, FALSE)
      }
      
      epsilon_arg <- epsilon.meta
      
      # Test gamma.s as a surrogate
      gamma.s.result.study <- test.surrogate.extension(
        yone = yone.study,
        yzero = yzero.study,
        sone = gamma.s.one.study,
        szero = gamma.s.zero.study,
        alpha,
        epsilon = epsilon_arg,
        power.want.s = NULL,
        u.y.hyp,
        alternative,
        paired = paired.study
      )
      
      # Extract relevant screening results
      rise.evaluate.results.allstudies[[ix]] <- data.frame(
        "study" = study,
        "epsilon" = gamma.s.result.study$epsilon.used,
        "marker" = "gamma",
        "n" = length(yone.study) + length(yzero.study),
        "u.y" = gamma.s.result.study$u.y,
        "u.s" = gamma.s.result.study$u.s,
        "delta" = gamma.s.result.study$delta.estimate,
        "ci.lower" = gamma.s.result.study$ci.delta[1],
        "ci.upper" = gamma.s.result.study$ci.delta[2],
        "sd" = gamma.s.result.study$sd.delta,
        "p.unadjusted" = gamma.s.result.study$p.delta
      )
      
      # Increase index
      ix <- ix + 1L
      
    }
    
    names(gamma.s) = all.study.names
    
    # Bind the per-study results into a dataframe
    evaluation.metrics.study <- bind_rows(rise.evaluate.results.allstudies) %>%
      mutate(p.adjusted = p.adjust(p.unadjusted, method = p.correction))
    
    # If some studies have exactly 0 standard error, report this
    if (any(evaluation.metrics.study$sd == 0)) {
      sd0.studies = evaluation.metrics.study %>%
        filter(sd == 0) %>%
        pull(study) %>%
        unique()
      
      message(
        paste0(
          "Note: studies '",
          paste(sd0.studies, collapse = ", "),
          "' have degenerate (exactly 0) estimation of standard error",
          " for the combined marker gamma. ",
          "This is likely due to small sample size and perfect pairwise",
          " concordance between gamma and the primary endpoint.",
          " These studies will be removed from the meta-analysis."
        )
      )
      
      evaluation.metrics.study = evaluation.metrics.study %>%
        filter(sd != 0)
    }
  }
  
  # Now do meta-analysis on these results
  # Extract the values of delta
  delta.marker = evaluation.metrics.study %>%
    pull(delta)
  
  # Extract the values of the standard deviation of delta
  sd.delta.marker = evaluation.metrics.study %>%
    pull(sd)
  
  if (nrow(evaluation.metrics.study) == 2) {
    message(
      paste0(
        "Note: gamma has only 2 studies available for estimation.",
        " Prediction intervals will not be available."
      )
    )
  }
  
  if (nrow(evaluation.metrics.study) == 1) {
    message(
      paste0(
        "Note: gamma has only 1 study available for estimation.",
        " Meta-analysis is not possible."
      )
    )
    
    return(
      list(
        "individual.metrics" = NULL,
        "evaluation.metrics.study" = evaluation.metrics.study,
        "evaluation.metrics.meta" = NULL,
        "gamma.s" = gamma.s,
        "gamma.s.plot" = list()
      )
    )
  }
  
  # Call the restricted maximum likelihood random-effects meta-analysis function
  delta.reml.marker = delta.reml.meta(
    delta = delta.marker,
    sd.delta = sd.delta.marker,
    epsilon = epsilon.meta,
    alpha = alpha,
    alternative = alternative,
    test = test,
    meta.analysis.method = meta.analysis.method
  )[["results"]]
  
  # Save the study weights in one of the results dataframes
  evaluation.metrics.study = evaluation.metrics.study %>%
    mutate(study.weight = delta.reml.marker[["weights.tau"]],
           study.weight.relative = delta.reml.marker[["weights.tau.relative"]]) %>%
    relocate(study.weight, .after = study) %>%
    relocate(study.weight.relative, .after = study.weight)
  
  # Bind per-marker results into a data frame
  evaluation.metrics.meta = bind_rows(delta.reml.marker) %>%
    mutate("marker" = "gamma") %>%
    relocate(marker, .before = n.studies)
  
  # Plots
  
  # If fit plot desired, return it
  # This plot shows the effect estimates for each study, with size proportional to sample size
  if (return.fit.plot) {
    # Set a minimum value for the axes
    plot.min.global <- evaluation.metrics.study %>%
      dplyr::select(u.y, u.s) %>%
      min() - 0.1
    
    # Compute CCC
    x <- evaluation.metrics.study$u.y
    y <- evaluation.metrics.study$u.s
    ccc <- (2 * cov(x, y)) / (var(x) + var(y) + (mean(x) - mean(y))^2)
    
    # Prepare legend sizing
    n_vals <- evaluation.metrics.study$n
    round_down_10 <- function(x)
      floor(x / 10) * 10
    round_up_10   <- function(x)
      ceiling(x / 10) * 10
    round_up_50   <- function(x)
      ceiling(x / 50) * 50
    
    if (length(unique(n_vals)) == 1) {
      # All sample sizes equal → single legend key
      legend_breaks <- unique(n_vals)
      legend_labels <- as.character(unique(n_vals))
    } else {
      min_label <- round_down_10(min(n_vals))
      max_label <- round_up_10(max(n_vals))
      mid_label <- round_up_50(median(n_vals))
      
      legend_breaks <- c(min(n_vals), mid_label, max(n_vals))
      legend_labels <- c(as.character(min_label),
                         as.character(mid_label),
                         as.character(max_label))
    }
    
    # Plot with CCC annotation and improved sizing
    # Plot with smallest n always visible (size_min = 5)
    fit.plot <- evaluation.metrics.study %>%
      ggplot(aes(y = u.y, x = u.s)) +
      geom_point(
        aes(size = n),
        shape = 21,
        alpha = 0.5,
        stroke = 1,
        fill = "#6FB1EF"
      ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        color = "#FF2128",
        linetype = "dashed",
        linewidth = 0.8,
        alpha = 0.5
      ) +
      # Annotate CCC in top-left
      annotate(
        "text",
        x = plot.min.global,
        y = 0.95,
        label = paste0("CCC = ", round(ccc, 2)),
        hjust = 0,
        vjust = 1,
        color = "red",
        size = 12
      ) +
      # Scale sizes relative to the smallest n
      scale_size_continuous(
        range = c(5, 20),
        # smallest n = size 5, largest = 14
        breaks = legend_breaks,
        labels = legend_labels
      ) +
      scale_x_continuous(limits = c(plot.min.global, 1.01),
                         expand = c(0, 0)) +
      scale_y_continuous(limits = c(plot.min.global, 1.01),
                         expand = c(0, 0)) +
      coord_fixed(ratio = 1) +
      labs(
        title = "Treatment effects on primary response and combined marker on evaluation data",
        x = "Treatment effect on combined marker",
        y = "Treatment effect on primary endpoint",
        size = "Study N",
        fill = NULL
      ) +
      theme_minimal(base_size = 25) +
      theme(
        plot.title = element_text(
          size = 25,
          hjust = 0.5,
          face = "bold"
        ),
        axis.title = element_text(size = 25),
        legend.position = "right"
      )
    
  } else {
    fit.plot = NULL
  }
  
  # If forest plot desired, return it
  if (return.forest.plot) {
    # Extract results for all studies with positive standard error
    evaluation.metrics.study.temp = evaluation.metrics.study
    
    df.gamma.temp = evaluation.metrics.meta %>%
      mutate(
        study = paste0("Pooled effect (", 100 * (1 - 2 * alpha), "% C.I.)"),
        epsilon = epsilon.meta,
        delta = mu.delta,
        ci.lower = ci.delta.lower,
        ci.upper = ci.delta.upper,
        pi.lower = pi.lower,
        pi.upper = pi.upper,
        u.y = NA,
        u.s = NA,
        study.weight = NA,
        study.weight.relative = NA,
        marker = "gamma",
        n = length(yone) + length(yzero),
        sd = se.delta,
        p.unadjusted = p,
        p.adjusted = NA
      ) %>%
      select(all_of(colnames(evaluation.metrics.study.temp)), pi.lower, pi.upper) %>%
      distinct()
    
    if (show.pooled.effect) {
      evaluation.metrics.study2 = bind_rows(evaluation.metrics.study.temp, df.gamma.temp)
    } else {
      evaluation.metrics.study2 = evaluation.metrics.study.temp
    }
    
    I2 = evaluation.metrics.meta$I2
    tau2 = evaluation.metrics.meta$tau2
    
    # Compute CCC
    x <- evaluation.metrics.study.temp$u.y
    y <- evaluation.metrics.study.temp$u.s
    ccc <- (2 * cov(x, y)) / (var(x) + var(y) + (mean(x) - mean(y))^2)
    
    # Separate studies vs summary
    studies.df <- evaluation.metrics.study2 %>% filter(study != paste0("Pooled effect (", 100 *
                                                                         (1 - 2 * alpha), "% C.I.)"))
    summary.row <- evaluation.metrics.study2 %>% filter(study == paste0("Pooled effect (", 100 *
                                                                          (1 - 2 * alpha), "% C.I.)"))
    
    # sort studies by effect size from most negative to most positive
    studies.df <- studies.df %>% arrange(delta)
    
    # vertical positions: top (k) down to 1; summary at y = 0
    k <- nrow(studies.df)
    studies.df <- studies.df %>% mutate(y = seq(from = k, to = 1))
    summary.row <- summary.row %>% mutate(y = 0)
    
    # Prepare combined plotting df and standard display labels
    plot.df <- bind_rows(studies.df, summary.row) %>%
      mutate(
        is.summary = (study == paste0(
          "Pooled effect (", 100 * (1 - 2 * alpha), "% C.I.)"
        )),
        study.label = study,
        label.pval = ifelse(
          is.na(p.unadjusted),
          "",
          formatC(p.unadjusted, format = "f", digits = 3)
        ),
        label.n = ifelse(is.na(n), "", as.character(n)),
        ci_width = ci.upper - ci.lower,
        summary_size = pmin(ci_width / 2, 6)
      )
    
    # clip confidence intervals to [-1, 1] before plotting
    plot.df <- plot.df %>%
      mutate(
        summary.ci.lower = ifelse(is.summary, ci.lower, NA_real_),
        summary.ci.upper = ifelse(is.summary, ci.upper, NA_real_),
        ci.lower = ifelse(is.summary, NA_real_, ci.lower),
        ci.upper = ifelse(is.summary, NA_real_, ci.upper),
        summary.ci.lower = ifelse(is.na(summary.ci.lower), NA_real_, pmax(pmin(
          summary.ci.lower, 1
        ), -1)),
        summary.ci.upper = ifelse(is.na(summary.ci.upper), NA_real_, pmax(pmin(
          summary.ci.upper, 1
        ), -1)),
        ci.lower = ifelse(is.na(ci.lower), NA_real_, pmax(pmin(ci.lower, 1), -1)),
        ci.upper = ifelse(is.na(ci.upper), NA_real_, pmax(pmin(ci.upper, 1), -1))
      )
    
    # add prediction-interval row (keeps its pi.lower/pi.upper values from summary.row),
    # placed below the pooled summary (y = -1), and prevent generic CI drawing
    if (show.pooled.effect) {
      pi.row <- summary.row %>%
        mutate(
          study = paste0(100 * (1 - alpha), "% Prediction interval"),
          study.label = paste0(100 * (1 - alpha), "% Prediction interval"),
          ci.lower = NA_real_,
          ci.upper = NA_real_,
          # keep pi.lower / pi.upper as they come from summary.row (do not overwrite)
          p.unadjusted = NA_real_,
          n = NA_integer_,
          study.weight = NA_real_,
          study.weight.relative = NA_real_,
          y = -1,
          is.summary = FALSE,
          label.pval = "",
          label.n = "",
          pi.lower = ifelse(is.na(pi.lower), NA_real_, pmax(pmin(pi.lower, 1), -1)),
          pi.upper = ifelse(is.na(pi.upper), NA_real_, pmax(pmin(pi.upper, 1), -1))
        )
      plot.df <- bind_rows(plot.df, pi.row)
    }
    
    # Labels for weighting
    weights.vec <- studies.df$study.weight.relative
    pct.vec <- weights.vec
    label.wgt.vec <- formatC(pct.vec, format = "f", digits = 1)
    
    # attach weight labels to plot.df (empty for summary / PI row)
    plot.df <- plot.df %>%
      left_join(tibble(study = studies.df$study, label.wgt = label.wgt.vec),
                by = "study") %>%
      mutate(label.wgt = ifelse(is.na(label.wgt), "", label.wgt))
    
    # prepare diamond polygon data for pooled summary (single diamond)
    diamond.df <- data.frame(x = numeric(0), y = numeric(0))
    if (any(plot.df$is.summary)) {
      s <- plot.df %>% filter(is.summary) %>% slice(1)
      if (!is.na(s$summary.ci.lower) &&
          !is.na(s$summary.ci.upper)) {
        h <- 0.10   # slim diamond height
        diamond.df <- data.frame(
          x = c(
            s$delta,
            s$summary.ci.upper,
            s$delta,
            s$summary.ci.lower
          ),
          y = c(s$y + h, s$y, s$y - h, s$y)
        )
      }
    }
    
    # Heterogeneity text from provided object evaluation.metrics.meta
    tau2.txt <- sub(
      "e-0?",
      "e-",
      # remove leading zero in exponent
      format(signif(evaluation.metrics.meta$tau2, 1), scientific = TRUE)
    )
    I2.txt   <- formatC(evaluation.metrics.meta$I2,
                        digits = 1,
                        format = "f")
    ccc.txt = formatC(ccc, digits = 2, format = "f")
    
    # Plot parameters
    base.text.size <- 14
    y.min <- if (show.pooled.effect)
      - 1.5
    else
      0
    y.max <- k + 1
    rel.w.left  <- 0.45
    rel.w.mid   <- 1.10
    rel.w.right <- 0.55
    x.min <- -1
    x.max <- 1
    
    # Left panel: study labels, bold the summary
    left.labels <- ggplot(plot.df, aes(y = y)) +
      geom_text(
        data = filter(plot.df, !is.summary),
        aes(x = 0, label = study.label),
        hjust = 0,
        size = 5
      ) +
      geom_text(
        data = filter(plot.df, is.summary),
        aes(x = 0, label = study.label),
        hjust = 0,
        size = 5,
        fontface = "bold"
      ) +
      scale_x_continuous(limits = c(-0.1, 0.95), expand = c(0, 0)) +
      scale_y_continuous(
        breaks = plot.df$y,
        limits = c(y.min, y.max),
        expand = c(0, 0)
      ) +
      coord_cartesian(clip = "off") +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = unit(c(1.2, 0.6, 1, 1.2), "lines")
      )
    
    epsilon.meta.rounded <- round(epsilon.meta, 3)
    
    shade_df <- data.frame(
      xmin = -epsilon.meta,
      xmax = epsilon.meta,
      ymin = y.min,
      ymax = y.max
    )
    
    forest.mid <- ggplot(plot.df, aes(x = delta, y = y)) +
      geom_rect(
        data = shade_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = "#B4B4B4",
        alpha = 0.3,
        color = NA
      ) +
      geom_errorbar(
        aes(xmin = ci.lower, xmax = ci.upper),
        width = 0.15,
        linewidth = 0.8,
        orientation = "y"
      ) +
      geom_point(
        data = filter(
          plot.df,
          !is.summary &
            study != paste0(100 * (1 - alpha), "% Prediction interval")
        ),
        shape = 16,
        size = 3.5
      ) +
      geom_polygon(
        data = diamond.df,
        aes(x = x, y = y),
        inherit.aes = FALSE,
        fill = "#CCCCCC",
        color = "black"
      ) +
      geom_segment(
        data = filter(plot.df, is.summary),
        aes(
          x = delta,
          xend = delta,
          y = y,
          yend = y.max
        ),
        inherit.aes = FALSE,
        linetype = "dotted",
        color = "#4C78A8",
        linewidth = 0.8
      ) +
      geom_errorbar(
        data = filter(plot.df, study == paste0(100 * (1 - alpha), "% Prediction interval")),
        aes(xmin = pi.lower, xmax = pi.upper),
        width = 0.15,
        linewidth = 1.5,
        color = "#EF5C52",
        orientation = "y"
      ) +
      scale_x_continuous(
        limits = c(x.min, x.max),
        breaks = seq(x.min, x.max, by = 0.5),
        expand = c(0, 0)
      ) +
      scale_y_continuous(
        breaks = plot.df$y,
        limits = c(y.min, y.max),
        expand = c(0, 0)
      ) +
      coord_cartesian(
        ylim = c(y.min, y.max),
        xlim = c(x.min, x.max),
        clip = "off"
      ) +
      labs(
        x = expression(delta),
        y = NULL,
        title = if (show.pooled.effect)
          "Random-effects meta-analysis of combined marker in evaluation data"
        else
          "Effect sizes of combined marker across studies in evaluation data"
      ) +
      theme_minimal(base_size = base.text.size) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = base.text.size + 4
        ),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = base.text.size + 4),
        axis.text.x = element_text(size = base.text.size),
        plot.margin = unit(c(1.2, 0.6, 1, 0.6), "lines")
      ) +
      geom_vline(
        xintercept = c(-epsilon.meta, epsilon.meta),
        linetype = "solid",
        color = "#2E2E2E",
        linewidth = 1
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "solid",
        color = "grey60"
      )
    
    # Right panel: p.unadjusted-value / weight / N
    right.df <- plot.df %>%
      mutate(
        display.pval = label.pval,
        display.wgt = ifelse(label.wgt == "", "", paste0(label.wgt, "%")),
        display.n   = label.n
      )
    
    if (show.pooled.effect) {
      right.table <- ggplot(right.df, aes(y = y)) +
        annotate(
          "text",
          x = 1,
          y = k + 0.9,
          label = "p-value",
          fontface = "bold",
          hjust = 0,
          size = 5
        ) +
        annotate(
          "text",
          x = 2.1,
          y = k + 0.9,
          label = "Weight",
          fontface = "bold",
          hjust = 0,
          size = 5
        ) +
        annotate(
          "text",
          x = 3.3,
          y = k + 0.9,
          label = "N",
          fontface = "bold",
          hjust = 0,
          size = 5
        ) +
        geom_text(aes(x = 1, label = display.pval),
                  hjust = 0,
                  size = 4.5) +
        geom_text(aes(x = 2.1, label = display.wgt),
                  hjust = 0,
                  size = 4.5) +
        geom_text(aes(x = 3.3, label = display.n),
                  hjust = 0,
                  size = 4.5) +
        scale_y_continuous(
          breaks = right.df$y,
          limits = c(y.min, y.max),
          expand = c(0, 0)
        ) +
        scale_x_continuous(limits = c(0.9, 4.0), expand = c(0, 0)) +
        coord_cartesian(clip = "off") +
        theme_void() +
        theme(plot.margin = unit(c(1.2, 1.2, 1, 0.6), "lines"))
    } else {
      right.table <- ggplot(right.df, aes(y = y)) +
        annotate(
          "text",
          x = 1,
          y = k + 0.9,
          label = "p-value",
          fontface = "bold",
          hjust = 0,
          size = 5
        ) +
        annotate(
          "text",
          x = 2.2,
          y = k + 0.9,
          label = "N",
          fontface = "bold",
          hjust = 0,
          size = 5
        ) +
        geom_text(aes(x = 1, label = display.pval),
                  hjust = 0,
                  size = 4.5) +
        geom_text(aes(x = 2.2, label = display.n),
                  hjust = 0,
                  size = 4.5) +
        scale_y_continuous(
          breaks = right.df$y,
          limits = c(y.min, y.max),
          expand = c(0, 0)
        ) +
        scale_x_continuous(limits = c(0.9, 2.8), expand = c(0, 0)) +
        coord_cartesian(clip = "off") +
        theme_void() +
        theme(plot.margin = unit(c(1.2, 1.2, 1, 0.6), "lines"))
    }
    
    # Combine panels and add bottom info row
    combined <- plot_grid(
      left.labels,
      forest.mid,
      right.table,
      nrow = 1,
      rel_widths = c(rel.w.left, rel.w.mid, rel.w.right),
      align = "h"
    )
    
    info.text <- bquote(
      "Tau-squared =" ~ .(tau2.txt) ~
        "|" ~ "I-Squared =" ~ .(I2.txt) ~ "%" ~
        "|" ~ "CCC =" ~ .(ccc.txt) ~
        "|" ~ epsilon == .(epsilon.meta.rounded)
    )
    
    info.grob <- ggdraw() + draw_label(
      info.text,
      x = 0.5,
      y = 0.5,
      hjust = 0.5,
      size = base.text.size + 2
    )
    
    if (show.pooled.effect) {
      forest.plot <- plot_grid(combined,
                               info.grob,
                               ncol = 1,
                               rel_heights = c(0.95, 0.05))
    } else {
      forest.plot <- combined
    }
    
  } else {
    forest.plot = NULL
  }
  
  if (return.all.evaluate) {
    rise.screen.results = rise.screen.meta(
      yone,
      yzero,
      sone,
      szero,
      studyone,
      studyzero,
      alpha,
      epsilon.study,
      epsilon.meta = epsilon.meta,
      power.want.s.study,
      u.y.hyp,
      p.correction,
      n.cores,
      alternative,
      paired.all,
      paired.studies,
      return.all.screen = T,
      return.all.weights = T,
      weight.mode,
      normalise.weights = T,
      return.forest.plot = F,
      return.fit.plot = F
    )
    
    individual.metrics = list(
      "individual.metrics.study" = rise.screen.results[["screening.metrics.study"]],
      "individual.metrics.meta" = rise.screen.results[["screening.metrics.meta"]]
    )
  } else {
    individual.metrics = NULL
  }
  
  
  return(
    list(
      "individual.metrics" = individual.metrics,
      "evaluation.metrics.study" = evaluation.metrics.study,
      "evaluation.metrics.meta" = evaluation.metrics.meta,
      "gamma.s" = gamma.s,
      "gamma.s.plot" = list("fit.plot" = fit.plot, "forest.plot" = forest.plot)
    )
  )
  
}
