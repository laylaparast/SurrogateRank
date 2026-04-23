#' Function to perform the screening stage of RISE-meta : Meta-Analysis of High-Dimensional Surrogate Markers
#'
#' @description The RISE screening algorithm is applied to each study using a rank-based measure of
#' treatment effect similarity. In the second stage, these effect estimates are combined using a
#' random-effects meta-analysis and the retained markers are those for which there is strong evidence
#' of surrogacy across many studies.
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
#' @param alpha significance level for determining surrogate candidates in both stages. Default is \code{0.05}.
#' @param power.want.s.study numeric in (0,1) - power desired for a test of treatment effect based on the
#'   surrogate candidate. Either this or \code{epsilon.study} argument must be specified.
#' @param epsilon.study numeric in (0,1) - non-inferiority margin for determining surrogate validity in the
#'   within-study screening phase. Either this or \code{power.want.s.study} argument must be specified.
#' @param epsilon.meta.mode character string specifying the mode to choose the value of the acceptable margin defined
#' by epsilon. By default, this is set to "user", where the value of epsilon is fixed by the user, defined by the
#' value of the argument \code{epsilon.meta}. The alternative is to set this as "mean.power", which corresponds to
#' taking the mean value of epsilon across studies such that the power to detect departures from the null within
#' each study is defined by the \code{power.want.s.study} argument.
#' @param epsilon.meta numeric in (0,1) - fixed non-inferiority margin for determining surrogate validity
#'   in the meta-analysis stage.
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
#' @param return.all.screen logical flag. If \code{TRUE} (default), a dataframe will be returned giving
#'   the screening results for all candidates. Else, only the significant candidates will be returned.
#' @param return.all.weights logical flag. If \code{FALSE} (default), a dataframe will be returned giving
#'   weights for significant markers screened. If \code{TRUE}, weights for all markers will be returned. Note
#'   that, if normalised weights are required, these will only be returned for significant markers, and raw
#'   weights will be returned in a second column.
#' @param weight.mode character giving the type of weighting to return. One of
#'   \code{c("diff.epsilon", "inverse.delta", or "none")}. The default is \code{"diff.epsilon"}, which calculates
#'   weights as the proportion of the interval between 0 and epsilon.study cut by the absolute value of delta,
#'   therefore giving delta = 0 a weight of 1 and delta = epsilon.study a weight of 0. Another option is
#'   \code{"inverse.delta"} where the weights are determined by taking the inverse of the absolute values of delta.
#'   When \code{"none"}, the weights are set to 1 for every marker.
#' @param normalise.weights logical flag. If \code{TRUE} (default), the weights are normalised by the sum of
#'   all the weights such that the maximum weight is 1, which can help with interpretability.
#' @param return.screen.plot logical flag. If \code{TRUE} (default), returns a forest plot of the top predictors, sorted by p-value,
#'   from the screening stage. The number of predictors to display is given by the \code{screen.plot.topN} argument, which has default
#'   value 15.
#' @param screen.plot.topN number of predictors to display in the screening results figure, default value is 15.
#' @param screen.plot.point.estimate logical flag. If \code{FALSE} (default), uses the \code{screen.plot.topN} argument to determine how many
#' markers to display on the screen plot. Otherwise, plots all the markers with a point estimate within the equivalence region.
#' @param return.forest.plot logical flag. If \code{TRUE} (default), a forest plot of the effect sizes for the
#' combined signature across studies, with its meta-analysis summary measure and prediction interval, will be included in the output.
#' @param return.fit.plot logical flag. If \code{TRUE} (default), a plot of the effects on the primary response
#' versus the effects on the combined surrogate signature for each study will be included in the output.
#' @param show.pooled.effect logical flag. If \code{TRUE} (default), the forest plot will show the pooled effect
#' estimate. Otherwise, it will just show the individual trial estimates.
#' @param return.study.similarity.plot logical flag. If \code{TRUE} (default), will return two plots showing the similarity
#' between study-wise marker signatures (i.e., the application of RISE to each study individually, with p-value correction within-study).
#' @param return.evaluate.results logical flag. If \code{TRUE} (default), returns results for combined marker gamma, evaluated on the same
#' data. Can be useful to set this as \code{FALSE} to save computational time if this is not of interest.
#' @param meta.analysis.method character giving the meta-analysis method to be used. The default is \code{RE}, corresponding to
#' random-effects meta-analysis, whereas setting this argument to \code{FE} uses fixed-effects meta-analysis.
#'
#' @return a list with elements \itemize{
#'   \item \code{screening.metrics.study} : dataframe of per-study results from RISE screening.
#'   For each candidate marker - study name, study sample size, estimate of delta, standard error of delta.
#'   \item \code{screening.metrics.meta} : dataframe of meta-analysis screening results.
#'   For each candidate marker - number of studies \code{n.studies},
#'   estimate of mean delta value \code{mu.delta},
#'   its standard error \code{se.delta}, confidence interval and prediction interval,
#'   estimate of tau-squared \code{tau2}, Cochran's Q-statistic and Higgins-Thompson I-Squared,
#'   unadjusted and adjusted meta-analysis p-values, and standardised weights.
#'   Note : if using the non-inferiority test (i.e. \code{alternative = "less"}),
#'          the intervals have width (1-\code{alpha})*100%,
#'          whereas the two-one-sided test (i.e. \code{alternative = "two.sided"})
#'          corresponds to a (1-2\code{alpha})*100% width.
#'   \item \code{significant.markers}: character vector of markers with meta-analysis p-values \code{< alpha}
#'   \item \code{screening.weights}: dataframe giving marker names and the standardised meta-analysis weights
#'   \item \code{evaluation.metrics.study} : dataframe of per-study results for the combined marker gamma, evaluated on the same data
#'   \item \code{evaluation.metrics.meta} : dataframe of meta-analysis results for the combined marker gamma, evaluated on the same data
#'   \item \code{gamma.s.plot}: if \code{return.forest.plot}, \code{return.fit.plot}, and/or  \code{return.study.similarity.plot}
#'   are \code{TRUE}, returns fitted evaluation plots on training data as a list.
#' }
#'
#' @import dplyr pbmcapply ggplot2 cowplot ComplexUpset tidyr
#' @export
#' @author Arthur Hughes
#'
#' @examples
#' data("example.data.highdim.multistudy.ipd")
#' yone <- example.data.highdim.multistudy.ipd$y1
#' yzero <- example.data.highdim.multistudy.ipd$y0
#' sone <- example.data.highdim.multistudy.ipd$s1
#' szero <- example.data.highdim.multistudy.ipd$s0
#' studyone <- example.data.highdim.multistudy.ipd$study1
#' studyzero <- example.data.highdim.multistudy.ipd$study0
#' rise.meta.screen.result <- rise.screen.meta(
#' yone, yzero, 
#' sone, szero, 
#' studyone, studyzero, 
#' epsilon.study = 0.2, epsilon.meta = 0.2
#' )
rise.screen.meta = function(yone,
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
                            return.all.screen = TRUE,
                            return.all.weights = FALSE,
                            weight.mode = "diff.epsilon",
                            return.screen.plot = TRUE,
                            screen.plot.topN = 15,
                            screen.plot.point.estimate = FALSE,
                            normalise.weights = TRUE,
                            return.forest.plot = TRUE,
                            return.fit.plot = TRUE,
                            show.pooled.effect = TRUE,
                            return.study.similarity.plot = TRUE,
                            return.evaluate.results = TRUE,
                            meta.analysis.method = "RE") {
  # DATA FORMATTING #
  ## Convert dataframes / vectors to numeric matrices
  to_numeric_matrix <- function(x) {
    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }
    if (is.vector(x) && !is.matrix(x)) {
      x <- matrix(x, ncol = 1)
    }
    x
  }
  
  sone  <- to_numeric_matrix(sone)
  szero <- to_numeric_matrix(szero)
  
  # If no column names on surrogate candidates, set them as the column indices
  if (is.null(colnames(sone))) {
    colnames(sone) <- paste0("marker", seq_len(ncol(sone)))
  }
  if (is.null(colnames(szero))) {
    colnames(szero) <- paste0("marker", seq_len(ncol(szero)))
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
  
  ## Retrieve the study names
  all.study.names = unique(studyone)
  
  ## if in paired mode, specified studies must have the same number of treated and untreated samples
  
  if (paired.all) {
    paired.studies = unique(studyone) # if in all-paired mode, set this argument to all study names
  }
  
  # Check all specified paired studies have the same number of treated and untreated samples
  for (study in paired.studies) {
    if (length(yone[studyone == study]) != length(yzero[studyzero == study])) {
      stop(
        paste0(
          "Paired mode is requested but the number of samples in yone does not match that of yzero in study '",
          study,
          "'"
        )
      )
    }
  }
  
  
  # APPLICATION OF RISE WITHIN EACH STUDY #
  rise.screen.results.allstudies <- vector("list", length(all.study.names)) # prepare storage for study results
  ix <- 1L # initialise index
  
  for (study in all.study.names) {
    # For each study
    message(paste0("Screening study '", study, "'")) # Print progress message
    
    yzero.study = yzero[studyzero == study] # Extract study samples
    yone.study = yone[studyone == study]
    
    szero.study = szero[studyzero == study, , drop = FALSE]
    sone.study = sone[studyone == study, , drop = FALSE]
    
    # Check if this study is paired
    if (is.null(paired.studies)) {
      paired.study = FALSE
    } else {
      paired.study = ifelse(study %in% paired.studies, TRUE, FALSE)
    }
    
    epsilon_arg <- if (is.null(power.want.s.study)) {
      epsilon.study
    } else {
      NULL
    }
    
    # Apply RISE screen function
    screen.results.study = rise.screen(
      yone = yone.study,
      yzero = yzero.study,
      sone = sone.study,
      szero = szero.study,
      epsilon = epsilon_arg,
      alpha = alpha,
      power.want.s = power.want.s.study,
      u.y.hyp = u.y.hyp,
      p.correction = p.correction,
      n.cores = n.cores,
      alternative = alternative,
      paired = paired.study,
      return.all.screen = TRUE,
      return.all.weights = F,
      weight.mode = weight.mode,
      normalise.weights = F,
      verbose = F
    )
    
    # Extract relevant screening results
    study_metrics <- screen.results.study[["screening.metrics"]]
    
    if (!"epsilon" %in% names(study_metrics)) {
      study_metrics$epsilon <- if (!is.null(epsilon_arg)) {
        epsilon_arg
      } else {
        NA_real_
      }
    }
    
    rise.screen.results.allstudies[[ix]] <- study_metrics %>%
      mutate(study = study) %>%
      dplyr::select(study,
                    epsilon,
                    marker,
                    n,
                    u.y,
                    u.s,
                    delta,
                    sd,
                    p_unadjusted,
                    p_adjusted)
    
    # Increase index
    ix <- ix + 1L
    
  }
  
  # Bind the per-study results into a dataframe
  rise.screen.results.allstudies.df <- bind_rows(rise.screen.results.allstudies)
  
  # Calculate average epsilon for meta-analysis
  if (epsilon.meta.mode == "mean.power") {
    epsilon.df = rise.screen.results.allstudies.df %>%
      dplyr::select(study, epsilon) %>%
      distinct()
    
    epsilon.meta = mean(epsilon.df$epsilon)
    
    message(
      "Using value ",
      round(epsilon.meta, 3),
      " for epsilon.meta, based on the mean of
            epsilon values across studies requiring specified power."
    )
  }
  
  if (return.study.similarity.plot) {
    sig_list <- rise.screen.results.allstudies.df %>%
      filter(p_adjusted < alpha) %>%
      group_by(study, n) %>%
      summarise(markers = list(unique(marker)), .groups = "drop") %>%
      mutate(study_label = paste0(study, " (N = ", n, ")")) %>%
      {
        setNames(.$markers, .$study_label)
      }
    
    # Plot Venn diagram
    venn.plot <- suppressWarnings(
      ggVennDiagram::ggVennDiagram(sig_list) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        theme(
          legend.position = "none",
          plot.margin = margin(10, 10, 10, 10)
        ) +
        coord_cartesian(clip = "off") +
        theme(text = element_text(size = 10))
    )
    
    # Prepare data frame for ComplexUpset
    upset_df <- sig_list %>%
      # Convert named list to long format
      tibble::enframe(name = "study", value = "marker") %>%
      unnest(marker) %>%
      mutate(present = TRUE) %>%
      pivot_wider(
        names_from = study,
        values_from = present,
        values_fill = FALSE
      )
    
    # Create the UpSet plot
    upset.plot = upset(
      upset_df,
      colnames(upset_df)[-1],
      # exclude the 'marker' column
      name = "Markers",
      base_annotations = list('Intersection size' = intersection_size(counts =
                                                                        TRUE))
    )
    
    similarity.plots = list("venn.plot" = venn.plot, "upset.plot" = upset.plot)
    
  }
  
  # Get names of all markers screened
  all.markers = rise.screen.results.allstudies.df %>%
    pull(marker) %>%
    unique()
  
  # Initialise list for per-marker random-effects restricted maximum likelihood meta-analysis
  delta.reml.marker = vector(mode = "list", length = length(all.markers))
  names(delta.reml.marker) = all.markers
  
  
  # If some markers in some studies have exactly 0 standard error, report this
  if (any(rise.screen.results.allstudies.df$sd == 0)) {
    sd0.studies = rise.screen.results.allstudies.df %>%
      filter(sd == 0) %>%
      pull(study) %>%
      unique()
    
    message(
      paste0(
        "Note: studies '",
        paste(sd0.studies, collapse = ", "),
        "' have degenerate (exactly 0) estimation of standard error",
        " for some markers. ",
        "This is likely due to small sample size and perfect pairwise",
        " concordance between gamma and the primary endpoint.",
        " These studies will be removed from the meta-analysis for those markers only."
      )
    )
    
    rise.screen.results.allstudies.df = rise.screen.results.allstudies.df %>%
      filter(sd != 0)
  }
  
  df.nstudies.marker = rise.screen.results.allstudies.df %>%
    group_by(marker) %>%
    summarise(n_studies = n(), .groups = "drop")
  
  if (any(df.nstudies.marker$n_studies == 1)) {
    ngenes = df.nstudies.marker %>%
      filter(n_studies == 1) %>%
      nrow()
    
    message(
      paste0(
        "Note: ",
        ngenes,
        " markers have only 1 studies available for estimation.",
        " These markers will be skipped in the meta-analysis."
      )
    )
  }
  
  # Pre-split data by marker once to avoid repeated filtering in the loop
  marker.data.list <- split(rise.screen.results.allstudies.df,
                            rise.screen.results.allstudies.df$marker)
  markers.to.run  <- intersect(all.markers, names(marker.data.list))
  markers.to.run  <- markers.to.run[sapply(marker.data.list[markers.to.run], nrow) > 1]
  
  run.one.marker <- function(df.summary.marker) {
    x   <- df.summary.marker$u.y
    y   <- df.summary.marker$u.s
    ccc <- (2 * cov(x, y)) / (var(x) + var(y) + (mean(x) - mean(y))^2)
    
    delta.marker        <- df.summary.marker$delta
    sd.delta.marker     <- df.summary.marker$sd
    
    res <- delta.reml.meta(
      delta        = delta.marker,
      sd.delta     = sd.delta.marker,
      epsilon      = epsilon.meta,
      alpha        = alpha,
      alternative  = alternative,
      test = test,
      meta.analysis.method = meta.analysis.method
    )[["results"]]
    
    res$weights.tau          <- NULL
    res$weights.tau.relative <- NULL
    res$ccc                  <- ccc
    res
  }
  
  delta.reml.marker <- parallel::mclapply(marker.data.list[markers.to.run], run.one.marker, mc.cores = n.cores)
  
  # Bind per-marker results into a data frame
  delta.reml.df = bind_rows(delta.reml.marker, .id = "marker")
  
  # Calculate surrogate strength weights according to user specification
  if (weight.mode == "diff.epsilon") {
    delta.reml.df = delta.reml.df %>%
      mutate(weight.strength = (epsilon.meta - abs(mu.delta)) / epsilon.meta)
  } else if (weight.mode == "inverse.delta") {
    delta.reml.df = delta.reml.df %>%
      mutate(weight.strength = ifelse (mu.delta != 0, 1 / abs(mu.delta), 1 / abs(mu.delta + 0.00001)))
  } else if (weight.mode == "none") {
    delta.reml.df = delta.reml.df %>%
      mutate(weight.strength = 1)
  }
  
  
  # Now calculate combined weights using heterogeneity and surrogate strength
  delta.reml.df = delta.reml.df %>%
    mutate(
      p.adjusted = p.adjust(p, method = p.correction),
      weight.heterogeneity = 1 / weights.tau.sum,
      weight = weight.strength * weight.heterogeneity,
      p.unadjusted = p,
      pi.delta.lower = pi.lower,
      pi.delta.upper = pi.upper
    ) %>%
    dplyr::select(
      marker,
      n.studies,
      mu.delta,
      se.delta,
      ci.delta.lower,
      ci.delta.upper,
      pi.delta.lower,
      pi.delta.upper,
      tau2,
      Q,
      I2,
      ccc,
      p.unadjusted,
      p.adjusted,
      weight
    )
  
  # Extract significant markers from meta-analysis
  significant.markers = delta.reml.df %>%
    filter(p.adjusted < alpha) %>%
    pull(marker)
  
  # If screening plot desired, return it
  if (return.screen.plot) {
    p_floor <- 1e-2   # practical lower bound for the colour scale
    
    if(screen.plot.point.estimate){
      screen.plot.topN = delta.reml.df %>% 
        filter(abs(mu.delta) <= epsilon.meta) %>% 
        nrow()
    }
    
    df_plot <- delta.reml.df %>%
      arrange(p.unadjusted) %>%
      slice_head(n = screen.plot.topN) %>%
      mutate(
        marker = factor(marker, levels = rev(unique(marker))),
        logp   = -log10(p.unadjusted),
        logp   = pmin(logp, -log10(p_floor)),
        sig    = p.adjusted < alpha,
        # Clip the CI to the plotting range
        ci.delta.lower = pmax(ci.delta.lower, -1),
        ci.delta.upper = pmin(ci.delta.upper, 1),
      )
    
    # Legend breaks on the natural p-value scale
    p_breaks <- c(1, 0.1, 0.05, p_floor)
    logp_breaks <- -log10(p_breaks)
    
    # Colour positions corresponding to the log scale
    colour_values <- scales::rescale(c(0, -log10(0.05), -log10(0.01), -log10(p_floor)), from = c(0, -log10(p_floor)))
    
    epsilon.meta.rounded = round(epsilon.meta, 3)
    
    # Light shading for equivalence region
    if (alternative == "two.sided") {
      lower_bound = -epsilon.meta
      vline_df <- data.frame(
        xintercept = c(-epsilon.meta, epsilon.meta),
        label = paste0("Equivalence margin = +/-", epsilon.meta.rounded)
      )
    } else {
      lower_bound = -1
      vline_df <- data.frame(
        xintercept = c(-2, epsilon.meta),
        label = paste0("Equivalence bound = ", epsilon.meta.rounded)
      )
    }
    
    shade_df <- data.frame(
      xmin = lower_bound,
      xmax = epsilon.meta,
      ymin = 0.5,
      ymax = nrow(df_plot) + 0.5
    )
    
    # Build the plot
    screen.plot <- ggplot(df_plot, aes(x = mu.delta, y = marker)) +
      
      # Shaded equivalence interval (behind points and CIs)
      geom_rect(
        data = shade_df,
        aes(
          xmin = lower_bound,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax
        ),
        fill = "#B4B4B4",
        alpha = 0.3,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      
      # Horizontal CI segments
      geom_segment(
        aes(
          x = ci.delta.lower,
          xend = ci.delta.upper,
          y = marker,
          yend = marker,
          color = logp
        ),
        linewidth = 1.1,
        lineend = "round"
      ) +
      
      # Points for estimates
      geom_point(aes(
        color = logp,
        shape = sig,
        size = sig
      )) +
      
      # Equivalence margin lines (with legend)
      geom_vline(
        data = vline_df,
        aes(xintercept = xintercept, linetype = label),
        color = "#2E2E2E",
        linewidth = 1,
        alpha = 0.8,
        show.legend = c(
          linetype = TRUE,
          color = FALSE,
          shape = FALSE
        )
      ) +
      
      # Vertical zero reference line
      geom_vline(
        xintercept = 0,
        color = "#B4B4B4",
        linewidth = 0.5,
        alpha = 0.5
      ) +
      
      # Colour scale for raw p-values
      scale_color_gradientn(
        colors = c("#2C7BB6", "grey80", "#D7191C", "#8B0000"),
        values = colour_values,
        limits = c(0, -log10(p_floor)),
        breaks = logp_breaks,
        labels = c("1", "0.1", "0.05", paste0(
          "<", format(p_floor, scientific = TRUE)
        )),
        name = "Raw p-value",
        oob = scales::squish
      ) +
      
      # Shape scale for adjusted significance
      scale_shape_manual(
        values = c(`TRUE` = 18, `FALSE` = 1),
        limits = c(TRUE, FALSE),
        drop = FALSE,
        labels = c(
          `TRUE` = bquote("Adjusted p" <= .(alpha)),
          `FALSE` = bquote("Adjusted p" > .(alpha))
        ),
        name = "Multiplicity-corrected \nsignificance",
        guide = guide_legend(
          override.aes = list(size = c(5, 4))
        )
      ) +
      scale_size_manual(
        values = c(`TRUE` = 5, `FALSE` = 4),
        guide = "none"  # no separate size legend
      ) +
      # Linetype scale for equivalence margins
      scale_linetype_manual(name = NULL, values = 1) +
      
      # Labels and title
      labs(
        x = expression("Pooled effect " ~ mu[delta]),
        y = NULL,
        title = glue::glue("Meta-analysis screening results: Top {screen.plot.topN} markers")
      ) +
      
      # Plot limits
      coord_cartesian(xlim = c(-1, 1)) +
      
      # Theme
      theme_minimal(base_size = 20) +
      theme(
        plot.title         = element_text(
          size = 25,
          face = "bold",
          hjust = 0.5
        ),
        axis.text.y        = element_text(size = 13),
        axis.text.x        = element_text(size = 15),
        axis.title.x       = element_text(size = 30),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.title       = element_text(size = 15),
        legend.text        = element_text(size = 13),
        plot.caption       = element_text(size = 13, hjust = 0)
      )
  }
  
  # If no significant markers, stop
  if (length(significant.markers) == 0) {
    message(
      "No significant markers found in meta-analysis.
          You could try to relax your significant criteria."
    )
    
    gamma.s.plot <- list("similarity.plots" = if (return.study.similarity.plot) {
      similarity.plots
    } else {
      NULL
    },
    "screen.plot" = if (return.screen.plot) {
      screen.plot
    } else {
      NULL
    })
    
    return(
      list(
        "screening.metrics.study" = rise.screen.results.allstudies.df,
        "screening.metrics.meta" = delta.reml.df,
        "significant.markers" = NULL,
        "screening.weights" = NULL,
        "gamma.s.plot" = gamma.s.plot
      )
    )
  }
  
  # Extract the weights of the significant markers
  weights.significant = delta.reml.df %>%
    filter(p.adjusted < alpha) %>%
    dplyr::select(marker, weight)
  
  # Normalise weights if requested
  if (normalise.weights) {
    weights.significant = weights.significant %>%
      mutate(weight = weight / max(weight))
  }
  
  # Now recompute per-study effects and RE summary of the combined marker
  # which is a weighted sum of the significant markers
  
  sone.significant = sone %>%
    as.data.frame() %>%
    dplyr::select(all_of(significant.markers))
  
  szero.significant = szero %>%
    as.data.frame() %>%
    dplyr::select(all_of(significant.markers))
  
  weights.vec <- weights.significant$weight
  names(weights.vec) <- weights.significant$marker
  
  # ensure order matches sone/szero column order
  stopifnot(all(colnames(sone.significant) %in% names(weights.vec)))
  weights.vec <- weights.vec[colnames(sone.significant)]  # re-order to match columns
  
  if (return.evaluate.results == FALSE) {
    return(
      list(
        "screening.metrics.study" = rise.screen.results.allstudies.df,
        "screening.metrics.meta" = delta.reml.df,
        "significant.markers" = significant.markers,
        "screening.weights" = weights.significant,
        "gamma.s.plot" = NULL
      )
    )
  }
  
  # Weighted combination of significant markers, called gamma
  gamma.one  <- as.numeric(as.matrix(sone.significant) %*% weights.vec)
  gamma.zero <- as.numeric(as.matrix(szero.significant) %*% weights.vec)
  
  # Initialise results across studies for gamma
  gamma.results.allstudies <- vector("list", length(all.study.names))
  ix <- 1L
  
  for (study in all.study.names) {
    # For each study
    yzero.study = yzero[studyzero == study]
    yone.study = yone[studyone == study]
    
    szero.study = gamma.zero[studyzero == study]
    sone.study = gamma.one[studyone == study]
    
    # Check if this study is paired
    if (is.null(paired.studies)) {
      paired.study = FALSE
    } else {
      paired.study = ifelse(study %in% paired.studies, TRUE, FALSE)
    }
    
    # epsilon_arg <- if (is.null(power.want.s.study)) {
    #   epsilon.meta
    # } else {
    #   NULL
    # }
    
    epsilon_arg = epsilon.meta
    
    # Apply surrogate test to gamma
    gamma.results.study = test.surrogate.extension(
      yone = yone.study,
      yzero = yzero.study,
      sone = sone.study,
      szero = szero.study,
      epsilon = epsilon_arg,
      alpha = alpha,
      power.want.s = NULL,
      u.y.hyp = u.y.hyp,
      alternative = alternative,
      paired = paired.study
    )
    
    # Extract relevant study results
    gamma.results.allstudies[[ix]] <- gamma.results.study %>%
      unlist() %>%
      t() %>%
      as.data.frame() %>%
      mutate(
        study = study,
        n = (length(yone.study) + length(yzero.study)),
        ci.delta.lower = ci.delta1,
        ci.delta.upper = ci.delta2,
        p = p.delta
      ) %>%
      dplyr::select(
        study,
        epsilon.used,
        n,
        delta.estimate,
        sd.delta,
        ci.delta.lower,
        ci.delta.upper,
        u.y,
        u.s,
        p
      )
    # Increase index
    ix <- ix + 1L
    
  }
  
  # Bind study results for gamma
  gamma.results.allstudies.df = bind_rows(gamma.results.allstudies)
  
  # If some studies have exactly 0 standard error, report this
  if (any(gamma.results.allstudies.df$sd.delta == 0)) {
    sd0.studies = gamma.results.allstudies.df %>%
      filter(sd.delta == 0) %>%
      pull(study) %>%
      unique()
    
    message(
      paste0(
        "Note: studies '",
        paste(sd0.studies, collapse = ", "),
        "' have degenerate (exactly 0) estimation of standard error",
        " for the combined marker delta. ",
        "This is likely due to small sample size and perfect pairwise",
        " concordance between gamma and the primary endpoint.",
        " These studies will be removed from the meta-analysis for delta."
      )
    )
    
    gamma.results.allstudies.df = gamma.results.allstudies.df %>%
      filter(sd.delta != 0)
  }
  
  # Now do random effects summary for gamma
  df.summary.marker = gamma.results.allstudies.df
  
  # Extract delta
  delta.gamma = df.summary.marker %>%
    pull(delta.estimate)
  
  # Extract standard deviation delta
  sd.delta.gamma = df.summary.marker %>%
    pull(sd.delta)
  
  # Random effects meta analysis
  delta.reml.gamma = delta.reml.meta(
    delta.gamma,
    sd.delta.gamma,
    epsilon = epsilon.meta,
    alpha = alpha,
    alternative = alternative,
    test = test,
    meta.analysis.method = meta.analysis.method
  )[["results"]]
  
  study.weights = data.frame(delta.reml.gamma$weights.tau.relative)
  
  # Initialise temporary dataframe for results
  study_label = if (alternative == "two.sided"){
    paste0("Pooled effect (", 100 * (1 - 2 * alpha), "% C.I.)")
  } else {
    paste0("Pooled effect (", 100 * (1 - alpha), "% C.I.)")
  }
  
  df.gamma.temp = data.frame(
    "study" = study_label,
    "n" = length(yone) + length(yzero),
    "delta.estimate" = delta.reml.gamma$mu.delta,
    "sd.delta" = delta.reml.gamma$se.delta,
    "ci.delta.lower" = delta.reml.gamma$ci.delta.lower,
    "ci.delta.upper" = delta.reml.gamma$ci.delta.upper,
    "pi.delta.lower" = delta.reml.gamma$pi.lower,
    "pi.delta.upper" = delta.reml.gamma$pi.upper,
    "u.y" = NA,
    "u.s" = NA,
    "p" = delta.reml.gamma$p,
    "study.weights" = NA
  )
  
  # Compute CCC
  x <- gamma.results.allstudies.df$u.y
  y <- gamma.results.allstudies.df$u.s
  ccc <- (2 * cov(x, y)) / (var(x) + var(y) + (mean(x) - mean(y))^2)
  
  # Initialise dataframe for evaluation output
  evaluation.metrics.meta = data.frame(
    "marker" = "gamma",
    "epsilon" = epsilon.meta,
    "n.studies" = delta.reml.gamma$n.studies,
    "mu.delta" = delta.reml.gamma$mu.delta,
    "se.delta" = delta.reml.gamma$se.delta,
    "ci.delta.upper" = delta.reml.gamma$ci.delta.upper,
    "ci.delta.lower" = delta.reml.gamma$ci.delta.lower,
    "pi.delta.upper" = delta.reml.gamma$pi.upper,
    "pi.delta.lower" = delta.reml.gamma$pi.lower,
    "tau2" = delta.reml.gamma$tau2,
    "Q" = delta.reml.gamma$Q,
    "I2" = delta.reml.gamma$I2,
    "ccc" = ccc,
    "p.unadjusted" = delta.reml.gamma$p,
    "p.adjusted" = delta.reml.gamma$p,
    "weight" = 1
  )
  
  # If fit plot desired, return it
  # This plot shows the effect estimates for each study, with size proportional to sample size
  if (return.fit.plot) {
    # Set a minimum value for the axes
    plot.min.global <- gamma.results.allstudies.df %>%
      dplyr::select(u.y, u.s) %>%
      min() - 0.1
    
    n_vals <- gamma.results.allstudies.df$n
    
    # Helper functions
    round_down <- function(x, step) {
      floor(x / step) * step
    }
    round_up   <- function(x, step) {
      ceiling(x / step) * step
    }
    
    if (length(unique(n_vals)) == 1) {
      # All sample sizes equal → single legend key
      {
        legend_breaks <- unique(n_vals)
        legend_labels <- as.character(unique(n_vals))
      }
    } else {
      {
        min_val <- min(n_vals)
        max_val <- max(n_vals)
        median_val <- median(n_vals)
        
        # Make mid value the closest integer to median but strictly between min and max
        mid_val <- round(median_val)
        if (mid_val <= min_val) {
          mid_val <- min_val + 1
        }
        if (mid_val >= max_val) {
          mid_val <- max_val - 1
        }
        
        legend_breaks <- c(min_val, mid_val, max_val)
        legend_labels <- as.character(legend_breaks)
      }
    }
    
    # Plot with CCC annotation and improved sizing
    # Plot with smallest n always visible (size_min = 5)
    fit.plot <- gamma.results.allstudies.df %>%
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
        title = "Treatment effects on primary response and combined marker on training data",
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
    # Extract results for all studies
    gamma.results.allstudies.df = gamma.results.allstudies.df %>%
      mutate(study.weights = delta.reml.gamma$weights.tau.relative)
    
    if (show.pooled.effect) {
      gamma.results.allstudies.df2 = bind_rows(gamma.results.allstudies.df, df.gamma.temp)
    } else {
      gamma.results.allstudies.df2 = gamma.results.allstudies.df
    }
    
    I2 = delta.reml.gamma$I2
    tau2 = delta.reml.gamma$tau2
    
    # Compute CCC
    x <- gamma.results.allstudies.df$u.y
    y <- gamma.results.allstudies.df$u.s
    ccc <- (2 * cov(x, y)) / (var(x) + var(y) + (mean(x) - mean(y))^2)
    
    # Separate studies vs summary
    studies.df <- gamma.results.allstudies.df2 %>% filter(study != study_label)
    summary.row <- gamma.results.allstudies.df2 %>% filter(study == study_label)
    
    # sort studies by effect size from most negative to most positive
    studies.df <- studies.df %>% 
      arrange(delta.estimate)
    
    # vertical positions: top (k) down to 1; summary at y = 0
    k <- nrow(studies.df)
    studies.df <- studies.df %>% mutate(y = seq(from = k, to = 1))
    summary.row <- summary.row %>% mutate(y = 0)
    
    # Prepare combined plotting df and standard display labels
    plot.df <- bind_rows(studies.df, summary.row) %>%
      mutate(
        is.summary = (study == study_label),
        study.label = study,
        label.pval = ifelse(is.na(p), "", formatC(
          p, format = "f", digits = 3
        )),
        label.n = ifelse(is.na(n), "", as.character(n))
      )
    
    # preserve summary CI for diamond and prevent generic errorbar from drawing for summary
    plot.df <- plot.df %>%
      mutate(
        summary.ci.lower = ifelse(is.summary, ci.delta.lower, NA_real_),
        summary.ci.upper = ifelse(is.summary, ci.delta.upper, NA_real_),
        ci.delta.lower = ifelse(is.summary, NA_real_, ci.delta.lower),
        ci.delta.upper = ifelse(is.summary, NA_real_, ci.delta.upper)
      )
    
    # add prediction-interval row (keeps its pi.delta.* values from summary.row),
    # placed below the pooled summary (y = -1), and prevent generic CI drawing
    if (show.pooled.effect) {
      pi.row <- summary.row %>%
        mutate(
          study = paste0(100 * (1 - alpha), "% Prediction interval"),
          study.label = paste0(100 * (1 - alpha), "% Prediction interval"),
          ci.delta.lower = NA_real_,
          ci.delta.upper = NA_real_,
          p = NA_real_,
          n = NA_integer_,
          study.weights = NA_real_,
          y = -1,
          is.summary = FALSE,
          label.pval = "",
          label.n = ""
        )
      plot.df <- bind_rows(plot.df, pi.row)
    }
    
    # clip CIs to the plotting interval prior to plotting
    plot.df <- plot.df %>%
      mutate(
        ci.delta.lower = ifelse(is.na(ci.delta.lower), NA_real_, pmax(pmin(
          ci.delta.lower, 1
        ), -1)),
        ci.delta.upper = ifelse(is.na(ci.delta.upper), NA_real_, pmax(pmin(
          ci.delta.upper, 1
        ), -1)),
        summary.ci.lower = ifelse(is.na(summary.ci.lower), NA_real_, pmax(pmin(
          summary.ci.lower, 1
        ), -1)),
        summary.ci.upper = ifelse(is.na(summary.ci.upper), NA_real_, pmax(pmin(
          summary.ci.upper, 1
        ), -1))
      )
    
    if (show.pooled.effect) {
      plot.df <- plot.df %>%
        mutate(
          pi.delta.lower = ifelse(is.na(pi.delta.lower), NA_real_, pmax(pmin(
            pi.delta.lower, 1
          ), -1)),
          pi.delta.upper = ifelse(is.na(pi.delta.upper), NA_real_, pmax(pmin(
            pi.delta.upper, 1
          ), -1))
        )
    }
    
    # Labels for weighting
    weights.vec <- studies.df$study.weights
    pct.vec <- weights.vec
    label.wgt.vec <- formatC(pct.vec, format = "f", digits = 1)
    
    # attach weight labels to plot.df (empty for summary)
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
        h <- 0.1   # reduced height (was 0.25)
        diamond.df <- data.frame(
          x = c(
            s$delta.estimate,
            s$summary.ci.upper,
            s$delta.estimate,
            s$summary.ci.lower
          ),
          y = c(s$y + h, s$y, s$y - h, s$y)
        )
      }
    }
    
    # Heterogeneity text from provided object delta.reml.gamma
    tau2.txt <- sub(
      "e-0?",
      "e-",
      # remove leading zero in exponent
      format(signif(delta.reml.gamma$tau2, 1), scientific = TRUE)
    )
    I2.txt   <- formatC(delta.reml.gamma$I2,
                        digits = 1,
                        format = "f")
    ccc.txt = formatC(ccc, digits = 2, format = "f")
    
    # Plot parameters
    base.text.size <- 14
    y.min <- if (show.pooled.effect) {
      -1.5
    } else {
      0
    }
    y.max <- k + 1
    rel.w.left  <- 0.45
    rel.w.mid   <- 1.10
    rel.w.right <- 0.55
    x.min <- -1
    x.max <- 1
    
    # Add shaded acceptable interval with legend
    # Light shading for equivalence region
    if (alternative == "two.sided") {
      lower_bound = -epsilon.meta
    } else {
      lower_bound = -1
    }
    
    shade_df <- data.frame(
      xmin = lower_bound,
      xmax = epsilon.meta,
      ymin = y.min,
      # bottom of plot
      ymax = y.max    # top of plot
    )
    
    # Define a factor for legend
    if (alternative == "two.sided") {
      shade_df$label <- paste0("Equivalence margin = +/-", epsilon.meta.rounded)
    } else {
      shade_df$label <- paste0("Equivalence margin = ", lower_bound)
    }
    
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
    
    # Middle panel: forest plot (uses delta.estimate and ci.delta.* columns)
    forest.mid <- ggplot(plot.df, aes(x = delta.estimate, y = y)) +
      geom_errorbar(
        aes(xmin = ci.delta.lower, xmax = ci.delta.upper),
        width = 0.15,
        linewidth = 0.8,
        orientation = "y"
      ) +
      geom_point(
        data = filter(
          plot.df,!is.summary &
            study != paste0(100 * (1 - alpha), "% Prediction interval")
        ),
        shape = 16,
        size = 3.5
      ) +
      # pooled summary drawn as diamond polygon (width = CI)
      geom_polygon(
        data = diamond.df,
        aes(x = x, y = y),
        inherit.aes = FALSE,
        fill = "#CCCCCC",
        color = "black"
      ) +
      # dashed vertical line from pooled estimate upward
      geom_segment(
        data = filter(plot.df, is.summary),
        aes(
          x = delta.estimate,
          xend = delta.estimate,
          y = y,
          yend = y.max
        ),
        inherit.aes = FALSE,
        linetype = "dotted",
        color = "#4C78A8",
        linewidth = 0.8
      ) +
      # prediction interval row (red horizontal line, no central point)
      geom_errorbar(
        data = filter(plot.df, study == paste0(
          100 * (1 - alpha), "% Prediction interval"
        )),
        aes(xmin = pi.delta.lower, xmax = pi.delta.upper),
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
          "Random-effects meta-analysis of combined marker in training data"
        else
          "Effect sizes of combined marker across studies in training data"
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
        xintercept = 0,
        linetype = "solid",
        color = "grey60"
      ) +
      # Shaded region for acceptable interval
      geom_rect(
        data = shade_df,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax
        ),
        alpha = 0.3,
        inherit.aes = FALSE,
        fill = "#B4B4B4",
        color = NA
      ) +
      # Solid vertical lines for the edges of acceptable interval
      geom_vline(
        xintercept = c(ifelse(alternative == "two.sided", -epsilon.meta, -2), epsilon.meta),
        linetype = "solid",
        color = "#2E2E2E",
        linewidth = 1
      )
    
    # Right panel: p-value / weight / N
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
      # no pooled effect: omit weight column entirely (only p-value and N)
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
    
    # Combine panels and add bottom info row (only when pooled effect shown)
    combined <- plot_grid(
      left.labels,
      forest.mid,
      right.table,
      nrow = 1,
      rel_widths = c(rel.w.left, rel.w.mid, rel.w.right),
      align = "h"
    )
    
    epsilon.meta.rounded = round(epsilon.meta, 3)
    
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
  
  gamma.s.plot <- list(
    "screen.plot" = if (return.screen.plot) {
      screen.plot
    } else {
      NULL
    },
    "fit.plot" = if (return.fit.plot) {
      fit.plot
    } else {
      NULL
    },
    "forest.plot" = if (return.forest.plot) {
      forest.plot
    } else {
      NULL
    },
    "similarity.plots" = if (return.study.similarity.plot) {
      similarity.plots
    } else {
      NULL
    }
  )
  
  return(
    list(
      "screening.metrics.study" = rise.screen.results.allstudies.df,
      "screening.metrics.meta" = delta.reml.df,
      "significant.markers" = significant.markers,
      "screening.weights" = weights.significant,
      "evaluation.metrics.study" = gamma.results.allstudies.df,
      "evaluation.metrics.meta" = evaluation.metrics.meta,
      "gamma.s.plot" = gamma.s.plot
    )
  )
  
}
