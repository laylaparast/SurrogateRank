#' Function to perform the screening stage of RISE : Two-Stage Rank-Based Identification of
#' High-Dimensional Surrogate Markers
#'
#' @description
#' A set of high-dimensional surrogate candidates are screened one-by-one to identify strong
#' candidates. Strength of surrogacy is assessed through a rank-based measure of the similarity in
#' treatment effects on a candidate surrogate and the primary response. P-values corresponding to
#' hypothesis testing on this measure are corrected for the high number of statistical tests
#' performed.
#'
#' @param yone numeric vector of primary response values in the treated group.
#' @param yzero numeric vector of primary response values in the untreated group.
#' @param sone matrix or dataframe of surrogate candidates in the treated group with dimension
#'   \code{n1 x p} where n1 is the number of treated samples and p the number of candidates. Sample
#'   ordering must match exactly yone.
#' @param szero matrix or dataframe of surrogate candidates in the untreated group with dimension
#'   \code{n0 x p} where n0 is the number of untreated samples and p the number of candidates. Sample
#'   ordering must match exactly yzero.
#' @param alpha significance level for determining surrogate candidates. Default is \code{0.05}.
#' @param power.want.s numeric in (0,1) - power desired for a test of treatment effect based on the
#'   surrogate candidate. Either this or \code{epsilon} argument must be specified.
#' @param epsilon numeric in (0,1) - non-inferiority margin for determining surrogate validity. Either
#'   this or \code{power.want.s} argument must be specified.
#' @param u.y.hyp hypothesised value of the treatment effect on the primary response on the probability
#'   scale. If not given, it will be estimated based on the observations.
#' @param p.correction character. Method for p-value adjustment (see \code{p.adjust()} function).
#'   Defaults to the Benjamini-Hochberg method (\code{"BH"}).
#' @param n.cores numeric giving the number of cores to commit to parallel computation in order to
#'   improve computational time through the \code{pbmcapply()} function. Defaults to \code{1}.
#' @param alternative character giving the alternative hypothesis type. One of
#'   \code{c("less","two.sided")}, where "less" corresponds to a non-inferiority test and "two.sided"
#'   corresponds to a two one-sided test procedure. Default is "two.sided".
#' @param paired logical flag giving if the data is independent or paired. If \code{FALSE} (default),
#'   samples are assumed independent. If \code{TRUE}, samples are assumed to be from a paired design.
#'   The pairs are specified by matching the rows of \code{yone} and \code{sone} to the rows of
#'   \code{yzero} and \code{szero}.
#' @param return.all.screen logical flag. If \code{TRUE} (default), a dataframe will be returned giving
#'   the screening results for all candidates. Else, only the significant candidates will be returned.
#' @param return.all.weights logical flag. If \code{FALSE} (default), a dataframe will be returned giving
#'   weights for significant markers screened. If \code{TRUE}, weights for all markers will be returned. Note
#'   that, if normalised weights are required, these will only be returned for significant markers, and raw
#'   weights will be returned in a second column.
#' @param weight.mode character giving the type of weighting to return. One of
#'   \code{c("inverse.delta","diff.epsilon", or "none")}. The default is \code{"inverse.delta"}, which means
#'   the weights are determined by taking the inverse of the absolute values of delta. If delta is exactly 0,
#'   this is uncomputable and the weight defaults to the inverse of the next closest absolute delta value. If
#'   delta is very close to 0, these estimates can be unstable and extreme. The \code{"diff.epsilon"} option
#'   seeks to aid this by calculating weights as the proportion of the interval between 0 and epsilon cut by
#'   the absolute value of delta, therefore giving delta = 0 a weight of 1 and delta = epsilon a weight of 0.
#'   When \code{"none"}, the weights are set to 1 for every marker.
#' @param normalise.weights logical flag. If \code{TRUE} (default), the weights are normalised by the sum of
#'   all the weights such that the maximum weight is 1, which can help with interpretability.
#' @param return.screen.plot logical flag. If \code{TRUE} (default), returns a plot of the top predictors, sorted by p-value,
#'   from the screening stage. The number of predictors to display is given by the \code{screen.plot.topN} argument, which has default
#'   value 15.
#' @param screen.plot.topN number of predictors to display in the screening results figure, default value is 15.
#' @param screen.plot.point.estimate logical flag. If \code{FALSE} (default), uses the \code{screen.plot.topN} argument to determine how many
#' markers to display on the screen plot. Otherwise, plots all the markers with a point estimate within the equivalence region.   
#' @param verbose logical flag. If \code{TRUE}, prints warning messages. 
#'   
#'
#' @return a list with elements \itemize{
#'   \item \code{screening.metrics} : dataframe of screening results (for each candidate marker - number of observations n,
#'     u.y, u.s, delta, CI, sd, epsilon, p-values).
#'   \item \code{significant.markers}: character vector of markers with \code{p_adjusted < alpha}
#'   \item \code{screening.weights}: dataframe giving marker names and the inverse absolute value of the
#'     associated deltas.
#' }
#'
#' @import dplyr pbmcapply ggplot2
#' @export
#' @author Arthur Hughes
#'
#' @examples
#' # Load high-dimensional example data
# data("example.data.highdim")
# yone <- example.data.highdim$y1
# yzero <- example.data.highdim$y0
# sone <- example.data.highdim$s1
# szero <- example.data.highdim$s0
# rise.screen.result <- rise.screen(yone, yzero, sone, szero, power.want.s = 0.8)
rise.screen <- function(yone,
                        yzero,
                        sone,
                        szero,
                        alpha = 0.05,
                        power.want.s = NULL,
                        epsilon = NULL,
                        u.y.hyp = NULL,
                        p.correction = "BH",
                        n.cores = 1,
                        alternative = "two.sided",
                        paired = FALSE,
                        return.all.screen = TRUE,
                        return.all.weights = FALSE,
                        weight.mode = "inverse.delta",
                        normalise.weights = TRUE,
                        return.screen.plot = TRUE,
                        screen.plot.topN = 15,
                        screen.plot.point.estimate = FALSE,
                        verbose = T) {
  # Data formatting
  ## Convert dataframes to numeric matrices
  if (is.data.frame(sone) | is.data.frame(szero)) {
    sone <- as.matrix(sone)
    szero <- as.matrix(szero)
  }
  
  # If no column names on surrogate candidates, set them as the column indices
  if (is.null(colnames(sone))) {
    colnames(sone) <- paste0("marker", 1:ncol(sone))
    colnames(szero) <- paste0("marker", 1:ncol(szero))
  }
  
  # Validity checks
  
  ## Check same number of samples in primary response and surrogates
  n0 <- length(yzero)
  n1 <- length(yone)
  
  if (nrow(szero) != n0) {
    stop("szero does not have the same number of samples as yzero.")
  }
  
  if (nrow(sone) != n1) {
    stop("sone does not have the same number of samples as yone.")
  }
  
  ## if in paired mode, yone/sone must have exactly the same number of samples as yzero/szero
  if (paired) {
    if (length(yone) != length(yzero)) {
      stop(
        "Paired mode is requested but the number of samples in yone does not match that of yzero."
      )
    } else if (length(sone) != length(szero)) {
      stop(
        "Paired mode is requested but the number of samples in sone does not match that of szero."
      )
    }
  }
  
  ## Check that either epsilon or power.want.s is specified
  if (is.null(epsilon) & is.null(power.want.s)) {
    stop("Must specify either epsilon or power.want.s.")
  }
  
  # Screen markers by applying surrogate test in parallel
  ## First define a function that we can then apply in parallel
  .test_marker <- function(idx) {
    args <- list(
      yone = yone,
      yzero = yzero,
      sone = sone[, idx],
      szero = szero[, idx],
      alternative = alternative,
      paired = paired,
      power.want.s = power.want.s,
      epsilon = epsilon,
      alpha = alpha
    )
    
    res <- do.call(test.surrogate.extension, args)
    
    c(
      delta = res$delta.estimate,
      ci_lower = res$ci.delta[1],
      ci_upper = res$ci.delta[2],
      sd = res$sd.delta,
      epsilon = res$epsilon.used,
      p_unadjusted = res$p.delta
    )
  }
  
  P <- ncol(sone)
  raw_list <- pbmclapply(1:P, .test_marker, mc.cores = n.cores)
  
  # Process the list results into a dataframe
  results <- do.call(rbind, raw_list) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(
      marker        = colnames(sone),
      p_adjusted    = p.adjust(p_unadjusted, method = p.correction)
    ) %>%
    dplyr::select(marker,
                  epsilon,
                  delta,
                  sd,
                  ci_lower,
                  ci_upper,
                  p_unadjusted,
                  p_adjusted)

  # Add a message to warn users about degenerate standard error estimation
  if(any(results$sd == 0)){
    n_zero <- sum(results$sd == 0)
    warning(sprintf(
      "For %d markers, the estimated standard error is 0. This occurs when the candidate marker and the primary endpoint show complete pairwise concordance in the observed sample, leading to a degenerate variance estimate. While this may reflect strong agreement between the marker and the endpoint, the sampling variability cannot be estimated, and standard inference (confidence intervals and p-values) is therefore not reliable. This situation is more likely with small sample sizes or when no discordant pairs are observed. Interpret the results for these markers with caution.",
      n_zero
    ))
  }
  
  # Calculate the treatment effect on the primary response
  u.y = SurrogateRank::delta.calculate.extension(
    yone = yone,
    yzero = yzero,
    sone = yone,
    szero = yzero,
    paired = paired
  )$u.y
  
  u.s = u.y - results$delta
  results$u.y = u.y
  results$u.s = u.s
  results$n = length(yone) + length(yzero)
  
  results = results %>% 
    dplyr::select(marker,
                  epsilon,
                  n,
                  u.y,
                  u.s,
                  delta,
                  sd,
                  ci_lower,
                  ci_upper,
                  p_unadjusted,
                  p_adjusted)
  
  
  # Output screen plot if desired
  
  if (return.screen.plot) {
    epsilon.val = unique(results$epsilon)
    p_floor <- 1e-2   # practical lower bound for the colour scale
    if (screen.plot.point.estimate) {
      screen.plot.topN = results  %>%
        filter(abs(delta) <= epsilon.val) %>%
        nrow()
    }
    
    df_plot <- results  %>%
      arrange(p_unadjusted) %>%
      slice_head(n = screen.plot.topN) %>%
      mutate(
        marker = factor(marker, levels = rev(unique(marker))),
        logp   = -log10(p_unadjusted),
        logp   = pmin(logp, -log10(p_floor)),
        sig    = p_adjusted < alpha,
        # Clip the CI to the plotting range
        ci_lower = pmax(ci_lower, -1),
        ci_upper = pmin(ci_upper, 1),
      )
    
    # Legend breaks on the natural p-value scale
    p_breaks <- c(1, 0.1, 0.05, p_floor)
    logp_breaks <- -log10(p_breaks)
    
    # Colour positions corresponding to the log scale
    colour_values <- scales::rescale(c(0, -log10(0.05), -log10(0.01), -log10(p_floor)), from = c(0, -log10(p_floor)))
    
    epsilon.val.rounded = round(epsilon.val, 3)
    # Light shading for equivalence region
    if (alternative == "two.sided") {
      lower_bound = -epsilon.val
      vline_df <- data.frame(
        xintercept = c(-epsilon.val, epsilon.val),
        label = paste0("Equivalence margin = +/-", epsilon.val.rounded)
      )
    } else {
      lower_bound = -1
      vline_df <- data.frame(
        xintercept = c(-2, epsilon.val),
        label = paste0("Equivalence bound = ", epsilon.val.rounded)
      )
    }
    
    shade_df <- data.frame(
      xmin = lower_bound,
      xmax = epsilon.val,
      ymin = 0.5,
      ymax = nrow(df_plot) + 0.5
    )
    
    # Build the plot
    screen.plot <- ggplot(df_plot, aes(x = delta, y = marker)) +
      
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
          x = ci_lower,
          xend = ci_upper,
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
        labels = c("1", "0.1", "0.05", paste0("<", format(
          p_floor, scientific = TRUE
        ))),
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
        guide = guide_legend(override.aes = list(size = c(5, 4)))
      ) +
      scale_size_manual(
        values = c(`TRUE` = 5, `FALSE` = 4),
        guide = "none"  # no separate size legend
      ) +
      # Linetype scale for equivalence margins
      scale_linetype_manual(name = NULL, values = 1) +
      
      # Labels and title
      labs(
        x = expression("Surrogacy parameter " ~ delta),
        y = NULL,
        title = glue::glue("RISE screening results: Top {screen.plot.topN} markers")
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
      ) + 
      guides(
        color = guide_colorbar(order = 1),
        shape = guide_legend(order = 2),
        linetype = guide_legend(order = 3)
      )
  }
  
  
  
  # Add a message to warn users about degenerate standard error estimation
  if(any(results$sd == 0)) {
    n_zero <- sum(results$sd == 0)
    if (verbose) {
      warning(
        sprintf(
          "For %d markers, the estimated standard error is 0. This occurs when the candidate marker and the primary endpoint show complete pairwise concordance in the observed sample, leading to a degenerate variance estimate. While this may reflect strong agreement between the marker and the endpoint, the sampling variability cannot be estimated, and standard inference (confidence intervals and p-values) is therefore not reliable. This situation is more likely with small sample sizes or when no discordant pairs are observed. Interpret the results for these markers with caution.",
          n_zero
        )
      )
    }
  }

  
  # Retreive names of significant markers
  significant_markers <- results %>%
    filter(p_adjusted < alpha) %>%
    pull(marker)
  
  # Weight calculation
  if (weight.mode == "inverse.delta") {
    # If in inverse delta weight mode, we must account for the fact that we cannot invert 0
    # Instead, we search for the next closest value and set the weight to this value
    
    min.nonzero.delta = min(results %>%
                              filter(delta != 0) %>%
                              pull(delta) %>%
                              abs())
    
    if (min.nonzero.delta == 0) {
      min.nonzero.delta = 1
    }
    
    screening.weights <- results %>%
      mutate(nonzero.delta = ifelse(delta != 0, delta, min.nonzero.delta)) %>%
      mutate(weight = 1 / (abs(nonzero.delta))) %>%
      dplyr::select(marker, weight)
    
  } else if (weight.mode == "diff.epsilon") {
    # In the diff.epsilon mode, we calculate the proportion of the distance between epsilon and 0 cut by delta
    screening.weights <- results %>%
      mutate(weight = (epsilon - abs(delta)) / epsilon) %>%
      dplyr::select(marker, weight)
    
  } else if (weight.mode == "none") {
    # If no weights desired, set all to 1
    screening.weights <- results %>%
      mutate(weight = 1) %>%
      dplyr::select(marker, weight)
  }
  
  if (normalise.weights) {
    # If desired to normalise weights, divide by the max weight
    max.weight = max(screening.weights$weight)
    
    screening.weights = screening.weights %>%
      mutate(weight_unstandardised = weight) %>%
      mutate(weight = ifelse(marker %in% significant_markers, weight /
                               max.weight, NA))
  } else if (!normalise.weights) {
    screening.weights = screening.weights %>%
      mutate(weight_unstandardised = weight)
  }
  
  # If desired that only significant weights be returned
  if (!return.all.weights) {
    screening.weights = screening.weights %>%
      filter(marker %in% significant_markers)
  }
  
  # If desired that only significant screening results be returned
  if (!return.all.screen) {
    results <- results %>%
      filter(marker %in% significant_markers)
  } else {
    results = results
  }
  
  plot <- list(
    "screen.plot" = if (return.screen.plot) {
      screen.plot
    } else {
      NULL
    }
  )
  
  return(
    list(
      screening.metrics   = results,
      significant.markers = significant_markers,
      screening.weights   = screening.weights,
      plot                = plot
    )
  )
}

