## quiets concerns of R CMD check about NSE data pronouns / pipeline variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "ci_lower", "ci_upper", "delta", "gamma.rank", "marker",
      "nonzero.delta", "p_adjusted", "p_unadjusted",
      "response.rank", "weight", "where", "u.s",
      "Q", "ci.delta.lower", "ci.delta.upper", "ci.delta1", "ci.delta2", "ci.lower", "ci.upper",
     "delta.estimate", "display.n", "display.pval", "display.wgt", "is.summary", "label.n",
      "label.pval", "label.wgt", "mu.delta", "n.studies", "p", "p.adjusted", "p.delta",
      "p.unadjusted", "sd.delta", "se.delta", "study.label", "study.weight",
      "study.weight.relative", "u.y", "weight.heterogeneity", "weight.strength",
      "weights.tau.sum", "y")
  )
}
