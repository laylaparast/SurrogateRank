## quiets concerns of R CMD check about NSE data pronouns / pipeline variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("ci_lower", "ci_upper", "delta", "gamma.rank", "marker",
      "nonzero.delta", "p_adjusted", "p_unadjusted",
      "response.rank", "weight", "where")
  )
}