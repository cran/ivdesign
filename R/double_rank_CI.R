#' Construct a two-sided confidence interval for the proportional
#' treatment effect in a cluster-level proportional treatment effect model
#'
#' \code{double_rank_CI} returns the two-sided level-alpha confidence
#' interval of the proportional treatment effect in a cluster-level
#' proportional treatment effect model.
#'
#' \code{double_rank_CI} constructs a two-sided level-alpha confidence interval
#' by interting the hypothesis test using a double_rank test. Function
#' \code{double_rank_CI} conducts a grid search with user-specified endpoints
#' and meshsize in order to construct the confidence interval. For more
#' details on the double_rank test, see \code{\link{double_rank}}.
#
#' @param R_t A length-K vector where K is equal to the number of
#'            clusters and the kth entry equal to
#'            the sum of unit-level outcomes in the encouraged cluster
#'            of the kth matched pair of two clusters.
#' @param R_c A length-K vector where K is equal to the number of
#'            clusters and the kth entry equal to
#'            the sum of unit-level outcomes in the control cluster
#'            of the kth matched pair of two clusters.
#' @param d_t A length-K vector where K is equal to the number of
#'            clusters and the kth entry equal to the sum of unit-level
#'            treatment received in the encouraged cluster
#'            of the kth matched pair of two clusters.
#' @param d_c A length-K vector where K is equal to the number of
#'            clusters and the kth entry equal to the sum of unit-level
#'            treatment received in the control cluster
#'            of the kth matched pair of two clusters.
#' @param Z_t A length-K vector where K is equal to the number of
#'            clusters and the kth entry equal to the encoruagement dose,
#'            i.e., the magnitude of the instrumental variable, of the
#'            encouraged cluster in the kth matched pair of two clusters.
#' @param Z_c A length-K vector where K is equal to the number of
#'            clusters and the kth entry equal to the encoruagement dose,
#'            i.e., the magnitude of the instrumental variable, of the
#'            control cluster in the kth matched pair of two clusters.
#' @param lower,upper The lower and upper endpoints of the interval to be searched.
#' @param meshsize The meshsize of the grid search.
#' @param psi A function specifying the score used in the test statistic.
#'            See Details of \code{\link{double_rank}}.
#' @param alpha The level of the confidence interval.
#'
#' @return A length-2 vector of two endpoints of the confidence interval.
#' @examples
#'
#' R_t = encouraged_clusters$aggregated_outcome
#' R_c = control_clusters$aggregated_outcome
#' d_t = encouraged_clusters$aggregated_treatment
#' d_c = control_clusters$aggregated_treatment
#' Z_t = encouraged_clusters$IV
#' Z_c = control_clusters$IV
#'
#'
#'# Construct a level 0.05 CI for the constant proportional
#'# treatment effect with the help of the double rank test using
#'# default psi(d_k, q_k) = d_k * q_k. Search from -0.1 to 0.1:
#'CI = double_rank_CI(R_t, R_c, d_t, d_c, Z_t, Z_c,
#'                     lower = -0.1, upper = 0.1)
#'
#'
#' @importFrom graphics plot abline
#' @export
#'


double_rank_CI <- function(R_t, R_c,
                           d_t, d_c,
                           Z_t, Z_c,
                           lower, upper,
                           meshsize = 0.001,
                           psi = NULL,
                           alpha = 0.05){
  p_value_vec = NULL
  beta_vec = seq(lower, upper, meshsize)
  for (beta_0 in beta_vec){
    p_value_temp = double_rank(beta_0, R_t, R_c,
                               d_t, d_c, Z_t, Z_c,
                               psi)$pval
    p_value_vec = c(p_value_vec, p_value_temp)
  }
  CI = beta_vec[which(p_value_vec > alpha)]
  plot(beta_vec, p_value_vec, xlab = 'beta_0', ylab = 'P Value')
  abline(h = alpha, col = 'red', lty = 'dashed')
  return(c(min(CI), max(CI)))
}


