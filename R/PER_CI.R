#' Construct a two-sided confidence interval for the pooled
#' effect ratio
#'
#' \code{PER_CI} returns the two-sided level-alpha confidence
#' interval of the pooled effect ratio in a cluster-randomized
#' encouragement experiment.
#'
#' \code{PER_CI} constructs a two-sided level-alpha confidence interval
#' by interting the corresponding hypothesis test for the pooled effect
#' ratio. See \code{\link{PER}} for details on the hypothesis tesing.
#' \code{PER_CI} conducts a grid search with user-specified endpoints
#' and meshsize in order to construct the confidence interval.
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
#' @param lower,upper The lower and upper endpoints of the interval to be searched.
#' @param Q A K times p design matrix containing the covariate information.
#'            See Details of the function \code{\link{PER}}.
#' @param meshsize The meshsize of the grid search.
#' @param alpha The level of the confidence interval.
#'
#' @return A length-2 vector of two endpoints of the confidence interval.
#' @examples
#' R_t = encouraged_clusters$aggregated_outcome
#' R_c = control_clusters$aggregated_outcome
#' d_t = encouraged_clusters$aggregated_treatment
#' d_c = control_clusters$aggregated_treatment
#'
#'# Construct 95% CI for the pooled effect ratio estimand
#'# using the default sample variance estimator, i.e.,
#'# setting Q = NULL.
#'CI = PER_CI(R_t, R_c, d_t, d_c, lower = -0.1, upper = 0.1,
#'            alpha = 0.05)
#'
#' @importFrom graphics plot abline
#' @export
#'

PER_CI <- function(R_t, R_c, d_t, d_c,
                   lower, upper, Q = NULL,
                   meshsize = 0.001,
                   alpha = 0.05){
  K = length(R_t)
  if (is.null(Q)) Q = as.matrix(rep(1,K), ncol = K)

  p_value_vec = NULL
  lambda_vec = seq(lower, upper, meshsize)
  for (lambda_0 in lambda_vec){
    p_value_temp = PER(lambda_0, R_t, R_c,
                       d_t, d_c, Q)$pval
    p_value_vec = c(p_value_vec, p_value_temp)
  }
  CI = lambda_vec[which(p_value_vec > alpha)]
  plot(lambda_vec, p_value_vec, xlab = 'lambda_0', ylab = 'P Value')
  abline(h = alpha, col = 'red', lty = 'dashed')
  return(c(min(CI), max(CI)))
}


