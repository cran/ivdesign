#' Two-sided test for the pooled effect ratio estimand
#'
#' \code{PER} returns the two-sided p-value testing
#' the pooled effect ratio equal to lambda_0 in a cluster-randomized
#' encouragement experiment.
#'
#' Q is used to construct a regression-assisted variance estimator.
#' Q is can in principle be any K times p design matrix such that p < K.
#' When Q is a column vector of 1's, the variance estimator is the classical
#' sample variance estimator. More generally, Q may contain any cluster-level
#' or even unit-level covariate information that are predictive of the
#' encouraged-minus-control difference in the observed aggregated outcomes.
#'
#
#' @param lambda_0 The magnitude of the pooled effect ratio estimand to
#'               be tested.
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
#' @param Q A K times p design matrix containing the covariate information.
#'            See Details.
#'
#' @return A list of five elements: two-sided p-value, deviate, test statistics,
#'         expectation of the test statistic under the null hypothesis,
#'         and variance of the test statistic under the null hypothesis.
#' @examples
#' R_t = encouraged_clusters$aggregated_outcome
#' R_c = control_clusters$aggregated_outcome
#' d_t = encouraged_clusters$aggregated_treatment
#' d_c = control_clusters$aggregated_treatment
#'
#'
#'# Test the pooled effect ratio estimand lambda = 0 using
#'# the default sample variance estimator, i.e., setting Q = NULL.
#'res = PER(0, R_t, R_c, d_t, d_c)
#'
#'# We may leverage observed covariates from both the encouraged
#'# and control clusters to construct less conservative variance
#'# estimator. The variance estimator will be less conservative if
#'# these covariate predict the treated-minus-control difference
#'# in the outcome. In this illustrated dataset, V1-V10 are simulated
#'# white noise; it is not surprising that they do not help
#'# reduce the variance.
#'Q = cbind(encouraged_clusters[,1:10], control_clusters[,1:10])
#'res_2 = PER(0, R_t, R_c, d_t, d_c, Q)
#'
#' @importFrom stats pchisq
#' @export
#'
PER <- function(lambda_0, R_t, R_c, d_t, d_c, Q = NULL){
  K = length(R_t)

  # Default Q is a column vector of 1s
  if (is.null(Q)) Q = as.matrix(rep(1,K), ncol = K)
  else Q = as.matrix(Q)

  treated_adj_outcome = R_t - lambda_0 * d_t
  control_adj_outcome = R_c - lambda_0 * d_c

  # treated-minus-control difference in adjusted outcome
  Y_k = treated_adj_outcome - control_adj_outcome

  hat_Q = Q %*% solve(t(Q) %*% Q) %*% t(Q)
  Y_Q = as.matrix(Y_k/sqrt(1 - diag(hat_Q)), ncol = K)

  # Construct the variance estimator
  var_est = as.numeric(1/K * t(Y_Q) %*% (diag(K) - hat_Q) %*% Y_Q)
  test_statistic = sum(Y_k)/K
  dev = sqrt(K)*test_statistic/sqrt(var_est)
  p_value = 1-pchisq(dev^2, df = 1)

  return(list(pval = p_value,
              deviate = dev,
              statistic = test_statistic,
              expectation = 0,
              variance = var_est))

}
