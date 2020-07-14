#' Two-sided double-rank test for Fisher's sharp null
#' hypothesis in a cluster-level proportional treatment effect model
#'
#' \code{double_rank} returns the two-sided p-value testing
#' Fisher's sharp null hypothesis in a cluster-level proportional
#' treatment effect model.
#'
#' Double-rank test statistics is a flexible family of nonparametric
#' test statistics. Function \code{psi} is a function that specifies
#' the relationship between d_k, the normalized rank of the absolute
#' treated-minus-control dose difference in the instrumental variable,
#' and q_k, the normalized rank of the absoluve treated-minus-control
#' dose difference in the observed outcome. For instance, psi(d_k, q_k)
#' = 1 yields the sign test, psi(d_k, q_k) = q_k yields the Wilcoxon
#' signed rank test. The default setting, psi(d_k, q_k) = d_k * q_k, yields
#' the dose-weighted signed rank test.
#'
#
#' @param beta_0 The magnitude of the proportional treatment effect to
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
#' @param Z_t A length-K vector where K is equal to the number of
#'            clusters and the kth entry equal to the encoruagement dose,
#'            i.e., the magnitude of the instrumental variable, of the
#'            encouraged cluster in the kth matched pair of two clusters.
#' @param Z_c A length-K vector where K is equal to the number of
#'            clusters and the kth entry equal to the encoruagement dose,
#'            i.e., the magnitude of the instrumental variable, of the
#'            control cluster in the kth matched pair of two clusters.
#' @param psi A function specifying the score used in the test statistic.
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
#' Z_t = encouraged_clusters$IV
#' Z_c = control_clusters$IV
#'
#'
#'# Test beta = 0 in the proportional treatment effect
#'# model with the help of the double rank test using
#'# default psi(d_k, q_k) = d_k * q_k:
#'res = double_rank(0, R_t, R_c, d_t, d_c, Z_t, Z_c)
#'
#'# Define a new psi function: psi(d_k, q_k) = q_k
#'psi_2 <- function(x, y) y
#'
#'# Using psi_2 and the double rank test is reduced to the
#'#Wilcoxon signed rank test.
#'res_2 = double_rank(0, R_t, R_c, d_t, d_c,
#'      Z_t, Z_c, psi = psi_2)
#'
#'
#' @importFrom stats pchisq
#'
#' @export
#'

double_rank <- function(beta_0, R_t, R_c,
                        d_t, d_c, Z_t, Z_c,
                        psi = NULL){

  K = length(R_t)
  # Treated-minus-control difference in the adjusted outcome
  treated_adj_outcome = R_t - beta_0 * d_t
  control_adj_outcome = R_c - beta_0 * d_c
  A_k = treated_adj_outcome - control_adj_outcome
  sign_A_k = (A_k > 0) + 0 # sign of A_k
  s_k = (A_k != 0) + 0 # A_k not equal to 0

  d_k = rank(abs(Z_t - Z_c))/(K + 1)
  q_k = rank(abs(A_k))/(K+1)

  if (is.null(psi))
    psi_d_q = d_k*q_k # Default is psi(d_k, q_k) = d_k * q_k
  else
    psi_d_q = psi(d_k, q_k)

  test_statistic_dr = sum(sign_A_k*psi_d_q)
  exp_dr = sum(s_k*psi_d_q/2)
  sd_dr = sqrt(sum(s_k * psi_d_q^2)/4)
  dev = (test_statistic_dr - exp_dr)/sd_dr
  # Compute two-sided p-value
  p_value = 1-pchisq(dev^2, df = 1)

  return(list(pval = p_value,
              deviate = dev,
              statistic = test_statistic_dr,
              expectation = exp_dr,
              variance = sd_dr^2))
}




