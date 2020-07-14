#' Two-sided test for the average cluster effect ratio estimand
#'
#' \code{ACER} tests (two-sided) if the average cluster effect
#' ratio (ACER) is equal to lambda.
#'
#'
#' @param num_t A length-K vector where K is equal to the number of
#'             clusters and the kth entry equal to the number of
#'             units in the encouraged cluster of the kth matched
#'             pair of two clusters.
#' @param num_c A length-K vector with the kth entry equal to the number of
#'             units in the control cluster of the kth matched
#'             pair of two clusters.
#' @param R_t A length-K vector with kth entry equal to
#'            the sum of unit-level outcomes in the encouraged cluster
#'            of the kth matched pair of two clusters.
#' @param R_c A length-K vector with the kth entry equal to
#'            the sum of unit-level outcomes in the control cluster
#'            of the kth matched pair of two clusters.
#' @param d_t A length-K vector with the kth entry equal to the sum
#'            of unit-level treatment received in the encouraged cluster
#'            of the kth matched pair of two clusters.
#' @param d_c A length-K vector with the kth entry equal to the sum
#'            of unit-level treatment received in the control cluster
#'            of the kth matched pair of two clusters.
#' @param lambda The magnitude of the average cluster effect ratio (ACER)
#'               to be tested.
#' @param alpha The level of the test.
#' @param kappa Minimum compliance rate.
#' @param gap Relative MIP optimality gap.
#' @param verbose If true, the solver output is enabled; otherwise,
#'             the solver output is disabled.
#'
#' @return A list of three elements: the optimal solution, the optimal
#'         objective value, and an indicator of whether or not the test
#'         is rejected.
#' @examples
#' \dontrun{
#' # To run the following example, Gurobi must be installed.
#'
#' R_t = encouraged_clusters$aggregated_outcome
#' R_c = control_clusters$aggregated_outcome
#' d_t = encouraged_clusters$aggregated_treatment
#' d_c = control_clusters$aggregated_treatment
#' num_t = encouraged_clusters$number_units
#' num_c = control_clusters$number_units
#'
#'# Test at level 0.05 if the ACER is equal
#'# to 0.2. Assume the minimum compliance rate across
#'# K clusters is at least 0.2. Set verbose = FALSE
#'# to suppress the output.
#'res = ACER(num_t, num_c, R_t, R_c, d_t, d_c,
#'           lambda = 0.2, alpha = 0.05, kappa = 0.2,
#'           verbose = FALSE)
#'
#'# The test is rejected
#'res$Reject
#'}
#' @importFrom stats sd qnorm qchisq
#' @export
#'

ACER <- function(num_t, num_c, R_t, R_c,
                 d_t, d_c, lambda,
                 alpha = 0.05, kappa = 0.1,
                 gap = 0.05, verbose = TRUE){

  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package \"gurobi\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package \"Matrix\" needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  # Compute the objective function Q matrix
  compute_Q_obj <- function(R, C_1, C_2){
    K = nrow(R)

    Q_obj = matrix(0, nrow = 2*K, ncol = 2*K)
    for (k in 1:K){
      Q_obj[k, k] = C_1 * R[k, 1]^2
    }

    for (k in (K+1):(2*K)){
      Q_obj[k, k] = C_1 * R[k - K, 2]^2
    }

    for (k in 1:K){
      Q_obj[k, k + K] = (-C_1)*R[k, 1]*R[k, 2]
      Q_obj[k + K, k] = (-C_1)*R[k, 1]*R[k, 2]
    }

    # C_2 = (K - 1 + chi^2)/(K^2*(K-1))

    for (i in 1:K){
      for (j in 1:K){
        if (i != j) Q_obj[i, j] = C_2 * R[i, 1]*R[j, 1]
      }
    }

    for (i in 1:K){
      for (j in 1:K){
        if (i != j) Q_obj[i + K, j + K] = C_2 * R[i, 2]*R[j, 2]
      }
    }

    for (i in 1:K){
      for (j in 1:K){
        if (i != j) Q_obj[i, j + K] = (-C_2) * R[i, 1]*R[j, 2]
      }
    }

    for (i in 1:K){
      for (j in 1:K){
        if(i != j) Q_obj[i + K, j] = (-C_2) * R[i, 2]*R[j, 1]
      }
    }
    return(Q_obj)
  }

  # Compute constant C1
  compute_C_1 <- function(K, alpha) return((1 - qchisq(1 - alpha/2, df=1))/K^2)

  # Compute constant C2
  compute_C_2 <- function(K, alpha) return((K - 1 + qchisq(1-alpha/2, df=1))/(K^2*(K-1)))

  # Compute vector c in the objective function
  compute_c_vec <- function(R, lambda){
    K = nrow(R)
    c_vec = c(R[,1], -R[,2])
    return(((-2*lambda)/K)*c_vec)
  }

  # Compute the quadratic constraint matrix Qc
  compute_Q_c <- function(i, K){
    Q_c = Matrix::Matrix(0, nrow = 4*K, ncol = 4*K, sparse = TRUE)
    Q_c[i, i + 2*K] = 1/2
    Q_c[i + 2*K, i] = 1/2
    return(Q_c)
  }

  #########################################################

  R = cbind(R_t, R_c)
  D = cbind(d_t, d_c)
  Num = cbind(num_t, num_c)
  K = dim(R)[1]

  # Upper bound on CO_kj
  ub_data = c(D[,1], Num[,2] - D[,2])

  # Compute the level-alpha/2 CI of the average compliance
  # rate using sample variance as a conservative variance
  # estimator

  D_k = D[,1]/Num[,1] - D[,2]/Num[,2]
  c_alpha_4 = qnorm(1 - alpha/4)
  L = mean(D_k) - c_alpha_4*sd(D_k)/sqrt(K)
  U = mean(D_k) + c_alpha_4*sd(D_k)/sqrt(K)

  ##################################################
  # Initiate a model
  model = list()

  ####################################################
  # Q is the quadratic part in the objective function
  C_1 = compute_C_1(K, alpha)
  C_2 = compute_C_2(K, alpha)
  Q_extract = cbind(Matrix::Matrix(0, 2*K, 2*K, sparse = TRUE),
                    Matrix::.sparseDiagonal(2*K, 1))
  Q_obj = Matrix::Matrix(compute_Q_obj(R, C_1, C_2), sparse = TRUE)
  Q_final_obj = t(Q_extract)%*%Q_obj%*%Q_extract
  model$Q = Q_final_obj

  ###################################################
  # c is the linear part in the objective function
  c_obj = compute_c_vec(R, lambda)
  c_final_obj = c(rep(0, 2*K), c_obj)
  model$obj = c_final_obj

  ########################################################
  # lb and ub are upper and lower bound of decision variables
  lb = c(pmin(rep(1, 2*K), ub_data), rep(0, 2*K))
  ub = c(ub_data, rep(1, 2*K))

  model$lb = lb
  model$ub = ub

  ###########################################################
  # Q_i in the quadratic constraint that enforces the
  # inverse condition

  for (i in 1:(2*K)){
    quad_constraint = list(Qc = compute_Q_c(i, K), rhs = 1, sense = '=')
    model$quadcon[[i]] = quad_constraint
  }

  #############################################################
  # Linear constraint: norm upper bound
  # N_vec = (n_11, n_21, ..., n_K1, n_12, n_22, ..., n_K2)
  N_vec = c(Num[,1], Num[,2])
  U_constraint = c(1/N_vec, rep(0, 2*K))
  L_constraint = c(-1/N_vec, rep(0, 2*K))

  A = matrix(c(U_constraint, L_constraint), nrow = 2, byrow = TRUE)
  A = Matrix::Matrix(A, sparse = TRUE)
  A2 = Matrix::.sparseDiagonal(2*K, -1)
  A2 = cbind(A2, Matrix::Matrix(0, nrow = 2*K, ncol = 2*K, sparse = TRUE))

  model$A = rbind(A, A2)
  RHS = pmin(ub_data, kappa*cbind(Num[,1], Num[,2]))
  model$rhs = c(2*K*U, -2*K*L, -RHS)

  model$vtype = c(rep('I', 2*K), rep('C', 2*K))
  model$modelsense = 'min'

  if (verbose) flag = 1
  else flag = 0

  res = gurobi::gurobi(model, params = list(NonConvex = 2, BestObjStop = -lambda^2,
                                    BestBdStop = -lambda^2,
                                    MIPGap = gap,
                                    OutputFlag = flag))

  if (is.null(res$objval))
    reject_status = TRUE
  else if (res$objval >= -lambda^2)
    reject_status = TRUE
  else reject_status = FALSE

  solution = list(ObjectiveValue = res$objval,
                  Solution = res$x,
                  Reject = reject_status)

  return(solution)
}

