#' @title Two-Way Likelihood Ratio Test under Tree Order Restriction
#' @description
#' Performs the likelihood ratio test (LRT) for testing
#' \eqn{H0: \theta_{0j} = \theta_{1j} = ... = \theta_{pj}} versus the
#' tree-ordered alternative \eqn{H1: \theta_{0j} \le \theta_{ij}} for all i,
#' with at least one strict inequality for all i.
#'
#' @param sample_data A list-of-lists: sample_data[[i]][[j]] = numeric vector for cell (i,j).
#' @param significance_level Significance level (default = 0.05).
#' @importFrom stats quantile rnorm var
#' @import Matrix
#' @import MASS
#' @export
#' @export
#' @export
TreeLRT <- function(sample_data, significance_level){
  set.seed(123)
  p <- length(sample_data)-1
  q <- length(sample_data[[1]])-1

  R_MLE <- function(X, n) {
    X1 <- X[-1]
    n1 <- n[-1]
    sorted_indices <- order(X1)
    X1_sorted <- X1[sorted_indices]
    n1_sorted <- n1[sorted_indices]
    if (length(X)==2){
      if (all(X1 < X[1])){
        new_X <- c((X[1]*n[1]+X[2]*n[2])/sum(n),(X[1]*n[1]+X[2]*n[2])/sum(n))
      }
      else{
        new_X <- X
      }
      return(new_X)
    }
    else {
      A <- numeric(length(X1_sorted))
      for (j in 2:length(X)) {
        A[j-1] <- (n[1] * X[1] + sum(n1_sorted[1:(j - 1)] * X1_sorted[1:(j - 1)])) /
          (n[1] + sum(n1_sorted[1:(j - 1)]))
      }
      if (all(X1 >= X[1])) {
        new_X <- X
      } else if (A[length(X)-2] >= X1_sorted[length(X)-1]) {
        X <- rep(A[length(X)-1], length(X))
        new_X <- X
      } else {
        comparisons <- logical(length(X1_sorted) - 1)
        comparisons1 <- logical(length(X1_sorted) - 1)
        stored_values <- numeric(0)
        for (k in 1:(length(X1_sorted) - 1)) {
          comparisons[k] <- A[k] < X1_sorted[k + 1]
          if(comparisons1[k] <- A[k] < X1_sorted[k + 1]) {
            for (s in 1:k) {
              stored_values[s] <- X1_sorted[s]
            }
            break
          }
        }
        selected_A_values <- A[comparisons]
        X[1] <- selected_A_values[1]
        for (l in 2:length(X)) {
          if (X[l] %in% stored_values) {
            X[l] <- selected_A_values[1]
          }
        }
        new_X <- X
      }
      return(new_X)
    }
  }

  ### MLE under Null parameter space
  Func_MLE_0<-function(sample_data){## Cell means
    cell_means <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j)
        mean(sample_data[[i]][[j]])
      )
    )
    cell_means <- matrix(cell_means,
                         nrow = p+1,
                         ncol = q+1,
                         byrow = TRUE)

    ## Cell sample sizes
    cell_ns <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j)
        length(sample_data[[i]][[j]])
      )
    )
    cell_ns <- matrix(cell_ns,
                      nrow = p+1,
                      ncol = q+1,
                      byrow = TRUE)

    ## Initial cell variances
    cell_vars <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j)
        ((cell_ns[i,j]-1)/cell_ns[i,j]) *
          var(sample_data[[i]][[j]])
      )
    )
    cell_vars <- matrix(cell_vars,
                        nrow = p+1,
                        ncol = q+1,
                        byrow = TRUE)

    ## Initial beta
    beta_0_0 <- colMeans(cell_means)

    ## Initial cell variances
    var_0_0 <- cell_vars

    repeat{

      ## Update beta
      beta_1_0 <- sapply(1:(q+1), function(j){
        w <- cell_ns[,j] / var_0_0[,j]
        sum(w * cell_means[,j]) / sum(w)
      })

      ## Update cell variances
      var_1_0 <- sapply(1:(p+1), function(i)
        sapply(1:(q+1), function(j)
          sum((sample_data[[i]][[j]] - beta_1_0[j])^2) /
            cell_ns[i,j]
        )
      )

      var_1_0 <- matrix(var_1_0,
                        nrow = p+1,
                        ncol = q+1,
                        byrow = TRUE)

      ## Stopping criterion
      if (max(abs(beta_1_0 - beta_0_0)) < 1e-5)
        break

      beta_0_0 <- beta_1_0
      var_0_0  <- var_1_0
    }

    return(var_1_0)}

  ### MLE under Full parameter space
  Func_MLE_1<-function(sample_data){
    theta_0 <- matrix(0, nrow = p+1, ncol = q+1)
    beta_0  <- numeric(q+1)
    var_0   <- matrix(0, nrow = p+1, ncol = q+1)

    theta_1 <- matrix(0, nrow = p+1, ncol = q+1)
    beta_1  <- numeric(q+1)
    var_1   <- matrix(0, nrow = p+1, ncol = q+1)

    ## Cell means
    cell_means <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j)
        mean(sample_data[[i]][[j]])
      )
    )
    cell_means <- matrix(cell_means,
                         nrow = p+1,
                         ncol = q+1,
                         byrow = TRUE)

    ## Cell sample sizes
    cell_ns <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j)
        length(sample_data[[i]][[j]])
      )
    )
    cell_ns <- matrix(cell_ns,
                      nrow = p+1,
                      ncol = q+1,
                      byrow = TRUE)

    ## Initial cell variances
    cell_vars <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j)
        ((cell_ns[i,j]-1)/cell_ns[i,j]) *
          var(sample_data[[i]][[j]])
      )
    )
    cell_vars <- matrix(cell_vars,
                        nrow = p+1,
                        ncol = q+1,
                        byrow = TRUE)

    ## Initial theta
    grand_mean <- mean(cell_means)
    row_means  <- rowMeans(cell_means)

    theta_0 <- sweep(cell_means, 1, row_means, "-") + grand_mean

    ## Initial beta
    beta_0 <- colMeans(cell_means) - grand_mean

    ## Initial variances
    var_0 <- cell_vars

    repeat{

      ## h matrix
      h <- cell_ns / var_0

      ## Y matrix
      Y <- sweep(cell_means, 2, beta_0, "-")

      ## Update theta
      theta_1 <- sapply(1:(q+1), function(j)
        R_MLE(Y[,j], h[,j])
      )

      theta_1 <- matrix(theta_1,
                        nrow = p+1,
                        ncol = q+1,
                        byrow = FALSE)

      ## Update beta
      t_b <- sapply(1:q, function(j)
        sum(h[,j] * (cell_means[,j] - theta_1[,j]) -
              h[,q+1] * (cell_means[,q+1] - theta_1[,q+1]))
      )

      Q_diag <- sapply(1:q, function(j)
        (p+1) * mean(h[,j])
      )

      if (q == 1) {
        Q <- matrix(Q_diag, nrow = 1, ncol = 1)
      } else {
        Q <- diag(Q_diag)
      }

      Q <- Q + (p+1) * mean(h[,q+1]) * matrix(1, nrow = q, ncol = q)

      beta_temp <- solve(Q, t_b)

      beta_1[1:q] <- beta_temp
      beta_1[q+1] <- -sum(beta_temp)

      ## Update variances
      var_1 <- sapply(1:(p+1), function(i)
        sapply(1:(q+1), function(j)
          sum((sample_data[[i]][[j]] -
                 theta_1[i,j] -
                 beta_1[j])^2) /
            cell_ns[i,j]
        )
      )

      var_1 <- matrix(var_1,
                      nrow = p+1,
                      ncol = q+1,
                      byrow = TRUE)

      ## Stopping criterion
      theta_diff <- apply(abs(theta_1 - theta_0), 2, max)

      if (all(theta_diff < 1e-5) &&
          max(abs(beta_1 - beta_0)) < 1e-5)
        break

      theta_0 <- theta_1
      beta_0  <- beta_1
      var_0   <- var_1
    }

    return(var_1)}### Func_MLE_1
  Func_MLE_1(sample_data)
  Calculate_lambda <- function(sample_data){

    ## Cell sample sizes
    cell_ns <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j)
        length(sample_data[[i]][[j]])
      )
    )

    cell_ns <- matrix(cell_ns,
                      nrow = p+1,
                      ncol = q+1,
                      byrow = TRUE)

    ## MLEs under H1 and H0
    mle1 <- Func_MLE_1(sample_data)
    mle0 <- Func_MLE_0(sample_data)

    ## Likelihood ratio
    lambda <- prod((mle1/mle0)^(cell_ns/2))

    return(lambda)
  }
  sample_data <- lapply(sample_data, function(row) {
    lapply(row, function(cell) cell[!is.na(cell)])
  })

  n_sample <- 10000

  ## Cell sample sizes
  cell_ns <- sapply(1:(p+1), function(i)
    sapply(1:(q+1), function(j)
      length(sample_data[[i]][[j]])
    )
  )

  cell_ns <- matrix(cell_ns,
                    nrow = p+1,
                    ncol = q+1,
                    byrow = TRUE)

  lambda_values_star <- numeric(n_sample)

  for (k in seq_len(n_sample)) {

    bootstrap_samples <- lapply(1:(p+1), function(i) {
      lapply(1:(q+1), function(j) {
        rnorm(cell_ns[i,j],
              mean = 0,
              sd = sqrt(var(sample_data[[i]][[j]])))
      })
    })

    lambda_values_star[k] <- Calculate_lambda(bootstrap_samples)
  }

  critical_value <- quantile(lambda_values_star,
                             probs = significance_level)

  test_statistic <- Calculate_lambda(sample_data)

  decision <- if (test_statistic < critical_value) {
    "Reject null hypothesis"
  } else {
    "Do not reject null hypothesis"
  }

  return(list(
    critical_value = critical_value,
    test_statistic = test_statistic,
    decision = decision
  ))
}



#' @title MinMin test for Two-Way Tree-Ordered Data
#' @description Performs the MinMin test for two-way ANOVA with tree-ordered treatment effects.
#' @param sample_data List-of-lists of numeric vectors: sample_data[[i]][[j]] = data in cell (i,j)
#' @param significance_level Significance level (default = 0.05)
#' @return A list with critical_value, test_statistic, and decision
#' @examples
#' TreeMinMin(sample_data, significance_level)
#' @export

TreeMinMin <- function(sample_data, significance_level = 0.05) {

  ## -----------------------------------
  ## DIMENSIONS
  ## -----------------------------------
  p <- length(sample_data)-1
  q <- length(sample_data[[1]])-1

  ## Remove NA values in every cell
  sample_data <- lapply(sample_data, function(row) {
    lapply(row, function(cell) cell[!is.na(cell)])
  })

  ## -----------------------------------
  ## Max Test Function  (Recomputed inside)
  ## -----------------------------------
  MinMin_Func <- function(sample_data) {

    # sample sizes (recomputed inside)
    cell_ns <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j) length(sample_data[[i]][[j]]))
    )
    cell_ns<-matrix(cell_ns,nrow = (p+1),ncol = (q+1),byrow = TRUE)
    # cell means
    cell_mean_data <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j) mean(sample_data[[i]][[j]]))
    )
    cell_mean_data<-matrix(cell_mean_data,nrow = (p+1),ncol = (q+1),byrow = TRUE)
    # row means

    M <- matrix(0, p, q+1)

    for(i in 1:p){
      for(j in 1:(q+1)){

        denom <- sqrt(
          var(sample_data[[i+1]][[j]]) / cell_ns[i+1,j] +
            var(sample_data[[1]][[j]])   / cell_ns[1,j]
        )

        M[i,j] <-
          (cell_mean_data[i+1,j] - cell_mean_data[1,j]) / denom
      }
    }

    ## MinMin statistic
    min(apply(M, 2, min))
  }

  ## OBSERVED TEST STATISTIC
  ## -----------------------------------
  obs_minmin <- MinMin_Func(sample_data)
  ## BOOTSTRAP SAMPLING
  boots_sample <- 10000

  # Precompute cell variances & sizes (based on original data)
  cell_vars <- sapply(1:(p+1), function(i)
    sapply(1:(q+1), function(j) var(sample_data[[i]][[j]]))
  )
  cell_vars<-matrix(cell_vars,nrow = (p+1),ncol = (q+1), byrow = TRUE)
  cell_ns <- sapply(1:(p+1), function(i)
    sapply(1:(q+1), function(j) length(sample_data[[i]][[j]]))
  )
  cell_ns<-matrix(cell_ns,nrow = (p+1),ncol = (q+1), byrow = TRUE)
  minmin_boot <- numeric(boots_sample)

  for (k in 1:boots_sample) {

    # generate bootstrap sample with same sizes & variances
    bootstrap_samples <- lapply(1:(p+1), function(i) {
      lapply(1:(q+1), function(j) {
        rnorm(
          n = cell_ns[i, j],
          mean = 0,
          sd = sqrt(cell_vars[i, j])
        )
      })
    })

    # compute max statistic in bootstrap
    minmin_boot[k] <- MinMin_Func(bootstrap_samples)
  }
  critical_val <- quantile(minmin_boot, probs = 1 - significance_level)
  decision <- ifelse(obs_minmin > critical_val,
                     "Reject null hypothesis",
                     "Do not reject null hypothesis")
  return(list(
    critical_value = critical_val,
    test_statistic = obs_minmin,
    decision = decision
  ))
}






#' @title MinMax test for Two-Way Tree-Ordered Data
#' @description Performs the MinMax test for two-way ANOVA with tree-ordered treatment effects.
#' @param sample_data List-of-lists of numeric vectors: sample_data[[i]][[j]] = data in cell (i,j)
#' @param significance_level Significance level (default = 0.05)
#' @return A list with critical_value, test_statistic, and decision
#' @examples
#' TreeMinMax(sample_data, significance_level)
#' @export

TreeMinMax <- function(sample_data, significance_level = 0.05) {

  ## -----------------------------------
  ## DIMENSIONS
  ## -----------------------------------
  p <- length(sample_data)-1
  q <- length(sample_data[[1]])-1

  ## Remove NA values in every cell
  sample_data <- lapply(sample_data, function(row) {
    lapply(row, function(cell) cell[!is.na(cell)])
  })

  ## -----------------------------------
  ## Max Test Function  (Recomputed inside)
  ## -----------------------------------
  MinMax_Func <- function(sample_data) {

    # sample sizes (recomputed inside)
    cell_ns <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j) length(sample_data[[i]][[j]]))
    )
    cell_ns<-matrix(cell_ns,nrow = (p+1),ncol = (q+1),byrow = TRUE)
    # cell means
    cell_mean_data <- sapply(1:(p+1), function(i)
      sapply(1:(q+1), function(j) mean(sample_data[[i]][[j]]))
    )
    cell_mean_data<-matrix(cell_mean_data,nrow = (p+1),ncol = (q+1),byrow = TRUE)
    # row means

    M <- matrix(0, p, q+1)

    for(i in 1:p){
      for(j in 1:(q+1)){

        denom <- sqrt(
          var(sample_data[[i+1]][[j]]) / cell_ns[i+1,j] +
            var(sample_data[[1]][[j]])   / cell_ns[1,j]
        )

        M[i,j] <-
          (cell_mean_data[i+1,j] - cell_mean_data[1,j]) / denom
      }
    }

    ## MinMin statistic
    min(apply(M, 2, max))
  }

  ## OBSERVED TEST STATISTIC
  ## -----------------------------------
  obs_minmax <- MinMax_Func(sample_data)
  ## BOOTSTRAP SAMPLING
  boots_sample <- 10000

  # Precompute cell variances & sizes (based on original data)
  cell_vars <- sapply(1:(p+1), function(i)
    sapply(1:(q+1), function(j) var(sample_data[[i]][[j]]))
  )
  cell_vars<-matrix(cell_vars,nrow = (p+1),ncol = (q+1), byrow = TRUE)
  cell_ns <- sapply(1:(p+1), function(i)
    sapply(1:(q+1), function(j) length(sample_data[[i]][[j]]))
  )
  cell_ns<-matrix(cell_ns,nrow = (p+1),ncol = (q+1), byrow = TRUE)
  minmax_boot <- numeric(boots_sample)

  for (k in 1:boots_sample) {

    # generate bootstrap sample with same sizes & variances
    bootstrap_samples <- lapply(1:(p+1), function(i) {
      lapply(1:(q+1), function(j) {
        rnorm(
          n = cell_ns[i, j],
          mean = 0,
          sd = sqrt(cell_vars[i, j])
        )
      })
    })

    # compute max statistic in bootstrap
    minmax_boot[k] <- MinMax_Func(bootstrap_samples)
  }
  critical_val <- quantile(minmax_boot, probs = 1 - significance_level)
  decision <- ifelse(obs_minmax > critical_val,
                     "Reject null hypothesis",
                     "Do not reject null hypothesis")
  return(list(
    critical_value = critical_val,
    test_statistic = obs_minmax,
    decision = decision
  ))
}




#' @title Run all three tree-ordered tests
#' @description
#' Convenience wrapper that runs LRT, Min, and Max tests and returns all results.
#' @param sample_data list-of-lists cells
#' @param significance_level numeric (default 0.05)
#' @return list(LRT=..., MinMin=..., MinMax=...)
#' @export
TreeTwoWay <- function(sample_data, significance_level) {
  lrt <- TreeLRT(sample_data, significance_level = significance_level)
  MinMin  <- TreeMinMin(sample_data, significance_level = significance_level)
  MinMax  <- TreeMinMax(sample_data, significance_level = significance_level)

  list(LRT = lrt, MinMin = MinMin, MinMax = MinMax)
}
