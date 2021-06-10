#' Title
#'
#' @param X Covariate matrix
#' @param beta A vector containing the beta coefficients
#' @param rho A vector containing the rho coefficients
#' @param lower_bound Probability lower boundary
#' @param upper_bound Probability upper boundary
#'
#' @return
#' @export
#'
#' @examples
gen <- function(X, beta, rho, lower_bound=0.000000001, upper_bound=0.999999999){
  n <- dim(X)[1]
  s <- length(rho)
  theta <- p <- NULL

  theta <- exp(X%*%beta)/(1 + exp(X%*%beta))
  theta[theta<=lower_bound] <- lower_bound
  theta[theta>=upper_bound] <- upper_bound
  p[1] <- theta[1]
  y <- stats::rbinom(1, 1, p[1])
  for (t in 2:n) {
    t_m <- ifelse(t%%s == 0, s, t%%s)
    c1  <-  theta[t]*(1-theta[t])
    c2  <-  theta[t-1]/(1-theta[t-1])
    if(y[t-1]==0) p[t]  <-  theta[t] - rho[t_m]*sqrt(c1*c2) else p[t]  <-  theta[t] + rho[t_m]*sqrt(c1/c2)
    if(p[t]>1) p[t] <- 1
    if(p[t]<0) p[t] <- 0
    y[t] <- stats::rbinom(n = 1, size = 1, prob = p[t])
  }
  return(y)
}
