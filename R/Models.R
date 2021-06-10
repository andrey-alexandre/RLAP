#' Calculates the analitical form of the log likelihood for Periodic Autocorrelated Logistic Regression
#'
#' @param init_par Initial estimative of the parameters
#' @param out List of parameters (response variable, covariable matrix and penalty factor)
#' @param weights Observation weights. Useful attribute for unbalanced data
#'
#' @return
#' @export
#'
#' @examples
loglik  <-  function(init_par, out, weights=NULL){
  lower_bound <- 1e-9; upper_bound <- 1-1e-9;
  y <- out$resp; x <- out$covar; penalty_factor <- out$penalty_factor;

  par <- init_par/penalty_factor;  n <- length(y); k <- dim(x)[2]

  beta  <-  par[1:k]
  rho  <-  par[(k+1):length(par)]

  s <- length(rho)
  per <- ifelse(2:n%%s == 0, s, 2:n%%s)
  if(is.null(weights)) weights <- rep(1, n-1) else weights <- weights[2:n]
  p <- vector("double", n)

  theta <- exp(x%*%beta)/(1+exp(x%*%beta)); theta <- as.vector(theta)
  theta[theta<=lower_bound] <- lower_bound;theta[theta>=upper_bound] <- upper_bound

  p[1]  <-  theta[1]

  c1  <-  theta[2:n]*(1-theta[2:n])
  c1aux1  <-  theta[2:n-1]*(1-theta[2:n-1]); c1aux <- c(0, c1aux1)
  c2  <-  theta[2:n-1]/(1-theta[2:n-1])
  c3  <-  (1-theta[2:n-1])/theta[2:n-1]
  c4  <-  (theta[2:n-1]^(-.5))*((1-theta[2:n-1])^(-1.5))
  c5  <-  (theta[2:n-1]^(-1.5))*((1-theta[2:n-1])^(-.5))

  A <- sqrt(c1*c1aux[-1])
  B <- sqrt(c1/c1aux[-1])

  p[2:n] <- theta[2:n] + rho[per]*sqrt(c1)*(y[2:n-1]-theta[2:n-1])/sqrt(c1aux1)

  p[p<=lower_bound] <- lower_bound;p[p>=upper_bound] <- upper_bound

  dldp <- ((y - p)/(p*(1-p)))[2:n]
  dldp_aux <- ((y[1] - p[1])/(p[1]*(1-p[1])))

  dBdv0 <- .5*(1-2*theta[2:n])/A
  dBdv1 <- -.5*(1-2*theta[2:n-1])*B/c1aux[-1]

  dpdv0 <- 1+.5*rho[per]*(1-2*theta[2:n])*(y[2:n-1]-theta[2:n-1])/sqrt(c1*c1aux1)
  dpdv1 <- -.5*rho[per]*sqrt(c1)*(y[2:n-1]-theta[2:n-1])*(1-2*theta[2:n-1])/(c1aux1^1.5)

  dv0db <- c1*x[2:n,]
  dv0db_aux <- theta[1]*(1-theta[1])*x[1,,drop=FALSE]

  dv1db <- c1aux1*x[2:n-1,]

  dpdb <- dpdv0*dv0db + dpdv1*dv1db #b = beta, B = B mesmo

  dldb <- dldp*dpdb

  mat.aux <- NULL
  for(v in 1:s) mat.aux <- cbind(mat.aux, ifelse(per==v, 1, 0))
  dpdr <- (y[2:n-1]-theta[2:n-1])*sqrt(c1)*mat.aux/sqrt(c1aux1)

  dldr <- dldp*dpdr

  grad  <-  c(colMeans(weights*dldb), colMeans(weights*dldr))/penalty_factor

  #hessian
  d2ldp2 <- (y[2:n]*(2*p[2:n] - 1) - p[2:n]^2)/(p[2:n]^2*(1 - p[2:n])^2)
  d2ldp2_aux <- (y[1]*(2*p[1] - 1) - p[1]^2)/(p[1]^2*(1 - p[1])^2)

  d2ldr2 <- array(data = 0, dim = c(s, s, n-1))
  for(i in 2:n-1){
    d2ldr2[, , i] <- weights[i] * d2ldp2[i] * t(dpdr[i, , drop = F]) %*% dpdr[i, , drop = F]
  }

  d2Bdv02 <- - A^(-1) * ( 1 + .25 * (c1^(-1)) * (1 - 2 * theta[2:n]) ^ 2 )
  d2pdv02 <- rho[per]*(y[2:n-1]-theta[2:n-1])*d2Bdv02

  d2Bdv12 <- B * (c1aux[-1] ^ -1) *  (1 + .75 * (c1aux[-1] ^ -1) * (1 - 2 * theta[2:n-1]) ^ 2 )
  d2pdv12 <- rho[per]*( d2Bdv12 * (y[2:n-1]-theta[2:n-1]) - 2 * dBdv1 )

  d2Bdv0dv1 <- -.25 * (A^(-1)) * (c1aux[-1]^(-1))*(1 - 2 * theta[2:n]) * (1 - 2 * theta[2:n-1])
  d2pdv0dv1 <- rho[per]*( d2Bdv0dv1 * (y[2:n-1]-theta[2:n-1]) - dBdv0 )


  d2ldb2 <- C <-  array(data = 0, dim = c(k, k, n-1))
  for(i in 2:n-1){
    d2ldb2[, , i] <-
      weights[i] * d2ldp2[i] * t(dpdb[i, , drop = F]) %*% dpdb[i, , drop = F] + #(I)
      weights[i] * dldp[i] * ( #(II)
        (t(dv0db[i, , drop = F]) %*% dv0db[i, , drop = F]) * d2pdv02[i] +
          (t(dv0db[i, , drop = F]) %*% dv1db[i, , drop = F]) * d2pdv0dv1[i] +
          (t(dv1db[i, , drop = F]) %*% dv0db[i, , drop = F]) * d2pdv0dv1[i] +
          (t(dv1db[i, , drop = F]) %*% dv1db[i, , drop = F]) * d2pdv12[i] +
          dpdv0[i] * theta[i+1]*(1-theta[i+1])*(1-2*theta[i+1])*(t(x[i+1,,drop=FALSE])%*%x[i+1,,drop=FALSE]) +
          dpdv1[i] * theta[i]*(1-theta[i])*(1-2*theta[i])*(t(x[i,,drop=FALSE])%*%x[i,,drop=FALSE])
      )
  }

  d2ldbdr <- array(data = 0, dim = c(k, s, n-1))
  for(i in 2:n-1){
    d2ldbdr[, , i] <- weights[i] * d2ldp2[i] *(t(dpdb[i, , drop = F]) %*% dpdr[i, , drop = F]) +
      weights[i] * dldp[i] * (
        (y[i] - theta[i])*dBdv0[i]*(t(dv0db[i,,drop=FALSE]) %*% mat.aux[i,,drop=FALSE]) +
          ((y[i] - theta[i])*dBdv1[i]-B[i])*(t(dv1db[i,,drop=FALSE]) %*% mat.aux[i,,drop=FALSE])
      )
  }

  hess  <-  cbind(rbind((rowSums(x = d2ldb2, dims = 2)), t(rowSums(x = d2ldbdr, dims = 2))),
                  rbind( rowSums(x = d2ldbdr, dims = 2), rowSums(x = d2ldr2, dims = 2)))

  logver <-  -mean(weights*y[2:n]*log(p[2:n])+weights*(1-y[2:n])*log(1-p[2:n]))

  attr(logver, 'gradient')  <-  -grad
  attr(logver, 'hessian')  <-  -hess
  return(logver)
}


#' Returns the gradient vector for log likelihood
#'
#' @param init_par Initial estimative of the parameters
#' @param out List of parameters (response variable, covariable matrix and penalty factor)
#' @param weights Observation weights. Useful attribute for unbalanced data
#'
#' @return
#' @export
#'
#' @examples
loglik_gradient <- function(init_par, out, weights=NULL){
  return(attr(loglik(init_par, out, weights=NULL), 'gradient'))
}


#' Returns the hessian matrix for log likelihood
#'
#' @param init_par Initial estimative of the parameters
#' @param out List of parameters (response variable, covariable matrix and penalty factor)
#' @param weights Observation weights. Useful attribute for unbalanced data
#'
#' @return
#' @export
#'
#' @examples
loglik_hessian <- function(init_par, out, weights=NULL){
  return(attr(loglik(init_par, out, weights=NULL), 'hessian'))
}


#' Fits the RLAP model
#'
#' @param resp The response variable
#' @param covar The covariate matrix
#' @param par0 Initial parameter
#' @param lamb some other parameter
#' @param penalty_factor also a parameter
#' @param weights  Observation weights. Useful attribute for unbalanced data
#' @param lower_bound Probability lower boundary
#' @param upper_bound Probability upper boundary
#' @param maxiter Maximum number of iterations
#' @param trace If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing.
#'
#' @return
#' @export
#'
#' @examples
rlap <- function(resp, covar, par0, penalty_factor=1, weights=NULL, maxiter=500, trace=0){
  out <- list(covar=covar, resp=resp, penalty_factor=penalty_factor)
  otim  <-  stats::optim(par = par0, fn = loglik, gr = loglik_gradient,
                         out=out, w=weights, method = "L-BFGS-B",
                         lower = c(rep(-Inf, dim(covar)[2]), rep(-0.999999, length(par0) - dim(covar)[2])),
                         upper = c(rep(Inf, dim(covar)[2]), rep(0.999999, length(par0) - dim(covar)[2])),
                         hessian = FALSE, control = list(trace=trace))
  seotim <-
    loglik_hessian(init_par = otim$par, out=out, w=weights) %>% solve() %>% diag()

  mat.res <- cbind(otim$par, seotim, 2*(1-stats::pnorm(abs(otim$par/seotim))))
  colnames(mat.res) <- c('Estimate','Std.Error','p-value')
  row.names(mat.res) <- c(colnames(covar), paste0('rho', seq(length(par0)-dim(covar)[2])))

  mat.res$coefficients <- fit_MCP[,1]
  mat.res$std.error <- fit_MCP[,1]
  mat.res$y <- resp
  mat.res$X <- covar

  class(mat.res) <- "rlap"
  return(mat.res)
}
