loglik  <-  function(x, y, li=1e-9, ls=1-1e-9, lamb=0, alpha=1, init_par, w=NULL){
  par <- init_par/alpha;  n <- length(y); k <- dim(x)[2]
  
  beta  <-  par[1:k]
  rho  <-  par[(k+1):length(par)]
  
  s <- length(rho)
  per <- ifelse(2:n%%s == 0, s, 2:n%%s)
  if(is.null(w)) w <- rep(1, n-1) else w <- w[2:n]
  p <- vector("double", n)
  
  theta <- exp(x%*%beta)/(1+exp(x%*%beta)); theta <- as.vector(theta)
  theta[theta<=li] <- li;theta[theta>=ls] <- ls
  
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
  
  p[p<=li] <- li;p[p>=ls] <- ls
  
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
  
  grad  <-  c(colMeans(w*dldb), colMeans(w*dldr))/alpha
  
  #hessian
  d2ldp2 <- (y[2:n]*(2*p[2:n] - 1) - p[2:n]^2)/(p[2:n]^2*(1 - p[2:n])^2)
  d2ldp2_aux <- (y[1]*(2*p[1] - 1) - p[1]^2)/(p[1]^2*(1 - p[1])^2)
  
  d2ldr2 <- array(data = 0, dim = c(s, s, n-1))
  for(i in 2:n-1){
    d2ldr2[, , i] <- w[i] * d2ldp2[i] * t(dpdr[i, , drop = F]) %*% dpdr[i, , drop = F]
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
      w[i] * d2ldp2[i] * t(dpdb[i, , drop = F]) %*% dpdb[i, , drop = F] + #(I)
      w[i] * dldp[i] * ( #(II)
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
    d2ldbdr[, , i] <- w[i] * d2ldp2[i] *(t(dpdb[i, , drop = F]) %*% dpdr[i, , drop = F]) +
      w[i] * dldp[i] * (
        (y[i] - theta[i])*dBdv0[i]*(t(dv0db[i,,drop=FALSE]) %*% mat.aux[i,,drop=FALSE]) +
          ((y[i] - theta[i])*dBdv1[i]-B[i])*(t(dv1db[i,,drop=FALSE]) %*% mat.aux[i,,drop=FALSE])
      )
  }
  
  hess  <-  cbind(rbind((rowSums(x = d2ldb2, dims = 2)), t(rowSums(x = d2ldbdr, dims = 2))),
                  rbind( rowSums(x = d2ldbdr, dims = 2), rowSums(x = d2ldr2, dims = 2)))
   
  logver <-  -mean(w*y[2:n]*log(p[2:n])+w*(1-y[2:n])*log(1-p[2:n]))

  attr(logver, 'gradient')  <-  -grad
  attr(logver, 'hessian')  <-  -hess
  return(logver)
}


gr.loglik <- function(x, y, li=1e-9, ls=1-1e-9, lamb=0, alpha=1, init_par, w=NULL){
  return(attr(loglik(x, y, li=1e-9, ls=1-1e-9, lamb=0, alpha=1, init_par, w=NULL), 'gradient'))
}

hess.loglik <- function(x, y, li=1e-9, ls=1-1e-9, lamb=0, alpha=1, init_par, w=NULL){
  return(attr(loglik(x, y, li=1e-9, ls=1-1e-9, lamb=0, alpha=1, init_par, w=NULL), 'hessian'))
}

logistic_MC <- function(resp,covar,par0, lamb=0, alpha=1, w0=NULL,trace=0,li=0.001,ls=0.999, maxiter=500){
  out <- list(covar=covar, resp=resp, li=li, ls=ls, lamb=lamb, alpha=alpha)
  if(lamb==0){
    otim  <-  stats::optim(par = par0, fn = loglik, gr = gr.loglik, out=out, w=w0, method = "L-BFGS-B",
                    lower = c(rep(-Inf, dim(covar)[2]), rep(-0.999999, length(par0) - dim(covar)[2])),
                    upper = c(rep(Inf, dim(covar)[2]), rep(0.999999, length(par0) - dim(covar)[2])),
                    hessian = FALSE, control = list(trace=trace))
    seotim <- diag(solve(hess.loglik(par = otim$par, out = out, w=w0)))
    mat.res <- cbind(otim$par, seotim, 2*(1-pnorm(abs(otim$par/seotim))))
    colnames(mat.res) <- c('Estimate','Std.Error','p-value')
  }else{
    par00 <- par0*alpha
    otim <- lbfgs(call_eval = loglik, call_grad = gr.loglik, vars = par00, out=out, w=w0,
                  invisible = as.numeric(trace==0),
                  orthantwise_c = lamb
    )
    mat.res <- cbind(otim$par/alpha)
    colnames(mat.res) <- c('Estimate')
  }
  row.names(mat.res) <- c(colnames(covar), paste0('rho', seq(length(par0)-dim(covar)[2])))
  return(mat.res)
}
