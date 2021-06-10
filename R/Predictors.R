predict.rlap <-  function(beta.vet,rho,covar,resp, li=0.000000001, ls=0.999999999){
  n <- length(resp)
  s <- length(rho)
  per <- ifelse(2:n%%s == 0, s, 2:n%%s)
  p <- theta <- NULL

  for(t in 1:n){
    theta[t] <- exp(sum(covar[t,]*beta.vet))/(1+exp(sum(covar[t,]*beta.vet)))
  }
  theta[theta<=li] <- li;theta[theta>=ls] <- ls
  p[1]  <-  theta[1]
  c1  <-  theta[2:n]*(1-theta[2:n])
  c2  <-  theta[2:n-1]/(1-theta[2:n-1])
  c3  <-  (1-theta[2:n-1])/theta[2:n-1]

  p[2:n]  <-  theta[2:n] - rho[per]*sqrt(c1)*(sqrt(c2)-resp[2:n-1]*(sqrt(c2)+sqrt(c3)))
  p[p<=li] <- li;p[p>=ls] <- ls


  return(p)
}

pfp.pfn.MC <- function(corte, out){
  y <- out$y; prob <- out$prob; s <- out$s
  N <- length(y); R <- N%/%s
  c0 <- corte[1:s]; c1 <- corte[s+1:s]
  yp <- NULL
  for(r in 1:R){
    for(v in 1:s) {
      t <- (r-1)*s+v
      if(t>1){
        if(y[t-1] == 0 & prob[t] <= c0[v]) yp[t] <- 0
        if(y[t-1] == 0 & prob[t] > c0[v]) yp[t] <- 1
        if(y[t-1] == 1 & prob[t] <= c1[v]) yp[t] <- 0
        if(y[t-1] == 1 & prob[t] > c1[v]) yp[t] <- 1
      }
    }
  }
  a <- length(which(y[2:N]==0 & yp[2:N]==0))
  b <- length(which(y[2:N]==0 & yp[2:N]==1))
  c <- length(which(y[2:N]==1 & yp[2:N]==0))
  d <- length(which(y[2:N]==1 & yp[2:N]==1))
  acc <- (a+d)/(a+b+c+d)
  esp <- ifelse(a+b==0, (a+c)/(a+b+c+d), a/(a+b))
  sen <- ifelse(c+d==0, (b+d)/(a+b+c+d), d/(c+d))
  pfp <- ifelse(b+d==0, (a+b)/(a+b+c+d), b/(b+d))
  pfn <- ifelse(a+c==0, (c+d)/(a+b+c+d), c/(a+c))
  return(list(pfp=pfp, pfn=pfn, sen=sen, esp=esp, acc=acc))
}

min.pfp.MC.func <- function(corte, out){
  out1 <- list(y=out$y, prob=out$prob, s=out$s)
  lamb <- out$lamb; esp.target <- out$esp.target
  res <- pfp.pfn.MC(corte = corte, out = out1)
  func <- (1-res$sen) + lamb*max(0, (1-res$esp-(1-esp.target)))
  return(-func)
}

min.pfp.MC <- function(corte.inic=NULL, lamb=30, esp.target=.95,
                       y, prob, s, maxiter=500){
  if(is.null(corte.inic) |
     (!is.null(corte.inic) & length(corte.inic)!= 2*s))
    corte.inic <- rep(.5, 2*s)
  out <- list(y=y, prob=prob, s=s, lamb=lamb, esp.target=esp.target)
  fit <- GA::ga(type = 'real-valued',
                fitness = min.pfp.MC.func,
                out=out,
                maxiter = maxiter,
                monitor = FALSE,
                lower = rep(0.01, 2*s), upper = rep(0.99, 2*s)
  )

  c0 <- attr(fit, 'solution')[1,1:s]; c1 <- attr(fit, 'solution')[1, s+1:s]
  names(c0) <- names(c1) <- paste('v=',1:s,sep='')
  val <- pfp.pfn.MC(corte = c(c0, c1),
                    out = list(y=y, prob=prob, s=s))

  if(1 - val$esp > 1 - esp.target)
    warning(paste('PFN <= ',esp.target,' (target) not founded!\n Choosing PFP based on minimum PFN!', sep=''))

  return(list(pfp=val$pfp,
              pfn=val$pfn,
              acc=val$acc,
              esp=val$esp,
              sen=val$sen,
              c0=c0,
              c1=c1,
              attain.esp.target=(1 - val$esp <= 1 - esp.target),
              esp.target=esp.target))
}

pfp.pfn.glm <- function(corte, out){
  y <- out$y; prob <- out$prob

  yp <- as.numeric(prob > corte)

  a <- length(which(y==0 & yp==0))
  b <- length(which(y==0 & yp==1))
  c <- length(which(y==1 & yp==0))
  d <- length(which(y==1 & yp==1))
  acc <- (a+d)/(a+b+c+d)
  esp <- ifelse(a+b==0, (a+c)/(a+b+c+d), a/(a+b))
  sen <- ifelse(c+d==0, (b+d)/(a+b+c+d), d/(c+d))
  pfp <- ifelse(b+d==0, (a+b)/(a+b+c+d), b/(b+d))
  pfn <- ifelse(a+c==0, (c+d)/(a+b+c+d), c/(a+c))


  return(list(pfp=pfp, pfn=pfn, sen=sen, esp=esp, acc=acc))
}

min.pfp.glm.func <- function(corte, out){
  out1 <- list(y=out$y, prob=out$prob)
  lamb <- out$lamb; esp.target <- out$esp.target
  res <- pfp.pfn.glm(corte = corte, out = out1)
  func <- (1-res$sen) + lamb*max(0, (1-res$esp-(1-esp.target)))
  return(-func)
}

min.pfp.glm <- function(corte.inic=NULL, lamb=30, esp.target=.95,
                        y, prob, maxiter=500){
  # INICIO FUNCAO
  if(missing(corte.inic)) corte.inic <- .5
  out <- list(y=y, prob=prob, lamb=lamb, esp.target=esp.target)
  fit <- GA::ga(type = 'real-valued',
                fitness = min.pfp.glm.func,
                out=out,
                maxiter = maxiter,
                monitor = FALSE,
                lower = 0.01, upper = 0.99
  )

  c0 <- attr(fit, 'solution')[1,1]; names(c0) <- "v=1"
  val <- pfp.pfn.glm(corte = c0,
                     out = list(y=y, prob=prob))

  if(1 - val$esp > 1 - esp.target)
    warning(paste('PFN <= ',esp.target,' (target) not founded!\n Choosing PFP based on minimum PFN!', sep=''))

  return(list(pfp=val$pfp,
              pfn=val$pfn,
              acc=val$acc,
              esp=val$esp,
              sen=val$sen,
              c0=c0,
              attain.esp.target=(1 - val$esp <= 1 - esp.target),
              esp.target=esp.target))
}


prev.glm <- function(fit, newY=NULL, newX, corte, li=0.000000001, ls=0.999999999){#Only new covariates!!
  prob.glm.p  <- stats::predict(object = fit, newx = newX, type=c("response"))
  n.desc <- dim(newX)[1]

  classif.glm <- as.numeric(prob.glm.p > corte)

  if(!is.null(newY)){
    a <- length(which(newY==0 & classif.glm==0))
    b <- length(which(newY==0 & classif.glm==1))
    c <- length(which(newY==1 & classif.glm==0))
    d <- length(which(newY==1 & classif.glm==1))
    acc <- (a+d)/(a+b+c+d)
    esp <- ifelse(a+b==0, (a+c)/(a+b+c+d), a/(a+b))
    sen <- ifelse(c+d==0, (b+d)/(a+b+c+d), d/(c+d))
    pfp <- ifelse(b+d==0, (a+b)/(a+b+c+d), b/(b+d))
    pfn <- ifelse(a+c==0, (c+d)/(a+b+c+d), c/(a+c))
    return(list(prev=classif.glm, pfp=pfp, pfn=pfn, esp=esp, sen=sen, acc = acc))
  }
  else{
    return(list(prev=classif.glm))
  }
}

prev.MC <- function(beta, rho, newY, newX, y, X, corte0, corte1){
  y.f  <-  c(y, newY)
  X.f  <-  rbind(X, newX)
  prob.f  <-  predict.rlap(beta.vet=beta, rho = rho, covar=X.f, resp=y.f)

  s <- length(rho)
  N <- length(y.f); R <- N%/%s

  N0 <- length(y)+1

  classif.MC <- NULL
  for(r in 1:R){
    for(v in 1:s) {
      t <- (r-1)*s+v
      if(t>1){
        if(y.f[t-1] == 0 & prob.f[t] <= corte0[v]) classif.MC[t] <- 0
        if(y.f[t-1] == 0 & prob.f[t] > corte0[v]) classif.MC[t] <- 1
        if(y.f[t-1] == 1 & prob.f[t] <= corte1[v]) classif.MC[t] <- 0
        if(y.f[t-1] == 1 & prob.f[t] > corte1[v]) classif.MC[t] <- 1
      }
    }
  }
  a <- length(which(y.f[N0:N]==0 & classif.MC[N0:N]==0))
  b <- length(which(y.f[N0:N]==0 & classif.MC[N0:N]==1))
  c <- length(which(y.f[N0:N]==1 & classif.MC[N0:N]==0))
  d <- length(which(y.f[N0:N]==1 & classif.MC[N0:N]==1))

  acc.MC <- (a+d)/(a+b+c+d)
  esp.MC <- ifelse(a+b==0, (a+c)/(a+b+c+d), a/(a+b))
  sen.MC <- ifelse(c+d==0, (b+d)/(a+b+c+d), d/(c+d))
  pfp.MC <- ifelse(b+d==0, (a+b)/(a+b+c+d), b/(b+d))
  pfn.MC <- ifelse(a+c==0, (c+d)/(a+b+c+d), c/(a+c))

  return(list(prev=classif.MC[N0:N], pfp=pfp.MC, pfn=pfn.MC, esp = esp.MC, sen = sen.MC, acc = acc.MC))
}
