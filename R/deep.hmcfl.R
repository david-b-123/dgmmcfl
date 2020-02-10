deepgmm_hmcfl <- function(y, layers, k, r,
                    it = 250, eps = 0.001, init = 'kmeans', init_est = 'factanal') {
  if (any(tolower(init) == c('kmeans', 'k-means', 'k')))
    init <- 'kmeans'
  if (any(tolower(init) == c('random', 'r')))
    init <- 'random'
  if (any(tolower(init) == c('hclass', 'h')))
    init <- 'hclass'
  if (any(tolower(init_est) == c('factanal', 'factana', 'fact', 'f')))
    init_est <- 'factanal'
  if (class(y) == 'data.frame')
    y <- as.matrix(scale(y))
  
  # check arguments
  # tmp <- valid_args(Y = y, layers = layers, k = k, r = r, it = it,
  # eps = eps, init = init)
  numobs <- nrow(y)
  p <- ncol(y)
  r <- c(p, r)
  
  
  # Initialing parameters
  lst <- list(w = list(), H = list(), mu = list(), psi = list(),
              psi.inv = list(), A = list(), xi = list(), omega = list(), D = list(), D.inv = list(), K = list())
  
  
  for (i in 1 : layers) {
    
    if (i == 1) {
      
      data <- y
    } else {
      
      data <- z[, 1 : r[i], drop = FALSE]
    }
    
    # provide initial parititioning of the observations
    s <- initial_clustering(data, k, i, init)
    
    
    # in case if one of the groups is small
    for  (j in 1 : k[i]) {
      if ((table(s)[j]) < 2) {
        s[sample(1 : numobs, 2, replace = FALSE)] <- j
      }
    }
    
    psi <- psi.inv <- array(0, c(r[i], r[i],k[i]))
    H <- array(0, c(r[i], r[i + 1],k[i]))
    mu <- matrix(0, r[i], k[i])
    
    if (init_est == "factanal") {
      # initialize parameters using factor analysis of covariance matrix
      
      i_lst <- factanal_para(data, s, k, r, i, numobs)
      
      lst$w[i] <- list(i_lst$w)
      lst$H[i] <- list(i_lst$H)
      lst$mu[i] <- list(i_lst$mu)
      lst$psi[i] <- list(i_lst$psi)
      lst$psi.inv[i] <- list(i_lst$psi.inv)
      
      z <- i_lst$z
      
      
      # i_lst<-EMMIXmfa::mcfa(Y = data,g = k[i],q = r[i+1],init_clust = as.factor(s),init_method = "gmf",itmax = 10,nkmeans = 35,nrandom = 35)
      xi <- array(NA, c(r[i+1], k[i]))
      omega <- array(NA, c(r[i+1], r[i+1], k[i]))
      D <- array(NA, c(r[i], r[i], k[i]))
      
      AB <- EMMIXmfa::gmf(scale(data), r[i+1], maxit = 200, lambda = 0.01, cor_rate = 0.9)
      
      svd_of_A <- svd(AB$A)
      A <- svd_of_A$u
      
      U <- data %*% A
      
      for (kki in 1 : k[i]) {
        
        indices  <- which(s == kki)
        xi[, kki]  <- apply(U[indices,, drop = FALSE], 2, mean)
        omega[,, kki] <- cov(U[indices,, drop = FALSE])
        D[,,kki] <- diag(apply((data - U %*% t(A))[indices,], 2, var))
      }
      
      
      if (i==1){
        lst$A[i] <- list(array(A,dim=c(r[i],r[i+1],1)))
      } else {
        lst$A[i] <- list(array(A,dim=c(r[i],r[i+1],k[i])))
      }
      
      lst$xi[i]<-list(xi)
      lst$omega[i]<-list(omega)
      
      lst$D[i] <- list(D)
      lst$D.inv[i] <- list(diag(apply(D,3,function(X){1/diag(X)})))
      
      
    } 
  }
  
  
  ##########################################################
  if (layers == 1) {
    out <- deep.sem.alg.1(y = y, numobs = numobs, p = p, r = r[2], k = k, A.list = lst$A, D.list = lst$D, D.list.inv = lst$D.inv, xi.list = lst$xi, omega.list = lst$omega, w.list = lst$w, it = it, eps = eps)

  }
  if (layers == 2) {
    
    out <- deep.sem.alg.2(y = y, numobs = numobs, p = p, r = r, k = k, A.list = lst$A, D.list = lst$D, D.list.inv = lst$D.inv, xi.list = lst$xi, omega.list = lst$omega, w.list = lst$w, it = it, eps = eps)

  }
  if (layers == 3) {
    out <- deep.sem.alg.3(y = y, numobs = numobs, p = p, r = r, k = k, A.list = lst$A, D.list = lst$D, D.list.inv = lst$D.inv, xi.list = lst$xi, omega.list = lst$omega, w.list = lst$w, it = it, eps = eps)
    
  }
  if (layers == 4) {
    out <- deep.sem.alg.4(y = y, numobs = numobs, p = p, r = r, k = k, A.list = lst$A, D.list = lst$D, D.list.inv = lst$D.inv, xi.list = lst$xi, omega.list = lst$omega, w.list = lst$w, it = it, eps = eps)
    
  }
  
  out$lik <- out$likelihood
  output <- out[c("A", "w", "xi", "omega","D", "lik", "bic", "aic", "clc",
                  "icl_bic", "s", "h")]
  output <- c(output, list(k = k, r = r[-1], numobs = numobs, layers = layers))
  output$call <- match.call()
  class(output) <- "dgmm"
  
  invisible(output)
}
