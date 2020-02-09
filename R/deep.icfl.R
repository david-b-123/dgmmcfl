deepgmm_icfl <- function(y, layers, k, r,
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
    
    psi <- psi.inv <- array(0, c(k[i], r[i], r[i]))
    H <- array(0, c(k[i], r[i], r[i + 1]))
    mu <- matrix(0, r[i], k[i])
    
    if (init_est == "factanal") {
      # initialize parameters using factor analysis of covariance matrix
      
      i_lst <- factanal_para(data, s, k, r, i, numobs)
      
      lst$w[i] <- list(i_lst$w)
      lst$H[i] <- list(i_lst$H)
      lst$mu[i] <- list(i_lst$mu)
      lst$psi[i] <- list(i_lst$psi)
      lst$psi.inv[i] <- list(i_lst$psi.inv)
      z <- array(i_lst$z,dim=c(numobs,r[i+1]))
      
      # i_lst<-EMMIXmfa::mcfa(Y = data,g = k[i],q = r[i+1],init_clust = as.factor(s),init_method = "gmf",itmax = 10,nkmeans = 25,nrandom = 25)
      lst$A[i] <- list(i_lst$H[,,1])
      xi_check <- array(0,dim=c(r[i+1],k[i]))
      omega_check <- array(0,dim=c(r[i+1],r[i+1],k[i]))
      
      for (j in 1:k[i]){
        indices <- which(s == j)
        if (r[i+1]==1){
          xi_check[,j] <- array(mean(z[indices,]),dim=c(r[i+1],1))
          omega_check[,,j] <- array(var(z[indices,]),dim=c(r[i+1],r[i+1],1))
        } else {
          xi_check[,j] <- array(colMeans(z[indices,]),dim=c(r[i+1],1))
          omega_check[,,j] <- array(cov(z[indices,]),dim=c(r[i+1],r[i+1],1))
        }
      }
      
      lst$xi[i]<-list(xi_check)
      lst$omega[i]<-list(omega_check)
      
      lst$D[i] <- list(i_lst$psi)
      lst$D.inv[i] <- list(i_lst$psi.inv)
      
      
      
    } else {
      # initialize parameters using probabilistic principal component analysis
      
      if (init_est != "ppca")
        stop("init_est has to be either `factanal` or `ppca`")
      
      i_lst <- ppca_para(data, s, k, r, i, numobs)
      
      lst$w[i] <- list(i_lst$w)
      lst$H[i] <- list(i_lst$H)
      lst$mu[i] <- list(i_lst$mu)
      lst$psi[i] <- list(i_lst$psi)
      lst$psi.inv[i] <- list(i_lst$psi.inv)
      z <- i_lst$z
      
      lst$A[i] <- list(i_lst$H[,,1])
      xi_check <- array(colMeans(z),dim=c(r[i+1],k[i]))
      omega_check <- array(cov(z),dim=c(r[i+1],r[i+1],k[i]))
      
      for (j in 1:k[i]){
        indices <- which(s == j)
        xi_check[,j] <- array(colMeans(z[indices,]),dim=c(r[i+1],1))
        omega_check[,,j] <- array(cov(z[indices,]),dim=c(r[i+1],r[i+1],1))
      }
      
      lst$xi[i]<-list(xi_check)
      lst$omega[i]<-list(omega_check)
      
      lst$D[i] <- list(i_lst$psi)
      lst$D.inv[i] <- list(i_lst$psi.inv)
      
    }
  }
  
  
  
  
  
  # init <- 'kmeans'
  # 
  # numobs <- nrow(y)
  # p <- ncol(y)
  # # r <- c(p, r)
  # 
  # # Initialing parameters
  # # lst <- list(w = list(), A = list(), xi = list(), omega = list(), D = list(), D.inv = list())
  # # print("done1")
  # 
  # 
  # 
  # for (i in 1 : layers) {
  #   print("done2")
  #   z_old <- matrix(0, nrow = numobs, ncol = r[i+1])
  #   
  #   if (i == 1) {
  #     data <- y
  #     
  #   } else {
  #     if (max(0,r[i+1]-dim(z)[2])>0){
  #       data<-cbind(z,y[,sample(c(1:dim(y)[2]),max(0,r[i+1]-dim(z)[2]+1),replace = F)])
  #     } else {
  #       data<-z
  #     }
  #   }
  #   s_mcfa <- initial_clustering(data, k, i, init)
  #   print(dim(data))
  #   
  #   i_lst<-EMMIXmfa::mcfa(Y = (data),g = k[i],q = r[i+1],itmax = 3,tol=10E-5,nkmeans = 2,nrandom = 2)
  #   
  #   lst$w[i] <- list(matrix(table(s_mcfa) / dim(data)[1]))
  #   print(r[i])
  #   print(r[i+1])
  #   lst$A[i] <- list(matrix(i_lst$A[1:r[i],],r[i],r[i+1]))
  # 
  #   lst$xi[i] <- list(i_lst$xi)
  #   lst$omega[i] <- list(i_lst$omega)
  #   
  #   lst$D[i] <- list(array(i_lst$D[1:r[i],1:r[i]],dim=c(r[i],r[i],k[i])))
  #   lst$D.inv[i] <- list(array(diag(1/diag(i_lst$D[1:r[i],1:r[i]])),dim=c(r[i],r[i],k[i])))
  #   
  #   
  #   for (j in 1:k[i]){
  #     if (r[i+1]==1){
  #       z_old[which(i_lst$clust==unique(i_lst$clust)[j]),] <- rnorm(length(which(i_lst$clust==unique(i_lst$clust)[j])),i_lst$xi[,j],i_lst$omega[,,j])
  #     }
  #     else if (length(which(i_lst$clust==unique(i_lst$clust)[j]))>0) {
  #       z_old[which(i_lst$clust==unique(i_lst$clust)[j]),] <- rmvnorm(length(which(i_lst$clust==unique(i_lst$clust)[j])),i_lst$xi[,j],make.positive.definite(makeSymm(i_lst$omega[,,j])))
  #     }
  #   }
  #   z <- z_old
  # }
  
  ##########################################################
  if (layers == 1) {
    out <- deep.sem.alg.1(y = y, numobs = numobs, p = p, r = r[2], k = k, A.list = lst$A, D.list = lst$D, D.list.inv = lst$D.inv, xi.list = lst$xi, omega.list = lst$omega, w.list = lst$w, it = it, eps = eps)
    # y = y; numobs = numobs; p = p; r = r[2]; k = k; A.list = lst$A; D.list = lst$D; D.list.inv = lst$D.inv; xi.list = lst$xi; omega.list = lst$omega; w.list = lst$w; it = it; eps = eps
    
  }
  if (layers == 2) {
    out <- deep.sem.alg.2(y = y, numobs = numobs, p = p, r = r, k = k, A.list = lst$A, D.list = lst$D, D.list.inv = lst$D.inv, xi.list = lst$xi, omega.list = lst$omega, w.list = lst$w, mu.list = lst$mu, H.list=lst$H, psi.list = lst$psi, it = it, eps = eps)
    # y = y; numobs = numobs; p = p; r = r; k = k; A.list = lst$A; D.list = lst$D; D.list.inv = lst$D.inv; xi.list = lst$xi; omega.list = lst$omega; w.list = lst$w; mu.list = lst$mu; H.list=lst$H; psi.list = lst$psi; it = it; eps = eps; column= column
  }
  if (layers == 3) {
    out <- deep.sem.alg.3(y = y, numobs = numobs, p = p, r = r, k = k, A.list = lst$A, D.list = lst$D, D.list.inv = lst$D.inv, xi.list = lst$xi, omega.list = lst$omega, w.list = lst$w, mu.list = lst$mu, H.list=lst$H, psi.list = lst$psi, it = it, eps = eps)
    # y = y; numobs = numobs; p = p; r = r; k = k; A.list = lst$A; D.list = lst$D; D.list.inv = lst$D.inv; xi.list = lst$xi; omega.list = lst$omega; w.list = lst$w; mu.list = lst$mu; H.list=lst$H; psi.list = lst$psi; it = it; eps = eps; column= column
    
  }
  if (layers == 4) {
    out <- deep.sem.alg.4(y = y, numobs = numobs, p = p, r = r, k = k, A.list = lst$A, D.list = lst$D, D.list.inv = lst$D.inv, xi.list = lst$xi, omega.list = lst$omega, w.list = lst$w, mu.list = lst$mu, H.list=lst$H, psi.list = lst$psi, it = it, eps = eps)
    # y = y; numobs = numobs; p = p; r = r; k = k; A.list = lst$A; D.list = lst$D; D.list.inv = lst$D.inv; xi.list = lst$xi; omega.list = lst$omega; w.list = lst$w; it = it; eps = eps
  }
  
  out$lik <- out$likelihood
  output <- out[c("A", "w", "xi", "omega","D", 'mu','psi', 'H',"lik", "bic", "aic", "clc",
                  "icl_bic", "s", "h")]
  output <- c(output, list(k = k, r = r[-1], numobs = numobs, layers = layers))
  output$call <- match.call()
  class(output) <- "dgmm"
  
  invisible(output)
}
