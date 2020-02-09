deepgmm_icfl_compute.lik <- function(y, numobs, k, xi.list, A.list, D.list, omega.list, w.list, mu.list, H.list, psi.list) {
  
  layers <- length(k)
  p <- ncol(y)
  py <- matrix(0, numobs)
  tot.k <- prod(k)
  py.s <- matrix(0, numobs, tot.k)
  pys <- matrix(0, numobs, tot.k)
  
  k.comb <- apply(t(k), 2, function(x) 1 : x)
  if (is.list(k.comb)) {
    k.comb <- expand.grid(k.comb)
  }
  if (is.matrix(k.comb)) {
    k.comb <- expand.grid(split(t(k.comb), 1 : ncol(k.comb)))
  }
  if (prod(k) == 1) {
    k.comb <- matrix(k.comb, nrow =1)
  }
  
  for (i in 1 : tot.k)  {
    
    mu.tot <- mu.list[[2]][, k.comb[i, 2]]
    
    var.tot <- psi.list[[2]][,, k.comb[i, 2]] + H.list[[2]][,,k.comb[i, 2] ] %*% t(H.list[[2]][,, k.comb[i, 2]])
    
    w.tot <- w.list[[1]][k.comb[i, 1]] * w.list[[2]][k.comb[i, 2]]
    
    if (layers > 2) {
      for (l in 3 : layers) {
        
        tot.H <- diag(dim(H.list[[2]][,,k.comb[i, 2] ])[1])
        for (m in 2 : (l-1)) {
          tot.H <- tot.H %*% H.list[[m]][,,k.comb[i, m] ]
        }
        
        mu.tot <- mu.tot + tot.H %*% mu.list[[l]][, k.comb[i, l]]
        
        var.tot <- var.tot + tot.H %*% (H.list[[l]][,,k.comb[i, l] ] %*%
                                          t(H.list[[l]][,,k.comb[i, l] ]) +
                                          psi.list[[l]][,,k.comb[i, l] ]) %*% t(tot.H)
        
        w.tot <- w.tot * w.list[[l]][k.comb[i, l]]
      }
    }
    
    var.tot<-as.matrix(var.tot)
    mu.tot<-as.matrix(mu.tot)
    
    
    if (length(k)==1){
      mu.tot <- A.list[[1]]  %*% mu.tot
      
      var.tot <- A.list[[1]] %*% var.tot %*%
        t(A.list[[1]]) +
        D.list[[1]][,,k.comb[i, 1]]
    }
    
    if (length(k)==2){
      mu.tot <- A.list[[1]]  %*% mu.tot
      
      var.tot <- A.list[[1]] %*% (var.tot) %*% t(A.list[[1]]) + D.list[[1]][,,k.comb[i,1]]
      
    }
    
    
    if (length(k)==3){
      mu.tot <- A.list[[1]]  %*% mu.tot
      
      var.tot <- A.list[[1]]  %*% (var.tot) %*% t(A.list[[1]]) + D.list[[1]][,,k.comb[i,1]]
      
    }
    
    
    if (length(k)==4){
      mu.tot <- A.list[[1]]  %*% mu.tot
      
      var.tot <- A.list[[1]]  %*% (var.tot) %*% t(A.list[[1]]) + D.list[[1]][,,k.comb[i,1]]
      
    }
    
    
    var.tot <- makeSymm(var.tot)
    var.tot <- make.positive.definite(var.tot)
    
    
    
    py.s[,i] <- dmvnorm(y, c(mu.tot), as.matrix(var.tot), log = TRUE)
    
    if (w.tot == 0) {
      w.tot <-  10^(-320)
    }
    pys[, i] <- log(w.tot) + py.s[, i]
  }
  
  # cc <- 705 - max(pys)
  # pys <- pys + cc
  # py <- rowSums(exp(pys))
  #
  # ps.y <- exp(pys) / matrix(py, numobs, tot.k)
  # ps.y <- ifelse(is.na(ps.y), 1/k, ps.y)
  # py <- exp(-cc) * py
  
  # followning lines were used instead of the lines above
  # to avod problmes with near zero pi_i's. 
  pys_max <- apply(pys, 1, max)
  pys <- sweep(pys, 1, pys_max, '-')
  pys <- exp(pys)
  ps.y <- sweep(pys, 1, rowSums(pys), '/')
  py <- exp(pys_max) * rowSums(pys)
  
  s <- matrix(0, nrow = numobs, ncol = layers)
  ps.y.list <- NULL
  for (l in 1 : layers)  {
    
    psi.y <- matrix(0, nrow = numobs, ncol = k[l])
    for (i in 1 : k[l]) {
      
      index <- (k.comb[, l] == i)
      psi.y[, i] <- rowSums(ps.y[, index, drop = FALSE])
    }
    s[, l] <- apply(psi.y, 1, which.max)
    ps.y.list[[l]] <- psi.y
  }
  
  y.tot<- rmvnorm(numobs, c(mu.tot), as.matrix(var.tot))
  
  return(list(py = py, py.s = py.s, ps.y = ps.y, k.comb = k.comb,
              s = s, ps.y.list = ps.y.list, y.tot=y.tot))
}


