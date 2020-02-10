deepgmm_hmcfl_compute.lik <- function(y, numobs, k, r, xi.list, A.list, D.list, omega.list, w.list) {
  
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
    

    
    if (length(k)==1){
      A1<-array(A.list[[1]][,,1],dim=c(r[1],r[2]))
      
      mu.tot <- A1 %*% xi.list[[1]][, k.comb[i, 1]]
      
      var.tot <- A1 %*% omega.list[[1]][,,k.comb[i, 1] ] %*%
                                    t(A1) +
                                    D.list[[1]][,,k.comb[i, 1]]
      
      w.tot <- w.list[[1]][k.comb[i, 1]] 
    }
    
    if (length(k)==2){
      A1<-array(A.list[[1]][,,1],dim=c(r[1],r[2]))
      A2<-array(A.list[[2]][,, k.comb[i,2]],dim=c(r[2],r[3]))
      
      mu.tot <- A1  %*% A2  %*% xi.list[[2]][, k.comb[i, 2]]
      
      var.tot <- A1 %*% (A2 %*% omega.list[[2]][,,k.comb[i, 2] ] %*%
                                                     t(A2) +
                                                     D.list[[2]][,,k.comb[i, 2]]) %*% t(A1) + D.list[[1]][,,k.comb[i,1]]
      
      w.tot <- w.list[[1]][k.comb[i, 1]] * w.list[[2]][k.comb[i, 2]] 
    }
    
    
    if (length(k)==3){
      A1<-array(A.list[[1]][,,1],dim=c(r[1],r[2]))
      A2<-array(A.list[[2]][,, k.comb[i,2]],dim=c(r[2],r[3]))
      A3<-array(A.list[[3]][,, k.comb[i,3]],dim=c(r[3],r[4]))
      
      mu.tot <- A1  %*% A2  %*% A3 %*% xi.list[[3]][, k.comb[i, 3]]
      
      var.tot <- A1  %*% (A2 %*% (A3 %*% omega.list[[3]][,,k.comb[i, 3] ] %*%
                                                     t(A3) +
                                                     D.list[[3]][,,k.comb[i, 3]])  %*% t(A2) + D.list[[2]][,,k.comb[i,2]] ) %*% t(A1) + D.list[[1]][,,k.comb[i,1]]
      
      w.tot <- w.list[[1]][k.comb[i, 1]] * w.list[[2]][k.comb[i, 2]] * w.list[[3]][k.comb[i, 3]]
    }

    
    if (length(k)==4){
      A1<-array(A.list[[1]][,,1],dim=c(r[1],r[2]))
      A2<-array(A.list[[2]][,,k.comb[i,2]],dim=c(r[2],r[3]))
      A3<-array(A.list[[3]][,,k.comb[i,3]],dim=c(r[3],r[4]))
      A4<-array(A.list[[4]][,,k.comb[i,4]],dim=c(r[4],r[5]))

      mu.tot <- A1  %*% A2  %*% A3  %*% A4  %*% xi.list[[4]][, k.comb[i, 4]]
      
      var.tot <- A1  %*% ( A2 %*% ( A3 %*% ( A4 %*% omega.list[[4]][,,k.comb[i, 4] ] %*%
                                                      t(A4) +
                                                      D.list[[4]][,,k.comb[i, 4]])  %*% t(A3) + D.list[[3]][,,k.comb[i,3]] ) %*% t(A2) + D.list[[2]][,,k.comb[i,2]] ) %*% t(A1) + D.list[[1]][,,k.comb[i,1]]
      
      w.tot <- w.list[[1]][k.comb[i, 1]] * w.list[[2]][k.comb[i, 2]] * w.list[[3]][k.comb[i, 3]] * w.list[[4]][k.comb[i, 4]]
    }
    
    var.tot <- make.positive.definite(var.tot)    
    var.tot <- makeSymm(var.tot)

    
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
  
  return(list(py = py, py.s = py.s, ps.y = ps.y, k.comb = k.comb,
              s = s, ps.y.list = ps.y.list))
}
