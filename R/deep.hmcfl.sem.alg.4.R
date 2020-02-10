deepgmm_hmcfl_sem.alg.4 <- function(y, numobs, p, r, k, A.list, D.list,
                           D.list.inv, xi.list, omega.list, w.list, it, eps) {
  HHHHH<-10
  likelihood <- NULL
  hh <- 0
  ratio <- Inf
  layers <- length(k)
  #################################
  #### compute the likelihood #####
  #################################
  out <- deepgmm_hmcfl_compute.lik(y, numobs, k,r=r, xi.list, A.list, D.list, omega.list, w.list)
  py <- out$py
  ps.y <- out$ps.y
  ps.y.list <- out$ps.y.list
  k.comb <- out$k.comb
  s <- out$s
  tot.k <- prod(k)
  temp <- sum(log(py))
  likelihood <- c(likelihood, temp)
  #################################
  z.list <- NULL
  while ((hh < it) & (ratio > eps )) {
    hh <- hh+1
    # print(hh)
    ###############################################################
    ###################  first layer ##############################
    ###############################################################
    l <- 1
    yy <- y
    z2 <- z2.one <- array(0, c(numobs, r[l + 1], k[l], k[l + 1], k[l + 2], k[l + 3]))
    zz2 <- array(0, c(numobs, r[l +1], r[l + 1], k[l], k[l + 1], k[l + 2], k[l + 3]))
    z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
    zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))
    
    for (p1 in 1 : k[l]) {
      for (p2 in 1 : k[l + 1]) {
        for (p3 in 1 : k[l + 2]) {
          for (p4 in 1 : k[l + 3]) {
            sigma.tilde.inv <- ginv(A.list[[l+1]][,,p2] %*% (ginv(A.list[[l+2]][,,p3] %*% (
              A.list[[l+3]][,,p4] %*% omega.list[[l+3]][,,p4]  %*% t(A.list[[l+3]][,,p4]) + D.list[[l+3]][,,p4]
            ) %*% t(A.list[[l+2]][,,p3])) )%*% t(A.list[[l+1]][,,p2]))
            
            A <- sigma.tilde.inv + t(A.list[[l]][,,1]) %*%
              (ginv(D.list[[l]][,, p1])) %*% (A.list[[l]][,,1])
            
            mu.tilde  <- c(A.list[[l+1]][,,p2] %*% A.list[[l+2]][,,p3] %*% A.list[[l+3]][,,p4] %*% xi.list[[l+3]][,p4])
            
            b <- c(sigma.tilde.inv %*% mu.tilde) +
              t(A.list[[l]][,,1]) %*% ginv(D.list[[l]][,,p1]) %*%
              (t(yy))
            
            chsi  <- ginv(A)
            
            roy <- chsi %*% b
            roy.quadro <- array(apply(roy, 2, function(x) x %*% t(x)),
                                c(r[l + 1], r[l + 1], numobs))
            zz2[,,, p1, p2, p3, p4] <- aperm(array(chsi,
                                               c(r[l + 1], r[l + 1], numobs)) +
                                           roy.quadro, c(3, 1, 2))
            z2.one[,, p1, p2, p3, p4] <- rmvnorm(numobs, rep(0, r[l + 1]), chsi) +
              t(roy)
            z2[,, p1, p2, p3, p4] <- t(roy)
          }
        }
      }
    }
    
    for (i1 in 1 : k[l + 1]) {
      for (i2 in 1 : k[l + 2]) {
          for (i3 in 1 : k[l + 3]) {
            prob <- ps.y.list[[l + 1]][, i1, drop = FALSE] *
              ps.y.list[[l + 2]][, i2, drop = FALSE] *
                ps.y.list[[l + 3]][, i3, drop = FALSE]
            z <- z + array(z2[,,, i1, i2, i3, drop = FALSE] *
                             array(rowSums(prob), c(numobs, r[l + 1], k[l], 1, 1, 1)),
                           c(numobs, r[l + 1], k[l]))
            z.one  <- z.one + array(z2.one[,,, i1, i2, i3, drop = FALSE] *
                                      array(rowSums(prob), c(numobs, r[l + 1], k[l], 1, 1, 1)),
                                    c(numobs, r[l + 1], k[l]))
            zz <- zz + array(zz2[,,,, i1, i2, i3, drop = FALSE] *
                               array(rowSums(prob), c(numobs, r[l + 1], r[l + 1], k[l], 1, 1, 1)),
                             c(numobs, r[l + 1], r[l + 1], k[l]))
          }
      }
    }
    
    z.list[[l]] <- aperm(z.one, c(3, 1, 2))
    
    
    
    
    
    ### M-step get A
    ### M-step get D
    ### M-step get xi
    Y<- yy
    
    # load("/home/davb/Documents/Samsung_T5/2019_11_21_Deep_Statistical_Models/Deep_Gaussian_Mixture_Model/scripts/deepgmm-branch_3_0_0/r_backup/backup_v1.RData")
    p <- ncol(Y)
    if (is.null(p))
      p <- 1
    n <- nrow(Y)
    
    
    A1 <- array(0, c(p, r[l+1]))
    A2 <- array(0, c(r[l+1], r[l+1]))
    # Di <- array(0,c(p,k[l]))
    
    tau<-(ps.y.list[[l]])
    
    omega<-array(omega.list[[l]],dim=c(r[l+1],r[l+1],k[l]))
    xi<-array(xi.list[[l]],dim=c(r[l+1],k[l]))
    AA<-array(A.list[[l]],dim=c(r[l],r[l+1],1))
    D<-D.list[[l]]
    
    w1<-rep(0,k[l])
    
    for (i in 1:k[l]) {
      
      inv_Di <- diag(1 / diag(D[,,i]))
      inv_Di_A <- sweep(AA[,,1], 1, diag(1/D[,,i]), '*')
      # gamma = (A Omega_i A^T + D) ^{-1} A Omega_i
      # using Woodbury Identity
      gamma <- (inv_Di - inv_Di_A %*% 
                  chol.inv(chol.inv(omega[,, i]) +  t(AA[,,1]) %*% inv_Di_A) %*%
                  t(inv_Di_A)) %*% AA[,,1] %*% omega[,, i]
      
      ti <- sum(tau[, i])
      xi_i <- xi[, i, drop = FALSE]
      
      # tau_ij * y_j  
      tY <- sweep(Y, 1, tau[, i], '*')
      
      # y_j - A xi_i
      Y_Axi_i <- sweep(Y, 2, AA[,,1] %*% xi_i, '-')
      
      # tau_ij * (y_j - A xi_i)
      tY_Axi_i <- sweep(Y_Axi_i, 1, tau[, i], '*')
      
      # xi_i = xi_i + gamma^T \sum{ tau_ij (y_j - A xi_i) } / sum {tau_ij}
      xi[, i] <- xi_i + t(gamma) %*% as.matrix(colSums(tY_Axi_i)) / ti
      
      # 
      omega[,, i] <- (diag(r[l+1]) - t(gamma) %*% AA[,,1]) %*% omega[,, i] +
        t(gamma) %*% t(Y_Axi_i) %*% tY_Axi_i %*% gamma / ti -
        (xi_i - xi[,i]) %*% t(xi_i - xi[, i])
      # tau_ij yj xi_i^T +   tau_ij yj (yj - A xi)^T gamma
      A1 <- A1 + colSums(tY) %*% t(xi_i) + t(Y) %*% tY_Axi_i %*% gamma
      #
      A2 <- A2 + (omega[,, i] + xi[, i] %*% t(xi[, i])) * ti
      
      #Di <- Di + diag(t(Y) %*% tY)
      # Di <- Di + colSums(Y * tY)
      
      beta<- ginv(AA[,,1] %*% omega[,,i] %*% t(AA[,,1]) + D[,,i])%*%D[,,i]
      
      D1 =  diag(colSums(tau[, i]%*%t(diag(D[,,i]%*%(diag(p)-beta))))) / ti
      D2 = diag(diag(t(beta)%*%t(tY_Axi_i)%*%Y_Axi_i%*%beta / ti))
      D[,,i] = D1+D2
      
      
      w1[i] <- ti / n
    }
    
    AA[,,1] <- try(A1 %*% chol.inv(A2))
    if (class(AA) == "try-error") {
      
      model <- "tried to invert an ill-conditioned or a singular matrix"
      class(model) <- "error"
      return(model)
    }
    
    
    # DD <- diag(Di - rowSums((AA %*% A2) * AA )) / n
    
    A.list[[l]] <- AA
    # D.list[[l]] <- DD
    D.list[[l]] <- D
    xi.list[[l]] <- xi
    omega.list[[l]] <- omega
    w.list[[l]] <- w1
    
    
    
    
    ###############################################################
    ###################  second layer #############################
    ###############################################################
    
    l <- 2
    yy <- matrix(0, numobs, r[l])
    zz <- z.list[[l - 1]]
    for (i in 1 : k[l - 1]) {
      yy <- yy + matrix(zz[i,,, drop = FALSE] *
                          array(rowSums(ps.y.list[[l - 1]][, i, drop = FALSE]),
                                c(1, numobs, r[l])), numobs, r[l])
    }
    z2 <- z2.one <- array(0, c(numobs, r[l + 1], k[l], k[l + 1], k[l + 2]))
    zz2 <- array(0, c(numobs, r[l + 1], r[l + 1], k[l], k[l + 1], k[l + 2]))
    z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
    zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))
    
    for (p1 in 1 : k[l]) {
      for (p2 in 1:k[l+1])  {
        for (p3 in 1:k[l+2])  {
          sigma.tilde.inv <- ginv(A.list[[l+1]][,,p2] %*% (
            A.list[[l+2]][,,p3] %*% omega.list[[l+2]][,,p3]  %*% t(A.list[[l+2]][,,p3]) + D.list[[l+2]][,,p3]
          ) %*% t(A.list[[l+1]][,,p2]))
          
          A <- sigma.tilde.inv + t(A.list[[l]][,,p1]) %*%
            (ginv(D.list[[l]][,, p1])) %*% (A.list[[l]][,,p1])
          
          mu.tilde  <- c(A.list[[l+1]][,,p2] %*% A.list[[l+2]][,,p3] %*% xi.list[[l+2]][,p3])
          
          b <- c(sigma.tilde.inv %*% mu.tilde) +
            t(A.list[[l]][,,p1]) %*% ginv(D.list[[l]][,,p1]) %*%
            (t(yy))
          
          
          chsi  <- ginv(A)
          
          roy  <- chsi %*% b
          roy.quadro  <- array(apply(roy, 2, function(x) x %*%
                                       t(x)), c(r[l + 1], r[l + 1], numobs))
          zz2[,,, p1, p2, p3]  <- aperm(array(chsi, c(r[l + 1], r[l + 1], numobs)) +
                                      roy.quadro, c(3, 1, 2))
          z2.one[,, p1, p2, p3]  <- rmvnorm(numobs, rep(0, r[l + 1]), chsi) + t(roy)
          z2[,, p1, p2, p3]  <- t(roy)
        }
      }
    }
    
    for (i1 in 1 : k[l + 1]) {
      for (i2 in 1 : k[l + 2]) {
        prob <- ps.y.list[[l + 1]][, i1, drop = FALSE] *
          ps.y.list[[l + 2]][, i2, drop = FALSE]
        z <- z + array(z2[,,, i1, i2, drop = FALSE] *
                         array(rowSums(prob), c(numobs, r[l + 1], k[l], 1, 1)),
                       c(numobs, r[l + 1], k[l]))
        z.one  <- z.one + array(z2.one[,,, i1, i2, drop = FALSE] *
                                  array(rowSums(prob), c(numobs, r[l + 1], k[l], 1, 1)),
                                c(numobs, r[l + 1], k[l]))
        zz <- zz + array(zz2[,,,, i1, i2, drop = FALSE] *
                           array(rowSums(prob), c(numobs, r[l + 1], r[l + 1], k[l], 1, 1)),
                         c(numobs, r[l + 1], r[l + 1], k[l]))
      }
    }
    
    
    z.list[[l]] <- aperm(z.one, c(3, 1, 2))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ### M-step get A
    ### M-step get D
    ### M-step get xi
    Y<- yy
    
    # load("/home/davb/Documents/Samsung_T5/2019_11_21_Deep_Statistical_Models/Deep_Gaussian_Mixture_Model/scripts/deepgmm-branch_3_0_0/r_backup/backup_v1.RData")
    p <- ncol(Y)
    if (is.null(p))
      p <- 1
    n <- nrow(Y)
    
    
    A1 <- array(0, c(p, r[l+1]))
    A2 <- array(0, c(r[l+1], r[l+1]))
    # Di <- array(0,c(p,k[l]))
    
    tau<-(ps.y.list[[l]])
    
    omega<-array(omega.list[[l]],dim=c(r[l+1],r[l+1],k[l]))
    xi<-array(xi.list[[l]],dim=c(r[l+1],k[l]))
    AA<-array(A.list[[l]],dim=c(r[l],r[l+1],k[l]))
    D<-D.list[[l]]
    
    w1<-rep(0,k[l])
    
    for (i in 1:k[l]) {
      
      inv_Di <- diag(1 / diag(D[,,i]))
      inv_Di_A <- sweep(AA[,,i], 1, diag(1/D[,,i]), '*')
      # gamma = (A Omega_i A^T + D) ^{-1} A Omega_i
      # using Woodbury Identity
      gamma <- (inv_Di - inv_Di_A %*% 
                  chol.inv(chol.inv(omega[,, i]) +  t(AA[,,i]) %*% inv_Di_A) %*%
                  t(inv_Di_A)) %*% AA[,,i] %*% omega[,, i]
      
      ti <- sum(tau[, i])
      xi_i <- xi[, i, drop = FALSE]
      
      # tau_ij * y_j  
      tY <- sweep(Y, 1, tau[, i], '*')
      
      # y_j - A xi_i
      Y_Axi_i <- sweep(Y, 2, AA[,,i] %*% xi_i, '-')
      
      # tau_ij * (y_j - A xi_i)
      tY_Axi_i <- sweep(Y_Axi_i, 1, tau[, i], '*')
      
      # xi_i = xi_i + gamma^T \sum{ tau_ij (y_j - A xi_i) } / sum {tau_ij}
      xi[, i] <- xi_i + t(gamma) %*% as.matrix(colSums(tY_Axi_i)) / ti
      
      # 
      omega[,, i] <- (diag(r[l+1]) - t(gamma) %*% AA[,,i]) %*% omega[,, i] +
        t(gamma) %*% t(Y_Axi_i) %*% tY_Axi_i %*% gamma / ti -
        (xi_i - xi[,i]) %*% t(xi_i - xi[, i])
      # tau_ij yj xi_i^T +   tau_ij yj (yj - A xi)^T gamma
      A1 <-  colSums(tY) %*% t(xi_i) + t(Y) %*% tY_Axi_i %*% gamma
      #
      A2 <-  (omega[,, i] + xi[, i] %*% t(xi[, i])) * ti
      
      AA[,,i] <- try(A1 %*% chol.inv(A2))
      
      #Di <- Di + diag(t(Y) %*% tY)
      # Di <- Di + colSums(Y * tY)
      
      beta<- ginv(AA[,,i] %*% omega[,,i] %*% t(AA[,,i]) + D[,,i])%*%D[,,i]
      
      D1 =  diag(colSums(tau[, i]%*%t(diag(D[,,i]%*%(diag(p)-beta))))) / ti
      D2 = diag(diag(t(beta)%*%t(tY_Axi_i)%*%Y_Axi_i%*%beta / ti))
      D[,,i] = D1+D2
      
      
      w1[i] <- ti / n
    }
    
    
    
    # DD <- diag(Di - rowSums((AA %*% A2) * AA )) / n
    
    A.list[[l]] <- AA
    # D.list[[l]] <- DD
    D.list[[l]] <- D
    xi.list[[l]] <- xi
    omega.list[[l]] <- omega
    w.list[[l]] <- w1
    
    
    
    
    
    
    
    
    
    
    
    
    ###############################################################
    ###################  third layer ##############################
    ###############################################################
    l <- 3
    yy <- matrix(0, numobs, r[l])
    zz <- z.list[[l-1]]
    for (i in 1 : k[l - 1]) {
      yy <- yy + matrix(zz[i,,, drop = FALSE] *
                          array(rowSums(ps.y.list[[l-1]][, i, drop = FALSE]),
                                c(1, numobs, r[l])), numobs, r[l])
    }
    
    z2 <- z2.one <- array(0, c(numobs, r[l + 1], k[l], k[l + 1]))
    zz2 <- array(0, c(numobs, r[l + 1], r[l + 1], k[l], k[l + 1]))
    z <- z.one <- array(0, c(numobs, r[l + 1], k[l]))
    zz <- array(0, c(numobs, r[l + 1], r[l + 1], k[l]))
    
    
    for (p1 in 1 : k[l]) {
      for (p2 in 1:k[l+1])  {
        # A is the inverse of \xi_{sl}
        # = (\Sigma_{sl}^(-1) + \Lambda_{sl}^T \Psi_{sl}^(-1) \Lambda_{sl} )
        A <- ginv((A.list[[l + 1]][,,p2] %*% as.matrix(array(omega.list[[ l + 1 ]],dim=c(r[l+2],r[l+2],k[l+1]))[,,p2]) %*% t(A.list[[l + 1]][,,p2])+D.list[[l+1]][,,p2])) +
          t(A.list[[l]][,,p1])%*%ginv(D.list[[l]][,,p1])%*%A.list[[l]][,,p1]
        
        # b contains part of p_{sl} (z^{l-1})
        # =  Sigma_{sl}^(-1) mu_{sl}
        b <- c(ginv((A.list[[l + 1]][,,p2] %*% as.matrix(array(omega.list[[ l + 1 ]],dim=c(r[l+2],r[l+2],k[l+1]))[,,p2]) %*% t(A.list[[l + 1]][,,p2])+D.list[[l+1]][,,p2])) %*% (c(A.list[[l + 1]][,,p2]%*%as.matrix(array(xi.list[[l+1]],dim=c(r[l+2],k[l+1]))[,p2]))))+ (t(A.list[[l]][,,p1])%*%(ginv(D.list[[l]][,,p1])) %*% (t(yy)))
        
        
        chsi  <- ginv(A)
        
        roy  <- chsi %*% b
        roy.quadro  <- array(apply(roy, 2, function(x) x %*%
                                     t(x)), c(r[l + 1], r[l + 1], numobs))
        zz2[,,, p1, p2]  <- aperm(array(chsi, c(r[l + 1], r[l + 1], numobs)) +
                                    roy.quadro, c(3, 1, 2))
        z2.one[,, p1, p2]  <- rmvnorm(numobs, rep(0, r[l + 1]), chsi) + t(roy)
        z2[,, p1, p2]  <- t(roy)
      }
    }
    
    for (i in 1 : k[l + 1]) {
      prob  <- ps.y.list[[l + 1]][, i, drop = FALSE]
      z  <- z + array(z2[,,, i, drop = FALSE] *
                        array(rowSums(prob), c(numobs, r[l + 1], k[l], 1)),
                      c(numobs, r[l + 1], k[l]))
      z.one  <- z.one + array(z2.one[,,, i, drop = FALSE] *
                                array(rowSums(prob), c(numobs, r[l + 1], k[l], 1)),
                              c(numobs, r[l + 1], k[l]))
      zz <- zz + array(zz2[,,,, i, drop = FALSE] * array(rowSums(prob),
                                                         c(numobs, r[l + 1], r[l + 1], k[l], 1)),
                       c(numobs, r[l + 1], r[l + 1], k[l]))
    }
    
    z.list[[l]] <- aperm(z.one, c(3, 1, 2))
    
    
    
    
    
    
    
    
    
    
    
    
    ### M-step get A
    ### M-step get D
    ### M-step get xi
    Y<- yy
    
    # load("/home/davb/Documents/Samsung_T5/2019_11_21_Deep_Statistical_Models/Deep_Gaussian_Mixture_Model/scripts/deepgmm-branch_3_0_0/r_backup/backup_v1.RData")
    p <- ncol(Y)
    if (is.null(p))
      p <- 1
    n <- nrow(Y)
    
    
    A1 <- array(0, c(p, r[l+1]))
    A2 <- array(0, c(r[l+1], r[l+1]))
    # Di <- array(0,c(p,k[l]))
    
    tau<-(ps.y.list[[l]])
    
    omega<-array(omega.list[[l]],dim=c(r[l+1],r[l+1],k[l]))
    xi<-array(xi.list[[l]],dim=c(r[l+1],k[l]))
    AA<-array(A.list[[l]],dim=c(r[l],r[l+1],k[l]))
    D<-D.list[[l]]
    
    w1<-rep(0,k[l])
    
    for (i in 1:k[l]) {
      
      inv_Di <- diag(1 / diag(D[,,i]))
      inv_Di_A <- sweep(AA[,,i], 1, diag(1/D[,,i]), '*')
      # gamma = (A Omega_i A^T + D) ^{-1} A Omega_i
      # using Woodbury Identity
      gamma <- (inv_Di - inv_Di_A %*% 
                  chol.inv(chol.inv(omega[,, i]) +  t(AA[,,i]) %*% inv_Di_A) %*%
                  t(inv_Di_A)) %*% AA[,,i] %*% omega[,, i]
      
      ti <- sum(tau[, i])
      xi_i <- xi[, i, drop = FALSE]
      
      # tau_ij * y_j  
      tY <- sweep(Y, 1, tau[, i], '*')
      
      # y_j - A xi_i
      Y_Axi_i <- sweep(Y, 2, AA[,,i] %*% xi_i, '-')
      
      # tau_ij * (y_j - A xi_i)
      tY_Axi_i <- sweep(Y_Axi_i, 1, tau[, i], '*')
      
      # xi_i = xi_i + gamma^T \sum{ tau_ij (y_j - A xi_i) } / sum {tau_ij}
      xi[, i] <- xi_i + t(gamma) %*% as.matrix(colSums(tY_Axi_i)) / ti
      
      # 
      omega[,, i] <- (diag(r[l+1]) - t(gamma) %*% AA[,,i]) %*% omega[,, i] +
        t(gamma) %*% t(Y_Axi_i) %*% tY_Axi_i %*% gamma / ti -
        (xi_i - xi[,i]) %*% t(xi_i - xi[, i])
      # tau_ij yj xi_i^T +   tau_ij yj (yj - A xi)^T gamma
      A1 <-  colSums(tY) %*% t(xi_i) + t(Y) %*% tY_Axi_i %*% gamma
      #
      A2 <-  (omega[,, i] + xi[, i] %*% t(xi[, i])) * ti
      
      AA[,,i] <- try(A1 %*% chol.inv(A2))
      
      #Di <- Di + diag(t(Y) %*% tY)
      # Di <- Di + colSums(Y * tY)
      
      beta<- ginv(AA[,,i] %*% omega[,,i] %*% t(AA[,,i]) + D[,,i])%*%D[,,i]
      
      D1 =  diag(colSums(tau[, i]%*%t(diag(D[,,i]%*%(diag(p)-beta))))) / ti
      D2 = diag(diag(t(beta)%*%t(tY_Axi_i)%*%Y_Axi_i%*%beta / ti))
      D[,,i] = D1+D2
      
      
      w1[i] <- ti / n
    }
    
    
    
    # DD <- diag(Di - rowSums((AA %*% A2) * AA )) / n
    
    A.list[[l]] <- AA
    # D.list[[l]] <- DD
    D.list[[l]] <- D
    xi.list[[l]] <- xi
    omega.list[[l]] <- omega
    w.list[[l]] <- w1
    
    
    
    
    
    
    
    
    
    ###############################################################
    ###################  fourth layer ##############################
    ###############################################################
    l <- 4
    yy <- matrix(0, numobs, r[l])
    zz <- z.list[[l-1]]
    for (i in 1 : k[l - 1]) {
      yy <- yy + matrix(zz[i,,, drop = FALSE] *
                          array(rowSums(ps.y.list[[l-1]][, i, drop = FALSE]),
                                c(1, numobs, r[l])), numobs, r[l])
    }
    
    z <- z.one<-array(0, c(numobs, r[l + 1], k[l]))
    zz <- array(0,c(numobs, r[l + 1], r[l + 1], k[l]))
    
    for (p1 in 1:k[l]) {
      
      # A is the inverse of \xi_{sl}
      # = (\Sigma_{sl}^(-1) + \Lambda_{sl}^T \Psi_{sl}^(-1) \Lambda_{sl} )
      A <- ginv(array(omega.list[[ l]],dim=c(r[l+1],r[l+1],k[l]))[,,p1]) +
        t(A.list[[l]][,,p1])%*%ginv(D.list[[l]][,,p1])%*%A.list[[l]][,,p1]
      
      # b contains part of p_{sl} (z^{l-1})
      # =  Sigma_{sl}^(-1) mu_{sl}
      b <- c((ginv(array(omega.list[[ l]],dim=c(r[l+1],r[l+1],k[l]))[,,p1]))%*% (c(array(xi.list[[l]],dim=c(r[l+1],k[l]))[,p1])))+ (t(A.list[[l]][,,p1])%*%(ginv(D.list[[l]][,,p1])) %*% (t(yy)))
      
      chsi <- ginv(A)
      if (!isSymmetric(chsi)) {
        chsi <- make.positive.definite(makeSymm(chsi))
      }
      
      roy <- chsi %*% b
      roy.quadro <- array(apply(roy, 2, function(x) x %*%
                                  t(x)), c(r[l + 1], r[l + 1], numobs))
      zz[,,, p1]  <- aperm(array(chsi, c(r[l + 1], r[l + 1], numobs)) +
                             roy.quadro, c(3, 1, 2))
      z.one[,, p1] <- t(roy) + rmvnorm(numobs, rep(0, r[l + 1]), chsi)
      z[,, p1]  <- t(roy)
    }
    
    z.list[[l]] <- aperm(z.one, c(3, 1, 2))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ### M-step get A
    ### M-step get D
    ### M-step get xi
    Y<- yy
    
    # load("/home/davb/Documents/Samsung_T5/2019_11_21_Deep_Statistical_Models/Deep_Gaussian_Mixture_Model/scripts/deepgmm-branch_3_0_0/r_backup/backup_v1.RData")
    p <- ncol(Y)
    if (is.null(p))
      p <- 1
    n <- nrow(Y)
    
    
    A1 <- array(0, c(p, r[l+1]))
    A2 <- array(0, c(r[l+1], r[l+1]))
    # Di <- array(0,c(p,k[l]))
    
    tau<-(ps.y.list[[l]])
    
    omega<-array(omega.list[[l]],dim=c(r[l+1],r[l+1],k[l]))
    xi<-array(xi.list[[l]],dim=c(r[l+1],k[l]))
    AA<-array(A.list[[l]],dim=c(r[l],r[l+1],k[l]))
    D<-D.list[[l]]
    
    w1<-rep(0,k[l])
    
    for (i in 1:k[l]) {
      
      inv_Di <- diag(1 / diag(D[,,i]))
      inv_Di_A <- sweep(AA[,,i], 1, diag(1/D[,,i]), '*')
      # gamma = (A Omega_i A^T + D) ^{-1} A Omega_i
      # using Woodbury Identity
      gamma <- (inv_Di - inv_Di_A %*% 
                  chol.inv(makeSymm(chol.inv(makeSymm(omega[,, i]))) +  t(AA[,,i]) %*% inv_Di_A) %*%
                  t(inv_Di_A)) %*% AA[,,i] %*% omega[,, i]
      
      ti <- sum(tau[, i])
      xi_i <- xi[, i, drop = FALSE]
      
      # tau_ij * y_j  
      tY <- sweep(Y, 1, tau[, i], '*')
      
      # y_j - A xi_i
      Y_Axi_i <- sweep(Y, 2, AA[,,i] %*% xi_i, '-')
      
      # tau_ij * (y_j - A xi_i)
      tY_Axi_i <- sweep(Y_Axi_i, 1, tau[, i], '*')
      
      # xi_i = xi_i + gamma^T \sum{ tau_ij (y_j - A xi_i) } / sum {tau_ij}
      xi[, i] <- xi_i + t(gamma) %*% as.matrix(colSums(tY_Axi_i)) / ti
      
      # 
      omega[,, i] <- (diag(r[l+1]) - t(gamma) %*% AA[,,i]) %*% omega[,, i] +
        t(gamma) %*% t(Y_Axi_i) %*% tY_Axi_i %*% gamma / ti -
        (xi_i - xi[,i]) %*% t(xi_i - xi[, i])
      # tau_ij yj xi_i^T +   tau_ij yj (yj - A xi)^T gamma
      A1 <-  colSums(tY) %*% t(xi_i) + t(Y) %*% tY_Axi_i %*% gamma
      #
      A2 <-  (omega[,, i] + xi[, i] %*% t(xi[, i])) * ti
      
      AA[,,i] <- try(A1 %*% chol.inv(A2))
      
      #Di <- Di + diag(t(Y) %*% tY)
      # Di <- Di + colSums(Y * tY)
      
      beta<- ginv(AA[,,i] %*% omega[,,i] %*% t(AA[,,i]) + D[,,i])%*%D[,,i]
      
      D1 =  diag(colSums(tau[, i]%*%t(diag(D[,,i]%*%(diag(p)-beta))))) / ti
      D2 = diag(diag(t(beta)%*%t(tY_Axi_i)%*%Y_Axi_i%*%beta / ti))
      D[,,i] = D1+D2
      
      
      w1[i] <- ti / n
    }
    
    
    
    # DD <- diag(Di - rowSums((AA %*% A2) * AA )) / n
    
    A.list[[l]] <- AA
    # D.list[[l]] <- DD
    D.list[[l]] <- D
    xi.list[[l]] <- xi
    omega.list[[l]] <- omega
    w.list[[l]] <- w1
    
    
    
    
    
    
    
    
    out <- deepgmm_hmcfl_compute.lik(y = y, numobs = numobs, k = k,r=r, xi.list = xi.list,A.list =  A.list, D.list = D.list,omega.list =  omega.list, w.list = w.list)
    py <- out$py
    ps.y <- out$ps.y
    ps.y.list <- out$ps.y.list
    k.comb <- out$k.comb
    s <- out$s
    
    lik <- sum(log(py))
    likelihood <- c(likelihood, lik)
    
    if (hh < HHHHH) {
      ratio <- 2 * eps
    }
    if (hh > HHHHH) {
      ratio <- (ma(likelihood, 5)[hh + 1] - ma(likelihood, 5)[hh]) /
        abs(ma(likelihood, 5)[hh])
    }
    # print(c(hh,lik,adjustedRandIndex(labels,s[,column])))
    
  }
  
  
  
  h <- 0
  for (j in 1 : layers) {
    if (j==1){
      h <- h + (k[j] - 1) + r[j] * k[j] + (r[j]*r[j+1]) - r[j+1]^2 + (r[j] + k[j])*r[j+1] + (r[j] * k[j] * (r[j] + 1) / 2) 
    } else {
      h <- h + (k[j] - 1) + r[j] * k[j] + k[j]*(r[j]*r[j+1] - r[j+1]^2 )+ (r[j] + k[j])*r[j+1]  + (r[j] * k[j] * (r[j] + 1) / 2) 
    }
  }
  
  
  
  

  bic <- -2 * lik + h * log(numobs)
  aic <- -2 * lik + 2 * h
  EN <- entr(ps.y.list[[4]])
  clc <- -2 * lik + 2 * EN
  icl.bic <- -2 * lik + 2 * EN + h * log(numobs)
  # print(c(hh,bic,adjustedRandIndex(labels,s[,column])))
  
  out <- list(A = list(A = A.list), w = list(w = w.list), xi = list(xi = xi.list), omega=list(omega = omega.list),
              D = list(D = D.list), likelihood = likelihood, bic = bic, aic = aic, clc = clc,
              s = s, icl_bic = icl.bic, h = h, ps.y = ps.y)
  return(out)
}
