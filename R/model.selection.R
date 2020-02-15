model_selection <- function (y, layers, g, seeds = 3, it = 50, eps = 0.001,
                             init = "kmeans", init_est = "factanal", criterion = "BIC", complete=F, method='dgmm', scale=F) {
  print("Remember to scale the data. Can be done through model_selection(..,scale=T)")
  
  bic.best <- Inf
  aic.best <- Inf
  p <- dim(y)[2]
  if (complete==F){
    pp <- 5
    ppp <- 4
    pppp <- 3
  } 
  # The following takes a long time for large p (high number of features, p >> 10)
  # ONLY RECOMMENDED FOR HIGH PERFORMANCE COMPUTERS
  if (complete==T){
    pp <- p-1
    ppp <- pp-1
    pppp <- ppp-1
  }
  
  if (layers == 1) {
    
    # There is no 1 layer dgmm.cfl, dgmm.hmcfl, dgmm.icfl - see EMMIXmfa::mcfa() for a shallow 1-layer MCFA instead
    method='dgmm'
    
    r <- c(1 : pp)
    bic <- array(NA, c(seeds, pp))
    bic.best <- Inf
    aic <- array(NA, c(seeds, pp))
    aic.best <- Inf
    
    for (i in 1 : seeds) for (rr in 1 : pp) {
      
      set.seed(i)
      out <- try(deepgmm_mcfa(y, 1, k = g, rr, it = it, eps = eps, init = init,
                              init_est = init_est, method = method, scale = scale))
      
      if (!is.character(out)) {
        if (criterion=="BIC") if (out$bic < bic.best) {
          out.best <- out
          bic.best <- out$bic
        }
        if (criterion=="AIC")
          if (out$aic<aic.best) {
            out.best <- out
            aic.best <- out$aic
          }
        
        bic[i, rr] <- out$bic
        aic[i, rr] <- out$aic}
    }
    
    if (criterion == "BIC")
      index <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    if (criterion == "AIC")
      index <- which(aic == min(aic, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    message(paste("Best Fit with init =", init, "and init_est =", init_est, "\n"))
    cat(paste0("Seed=", index[1], 
               ", r = ", index[2],
               ", BIC: ", round(out.best$bic, 2), 
               ", AIC: ", round(out.best$aic, 2)))
  }
  
  if (layers == 2) {
    
    r <- matrix(c(dim(y)[2]-1,dim(y)[2]-2),nrow=1)
    
    k <- 
      do.call('rbind',lapply(c(1:4),function(A){
        do.call('rbind',lapply(c(1:4),function(B){
          internal_k<-c(A,B)
          
          if (method=='dgmm'){
            internal_k[1] <- g
          }
          if (method=='dgmm.cfl'){
            internal_k[length(internal_k)] <- g
          }              
          if (method=='dgmm.hmcfl'){
            internal_k[length(internal_k)] <- g
          }              
          if (method=='dgmm.icfl'){
            internal_k[2] <- g
          }
          
          
          return(internal_k)
          
        }))
      }))
    
    k<-unique(k)
    
    bic <- array(NA, c(seeds, nrow(k), nrow(r)))
    bic.best <- Inf
    aic <- array(NA, c(seeds, nrow(k), nrow(r)))
    aic.best <- Inf
    
    for (i in 1 : seeds)
      for (kk in 1 : nrow(k))
        for (rr in 1 : nrow(r)) {
          set.seed(i)
          
          print(paste("Seed: ",i,"/",seeds, "   | | |   ","Mixture models searched: ",kk,"/",nrow(k),"   | | |   ","Latent dimensions searched: ",rr,"/",nrow(r),  sep=""))
          
          out <- try(deepgmm_mcfa(y, 2, k[kk, ], r[rr, ], it = it, eps = eps,
                                  init = init, init_est = init_est, method = method, scale = scale ))
          if (!is.character(out)) {
            if (criterion == "BIC")
              if (out$bic < bic.best) {
                out.best <- out
                bic.best <- out$bic
              }
            if (criterion == "AIC")
              if (out$aic < aic.best) {
                out.best <- out
                aic.best <- out$aic
              }
            
            bic[i, kk, rr] <- out$bic
            aic[i, kk, rr] <- out$aic
          }
        }
    
    if (criterion == "BIC")
      index <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    if (criterion == "AIC")
      index <- which(aic == min(aic, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    message(paste("Best Fit with init =", init, "and init_est =", init_est, "\n"))
    cat(paste0("Seed = ", index[1], 
               ", k =", paste(k[index[2],], collapse=" "),
               ", r = ", paste(r[index[3], ], collapse = " "),
               ", BIC:", round(out.best$bic, 2), 
               ", AIC:", round(out.best$aic, 2)))
  }
  
  if (layers == 3) {
    
    r <- as.matrix(expand.grid(1 : pp, 1 : ppp, 1 : ppp))
    r <- r[((r[, 1]) > (r[, 2])) & ((r[, 2]) > (r[, 3])), ]
    
    k <- 
      do.call('rbind',lapply(c(1:3),function(A){
        do.call('rbind',lapply(c(1:3),function(B){
          do.call('rbind',lapply(c(1:3),function(C){
            internal_k<-c(A,B,C)
            
            if (method=='dgmm'){
              internal_k[1] <- g
            }
            if (method=='dgmm.cfl'){
              internal_k[length(internal_k)] <- g
            }              
            if (method=='dgmm.hmcfl'){
              internal_k[length(internal_k)] <- g
            }              
            if (method=='dgmm.icfl'){
              internal_k[2] <- g
            }
            
            
            return(internal_k)
            
          }))
        }))
      }))
    
    k<-unique(k)
    
    bic <- array(NA, c(seeds, nrow(k), nrow(r)))
    bic.best <- Inf
    aic <- array(NA, c(seeds, nrow(k), nrow(r)))
    aic.best <- Inf
    
    for (i in 1 : seeds)
      for (kk in 1 : nrow(k))
        for (rr in 1 : nrow(r)) {
          
          print(paste("Seed: ",i,"/",seeds, "   | | |   ","Mixture models searched: ",kk,"/",nrow(k),"   | | |   ","Latent dimensions searched: ",rr,"/",nrow(r),  sep=""))
          
          set.seed(i)
          out <- try(deepgmm_mcfa(y, 3, k[kk, ], r[rr, ], it = it, eps = eps,
                                  init = init, init_est = init_est, method = method, scale = scale ))
          if (!is.character(out)) {
            if (criterion=="BIC")
              if (out$bic < bic.best) {
                out.best <- out
                bic.best <- out$bic
              }
            
            if (criterion=="AIC")
              if (out$aic<aic.best) {
                out.best <- out
                aic.best <- out$aic
              }
            
            bic[i,kk,rr] <- out$bic
            aic[i,kk,rr] <- out$aic
          }
        }
    
    if (criterion == "BIC")
      index <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    if (criterion == "AIC")
      index <- which(aic == min(aic, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    message(paste("Best Fit with init =", init, "and init_est =", init_est, "\n"))
    
    cat(paste0("Seed = ", index[1], 
               ", k = ", paste(k[index[2], ], collapse=" "),
               ", r = ", paste(r[index[3], ], collapse=" "),
               ", BIC: ", round(out.best$bic, 2), 
               ", AIC: ", round(out.best$aic, 2)))
  }
  if (layers == 4){
    
    r <- as.matrix(expand.grid(1 : pp, 1 : ppp, 1 : ppp, 1:pppp))
    r <- r[((r[, 1]) > (r[, 2])) & ((r[, 2]) > (r[, 3])) & ((r[, 3]) > (r[, 4])), ]
    
    k <- 
      do.call('rbind',lapply(c(1:3),function(A){
        do.call('rbind',lapply(c(1:3),function(B){
          do.call('rbind',lapply(c(1:3),function(C){
            do.call('rbind',lapply(c(1:3),function(D){
              internal_k<-c(A,B,C,D)
              
              if (method=='dgmm'){
                internal_k[1] <- g
              }
              if (method=='dgmm.cfl'){
                internal_k[length(internal_k)] <- g
              }              
              if (method=='dgmm.hmcfl'){
                internal_k[length(internal_k)] <- g
              }              
              if (method=='dgmm.icfl'){
                internal_k[2] <- g
              }
              
              
              return(internal_k)
              
            }))
          }))
        }))
      }))
    
    k<-unique(k)
    
    bic <- array(NA, c(seeds, nrow(k), nrow(r)))
    bic.best <- Inf
    aic <- array(NA, c(seeds, nrow(k), nrow(r)))
    aic.best <- Inf
    
    for (i in 1 : seeds)
      for (kk in 1 : nrow(k))
        for (rr in 1 : nrow(r)) {
          
          print(paste("Seed: ",i,"/",seeds, "   | | |   ","Mixture models searched: ",kk,"/",nrow(k),"   | | |   ","Latent dimensions searched: ",rr,"/",nrow(r),  sep=""))

          set.seed(i)
          out <- try(deepgmm_mcfa(y, 4, k[kk, ], r[rr, ], it = it, eps = eps,
                                  init = init, init_est = init_est, method = method, scale = scale ))
          if (!is.character(out)) {
            if (criterion=="BIC")
              if (out$bic < bic.best) {
                out.best <- out
                bic.best <- out$bic
              }
            
            if (criterion=="AIC")
              if (out$aic<aic.best) {
                out.best <- out
                aic.best <- out$aic
              }
            
            bic[i,kk,rr] <- out$bic
            aic[i,kk,rr] <- out$aic
          }
        }
    
    if (criterion == "BIC")
      index <- which(bic == min(bic, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    if (criterion == "AIC")
      index <- which(aic == min(aic, na.rm = TRUE), arr.ind = TRUE)[1, ]
    
    message(paste("Best Fit with init =", init, "and init_est =", init_est, "\n"))
    
    cat(paste0("Seed = ", index[1], 
               ", k = ", paste(k[index[2], ], collapse=" "),
               ", r = ", paste(r[index[3], ], collapse=" "),
               ", BIC: ", round(out.best$bic, 2), 
               ", AIC: ", round(out.best$aic, 2)))
    
  }
  
  
  
  out <- list(fit = out.best, bic = bic, aic = aic)
  invisible(out)
}
