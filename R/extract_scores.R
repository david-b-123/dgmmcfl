extract_scores <- function(data, model.fit,method="dgmm.cfl",plot=F){
  
  if (method=="dgmm"){
    print("DGMM can not be plotted")
  }
  if (method=="dgmm.cfl"){
    
    # PLOT DGMM.CFL
    dgmm_plotting<- model.fit
    
    k<-dgmm_plotting$k
    r<-dgmm_plotting$r
    layers<-length(k)
    
    A1<-dgmm_plotting$A$A[[1]]
    D1<-dgmm_plotting$D$D[[1]]
    xi1<-dgmm_plotting$xi$xi[[1]]
    omega1<-dgmm_plotting$omega$omega[[1]]
    w1<-dgmm_plotting$w$w[[1]]
    
    u_score=as.vector(c(),mode = "list")
    u_internal<-0
    for (i in 1:k[1]){
      gamma1=ginv(A1%*%omega1[,,i]%*%t(A1)+D1[,,i])%*%A1%*%omega1[,,i]
      u_internal<-u_internal+(w1[i]*(c(xi1[,i])+t(gamma1)%*%t(as.matrix(data)-c(A1%*%c(xi1[,i])))))
    }
    u_score[[1]]<-t(u_internal)
    
    for (l in 2:c(layers)){
      A<-dgmm_plotting$A$A[[l]]
      D<-dgmm_plotting$D$D[[l]]
      xi<-dgmm_plotting$xi$xi[[l]]
      omega<-dgmm_plotting$omega$omega[[l]]
      w<-dgmm_plotting$w$w[[l]]
      
      
      u_internal=0
      for (i in 1:k[l]){
        gamma=ginv(A%*%omega[,,i]%*%t(A)+D[,,i])%*%A%*%omega[,,i]
        u_internal<-u_internal+(w[i]*(c(xi[,i])+t(gamma)%*%(t(u_score[[l-1]])-c(A%*%c(xi[,i])))))
      }
      
      u_score[[l]]<-t(u_internal)
    }
    
  }
  
  if (method=="dgmm.icfl"){
    
    # PLOT DGMM.ICFL
    dgmm_plotting<- model.fit
    
    k<-dgmm_plotting$k
    r<-dgmm_plotting$r
    layers<-length(k)
    
    A1<-dgmm_plotting$A$A[[1]]
    D1<-dgmm_plotting$D$D[[1]]
    xi1<-dgmm_plotting$xi$xi[[1]]
    omega1<-dgmm_plotting$omega$omega[[1]]
    w1<-dgmm_plotting$w$w[[1]]
    
    u_score=as.vector(c(),mode = "list")
    u_internal<-0
    for (i in 1:k[1]){
      gamma1=ginv(A1%*%omega1[,,i]%*%t(A1)+D1[,,i])%*%A1%*%omega1[,,i]
      u_internal<-u_internal+(w1[i]*(c(xi1[,i])+t(gamma1)%*%t(as.matrix(data)-c(A1%*%c(xi1[,i])))))
    }
    u_score[[1]]<-t(u_internal)
    
  }
  
  
  if (method=="dgmm.hmcfl"){
    # PLOT DGMM.HMCFL
    dgmm_plotting<- model.fit
    
    k<-dgmm_plotting$k
    r<-dgmm_plotting$r
    
    A1<-dgmm_plotting$A$A[[1]][,,1]
    D1<-dgmm_plotting$D$D[[1]]
    xi1<-dgmm_plotting$xi$xi[[1]]
    omega1<-dgmm_plotting$omega$omega[[1]]
    w1<-dgmm_plotting$w$w[[1]]
    
    u_score=as.vector(c(),mode = "list")
    u_internal<-0
    for (i in 1:k[1]){
      gamma1=ginv(A1%*%omega1[,,i]%*%t(A1)+D1[,,i])%*%A1%*%omega1[,,i]
      u_internal<-u_internal+(w1[i]*(c(xi1[,i])+t(gamma1)%*%t(as.matrix(data)-c(A1%*%c(xi1[,i])))))
    }
    u_score[[1]]<-t(u_internal)
    
  }
  
  
  if (plot==T){
    for (i in 1:length(u_score)){
      # This will plot all the dimensions
      plot(data.frame(u_score[[l]]),col=as.factor(labels))
      
      # UMAP
      # Further visualisation of factor scores
      plot(umap::umap(u_score[[l]])$layout,col=as.factor(labels)) # DGMM-HMCFL
      # Compared to original data
      plot(umap::umap(data)$layout,col=as.factor(labels)) # tSNE
      
      
      
      # tSNE
      # Further visualisation of factor scores 
      plot(tsne::tsne(u_score[[l]],k=2),col=as.factor(labels)) # DGMM-HMCFL
      # Compared to original data
      plot(tsne::tsne(data,k=2),col=as.factor(labels)) # tSNE
    }
  }
  
  return(u_score)
  
}