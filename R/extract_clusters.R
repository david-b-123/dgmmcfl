extract_clusters <- function(model.fit,method="dgmm.cfl"){
  
  if (method=="dgmm"){
    dgmm_clusters<- model.fit
    clusters<-dgmm_clusters$s[,1]
  }
  if (method=="dgmm.cfl"){    
    dgmm_clusters<- model.fit
    layers<-length(dgmm_clusters$k)
    clusters<-dgmm_clusters$s[,layers]
  }
  
  if (method=="dgmm.icfl"){
    dgmm_clusters<- model.fit
    layers<-length(dgmm_clusters$k)
    clusters<-dgmm_clusters$s[,2]
  }
  
  if (method=="dgmm.hmcfl"){
    dgmm_clusters<- model.fit
    layers<-length(dgmm_clusters$k)
    clusters<-dgmm_clusters$s[,layers]
  }
  
  return(clusters)
  
}