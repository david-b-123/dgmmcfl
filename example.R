# install.packages("EMMIXmfa")
# install.packages("propr")
# install.packages("corpcor")
# install.packages("mvtnorm")

library(EMMIXmfa)
library(propr)
library(corpcor)
library(mvtnorm)

# path_to : path to the "R" directory in github download
path_to<-"~~~~~~~~~~~~~~~~~~~~~/dgmmcfl-master/R/"










files<-list.files(path_to,pattern = ".R")
lapply(files,function(X){
  source(paste(path_to,X,sep=""))
})





# Toy dataset (single cell rna-seq by Goolam et al)
# Goolam, M. et al. Heterogeneity in Oct4 and Sox2 Targets Biases Cell Fate in 4-Cell Mouse Embryos. Cell 165, 61–74 (2016)

name="goolam"
data_path=paste(path_to,"/../toy_data/goolam.rds",sep="")
dataset_rds<-readRDS(data_path)

# Extract true labels
labels<-dataset_rds@colData@listData$cell_type1

# Extract dataset
dataset_aim<-dataset_rds@assays@.xData$data$counts


# The following filtering steps takes the same process as: SC3, a state of the art R package for single-cell data

# Filter low counts 
raw_dataset_aim<-dataset_aim[apply(dataset_aim,1,function(X){
  c(sum(X!=0)>c(length(X)*0.06))
}),]

# Transform and calculate Princpal Components
dataset_aim<-log2((1+raw_dataset_aim))
phis<- dismay::dismay(data.frame(raw_dataset_aim),metric = "phi_s")
pearson<- cor(data.frame(dataset_aim))
spearman<- cor(data.frame(dataset_aim),method = "spearman")
euclidean<- dist(t(data.frame(dataset_aim)),diag = T,upper = T)

phi_results<-eigen(phis)$vectors[,c(1:10)]
pea_results<-eigen(pearson)$vectors[,c(1:10)]
spe_results<-eigen(spearman)$vectors[,c(1:10)]
euc_results<-eigen(euclidean)$vectors[,c(1:10)]

# DGMM methods requires a dataframe
y<-data.frame(phi_results,pea_results,spe_results,euc_results)

# These methods take a long time to grid search the parameter space to find the optimal model, based on BIC.
# If taking too long, set complete==F, this grid searches naively over a SMALLER model parameter space

# Also, I can parrallelise this, if needed, to work on multiple cores - this would speed it up by a certain factor depending on the number of cores used.
# I plan on adding in parallelisation later.

# Run DGMM - THIS MODEL CANNOT BE VISUALISED
# Possible Layers = {1,2,3}
dgmm_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm',scale = T)

# RUN DGMM-CFL - ALL LAYERS CAN BE VISUALISED 
# Possible Layers = {2,3,4}
dgmm.cfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.cfl',scale = T)

# RUN DGMM-ICFL - ONLY FIRST LAYER CAN BE VISUALISED
# Possible Layers = {2,3,4}
dgmm.icfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.icfl',scale = T)

# RUN DGMM-HMCFL - ONLY FIRST LAYER CAN BE VISUALISED
# Possible Layers = {2,3,4}
dgmm.hmcfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.hmcfl',scale = T)


# dgmm_model cannot be visualised
# dgmm.cfl can have all layers be visualised
# dgmm.icfl can only have first layer visualised
# dgmm.hmcfl can only have first layer visualised

# PLOT DGMM.CFL
dgmm_plotting<- dgmm.cfl_model$fit

k<-dgmm_plotting$k
r<-dgmm_plotting$r

A1<-dgmm_plotting$A$A[[1]]
D1<-dgmm_plotting$D$D[[1]]
xi1<-dgmm_plotting$xi$xi[[1]]
omega1<-dgmm_plotting$omega$omega[[1]]
w1<-dgmm_plotting$w$w[[1]]

u11=0
for (i in 1:k[1]){
  gamma1=ginv(A1%*%omega1[,,i]%*%t(A1)+D1[,,i])%*%A1%*%omega1[,,i]
  u11<-u11+(w1[i]*(c(xi1[,i])+t(gamma1)%*%t(as.matrix(y)-c(A1%*%c(xi1[,i])))))
}


A2<-dgmm_plotting$A$A[[2]]
D2<-dgmm_plotting$D$D[[2]]
xi2<-dgmm_plotting$xi$xi[[2]]
omega2<-dgmm_plotting$omega$omega[[2]]
w2<-dgmm_plotting$w$w[[2]]


u21=0
for (i in 1:k[2]){
  gamma2=ginv(A2%*%omega2[,,i]%*%t(A2)+D2[,,i])%*%A2%*%omega2[,,i]
  u21<-u21+(w2[i]*(c(xi2[,i])+t(gamma2)%*%(as.matrix(u11)-c(A2%*%c(xi2[,i])))))
}

# This will plot all the dimensions
plot(data.frame(t(u11)),col=as.factor(labels))
plot(data.frame(t(u21)),col=as.factor(labels))

# UMAP
# Further visualisation of factor scores
plot(umap::umap(t(u21))$layout,col=as.factor(labels)) # DGMM-HMCFL
# Compared to original data
plot(umap::umap(y)$layout,col=as.factor(labels)) # tSNE



# tSNE
# Further visualisation of factor scores 
plot(tsne::tsne(t(u21),k=2),col=as.factor(labels)) # DGMM-HMCFL
# Compared to original data
plot(tsne::tsne(y,k=2),col=as.factor(labels)) # tSNE






# PLOT DGMM.ICFL
dgmm_plotting<- dgmm.icfl_model$fit

k<-dgmm_plotting$k
r<-dgmm_plotting$r

A1<-dgmm_plotting$A$A[[1]]
D1<-dgmm_plotting$D$D[[1]]
xi1<-dgmm_plotting$xi$xi[[1]]
omega1<-dgmm_plotting$omega$omega[[1]]
w1<-dgmm_plotting$w$w[[1]]

u11=0
for (i in 1:k[1]){
  gamma1=ginv(A1%*%omega1[,,i]%*%t(A1)+D1[,,i])%*%A1%*%omega1[,,i]
  u11<-u11+(w1[i]*(c(xi1[,i])+t(gamma1)%*%t(as.matrix(y)-c(A1%*%c(xi1[,i])))))
}

# This will plot all the dimensions
plot(data.frame(t(u11)),col=as.factor(labels))


# UMAP
# Further visualisation of factor scores
plot(umap::umap(t(u21))$layout,col=as.factor(labels)) # DGMM-HMCFL
# Compared to original data
plot(umap::umap(y)$layout,col=as.factor(labels)) # tSNE



# tSNE
# Further visualisation of factor scores 
plot(tsne::tsne(t(u21),k=2),col=as.factor(labels)) # DGMM-HMCFL
# Compared to original data
plot(tsne::tsne(y,k=2),col=as.factor(labels)) # tSNE





# PLOT DGMM.HMCFL
dgmm_plotting<- dgmm.hmcfl_model$fit

k<-dgmm_plotting$k
r<-dgmm_plotting$r

A1<-dgmm_plotting$A$A[[1]][,,1]
D1<-dgmm_plotting$D$D[[1]]
xi1<-dgmm_plotting$xi$xi[[1]]
omega1<-dgmm_plotting$omega$omega[[1]]
w1<-dgmm_plotting$w$w[[1]]

u11=0
for (i in 1:k[1]){
  gamma1=ginv(A1%*%omega1[,,i]%*%t(A1)+D1[,,i])%*%A1%*%omega1[,,i]
  u11<-u11+(w1[i]*(c(xi1[,i])+t(gamma1)%*%t(as.matrix(y)-c(A1%*%c(xi1[,i])))))
}

# This will plot all the dimensions
plot(data.frame(t(u11)),col=as.factor(labels))


# UMAP
# Further visualisation of factor scores
plot(umap::umap(t(u21))$layout,col=as.factor(labels)) # DGMM-HMCFL
# Compared to original data
plot(umap::umap(y)$layout,col=as.factor(labels)) # tSNE



# tSNE
# Further visualisation of factor scores 
plot(tsne::tsne(t(u21),k=2),col=as.factor(labels)) # DGMM-HMCFL
# Compared to original data
plot(tsne::tsne(y,k=2),col=as.factor(labels)) # tSNE



