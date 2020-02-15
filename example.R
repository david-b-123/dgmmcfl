# install.packages("EMMIXmfa")
# install.packages("propr")
# install.packages("corpcor")
# install.packages("mvtnorm")
# install.packages("tsne")
# install.packages("umap")

library(EMMIXmfa)
library(propr)
library(corpcor)
library(mvtnorm)
library(tsne)
library(umap)

# path_to : path to the "R" directory in github download
path_to<-"/home/davb/Documents/Samsung_T5/2019_11_21_Deep_Statistical_Models/test_code_all/dgmmcfl-master_v2/R/"










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

phi_results<-eigen(phis)$vectors[,c(1:10)]
pea_results<-eigen(pearson)$vectors[,c(1:10)]
spe_results<-eigen(spearman)$vectors[,c(1:10)]

# DGMM methods requires a dataframe
y<-data.frame(phi_results,pea_results,spe_results)
# These methods take a long time to grid search the parameter space to find the optimal model, based on BIC.
# If taking too long, set complete==F, this grid searches naively over a SMALLER model parameter space

# Also, I can parrallelise this, if needed, to work on multiple cores - this would speed it up by a certain factor depending on the number of cores used.
# I plan on adding in parallelisation later.

# Run DGMM - THIS MODEL CANNOT BE VISUALISED
# Possible Layers = {1,2,3}
dgmm_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm',scale = T)
print(adjustedRandIndex(dgmm_model$fit$s[,1],labels))

# RUN DGMM-CFL - ALL LAYERS CAN BE VISUALISED 
# Possible Layers = {2,3,4}
dgmm.cfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.cfl',scale = T)
print(adjustedRandIndex(dgmm.cfl_model$fit$clusters,labels))

# RUN DGMM-ICFL - ONLY FIRST LAYER CAN BE VISUALISED
# Possible Layers = {2,3,4}
dgmm.icfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.icfl',scale = T)
print(adjustedRandIndex(dgmm.icfl_model$fit$clusters,labels))

# RUN DGMM-HMCFL - ONLY FIRST LAYER CAN BE VISUALISED
# Possible Layers = {2,3,4}
dgmm.hmcfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.hmcfl',scale = T)
print(adjustedRandIndex(dgmm.hmcfl_model$fit$clusters,labels))

