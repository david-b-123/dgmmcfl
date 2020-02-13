# install.packages("EMMIXmfa")
# install.packages("propr")
# install.packages("corpcor")
# install.packages("mvtnorm")

library(EMMIXmfa)
library(propr)
library(corpcor)
library(mvtnorm)

# path_to : path to the "R" directory in github download
path_to<-"/home/davb/Documents/Samsung_T5/2019_11_21_Deep_Statistical_Models/test_code/dgmmcfl-master/R/"










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

phi_results<-eigen(phis)$vectors[,c(1:5)]
pea_results<-eigen(pearson)$vectors[,c(1:5)]
spe_results<-eigen(spearman)$vectors[,c(1:5)]
euc_results<-eigen(euclidean)$vectors[,c(1:5)]

# DGMM methods requires a dataframe
y<-data.frame(phi_results,pea_results,spe_results,euc_results)

# These methods take a long time to grid search the parameter space to find the optimal model, based on BIC.
# If taking too long, set complete==F, this grid searches naively over a SMALLER model parameter space

# Also, I can parrallelise this, if needed, to work on multiple cores - this would speed it up by a certain factor depending on the number of cores used.
# I plan on adding in parallelisation later.

# Run DGMM - THIS MODEL CANNOT BE VISUALISED
# Possible Layers = {1,2,3}
dgmm_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 1,it = 1,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm',scale = T)

# RUN DGMM-CFL - ALL LAYERS CAN BE VISUALISED 
# Possible Layers = {2,3,4}
dgmm.cfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 1,it = 1,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.cfl',scale = T)

# RUN DGMM-ICFL - ONLY FIRST LAYER CAN BE VISUALISED
# Possible Layers = {2,3,4}
dgmm.icfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.icfl',scale = T)

# RUN DGMM-HMCFL - ONLY FIRST LAYER CAN BE VISUALISED
# Possible Layers = {2,3,4}
dgmm.hmcfl_model<-model_selection(y,layers = 2,g = length(unique(labels)),seeds = 10,it = 150,eps = 10E-4,criterion = "BIC",complete = F,method = 'dgmm.hmcfl',scale = T)


# These scores iterate over all dimensions of the latent factor scores to find the best BIC fit.
# Thus, it is possible to receive a plotting dimension of 1, which is of course unsuitable

# To plot the last layer onto two dimensions, set the following:
# r = c(R1, 2) for 2 layers
# r = c(R1, R2, 2) for 3 layers
# r = c(R1, R2, R3, 2) for 4 layers

bic.best <- Inf
aic.best <- Inf

seed_setting=50
K1=3
R1=6

# Run to completion
for (k1 in 1:K1){ # iterate from 1 to 3 mixtures
  for (r1 in 3:R1){ # iterate from 3 to 7 - cannot work on a dataset with 3 or less dimensions
    for (seed in 1:seed_setting){ # run EM algorithm - it does not guarantee global optima, and so the best BIC will be chosen
      
      print(paste("Seed: ",seed,"/",seed_setting, "   | | |   ","Mixture models searched: ",k1,"/",K1,"   | | |   ","Latent dimensions searched: ",r1-2,"/",R1-2,  sep=""))
      set.seed(seed)
      out <- try(deepgmm_mcfa(y,layers = 2,k = c(k1,length(unique(labels))),r = c(r1,2),it = 150,eps = 10E-5,method = "dgmm.cfl",scale = T),silent=F)
      
      if (!is.character(out)) {
        if (out$bic < bic.best) {
          out.best <- out
          bic.best <- out$bic
          
          par(mfcol=c(1,2))
          plot(data.frame(out.best$factor_scores[[2]]),col=as.factor(labels),xlab="Score 1", ylab="Score 2",main="True Labels")
          plot(data.frame(out.best$factor_scores[[2]]),col=as.factor(out.best$clusters),xlab="Score 1", ylab="Score 2",main="DGMM labels")

        }
        
      }
    }
  }
}
