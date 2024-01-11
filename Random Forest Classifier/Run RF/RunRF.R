# Rewrite the random forest classifier
# This one tests the feature selection method on the MGH Data against random genes

setwd('/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/Weekly Analyses/Liver Paper Figures/Random_Forest_Classifier_Tuning/')
source('Scripts/CrossValidationLoops_05-26-2020.R')
source('Scripts/FeatureSelection.R')
source('Scripts/RandomForestClassifier.R')
source('Scripts/ROC_plotter.R')
source('Scripts/HomogenizeGGPlots.R')

# Scripts for this batch running
source('GenAUC_MGH.R')
source('GenFeats_MGH.R')
source('GenAUC_KMU.R')

library(ranger)
library(caret)
library(ROCR)
library(pheatmap)
library(gridGraphics)
library(ggplot2)
library(RColorBrewer)
library(parallel)
date <- Sys.Date()
seed = 3070

# load data and set thresholds
data <- readRDS('Files/KMU_MGHDataHybrid_2021-06-25.RDS')

# Load feature subsetting
LM22Genes <- read.csv('Files/2019-12-11_LM22Panels.csv', stringsAsFactors = F)

#IMPORTANT RUN VARIBLES#
########################
pVal.threshs <- 0.25
l2fc.threshs <- 1.5
folds <- 10           
reps <- 10            
########################

set.seed(3070)

MGHCounts <- data[[1]][[1]]
MGHLookup <- data[[1]][[2]]
MGHLookup$sample_name <- MGHLookup$sample_id

KMUCounts <- data[[2]][[1]]
KMULookup <- data[[2]][[2]]

MGH_Res <- matrix(nrow = length(pVal.threshs), ncol=length(l2fc.threshs), dimnames = mtx.names)
MGH_Rand_Res <- matrix(nrow = length(pVal.threshs), ncol=length(l2fc.threshs), dimnames = mtx.names)
KMU_Res <- matrix(nrow = length(pVal.threshs), ncol=length(l2fc.threshs), dimnames = mtx.names)
KMU_Rand_Res <- matrix(nrow = length(pVal.threshs), ncol=length(l2fc.threshs), dimnames = mtx.names)

for(p.idx in 1:length(pVal.threshs)){
  for(fc.idx in 1:length(l2fc.threshs)){
    res <- GenAUC_MGH(MGHCounts = MGHCounts,
                      MGHLookup = MGHLookup,
                      pVal.thresh = pVal.threshs[p.idx],
                      l2fc.thresh = l2fc.threshs[fc.idx],
                      folds = folds,
                      reps = reps,
                      LM22Genes = LM22Genes)
    
    MGH_Res[p.idx,fc.idx] <- res[1]
    MGH_Rand_Res[p.idx,fc.idx] <- res[2]
    
    feats <- GenFeats_MGH(MGHCounts = MGHCounts,
                          MGHLookup = MGHLookup,
                          pVal.thresh = pVal.threshs[p.idx],
                          l2fc.thresh = l2fc.threshs[fc.idx],
                          LM22Genes = LM22Genes)
    res <- GenAUC_KMU(KMUCounts = KMUCounts,
                      KMULookup = KMULookup,
                      features = feats,
                      folds = folds,
                      reps = reps)
    KMU_Res[p.idx, fc.idx] <- res[1]
    KMU_Rand_Res[p.idx, fc.idx] <- res[2]
  }
}

# Any reversed classifiers should be flipped
MGH_Res[MGH_Res<0.5] <- 1-MGH_Res[MGH_Res<0.5]
KMU_Res[KMU_Res<0.5] <- 1-KMU_Res[KMU_Res<0.5]
MGH_Rand_Res[MGH_Rand_Res<0.5] <- 1-MGH_Rand_Res[MGH_Rand_Res<0.5]
KMU_Rand_Res[KMU_Rand_Res<0.5] <- 1-KMU_Rand_Res[KMU_Rand_Res<0.5]

# Turn Col and Row names into friendly formats
data <- list(KMU_Res=KMU_Res,
             KMU_Rand_Res=KMU_Rand_Res,
             MGH_Res=MGH_Res,
             MGH_Rand_Res=MGH_Rand_Res)

saveRDS(data, file = paste0('MGHKMUOutput', date, '.RDS'))

