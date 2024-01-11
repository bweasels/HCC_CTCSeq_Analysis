##Redo my dataloading for the 10X Data
setwd("/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/Weekly Analyses/12-10-21/")
source('scripts/SeuratUtilities.R')
date <- Sys.Date()
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)

#Load counts and perform sctransform and intersample normalizations
if(!file.exists('data.integrated_filtered.RDS')){
  
  # get the samples and make a holder
  sample.dirs <- list.dirs('/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/10X_LiverCancer/counts_holder/')
  sample.dirs <- grep('hg19', sample.dirs, value = T)
  samples <- toupper(gsub('.*//(.*)/hg19', '\\1', sample.dirs))
  
  # Filter out the high rbc cells
  samps_toKeep <- c("CLD_025", "CLD_073", "CLD_074", "CLD_085",
                    "CLD_134", "CTC_TX_090", "CTC_TX_092", "CTC_TX_101")
  sample.dirs <- sample.dirs[samples%in%samps_toKeep]
  samples <- samples[samples%in%samps_toKeep]
  
  #load the data and plot QC using my homebrew function
  data.integrated <- loadAndIntegrateSamples(sample.dirs = sample.dirs, samples = samples, filterRBC = T, minCells = 50)
  saveRDS(data.integrated, 'data.integrated_filtered.RDS')
  
}else{
  #if we've already did this, just load the saved object
  data.integrated <- readRDS('data.integrated_filtered.RDS')
}