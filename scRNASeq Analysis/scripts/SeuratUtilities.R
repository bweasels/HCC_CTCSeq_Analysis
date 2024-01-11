plotQCBasics <- function(data, label){
  metadata <- data@meta.data
  counts <- data@assays$RNA@data
  #Make the QC plots
  plot1 <- advnUMIPlots(counts = counts, metadata = metadata, title = paste('Histogram of nUMI for', label))
  plot2 <- VlnPlot(data, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, combine = T) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  plot3 <- FeatureScatter(data, feature1 = 'nCount_RNA', feature2 = 'percent.mt') + ggtitle('')
  plot4 <- FeatureScatter(data, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') + ggtitle('')
  
  combo <- plot_grid(plot1, plot2, plot3, plot4, ncol = 2, align='h')
  title <- ggdraw() + draw_label(label)
  
  return(plot_grid(title, combo, ncol=1, rel_heights=c(0.05,1)))
}

advnUMIPlots <- function(counts, metadata, title){
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(cowplot)
  
  #get the most common gene for reach cell
  maxGene <- rownames(counts)[apply(counts, 2, which.max)]
  metadata$mostCommonGene <- maxGene
  
  #define the bins and get the most common max gene from each
  bins <- seq(from = min(metadata$nFeature_RNA), to = max(metadata$nFeature_RNA), length = 50)
  colors <- bins
  freq <- bins
  for(j in 1:length(colors)){
    if(j !=length(colors)){
      mostCommonGene <- as.factor(maxGene[metadata$nFeature_RNA>bins[j]&metadata$nFeature_RNA<bins[j+1]])
      freq[j] <- sum(metadata$nFeature_RNA>bins[j]&metadata$nFeature_RNA<bins[j+1])
    }else{
      mostCommonGene <- as.factor(maxGene[metadata$nFeature_RNA>bins[j]])
      freq[j] <- sum(metadata$nFeature_RNA>bins[j])
    }
    if(length(mostCommonGene)>0){
      colors[j] <- names(sort(summary(mostCommonGene), decreasing = T)[1])
    }else{
      colors[j] <- "null"
      freq[j] <- 0
    }
  }
  colors <- as.factor(colors)
  
  data <- data.frame(nUMI = bins, Frequency = freq, color = colors)
  ggplot(data,aes(x=nUMI, y=Frequency, fill=colors))+geom_bar(color = 'white', stat = 'identity')+ggtitle(title)
}

loadAndIntegrateSamples <- function(sample.dirs, samples, filterRBC = FALSE, minCells=150){
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(cowplot)
  date <- Sys.Date()
  
  #make a holder
  data.list <- vector('list', length = length(samples))
  names(data.list) <- samples
  
  #read in data and output QC Plots
  pdf(paste(date, '_QCPlots.pdf'), width = 15, height = 10)
  for(i in 1:length(data.list)){
    #load data & Label it as the sample ID
    data <- Read10X(data.dir = sample.dirs[i])
    data <- CreateSeuratObject(data, min.cells = 3, min.features=200)
    data@meta.data$orig.ident <- samples[i]
    data <- PercentageFeatureSet(data, pattern = '^MT-', col.name = 'percent.mt')
    
    #plot the QC Basics
    plot(plotQCBasics(data, samples[i]))
    
    if(filterRBC){
      # Filter out high RBC Percentage genes
      RBCGenes <- c("HBB", "HBA1", "HBA2", "HBD", "HBM")
      counts <- as.matrix(data@assays$RNA@counts)
      data$RBCPercent <- colSums(counts[RBCGenes,])/colSums(counts)
      data <- data[,data$RBCPercent<0.5]
      
    }
    
    if(ncol(data)>minCells){
      #Samples with n UMI < 750 are platelets and RBCs and >2500 is possibly doublets, so filter all these out
      data <- SCTransform(data, vars.to.regress = 'percent.mt', verbose = F)
      
      #For the sake of the repeats, we wont regress to a variable
      #data <- SCTransform(data, verbose = F)
      data.list[[i]] <- data
      print(paste('Completed Loading Sample:', names(data.list)[i]))
      rm(data)
    }
  }
  dev.off()
  print(data.list)
  #data.list <- data.list[lapply(data.list, function(x) !is.na(x))]
  for (i in 1:length(data.list)){
    if (is.null(data.list[[i]])){
      data.list[[i]] <- NULL
    }
  }
    
  #select features for downstream integration and calculae all the necessary pearson residuals
  options(future.globals.maxSize=4000*1024^2)
  data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 2000)
  data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features, verbose = F)
  
  #Identify anchors and integrate the datasets
  data.list <- FindIntegrationAnchors(object.list = data.list, normalization.method = 'SCT', anchor.features = data.features, verbose = T)
  data.integrated <- IntegrateData(anchorset = data.list, normalization.method = 'SCT', verbose = T)
  #saveRDS(object = data.integrated, 'data.integrated.RDS')
  return(data.integrated)
}

sampleComposition <- function(dataset){
  plot <- ggplot(dataset@meta.data, aes(x=orig.ident, fill = orig.ident))
  plot <- plot + geom_bar()+NoLegend()+theme(axis.text.x=element_text(angle=90, hjust=1))
  return(plot)
}