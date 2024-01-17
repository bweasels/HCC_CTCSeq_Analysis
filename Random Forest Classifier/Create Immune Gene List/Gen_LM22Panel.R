#Make a wbc subtype panel from the CIBERSORT data
library(pheatmap)
date <- Sys.Date()
setwd('E:/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/Weekly Analyses/09-23-19/Regenerate Cell Panel LM22/')

#load Cibersort's 22 sample panel
cibersortLM22Counts <- read.table('E:/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/Weekly Analyses/08-20-19/CellTypePanel/LM22-ref-sample.txt', sep = '\t', row.names = 1)
cibersortLM22Ref <- read.table('E:/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/Weekly Analyses/08-20-19/CellTypePanel/LM22-classes.txt', sep = '\t', row.names = 1)

#look for samples which are not 0 in the ref file (0s are unused)
cibersortSamps <- colSums(cibersortLM22Ref)
cibersortLM22Counts <- cibersortLM22Counts[,cibersortSamps!=0]
cibersortLM22Ref <- cibersortLM22Ref[,cibersortSamps!=0]

#map the Counts Data to the samples indexed in the ref
sampleID <- c()
for(i in 1:ncol(cibersortLM22Counts)){
  temp <- cibersortLM22Ref[,i, drop = F]
  sample <- rownames(temp)[temp==1]
  sampleID <- c(sampleID, sample)
}

#make the output foldchange matrix
colnames(cibersortLM22Counts) <- sampleID
cellTypes <- unique(sampleID)
FoldChange <- matrix(nrow = nrow(cibersortLM22Counts), ncol = length(cellTypes), dimnames = list(rownames(cibersortLM22Counts), cellTypes))

#put the fold change (celltype vs all else) into the foldchange matrix
for(i in 1:length(cellTypes)){
  cellData <- cibersortLM22Counts[,grepl(cellTypes[i], colnames(cibersortLM22Counts))]
  otherData <- cibersortLM22Counts[,!grepl(cellTypes[i], colnames(cibersortLM22Counts))]
  FoldChange[,i] <- rowMeans(cellData)/rowMeans(otherData)
}

signature <- vector('list', length = length(cellTypes))
names(signature) <- cellTypes

#plot out the histogram of the log2 fold change
pdf(paste0('HistogramOfLog2FoldChange_', date, '.pdf'))
hist(log2(temp), 50, main = 'Histogram of L2FC Cibersort')
dev.off()

#select genes with a fold change > 8 for the signature
for(i in 1:length(signature)){
  selected <- FoldChange[FoldChange[,i]>3, i]
  signature[[i]] <- names(sort(selected, decreasing = T))
}

#make the heatmap to make sure that I've selected a good panel
outputMatrix <- matrix(nrow = 0, ncol = length(signature))

cibersortLM22Ordered <- cibersortLM22Counts[,1]
#reorder the counts so they go in order
for(i in 1:length(cellTypes)){
  cibersortLM22Ordered <- cbind(cibersortLM22Ordered, cibersortLM22Counts[,grep(cellTypes[i], colnames(cibersortLM22Counts))])
}
cibersortLM22Ordered <- cibersortLM22Ordered[,-1]
rownames(cibersortLM22Ordered) <- rownames(cibersortLM22Counts)
cibersortLM22Ordered <- log10(cibersortLM22Ordered+1)

#Look at the "upregulated genes" for each item in the panel
for(i in 1:length(signature)){
  data <- cibersortLM22Ordered[match(signature[[i]], rownames(cibersortLM22Ordered)),]
  outputMatrix <- rbind(outputMatrix, data)
}

outputMatrix <- t(scale(t(outputMatrix), center = T, scale = F))


#heatmap it out (we got a good pattern! so works within our data)
pdf(paste0(date, '_LM22DerivedPanelValidation.pdf'), width = 40, height = 20)
pheatmap(outputMatrix, cluster_rows = F, cluster_cols = F)
dev.off()

#output the panel
outputPanels <- matrix(nrow = max(lengths(signature)), ncol = length(signature))
colnames(outputPanels) <- names(signature)
for(i in 1:length(signature)){
  genes <- signature[[i]]
  for(j in 1:length(genes)){
    outputPanels[j,i] <- genes[j]
  }
}

write.csv(outputPanels, paste0(date,'_LM22Panels.csv'), row.names = F)
