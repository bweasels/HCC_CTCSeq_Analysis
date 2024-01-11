##Calculate Gene Variance for each gene
setwd('/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/Weekly Analyses/06-20-23/Recreate Single Cell Figures/')
source('/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/Weekly Analyses/Liver Paper Figures/Finalized Figures/HomogenizeGGPlots.R')
date <- Sys.Date()
library(Seurat)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(cowplot)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(magick)
library(DESeq2)
library(circlize)

# Load data and set metadata
data.integrated <- readRDS('../../12-10-21/data.integrated_filtered.RDS')

#trim doublets and over 5% mitochondria & plot QC
metadata <- data.integrated@meta.data
data.integrated <- data.integrated[,metadata$nFeature_RNA<(median(metadata$nFeature_RNA)*2) & metadata$percent.mt<5]
data.integrated$diagnosis <- ifelse(grepl("CLD", data.integrated$orig.ident), 'CLD', 'HCC')

# Filter out high RBC Percentage genes
RBCGenes <- c("HBB", "HBA1", "HBA2", "HBD", "HBM")
counts <- as.matrix(data.integrated@assays$RNA@counts)
data.integrated$RBCPercent <- colSums(counts[RBCGenes,])/colSums(counts)
data.integrated <- data.integrated[,data.integrated$RBCPercent<0.05]

# Find low percentage genes and remove them
counts <- as.matrix(data.integrated@assays$RNA@counts)
percentDropout <- apply(X = counts, MARGIN = 1, FUN = function(x) sum(x==0)/length(x))
data.integrated <- data.integrated[rownames(data.integrated)%in%names(percentDropout)[percentDropout<0.9],]

# Run PCA, find clusters and then run UMAP on filtered data
data.integrated <- RunPCA(data.integrated, dims = 30, verbose = F)
data.integrated <- FindNeighbors(data.integrated, dims = 1:6, reduction='pca')
data.integrated <- FindClusters(data.integrated, resolution = 0.5)
data.integrated <- RunUMAP(data.integrated, reduction = 'pca', dims = 1:6)
data.integrated$FilteredClusters <- Idents(data.integrated)

activationGenes <- c("KLRK1", "KLRD1", "KLRC2", "NCR3","NCR2", "NCR1",
                     "CD226", "FCGR3A", "GZMB", "RF1", "SH2D1B", "KIR2DS1",
                     "KIR2DS2", "KIR2DS3", "KIR2DS4", "IL2RA", "IL2RB",
                     "IL12", "IL12RB1", "IL12RB2", "IL15RA", "IL18R1", 
                     "IFGNR1", "IFGNR2")

activationCytokines <- c("IFNA1", "IFNA2", "IFNB1", "IL2", "IL12A",
                         "IL12B", "IL15", "IL18", "IFNG")

Inhibition <- c("NCR3LG1", "PDCD1", "KLRC1", "TIGIT", "KIR2DL1",
                "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B",
                "KIR3DL1", "KIR3DL2", "KIR3DL3", "IL10RA", "IL10RB")

inhibitionCytokines <- c("IGFB", "IL10")

data.integrated$activationGenes <- colSums(data.integrated[activationGenes,])
data.integrated$inhibition <- colSums(data.integrated[Inhibition,])

#add the score for LM22Extend
lm22Genes <- read.csv('/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/CTC-Seq/feature/2019-09-30_LM22Panels.csv', stringsAsFactors=F)
lm22Genes$EndothelialCells <- c('CLEC4G', 'CLEC4M', 'FLT1', 'PECAM1', 'VWF', 'CD34', rep('', nrow(lm22Genes)-6))
lm22Genes$Hepatocytes <- c('ALB', 'TF', 'CYP3A4', 'CYP2E1', 'APOB', 'ASGR1', 'PCK1', 'HPP', 'APOE', 'ASS1', rep('', nrow(lm22Genes)-10))

lm22Extend <- read.csv("/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/CTC-Seq/feature/2019-11-1_pathwayGenes.csv", stringsAsFactors = F)$LM22.Extended.Curated
lm22Extend <- lm22Extend[lm22Extend!='']

lm22ExtendScore <- data.integrated[rownames(data.integrated)%in%lm22Extend,]
lm22ExtendScore <- colSums(as.matrix(lm22ExtendScore@assays$integrated@data))
metadata <- metadata[match(names(lm22ExtendScore), rownames(metadata)),]
metadata$lm22ExtendedScore <- lm22ExtendScore

#Conduct DESeq on the clusters
if(!file.exists("AllMarkers_FilteredDataset_2023-06-27.csv")){
  markers <- FindAllMarkers(data.integrated)
  write.csv(markers, paste0('AllMarkers_FilteredDataset_2023-06-27.csv'), row.names = F)
}else{
  markers <- read.csv("AllMarkers_FilteredDataset_2023-06-27.csv")
}

#Extract to a per cluster enriched gene list
markerList <- vector(mode = 'list', length = length(unique(markers$cluster)))
TopMarkers <- matrix(nrow = 100, ncol = length(markerList), dimnames = list(seq(1, 100), paste0("Cluster ", seq(1,length(markerList)))))
for (i in 1:length(unique(markers$cluster))){
  markerList[[i]] <- markers[markers$cluster==unique(markers$cluster)[i],]
  TopMarkers[,i] <- markerList[[i]]$gene[1:100]
}
write.csv(TopMarkers, paste0('TopMarkers_ByCluster_', date, '.csv'))

TopMarkers <- gsub("\\.[0-9]$", "", TopMarkers)

#get human protein atlas annotation
hpa_annotate <- read.table('proteinatlas.tsv', header = T, stringsAsFactors = F, sep = "\t", na.strings = "NA")
cellTypes <- hpa_annotate$RNA.blood.cell.specific.NX
names(cellTypes) <- hpa_annotate$Gene

#match the annotation and write it out
clusterTable <- TopMarkers
for(i in 1:ncol(TopMarkers)){
  clustergenes <- TopMarkers[,i]
  clusterTypes <- cellTypes[match(clustergenes, names(cellTypes))]
  clusterTable[,i] <- clusterTypes
}

#subset the hpa_annotate for the expression data of interest
typesofInterest <- c('RNA...granulocytes..NX', 'RNA...monocytes..NX', 'RNA...dendritic.cells..NX', 'RNA...NK.cells..NX', 'RNA...T.cells..NX', 'RNA...B.cells..NX')
expressionTypes <- hpa_annotate[,grep(paste(typesofInterest, collapse = '|'), colnames(hpa_annotate), value = T)]
expressionTypes$genes <- hpa_annotate$Gene
colnames(expressionTypes) <- gsub('Tissue.RNA...(.*)..NX.', '\\1', colnames(expressionTypes))

for(i in 1:ncol(clusterTable)){
  tempGenes <- TopMarkers[,i]
  table <- expressionTypes[match(tempGenes, expressionTypes$genes),!grepl('genes', colnames(expressionTypes))]
  table <- as.matrix(table)
  for(j in 1:length(tempGenes)){
    if(!is.na(table[j,1])){
      cellTypeOrder <- table[j,order(table[j,],decreasing = T)]
      string <- paste0("EXPRESSION: ",names(cellTypeOrder)[1], ":", cellTypeOrder[1], " ; ", names(cellTypeOrder)[2], ':', cellTypeOrder[2], collapse = "")
      tempGenes[j] <- string
    }else{
      tempGenes[j] <- NA
    }
  }
  undefined <- is.na(clusterTable[,i])
  clusterTable[undefined,i] <- tempGenes[undefined]
}

write.csv(clusterTable, paste0('ProteinAtlasAnnotated_TopGenes_', date, '.csv'))

clusters <- c("T Cells", 'NK Cells', 'NK Cells', 'NK Cells', 
              "NK Cells", 'Other Myeloid Lineage Cells', "Other Myeloid Lineage Cells", "Platelets",
              "B Cells", "Platelets", 'B Cells')

# Make heatmap of the top markers to validate the clustering
topGenes <- unique(as.vector(TopMarkers[1:20,]))
expr_mtx <- data.integrated@assays$integrated@data
top_expr <- as.matrix(expr_mtx[topGenes,])
cell_labels <- paste0('Cluster:', Idents(data.integrated), '\n', clusters[Idents(data.integrated)])

pdf(paste0('ClusterAssignmentValidation_', date, '.pdf'), width = 30, height = 20)
Heatmap(top_expr, column_split = cell_labels, show_column_dend = F, show_row_dend = F,
        show_column_names = F)
dev.off()

idents <- sort(unique(Idents(data.integrated)))
data.integrated$IndependentID <- clusters[match(Idents(data.integrated), idents)]
data.integrated$ClusterNumber <- Idents(data.integrated)
Idents(data.integrated) <- data.integrated$IndependentID

# Run UMAP 
umap <- data.integrated@reductions$umap
dat = data.frame('UMAP1' = umap@cell.embeddings[,1],
                 'UMAP2' = umap@cell.embeddings[,2],
                 'Identity' = data.integrated$IndependentID,
                 'ClusterNum' = Idents(data.integrated))


p1 <- DimPlot(data.integrated, pt.size = 1)
p1 <- LabelClusters(p1, id = 'ident', repel = T)
p1 <- p1 + theme(legend.position = 'none', axis.line = element_blank())
p1 <- p1 + theme(axis.text = element_blank(), axis.ticks = element_blank())
p1 <- p1 + labs(x = 'UMAP 1', y = 'UMAP 2') + theme(axis.text = element_text(size = 10))

# go to the NK Cell Subset
NKCells <- data.integrated[,data.integrated$IndependentID=='NK Cells']

# Run PCA, find clusters and then run UMAP
NKCells <- RunPCA(NKCells, dims = 30, verbose = F) 
NKCells <- FindNeighbors(NKCells, dims = 1:6, reduction='pca')
NKCells <- FindClusters(NKCells, resolution = 0.5)
NKCells <- RunUMAP(NKCells, reduction = 'pca', dims = 1:6)
NKCells$NKClusters <- Idents(NKCells)

metadata <- NKCells@meta.data
metadata <- data.frame(sampleID = metadata$orig.ident,
                       Diagnosis = metadata$diagnosis,
                       InhibitedScore = metadata$inhibition,
                       ActivatedScore = metadata$activationGenes)

#metadata <- pivot_longer(metadata, c('InhibitedScore', 'ActivatedScore'), names_to = "Panel")
metadata_hcc <- metadata[metadata$Diagnosis=='HCC',]
metadata_cld <- metadata[metadata$Diagnosis=='CLD',]

# Make histograms
hist_hcc <- ggplot(data = metadata_hcc, aes(x = ActivatedScore)) + geom_histogram(position = 'identity', bins = 20)
#hist_hcc <- hist_hcc + theme_classic() + theme(axis.title = element_text(size = 20, face = 'bold'), plot.title = element_text(size = 20, face = 'bold'))
hist_hcc <- hist_hcc + theme_classic() + theme(axis.title = element_text(size = 20, face = 'bold'), plot.title = element_blank())
hist_hcc <- hist_hcc + theme(axis.title.x = element_blank(), axis.text = element_text(size = 15))
hist_hcc_active <- hist_hcc + labs(y = 'NK Cells (HCC)', title = 'Activation Scores') + xlim(-5, 20)

hist_hcc <- ggplot(data = metadata_hcc, aes(x = InhibitedScore)) + geom_histogram(position = 'identity', bins = 20)
#hist_hcc <- hist_hcc + theme_classic() + theme(axis.text = element_text(size = 15), plot.title = element_text(size = 20, face = 'bold'))
hist_hcc <- hist_hcc + theme_classic() + theme(axis.text = element_text(size = 15), plot.title = element_blank())
hist_hcc <- hist_hcc + theme(axis.title = element_blank())
hist_hcc_inhib <- hist_hcc + labs(title = 'Inhibition Scores') + xlim(-4, 8)

hist_cld <- ggplot(data = metadata_cld, aes(x = ActivatedScore)) + geom_histogram(position = 'identity', bins = 20)
hist_cld <- hist_cld + theme_classic() + theme(axis.text = element_text(size = 15), plot.title = element_blank())
hist_cld <- hist_cld + theme(axis.title = element_text(size = 20, face = 'bold'))
hist_cld_active <- hist_cld + labs(x='Activation Score', y = 'NK Cells (CLD)', title = 'Distribution of Activation Scores\nCLD') + xlim(-5, 20)

hist_cld <- ggplot(data = metadata_cld, aes(x = InhibitedScore)) + geom_histogram(position = 'identity', bins = 20)
hist_cld <- hist_cld + theme_classic() + theme(axis.text = element_text(size = 15), plot.title = element_blank())
hist_cld <- hist_cld + theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20, face = 'bold'))
hist_cld_inhib <- hist_cld + labs(x='Inhibition Score', title = 'Distribution of Inhibition Scores\nCLD') + xlim(-4, 8)


histograms <- plot_grid(hist_hcc_active, hist_hcc_inhib,
                        hist_cld_active, hist_cld_inhib,
                        nrow = 2, ncol = 2)




pdf(paste0("UMAP_Histograms_", date, '.pdf'), 10, 10)
plot(p1)
plot(histograms)
dev.off()

# Settled on score 8 and 1
NKCells$highActivation <- NKCells$activationGenes>8
NKCells$highInhibition <- NKCells$inhibition>1

# Get all Expression and COunts for the NK heatmap
Expr_NK <- as.matrix(NKCells@assays$integrated@data)
Expr_NK <- cbind(Expr_NK[,NKCells$highActivation], Expr_NK[,NKCells$highInhibition])
col_labels <- c(rep('Activated NK Cells', sum(NKCells$highActivation)), rep('Inhibited NK Cells', sum(NKCells$highInhibition)))
col_labels <- factor(col_labels, levels = c('Activated NK Cells', 'Inhibited NK Cells'))
active_genes <- c('GZMB', 'IL2RB', 'FCGR3A', 'KLRD1', 'NCR3', 'SH2D1B')
inhibit_genes <- c('KLRC1', 'TIGIT')
Expr_sub <- rbind(Expr_NK[active_genes,], Expr_NK[inhibit_genes,])
expr_split <- c(rep('Activation Genes', length(active_genes)), rep('Inhibition Genes', length(inhibit_genes)))
expr_split <- factor(expr_split, levels = c('Activation Genes', 'Inhibition Genes'))
overall_heatmap <- Heatmap(Expr_sub, column_split = col_labels, row_split = expr_split,
                           show_column_dend = F, show_row_dend = F, show_column_names = F,
                           cluster_column_slices = F, cluster_row_slices = F, name = 'Expression')
pdf(paste0('PanNKExpressionHistogram_', date, '.pdf'), width = 15, height = 10)
plot(overall_heatmap)
dev.off()

# Explore the TIGIT vs KLRC1 in CTC_TX_092
CTC_TX_092 <- NKCells[,grep('CTC_TX_092', NKCells$orig.ident)]
Cnts_092 <- as.matrix(CTC_TX_092@assays$RNA@counts)
Expr_092 <- as.matrix(CTC_TX_092@assays$integrated@data)

# Subset for high inhibition CTC_TX_092
CTC_TX_092_subset <- CTC_TX_092[,CTC_TX_092$highInhibition]
CTC_TX_092_subset$highTIGIT <- Expr_092['TIGIT',]>1
CTC_TX_092_subset$highKLRC1 <- Expr_092['KLRC1',]>1
CTC_TX_092_subset <- CTC_TX_092_subset[,!(CTC_TX_092_subset$highTIGIT&CTC_TX_092_subset$highKLRC1)]
 
# Find markers
Idents(CTC_TX_092_subset) <- ifelse(CTC_TX_092_subset$highTIGIT, 'High TIGIT', 'High KLRC1')
diff <- FindMarkers(CTC_TX_092_subset, ident.1 = 'High TIGIT', ident.2 = 'High KLRC1')
diff$negLogPAdj <- -log10(diff$p_val_adj)

# Heatmap of all NK Cells & Cluster seperate HCC & CLD
Expr <- as.data.frame(CTC_TX_092_subset@assays$integrated@data)
sigExpr <- Expr[rownames(Expr)%in%c(rownames(diff)[diff$p_val_adj<0.05], 'TIGIT', 'KLRC1'),]
high_sigExpr <- sigExpr[rownames(sigExpr)%in%c(rownames(diff)[diff$AvgExpr>0.5], 'TIGIT', 'KLRC1'),]

l2fcThreshs <- c(2, 4, 8, 10)


colAnno.all <- Idents(CTC_TX_092_subset)
genes<- c('KLRC1', 'TIGIT', 'GZMH', 'CLIC3', 'CD7', 'SELL', 'CD3G', 'CD3D', 'IL32') 
dat <- sigExpr[genes,]

write.table(colnames(CTC_TX_092_subset), paste0('HighTIGITKLRC1Cells_', date, '.txt'),
            row.names = F, col.names = F, quote = F)

library(dendextend)
dat1 <- dat[,colAnno.all=='High KLRC1']
dat2 <- dat[,colAnno.all=='High TIGIT']
clust1 <- as.dendrogram(hclust(dist(t(dat1)), method = 'ward.D'))
clust2 <- as.dendrogram(hclust(dist(t(dat2)), method = 'ward.D'))
group1 <- cutree(clust1, h = 250)
group1 <- c('A', 'B', 'C')[group1]
group2 <- cutree(clust2, h = 250)
group2 <- c('D', 'E')[group2]

colGrp1 <- HeatmapAnnotation(Clustering = group1, 
                             col = list(Clustering = c('A' = 'red', 'B' = 'yellow', 'C' = 'orange')),
                             show_annotation_name = F,
                             show_legend = F)
colGrp2 <- HeatmapAnnotation(Clustering = group2, 
                             col = list(Clustering = c('D' = 'purple', 'E' = 'blue')),
                             show_annotation_name = F,
                             show_legend = F)

p1 <- Heatmap(dat1,
              heatmap_legend_param = list(title = 'Expression', at = c(-2, 4, 10)),
              show_column_names = FALSE,
              cluster_rows = FALSE,
              cluster_columns = clust1,
              top_annotation = colGrp1,
              column_title = NULL,
              show_column_dend = TRUE,
              show_row_dend = FALSE,
              show_row_names = FALSE,
              name = 'ht1',
              border = F,
              show_heatmap_legend = FALSE)

p2 <- Heatmap(dat2,
              heatmap_legend_param = list(title = 'Expression', at = c(-2, 4, 10)),
              show_column_names = FALSE,
              cluster_rows = FALSE,
              cluster_columns = clust2,
              top_annotation = colGrp2,
              column_title = NULL,
              show_column_dend = TRUE,
              show_row_dend = FALSE,
              border = F)


title  <- ggdraw() + draw_text('TIGIT v KLRC1 | abs(l2fc) > 2 & -log10(P Adj) > 25', x = 0.5, y = 0.5)
p <- plot_grid(plotlist = list(title, grid.grabExpr({draw(p1+p2); grid.lines(x=c(0.02,0.85), y=.986, draw=T, gp=gpar(lty=2, col='maroon'))})),
               rel_heights = c(1, 20),
               nrow = 2)
pdf(paste0('HeatmapKLRC1vsTIGIT_', date, '.pdf'), height = 11, width = 10)
plot(p)
dev.off()


