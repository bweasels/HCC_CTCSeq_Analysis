##Look at individual LM22 WBC Scores betweeen KMU high and low and see if there's a difference.

#########BOILERPLATE CODE FROM LOAD COUNTS###################
source('/Dropbox (Partners HealthCare)/Templates/RPM.R')
source('/Dropbox (Partners HealthCare)/CTC-Seq/scripts/depleteRBC.R')
setwd('/Dropbox (Partners HealthCare)/CTC-Seq/input/')

threshold <- 0.5
#load the data
counts.df <- read.csv('2019-08-12_KMUCounts_star_genes_erc.csv', row.names = 1)
colnames(counts.df) <- gsub('^X(.*)', '\\1', colnames(counts.df))
key <- read.csv('2019-08-12_KMUMetadata.csv', row.names = 1)

#confirm that the metadata is in alignment with the counts table
if(sum(colnames(counts.df)!=key$Sample.ID)>0){
  stop('WARNING: Lookup table and RNASeq data tables are not similarly ordered.')
}

#create lookup df
patientStatus <- ifelse(grepl('CLD|LC', key$Patient.Status), 'CLD', 'HCC')
lookup.df <- data.frame(sample_name = key$Sample.ID, diagnosis = patientStatus)

#get the total number of counts from each sample and remove lowQ samples
colSums <- colSums(counts.df)

counts.df <- counts.df[,colSums>100000]
lookup.df <- lookup.df[colSums>100000,]

#deplete RBCs
#counts.df <- depleteRBC(counts.df)

#generate liverScore
liverPanel <- c('FGB', 'FGG', 'ALB', 'AFP', 'FABP1', 'GPC3', 'RBP4', 'RBP4', 'AHSG', 'APOH', 'TF')
rpm.df <- apply(counts.df, 2, function(x) x*1e6/sum(x))
rpm.df <- log10(rpm.df+1)
liverScore <- colSums(rpm.df[rownames(rpm.df) %in% liverPanel,])

################END BOILERPLATE CODE##########################
library(pheatmap)
library(ggplot2)
library(gridExtra)
setwd('/Dropbox (Partners HealthCare)/Weekly Analyses/10-15-19/Hand Deconvolution and Immune-Liver signature/')
date <- Sys.Date()

#order RPM.df to be lowq followed by high q
rpm.df <- cbind(rpm.df[,liverScore<0.5], rpm.df[,liverScore>0.5])
liverScore <- c(liverScore[liverScore<0.5], liverScore[liverScore>0.5])

#load data and make output matrix
cibersortPanel <- read.csv('/Dropbox (Partners HealthCare)/CTC-Seq/feature/2019-09-30_LM22Panels.csv')
cibersortScores <- matrix(ncol = ncol(rpm.df), nrow = ncol(cibersortPanel), dimnames = list(colnames(cibersortPanel), colnames(rpm.df)))

#Put the expression scores into the matrix & tell it to make boxplots to compare the sample scores
pdf(paste0('CibersortScores-lowvshighLiverScore_', date, '.pdf'))
for (i in 1:ncol(cibersortPanel)){
  score <- colSums(rpm.df[rownames(rpm.df) %in% cibersortPanel[,i],])
  cibersortScores[i,] <- score
  
  score <- data.frame(row.names = names(score), score = score, quality = c(rep('Low', sum(liverScore<0.5)), rep('High', sum(liverScore>0.5))))
  g <- ggplot(score, aes(x=quality, y=score, color=quality)) + geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.2))
  g <- g + ggtitle(rownames(cibersortScores)[i]) + xlab('Sample Quality') + ylab('Log10(Expression Score)')
  plot(g)
}

#Make heatmaps of the scores
pheatmap(cibersortScores,
         cluster_cols = F,
         show_colnames = F,
         gaps_col = sum(liverScore<0.5),
         main = 'Cibersort gene Signature Scores Unscaled')

cibersortScores <- t(scale(t(cibersortScores), scale=F, center = T))

pheatmap(cibersortScores,
         cluster_cols = F,
         show_colnames = F,
         gaps_col = sum(liverScore<0.5),
         main = 'Cibersort gene Signature Scores Centered')
dev.off()

##These didn't work very well :( So now brute force going through each gene and finding what are associated with each score
cibersortGenes <- unlist(cibersortPanel)
cibersortGenes <- unique(cibersortGenes)

#Make sure that liverscores and RPM are ordered the same way
if(sum(names(liverScore)!=colnames(rpm.df))>0){
  error('LiverScore and RPM.df are not aligned!!!')
}

geneCorrelation <- matrix(nrow = length(cibersortGenes), ncol = 2, dimnames = list(cibersortGenes, c('Correlation Coefficient', 'p.value')))

for(i in 1:length(cibersortGenes)){
  expression <- rpm.df[match(cibersortGenes[i], rownames(rpm.df)),]
  
  #some of these genes don't exist, so skip if that happens
  if(sum(is.na(expression))==0){
    correlation <- cor.test(liverScore, expression)
    geneCorrelation[i,1] <- correlation$estimate
    geneCorrelation[i,2] <- correlation$p.value
  }
}
#convert geneCorrelation to a data.frame, remove NAs and calculate p.adjusted value
geneCorrelation <- as.data.frame(geneCorrelation)
geneCorrelation <- geneCorrelation[!is.na(geneCorrelation$p.value),]
geneCorrelation$p.adj <- p.adjust(geneCorrelation$p.value, method = 'BH', n = length(cibersortGenes))

#Select for genes with a p.adj < 0.05 and order them
geneCorrelation <- geneCorrelation[geneCorrelation$p.adj<0.05,]
geneCorrelation <- geneCorrelation[order(geneCorrelation$p.adj),]

##See which cell type is producing a good number of useful genes and use that for a story
correlatedCellType <- vector('list', length = ncol(cibersortPanel))
names(correlatedCellType) <- colnames(cibersortPanel)

for(i in 1:length(correlatedCellType)){
  genes <- unlist(cibersortPanel[,i])
  correlations <- geneCorrelation[rownames(geneCorrelation)%in%genes, ]
  correlatedCellType[[i]] <- correlations
}

sapply(correlatedCellType, function(x) nrow(x))

#see how the depleted panel correlates
depletePanel <- read.csv('/Dropbox (Partners HealthCare)/CTC-Seq/input/Molecular_Depletion_Panel.csv')
depletePanel <- depletePanel$Gene.ID
depleteScore <- c()
for(i in 1:ncol(rpm.df)){
  score <- sum(rpm.df[rownames(rpm.df)%in%depletePanel,i])
  depleteScore <- c(depleteScore, score)
}

#it correlates positively???!??? Dude this makes literally no sense.
depleteCorrelation <- cor.test(liverScore, depleteScore)

##Make boxplots for each of the cell signature scores
pdf(paste0('SignatureScoresofCorrelatedGenes_',date, '.pdf'), width = 16)

nGenes <- data.frame(panel.names = c(names(correlatedCellType), 'RBC.Deplete.Score'), nGenes = c(sapply(correlatedCellType, function(x) nrow(x)), length(depleteScore)))
g <- ggplot(nGenes, aes(x=panel.names, y=nGenes, fill=panel.names)) + geom_bar(stat='identity')
g <- g + ggtitle('Number of Significantly Correlated Genes From Each Cell Type') + xlab('') + ylab('Number of Genes') 
g <- g + guides(fill=FALSE) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
plot(g)

for(i in 1:length(correlatedCellType)){
  genes <- rownames(correlatedCellType[[i]])
  
  #dealing with edge cases
  if(length(genes)!=0){
    
    if(length(genes)>2){
      score <- colSums(rpm.df[rownames(rpm.df)%in%genes,])
    }else if(length(genes)==1){
      score <- rpm.df[rownames(rpm.df)%in%genes,]
    }
    #make the ggplot dataframe
    data <- data.frame(rownames = names(score), CellTypeScore = score, Quality = c(rep('Low', sum(liverScore<0.5)), rep('High', sum(liverScore>0.5))), LiverScore = liverScore)
    
    correlation <- cor.test(data$LiverScore, data$CellTypeScore)
    p.val <- correlation$p.value
    corr <- correlation$estimate
    
    #plot boxplots
    box <- ggplot(data, aes(x=Quality, y=CellTypeScore, color=Quality)) + geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.2))
    box <- box + ggtitle(paste(names(correlatedCellType)[i], 'Correlated Genes Score\n', 'P.value:', round(p.val, digits = 3), '\nCorrelation:', round(corr, digits = 3))) + xlab('Sample Quality') + ylab('Log10(Expression Score)')
    
    #plot scatter plots
    scatter <- ggplot(data, aes(x=LiverScore, y=CellTypeScore, color=Quality)) + geom_point()
    scatter <- scatter + ggtitle(paste(names(correlatedCellType)[i], 'vs Liver Panel')) + xlab('Liver Score') + ylab('Cell Panel Score')
    
    #arrange them side by side
    grid.arrange(box, scatter, ncol = 2)
  }
}

#repeat the exercise for RBCDepletionGenes
depleteCorrelation <- cor.test(liverScore, depleteScore)
data <- data.frame(CellTypeScore = depleteScore, Quality = c(rep('Low', sum(liverScore<0.5)), rep('High', sum(liverScore>0.5))), LiverScore = liverScore)
box <- ggplot(data, aes(x=Quality, y=CellTypeScore, color=Quality)) + geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.2))
box <- box + ggtitle(paste('Depletion Panel Score\n', 'P.value:', round(depleteCorrelation$p.value, digits = 3), '\nCorrelation:', round(depleteCorrelation$estimate, digits = 3))) + xlab('Sample Quality') + ylab('Log10(Expression Score)')

#plot scatter plots
scatter <- ggplot(data, aes(x=LiverScore, y=CellTypeScore, color=Quality)) + geom_point()
scatter <- scatter + ggtitle('Depletion Panel vs Liver Panel') + xlab('Liver Score') + ylab('Depletion Panel Score')

grid.arrange(box, scatter, ncol = 2)
dev.off()

bestGenes <- c()
##Gather the highest correlation genes 
for(i in 1:length(correlatedCellType)){
  data <- correlatedCellType[[i]]
  bestGenes <- c(bestGenes, rownames(data)[data$`Correlation Coefficient`>0.8])
}

###Seems like the best genes are just the same as our panel - lets try depleting them and seeing how the panels perform
bestGenes <- unique(bestGenes)

for(i in 1:length(correlatedCellType)){
  data <- correlatedCellType[[i]]
  data <- data[!rownames(data)%in%liverPanel,]
  correlatedCellType[[i]] <- data
}

##Make boxplots for each of the cell signature scores
pdf(paste0('SignatureScores_LiverPanelDepleted_',date,'.pdf'), width = 16)

nGenes <- data.frame(panel.names = c(names(correlatedCellType), 'RBC.Deplete.Score'), nGenes = c(sapply(correlatedCellType, function(x) nrow(x)), length(depleteScore)))
g <- ggplot(nGenes, aes(x=panel.names, y=nGenes, fill=panel.names)) + geom_bar(stat='identity')
g <- g + ggtitle('Number of Significantly Correlated Genes From Each Cell Type') + xlab('') + ylab('Number of Genes') 
g <- g + guides(fill=FALSE) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
plot(g)

for(i in 1:length(correlatedCellType)){
  genes <- rownames(correlatedCellType[[i]])
  
  #dealing with edge cases
  if(length(genes)!=0){
    
    if(length(genes)>2){
      score <- colSums(rpm.df[rownames(rpm.df)%in%genes,])
    }else if(length(genes)==1){
      score <- rpm.df[rownames(rpm.df)%in%genes,]
    }
    #make the ggplot dataframe
    data <- data.frame(rownames = names(score), CellTypeScore = score, Quality = c(rep('Low', sum(liverScore<0.5)), rep('High', sum(liverScore>0.5))), LiverScore = liverScore)
    
    correlation <- cor.test(data$LiverScore, data$CellTypeScore)
    p.val <- correlation$p.value
    corr <- correlation$estimate
    
    #plot boxplots
    box <- ggplot(data, aes(x=Quality, y=CellTypeScore, color=Quality)) + geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.2))
    box <- box + ggtitle(paste(names(correlatedCellType)[i], 'Correlated Genes Score\n', 'P.value:', round(p.val, digits = 3), '\nCorrelation:', round(corr, digits = 3))) + xlab('Sample Quality') + ylab('Log10(Expression Score)')
    
    #plot scatter plots
    scatter <- ggplot(data, aes(x=LiverScore, y=CellTypeScore, color=Quality)) + geom_point()
    scatter <- scatter + ggtitle(paste(names(correlatedCellType)[i], 'vs Liver Panel')) + xlab('Liver Score') + ylab('Cell Panel Score')
    
    #arrange them side by side
    grid.arrange(box, scatter, ncol = 2)
  }
}

#repeat the exercise for RBCDepletionGenes
depleteCorrelation <- cor.test(liverScore, depleteScore)
data <- data.frame(CellTypeScore = depleteScore, Quality = c(rep('Low', sum(liverScore<0.5)), rep('High', sum(liverScore>0.5))), LiverScore = liverScore)
box <- ggplot(data, aes(x=Quality, y=CellTypeScore, color=Quality)) + geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.2))
box <- box + ggtitle(paste('Depletion Panel Score\n', 'P.value:', round(depleteCorrelation$p.value, digits = 3), '\nCorrelation:', round(depleteCorrelation$estimate, digits = 3))) + xlab('Sample Quality') + ylab('Log10(Expression Score)')

#plot scatter plots
scatter <- ggplot(data, aes(x=LiverScore, y=CellTypeScore, color=Quality)) + geom_point()
scatter <- scatter + ggtitle('Depletion Panel vs Liver Panel') + xlab('Liver Score') + ylab('Depletion Panel Score')

grid.arrange(box, scatter, ncol = 2)
dev.off()


##############################Winnowing out the Interesting Genes######################
#it seems like the following cell types have a decent correlateion without liver panel
#T.cells.gamma.delta 0.853
#B.cells.naive 0.86
#B.cells.memory 0.862
#T.cells.follicular.helper 0.842
#T.cells.regulatory 0.822
#Lets combine these gene sets and see how we can affect correlation

selected <- correlatedCellType[grep('gamma.delta|B.cells|follicular.helper|T.cells.regulatory|CD4.naive', names(correlatedCellType))]
collatedSelected <- selected[[1]]
genes <- rownames(selected[[1]])
for(i in 2:length(selected)){
  temp <- selected[[i]]
  genes <- c(genes, rownames(temp))
  collatedSelected <- rbind(collatedSelected, temp)
}
genes <- unique(genes)

pdf(paste0('HybridPanelStatistics_', date, '.pdf'), width = 14)
nGenes <- data.frame(panel.names = names(selected), nGenes = sapply(selected, function(x) nrow(x)))
g <- ggplot(nGenes, aes(x=panel.names, y=nGenes, fill=panel.names)) + geom_bar(stat='identity')
g <- g + ggtitle('Number of Significantly Correlated Genes From Each Cell Type') + xlab('') + ylab('Number of Genes') 
g <- g + guides(fill=FALSE) + theme(axis.text.x=element_text(angle = 90, hjust = 1))
plot(g)


hybridPanel <- data.frame(row.names = colnames(rpm.df), 
                          PanelScore = colSums(rpm.df[rownames(rpm.df)%in%genes,]),
                          liverScore = liverScore,
                          Quality = c(rep('Low', sum(liverScore<0.5)), rep('High', sum(liverScore>0.5))))
hybridCorr <- cor.test(hybridPanel$liverScore, hybridPanel$PanelScore)

box <- ggplot(hybridPanel, aes(x=Quality, y=PanelScore, color=Quality)) + geom_boxplot() + geom_jitter(shape = 16, position=position_jitter(0.2))
box <- box + ggtitle(paste('Hybrid Panel Score\n', 'P.value:', round(hybridCorr$p.value, digits = 3), '\nCorrelation:', round(hybridCorr$estimate, digits = 3))) + xlab('Sample Quality') + ylab('Log10(Expression Score)')

#plot scatter plots
scatter <- ggplot(hybridPanel, aes(x=liverScore, y=PanelScore, color=Quality)) + geom_point()
scatter <- scatter + ggtitle('Hybrid Panel vs Liver Panel') + xlab('Liver Score') + ylab('Hybrid Panel Score')

grid.arrange(box, scatter, ncol = 2)

hist(hybridPanel$PanelScore[liverScore>0.5], 
     main = 'Distribution of Hybrid Panel Scores\nRed: High Liver Panel Score\nBlue:Low Liver Panel Score',
     col = rgb(1,0,0,0.5),
     xlim=c(0,100),
     ylim=c(0,25), 
     breaks = seq(from=0, to=100, by = 2))
hist(hybridPanel$PanelScore[liverScore<0.5],
     col = rgb(0,0,1,0.5),
     add=T,
     breaks = seq(from=0, to=100, by = 2))

dev.off()

####Write out the Hybrid Panel Score
write.csv(genes, paste0('HybridPanelGenes_', date, '.csv'), row.names = F)
