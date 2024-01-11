featureSelectPValInd <- function(counts, lookup, pValThresh, l2fcThresh, subtitle, plot=T, LM22Genes){

  # Get LM22 Genes
  LM22Genes <- unique(unlist(LM22Genes))
  
  # RPM Normalize and trim to immune genes
  rpm <- apply(counts, 2, function(x) x*1e6/sum(x))
  rpm <- rpm[rownames(rpm)%in%LM22Genes,]
  rpm <- rpm[rowSums(rpm)>10,]
  
  # calculate l2fc and pval
  l2fc <- log2(rowMeans(rpm[,lookup$diagnosis=='HCC'])/rowMeans(rpm[,lookup$diagnosis=='CLD']))
  pval <- apply(rpm, 1, function(x) t.test(x[lookup$diagnosis=='HCC'], x[lookup$diagnosis=='CLD'])$p.value)
  
  # stick in a dataframe and plot out the distribution
  data <- data.frame(l2fc = l2fc,
                     pval = pval,
                     row.names = names(l2fc),
                     negLogPVal = -log10(pval),
                     col = ifelse(pval<pValThresh & abs(l2fc)>l2fcThresh, 'red', 'black'),
                     pch = ifelse(pval<pValThresh & abs(l2fc)>l2fcThresh, 1, 20))
  data <- data[!is.na(data$l2fc),]
  data <- data[!is.infinite(data$l2fc),]
  
  if(plot){
    plot(data$l2fc,
         data$negLogPVal,
         col = data$col,
         pch = data$pch,
         xlab = 'Log2(Fold Change)',
         ylab = '-log10(pValue)',
         main = paste0('Significant Genes\n', subtitle))
  }
  # Subset the data and reorder
  data <- data[data$pval<pValThresh & abs(data$l2fc)>l2fcThresh,]
  
  if (nrow(data) > 5){
    rpm <- rpm[rownames(rpm)%in%rownames(data),, drop = F]
    rpm <- cbind(rpm[,lookup$diagnosis=='HCC', drop = F], rpm[,lookup$diagnosis=='CLD', drop = F])
    lookup <- lookup[match(colnames(rpm), lookup$sample_name),]
    
    # heatmap of the top points
    if(plot){
      annotation <- data.frame(row.names = lookup$sample_name, diagnosis = lookup$diagnosis)
      pheatmap(log10(rpm+1), show_rownames = F, show_colnames = F, border_color = NA, annotation_col = annotation,
               main = paste0('Log10 Normalized Counts of Top Variable Genes\n', subtitle), cluster_cols = F, gaps_col = sum(annotation$diagnosis=='HCC'))
    }
  }

return(rownames(rpm))
}

expressionMatchRand <- function(genes, counts, lookup, subtitle, plot = T){
  
  # Make RPMs and get mean and sd for all genes
  rpm <- apply(counts, 2, function(x) x*1e6/sum(x))
  overallStats <- data.frame(sd = apply(rpm, 1, sd),
                             mean = apply(rpm, 1, mean),
                             row.names = rownames(rpm))
  
  featureStats <- overallStats[rownames(overallStats)%in%genes,]
  
  randGenes <- rownames(featureStats)
  for(i in 1:nrow(featureStats)){
    if(featureStats$sd[i]!=0){
      minVal <- featureStats$mean[i]-featureStats$sd[i]
      maxVal <- featureStats$mean[i]+featureStats$sd[i]
      potentialGenes <- overallStats[overallStats$mean>minVal & overallStats$mean<maxVal,]
      randGenes[i] <- sample(rownames(potentialGenes), size = 1)
    }
  }
  
  rpm <- rpm[rownames(rpm)%in%randGenes,, drop = F]
  rpm <- cbind(rpm[,lookup$diagnosis=='HCC', drop = F], rpm[,lookup$diagnosis=='CLD', drop = F])
  lookup <- lookup[match(colnames(rpm), lookup$sample_name),]
  
  # heatmap of the top points
  if(plot){
    annotation <- data.frame(row.names = lookup$sample_name, diagnosis = lookup$diagnosis)
    pheatmap(log10(rpm+1), show_rownames = F, show_colnames = F, border_color = NA, annotation_col = annotation,
             main = paste0('Expression Matched Random Genes\n', subtitle), cluster_cols = F, gaps_col = sum(annotation$diagnosis=='HCC'))
  }
  return(randGenes)
}