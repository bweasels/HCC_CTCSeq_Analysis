StratifiedKmeansSplit <- function(lookup, nFolds){
  require(caret)
  
  #break out into hcc and cld
  lookup.HCC <- lookup[lookup$diagnosis=='HCC',]
  lookup.CLD <- lookup[lookup$diagnosis=='CLD',]
  
  # make the folds
  output <- vector(mode='list', length = nFolds)
  folds.HCC <- createFolds(lookup.HCC$diagnosis, k = nFolds)
  folds.CLD <- createFolds(lookup.CLD$diagnosis, k = nFolds)
  
  # translate the hcc and cld indices to overall indices
  for(i in 1:length(output)){
    hccIndices <- match(lookup.HCC$sample_name[folds.HCC[[i]]], lookup$sample_name)
    cldIndices <- match(lookup.CLD$sample_name[folds.CLD[[i]]], lookup$sample_name)
    output[[i]] <- c(hccIndices, cldIndices)
  }
  return(output)
}

# function to plot CV disributions to make sure splits are random

plotCVDistributions <- function(lookup, splits, title){
  require(pheatmap)
  data <- data.frame(id = lookup$sample_name,
                     diagnosis = lookup$diagnosis)
  for(i in 1:length(splits)){
    selection <- vector('numeric', length = nrow(lookup))
    selection[splits[[i]]] <- 1
    data <- cbind(data, selection)
  }
  
  plotData <- t(data[,3:ncol(data)])
  colnames(plotData) <- data$id
  
  annotation <- data.frame(diagnosis = data$diagnosis,
                           row.names = data$id)
  plotData <- plotData[,c(match(rownames(annotation)[annotation$diagnosis=='CLD'], colnames(plotData)),
                          match(rownames(annotation)[annotation$diagnosis=='HCC'], colnames(plotData)))]
  
  graph <- pheatmap(plotData,
                   cluster_rows = F,
                   cluster_cols = F,
                   annotation_col = annotation,
                   gaps_col = sum(annotation$diagnosis=='CLD'),
                   show_colnames = F,
                   legend = F,
                   border_color = NA,
                   main = title)
  return(graph)
}