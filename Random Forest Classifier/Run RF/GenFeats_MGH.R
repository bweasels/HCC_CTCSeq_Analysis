GenFeats_MGH <- function(MGHCounts, MGHLookup, pVal.thresh, l2fc.thresh, LM22Genes){
  
  #perform feature selection
  features <- featureSelectPValInd(counts = MGHCounts,
                                   lookup = MGHLookup,
                                   pValThresh = pVal.thresh,
                                   l2fcThresh = l2fc.thresh,
                                   subtitle = 'Final MGH Feature Selection',
                                   plot = F,
                                   LM22Genes = LM22Genes)
  return(features)
}
