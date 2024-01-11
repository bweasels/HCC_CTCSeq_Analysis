GenAUC_MGH <- function(MGHCounts, MGHLookup, pVal.thresh, l2fc.thresh, folds, reps, LM22Genes){
  # Make the k means stratified split
  MGH_CV <- StratifiedKmeansSplit(MGHLookup, folds)
  
  # Make holder for the predictions and the importances
  RandVotes <- array(NA, dim = c(length(MGHLookup$sample_name), reps, folds), dimnames = list(MGHLookup$sample_name, paste0('rep',1:reps), paste0('fold', 1:folds)))
  MGHVotes <- array(NA, dim = c(length(MGHLookup$sample_name), reps, folds), dimnames = list(MGHLookup$sample_name, paste0('rep',1:reps), paste0('fold', 1:folds)))
  MGHImps <- vector(mode = 'list', length = folds)
  
  # Run folds and reps
  for(fold in 1:folds){
    featuresMGH <- featureSelectPValInd(counts = MGHCounts[,-MGH_CV[[fold]]],
                                        lookup = MGHLookup[-MGH_CV[[fold]], ],
                                        pValThresh = pVal.thresh,
                                        l2fcThresh = l2fc.thresh,
                                        subtitle = paste('MGH Fold', fold),
                                        plot = F,
                                        LM22Genes = LM22Genes)
    featuresRand <- expressionMatchRand(genes = featuresMGH, 
                                        counts = MGHCounts[,-MGH_CV[[fold]]],
                                        lookup = MGHLookup[-MGH_CV[[fold]], ],
                                        subtitle = paste('MGH Fold', fold),
                                        plot = F)
    
    #Break to training set and subset selected features
    MGHCounts.train <- data.frame(t(MGHCounts[rownames(MGHCounts)%in%featuresMGH, -MGH_CV[[fold]]]))
    MGHCounts.train$diagnosis <- MGHLookup$diagnosis[-MGH_CV[[fold]]]
    MGHRand.train <- data.frame(t(MGHCounts[rownames(MGHCounts)%in%featuresRand, -MGH_CV[[fold]]]))
    MGHRand.train$diagnosis <- MGHLookup$diagnosis[-MGH_CV[[fold]]]
    
    # Make Test set and similarly trim to feature set
    MGHCounts.test <- data.frame(t(MGHCounts[rownames(MGHCounts)%in%featuresMGH, MGH_CV[[fold]]]))
    MGHRand.test <- data.frame(t(MGHCounts[rownames(MGHCounts)%in%featuresRand, MGH_CV[[fold]]]))
    MGHLookup.test <- MGHLookup[MGH_CV[[fold]],]
    
    # make the importance matrices to hold values
    MGH.impval <- matrix(NA, nrow = length(featuresMGH), ncol = reps, dimnames = list(colnames(MGHCounts.test), paste0('rep', 1:reps)))
    
    Rand.predict <- matrix(NA, nrow = nrow(MGHRand.test), ncol = reps, dimnames = list(rownames(MGHRand.test), paste0('rep', 1:reps)))
    MGH.predict <- matrix(NA, nrow = nrow(MGHCounts.test), ncol = reps, dimnames = list(rownames(MGHCounts.test), paste0('rep', 1:reps)))
    
    # run RF reps (10) times to make sure that there isn't any seeding variability
    for(rep in 1:reps){
      MGHRes <- RandomForestClassifier(trainData = MGHCounts.train, 
                                       testCounts = MGHCounts.test,
                                       testLookup = MGHLookup.test)
      MGHVotes[rownames(MGHRes[[1]]), rep, fold] <- MGHRes[[1]][,2]
      MGH.impval[,rep] <- MGHRes[[2]][match(rownames(MGH.impval), names(MGHRes[[2]]))]
      
      RandRes <- RandomForestClassifier(trainData = MGHRand.train, 
                                        testCounts = MGHRand.test,
                                        testLookup = MGHLookup.test)
      
      # the second column of the first object in results holds the averaged votes for the subsetted samples so stick that into our mega votes matrix
      RandVotes[rownames(RandRes[[1]]), rep, fold] <- RandRes[[1]][,2]
      MGHVotes[rownames(MGHRes[[1]]), rep, fold] <- MGHRes[[1]][,2]
      
      # save the good model's importances
      MGH.impval[,rep] <- MGHRes[[2]][match(rownames(MGH.impval), names(MGHRes[[2]]))]
    }
    MGHImps[[fold]] <- MGH.impval
    print(paste('Finished MGH Fold', fold))
    graphics.off()
  }
  
  # Make the full ROC object and get the AUC from it to save.
  pdf(paste0("FinalAUCs_MGH_", date, '.pdf'))
  ROC_MGH <- ROC_plotter(votes = MGHVotes,
                         key = MGHLookup$diagnosis,
                         title = 'MGH RF Model AUC',
                         plot = T)
  ROC_Rand <- ROC_plotter(votes = RandVotes,
                          key = MGHLookup$diagnosis,
                          title = 'MGH RF Rand AUC',
                          plot = T)
  dev.off()
  saveRDS(MGHImps, file = paste0('MGH_Imps_', date, '.RDS'))
  saveRDS(list(Diagnostic = list(votes = MGHVotes, key = MGHLookup$diagnosis, features = featuresMGH),
               Random = list(votes = RandVotes, key = MGHLookup$diagnosis, features = featuresRand)),
          file = paste0('MGH_VotesSet_', date ,'.RDS'))
  
  output <- c(MGH_AUC = 0, Rand_AUC = 0) 
  output[1] <- ROC_MGH$modelROC$auc
  output[2] <- ROC_Rand$modelROC$auc
  return(output)
}