GenAUC_KMU <- function(KMUCounts, KMULookup, features, folds, reps){
  
  KMU_votes <- array(NA, dim = c(ncol(KMUCounts), folds, reps), dimnames = list(colnames(KMUCounts), paste0('fold', 1:folds), paste0('rep', 1:reps)))
  KMU_Randvotes <- array(NA, dim = c(ncol(KMUCounts), folds, reps), dimnames = list(colnames(KMUCounts), paste0('fold', 1:folds), paste0('rep', 1:reps)))
  KMU_imps <- vector(mode = 'list', length = folds)
  
  # get the K means splits and plot out
  KMU_CV <- StratifiedKmeansSplit(KMULookup, folds)
  
  #pdf(paste0('Plots/KMeansDistributionValidation_KMU_', date, '.pdf'))
  #plotCVDistributions(KMULookup, KMU_CV, 'KMU Crossvalidation Splits')
  #dev.off()
  
  # Make holder for the predictions and the importances
  KMUVotes <- array(NA, dim = c(length(KMULookup$sample_name), reps, folds), dimnames = list(KMULookup$sample_name, paste0('rep',1:reps), paste0('fold', 1:folds)))
  RandVotes <- array(NA, dim = c(length(KMULookup$sample_name), reps, folds), dimnames = list(KMULookup$sample_name, paste0('rep',1:reps), paste0('fold', 1:folds)))
  KMUImps <- vector(mode = 'list', length = folds)
  
  #pdf(paste0('Plots/FeatQCPlots_KMUFinal_', date, '.pdf'))
  for(fold in 1:folds){
    
    # conduct feature selection on the training set
    RandFeatures <- expressionMatchRand(genes = features,
                                        counts = KMUCounts,
                                        lookup = KMULookup,
                                        subtitle = paste('Fold', fold),
                                        plot = F)
    
    # break to training set and subset to selected features
    KMUCounts.train <- data.frame(t(KMUCounts[rownames(KMUCounts)%in%features, -KMU_CV[[fold]]]))
    KMUCounts.train$diagnosis <- KMULookup$diagnosis[-KMU_CV[[fold]]]
    RandCounts.train <- data.frame(t(KMUCounts[rownames(KMUCounts)%in%RandFeatures, -KMU_CV[[fold]]]))
    RandCounts.train$diagnosis <- KMULookup$diagnosis[-KMU_CV[[fold]]]
    
    # Make Test set and trim to feature set
    RandCounts.test <- data.frame(t(KMUCounts[rownames(KMUCounts)%in%RandFeatures, KMU_CV[[fold]]]))
    KMUCounts.test <- data.frame(t(KMUCounts[rownames(KMUCounts)%in%features, KMU_CV[[fold]]]))
    KMULookup.test <- KMULookup[KMU_CV[[fold]],]
    
    # make the importance matrices to hold values
    KMU.impval <- matrix(NA, nrow = length(features), ncol = reps, dimnames = list(colnames(KMUCounts.test), paste0('rep', 1:reps)))
    
    KMU.predict <- matrix(NA, nrow = nrow(KMUCounts.test), ncol = reps, dimnames = list(rownames(KMUCounts.test), paste0('rep', 1:reps)))
    Rand.predict <- matrix(NA, nrow = nrow(KMUCounts.test), ncol = reps, dimnames = list(rownames(KMUCounts.test), paste0('rep', 1:reps)))
    
    # run RF reps (10) times to make sure that there isn't any seeding variability
    for(rep in 1:reps){
      KMURes <- RandomForestClassifier(trainData = KMUCounts.train, 
                                       testCounts = KMUCounts.test,
                                       testLookup = KMULookup.test)
      KMUVotes[rownames(KMURes[[1]]), rep, fold] <- KMURes[[1]][,2]
      KMU.impval[,rep] <- KMURes[[2]][match(rownames(KMU.impval), names(KMURes[[2]]))]
      
      RandRes <- RandomForestClassifier(trainData = RandCounts.train, 
                                        testCounts = RandCounts.test,
                                        testLookup = KMULookup.test)
      RandVotes[rownames(KMURes[[1]]), rep, fold] <- RandRes[[1]][,2]
    }
    KMUImps[[fold]] <- KMU.impval
    print(paste('Finished KMU Fold', fold))
  }

  #Generate AUCs
  pdf(paste0("FinalAUCs_KMU_", date, '.pdf'))
  ROC_KMU <- ROC_plotter(votes = KMUVotes,
                         key = KMULookup$diagnosis,
                         title = paste0('KMU RF Model AUC'),
                         plot = T)
  ROC_Rand <- ROC_plotter(votes = RandVotes,
                          key = KMULookup$diagnosis,
                          title = 'KMU RF Random Control AUC',
                          plot = T)
  dev.off()
  
  saveRDS(KMUImps, file = paste0('KMU_Imps_', date, '.RDS'))
  saveRDS(list(Diagnostic = list(votes = KMUVotes, key = KMULookup$diagnosis, features = features),
               Random = list(votes = RandVotes, key = KMULookup$diagnosis, features = RandFeatures)),
          file = paste0('KMU_VotesSet_', date ,'.RDS'))
  # Save in output variable and return
  output <- c(KMU_AUC=0, Rand_AUC=0)
  output[1] <- ROC_KMU$modelROC$auc
  output[2] <- ROC_Rand$modelROC$auc
  return(output)
}