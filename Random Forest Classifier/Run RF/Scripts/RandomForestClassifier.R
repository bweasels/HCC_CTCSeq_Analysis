## Function to run random forest
RandomForestClassifier <- function(trainData, testCounts, testLookup) {
  require(ranger)
  rf.model <- ranger(diagnosis ~ . ,
                     data = trainData,
                     num.trees = 500,
                     write.forest = TRUE,
                     importance = 'impurity')
  pred.model <- predict(rf.model,
                      data = testCounts,
                      predict.all = TRUE)
  
  # convert factor levels to percent voted 1 or 2
  factor_levels <- rf.model$forest$levels
  votes <- pred.model$predictions
  votes <- apply(votes, 2, function(x) factor_levels[x]) #replaces the 1 and 2 votes with cancer/not cancer depending on the RF's order
  
  # I don't trust that 1 is 'not cancer' and 2 is always going to be 'cancer', so I rely on the forest's own ordering
  predictions <- array(NA, dim = c(length(testLookup$diagnosis), 2), 
                       dimnames = list(testLookup$sample_name, factor_levels))
  
  # sum up the votes in either direction and spit out the percent the model believed that it was cancer/not cancer
  predictions[,1] <- apply(votes, 1, function(x) sum(grepl(factor_levels[1], x))/length(x))      
  predictions[,2] <- apply(votes, 1, function(x) sum(grepl(factor_levels[2], x))/length(x))

  return(list(predictions, importance(rf.model)))
}