##Function to deplete RBC counts

depleteRBC <- function(counts){
  depletePanel <- read.csv('Molecular_Depletion_Panel.csv')
  counts <- counts[!rownames(counts)%in%depletePanel$Gene.ID,]
  return(counts)
}