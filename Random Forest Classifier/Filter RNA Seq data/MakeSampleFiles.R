# Rewrite CTC-Seq to validate our LFC method
# Process:
# Step 1: Subset for MGH and KMU data with p val < 0.1
# Step 2: Subset the data for immune genes
# Step 3: Select for genes that are expressed in the same direction
#####################LOAD KMU AND MGH COUNTS#####################
setwd('/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/CTC-Seq/input/')
source('depleteRBC.R')
lookup.df <- read.csv('HCC Sample Lookup.csv', stringsAsFactors = F)
counts.df <- read.csv('HCC CLD RawCounts copy.csv', header = TRUE, row.names = 1)

#Confirm that the RNASeq is in same order as lookup table
if(sum(colnames(counts.df)!=lookup.df$sample_id)>0){ 
  stop('WARNING: Lookup table and RNASeq data tables are not similarly ordered.')}

#########TRIM OUT EVERYTHING EXCEPT CANCER AND ALMOST CANCER (HCC/CLD)
#See which samples are CLD (chronic liver disease -or- pre cancer) and which samples are HCC (hepatic cell carcinoma -or- actually cancer)
HCC_CLD <- lookup.df$class == "CLD" | lookup.df$class == "HCC" & lookup.df$active == TRUE & lookup.df$repeat_draw == FALSE 
counts.df <- counts.df[,HCC_CLD]
lookup.df <- lookup.df[HCC_CLD,]

#get the total number of counts from each sample and remove lowQ samples
colSums <- colSums(counts.df)

counts.df <- counts.df[,colSums>100000]
MGHlookup.df <- lookup.df[colSums>100000,]

#deplete RBCs
MGHcounts.df <- depleteRBC(counts.df)
######################LOAD KMU COUNTS############################
threshold.ImmuneSig <- 1.5 # New immunePanel Threshold
threshold.Liver <- 0.5 # liver panel threshold
#load the data
counts.df <- read.csv('2019-08-12_KMUCounts_star_genes_erc.csv', row.names = 1)
colnames(counts.df) <- gsub('^X(.*)', '\\1', colnames(counts.df))
key <- read.csv('2019-08-12_KMUMetadata.csv', row.names = 1)

# Remove the sample that Irun asked for (KT00105K1 from HCC group)
key = key[!grepl('KT00105K1', key$Sample.ID), ]
counts.df = counts.df[,!grepl('KT00105K1', colnames(counts.df))]


#confirm that the metadata is in alignment with the counts table
if(sum(colnames(counts.df)!=key$Sample.ID)>0){
  stop('WARNING: Lookup table and RNASeq data tables are not similarly ordered.')
}

#create lookup df
patientStatus <- ifelse(grepl('CLD|LC', key$Patient.Status), 'CLD', 'HCC')
lookup.df <- data.frame(sample_name = key$Sample.ID, diagnosis = patientStatus)
lookup.df$diagnosis[grep('KT00076K1', lookup.df$sample_name)] <- 'CLD'
lookup.df$diagnosis[grep('KT00158K1', lookup.df$sample_name)] <- 'HCC'

key$diagnosisFinal = lookup.df$diagnosis


#get the total number of counts from each sample and remove lowQ samples
colSums <- colSums(counts.df)

counts.df <- counts.df[,colSums>100000]
lookup.df <- lookup.df[colSums>100000,]
key <- key[colSums>100000,]

#deplete RBCs
counts.df <- depleteRBC(counts.df)

#generate immunePanel Score
immunePanel <- unlist(read.csv('immunePanelGenes_2019-12-08.csv', stringsAsFactors = F))
liverPanel <- c('FGB', 'FGG', 'ALB', 'AFP', 'FABP1', 'GPC3', 'RBP4', 'RBP4', 'AHSG', 'APOH', 'TF')
rpm.df <- apply(counts.df, 2, function(x) x*1e6/sum(x))
#rpm.df <- log10(rpm.df+1)
hybridScore <- log10(colSums(rpm.df[rownames(rpm.df) %in% immunePanel,])+1)
liverScore <- log10(colSums(rpm.df[rownames(rpm.df) %in% liverPanel,])+1)

#remove all samples with a liver score < 0.5 (less than 3 RPM for the panel)
KMUhybridCounts.df <- counts.df[,hybridScore>threshold.ImmuneSig]
KMUhybridLookup.df <- lookup.df[hybridScore>threshold.ImmuneSig,]
KMUHybridKey <- key[hybridScore>threshold.ImmuneSig,]

KMUliverCounts.df <- counts.df[,liverScore>threshold.Liver]
KMUliverLookup.df <- lookup.df[liverScore>threshold.Liver,]

#############################END LOAD COUNTS######################

date <- Sys.Date()
setwd('/OneDrive/MGHWork/Dropbox (Partners HealthCare)/Dropbox (Partners HealthCare)/Weekly Analyses/05-26-20/Make Sample Files for validation/')

KMUhybridCounts.df <- KMUhybridCounts.df[rownames(KMUhybridCounts.df)%in%rownames(MGHcounts.df),]
KMUliverCounts.df <- KMUliverCounts.df[rownames(KMUliverCounts.df)%in%rownames(MGHcounts.df),]
MGHcounts.df <- MGHcounts.df[match(rownames(KMUhybridCounts.df), rownames(MGHcounts.df)),]

data.hybrid <- list(list(MGHcounts.df, MGHlookup.df), list(KMUhybridCounts.df, KMUhybridLookup.df))
data.liver <- list(list(MGHcounts.df, MGHlookup.df), list(KMUliverCounts.df, KMUliverLookup.df))
saveRDS(data.hybrid, paste0('KMU_MGHDataHybrid_', date,'.RDS'))
saveRDS(data.liver, paste0('KMU_MGHDataLiver_', date,'.RDS'))

write.csv(KMUhybridCounts.df, paste0('KMUCounts_filtered_', date, '.csv'))
write.csv(KMUHybridKey, paste0('KMUmetadata_filtered_', date, '.csv'))
write.csv(MGHcounts.df, paste0('MGHCounts_filtered_', date, '.csv'))
write.csv(MGHlookup.df, paste0('MGHmetadata_filtered_', date, '.csv'))
