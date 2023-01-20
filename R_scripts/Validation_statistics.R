# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part II: Validation_statistics


##########R Pipeline - Part II: Validation_statistics##########
print(Sys.time())
start_time <- Sys.time()
print("R pipeline - Part II: Validation_statistics - start!")
# Part II: Validation_statistics


## check projection1 present
if(AMOUNT_OF_PROJECTIONS < 1){
  stop("ERROR: Part II: validation statistics stopped because no Projection1 present on which calculations are performed.")
}


##adjustmetns config that are fixed => use DEFAULT_PC_AMOUNT to be sure correct hotellings ellipse 
DEFAULT_PC_AMOUNT <- "calculates PCs depending on sample amount (n-1)"
FIXED_PC_AMOUNT <- "set amount of PCs to calculate to certain integer, ! less than amount of samples present in projections"

#PCA amount of PCs in projection:
FIXED_PC_AMOUNT <- 20   #integer, less than amount of samples present in projections
PC_AMOUNT <- DEFAULT_PC_AMOUNT

#style
SIZE_POINTS <- 3   #size samplepoints on scoreplots PCA. default: 5, big dataset: 3


## data_loading
setwd(path_data_in)
COLLUMN_NR_COMPARISON1 <- COLLUMN_NR_LAST_BEFORE_COMPARISONS + 1        #column number of the column 'Comparison1' for pairwise comparison (OPLSDA)
COLLUMN_NR_MULTIPLE_COMPARISON1 <- COLLUMN_NR_COMPARISON1 + AMOUNT_OF_COMPARISONS         #column number of the column 'MultipleComparison1' for multiple comparison (PLSDA, LIMMA)
COLLUMN_NR_PROJECTION1 <- COLLUMN_NR_MULTIPLE_COMPARISON1 + AMOUNT_OF_MULTIPLE_COMPARISONS
COLLUMN_NR_START_VARIABLES <- COLLUMN_NR_PROJECTION1 + AMOUNT_OF_PROJECTIONS

sampleMetadata <- load_sampleMetadata(INPUT_SAMPLES)
#Add "X" to sampleNames
sampleMetadata[,1] <- paste(replicate(nrow(sampleMetadata),"X"), sampleMetadata[,1], sep="")
check_nrow_SM <- nrow(sampleMetadata)

if (INPUT_VARIABLES == VARIABLEMETADATA_FROM_PIPELINE){
  INPUT_VARIABLES <- paste(name_project, '_variableMetadata.txt', sep="")
  if(exists("COLLUMN_NR_START_SAMPLES") == FALSE){ 		#if after merge, will be given value 21, so do not change
    COLLUMN_NR_START_SAMPLES <- 20  #always 20 (auto and manual must be same format)
  }
}
variableMetadata <- load_variableMetadata(INPUT_VARIABLES)
#note: automatically adds X to samples at this point


## set directory to output
setwd(path_data_out)


## merge variables from variableMetadata into the sampleMetadata
### make sampleMatrix
suppressMessages(library(data.table))
variableMetadata_from_start_samples <- subset(variableMetadata, select = -c(2:(COLLUMN_NR_START_SAMPLES-1))) #remove info exept CompID
sampleMatrix <- transpose(variableMetadata_from_start_samples)
colnames(sampleMatrix) <- variableMetadata_from_start_samples[ ,1]
sampleMatrix <- as.data.frame(sapply(sampleMatrix, as.numeric))
sampleMatrix$SampleName <- colnames(variableMetadata_from_start_samples)
sampleMatrix <- sampleMatrix[-1,] #remove first row (colnames compids + samplename; so 1 extra col than no of variables), is captured in colnames
#write_dataframe_as_txt_file(sampleMatrix, 'sampleMatrix.txt')

#merge variable intensities from compIDs to correct sample (if order not same), output order is sorted by samplename
sampleMetadata <- merge_accoding_to_SampleName(sampleMetadata, sampleMatrix)
write_dataframe_as_txt_file(sampleMetadata, 'sampleMetadata_variableMatrix_merged.txt')


##Check sampleMatrx 
#see if same number of rows (amount of sampleNames) as when loaded prev, so merge was succesfull
if(nrow(sampleMetadata) != check_nrow_SM){
  stop("ERROR: Part II: multivariate analysis stopped because SampleMetadata and VariableMetadata are incompatible to merge into correct sampleMatrix.")
} 




###### MUTUAL PARTS (projection1) ###### 
#projection1 contains groups listed in comps
projection <- 1

print(paste("start validation statistics for Projection ", projection, ".", sep=""))

### remove non-essential data for comparison
#remove samples that are not included in comp (Na instead of 0/1)
nsamples_metadata <- sampleMetadata #also QC can be in PCA
comp_ <- nsamples_metadata[,COLLUMN_NR_PROJECTION1+projection-1]
stopifnot(class(comp_)=='integer') #need {0,1,2, ...} arguments, rest is NA and not ""!
samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    

#remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
#samples_matrix_comp_no0 <- samples_matrix_comp[, !apply(samples_matrix_comp == 0, 2, all)]  #remove 0
samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x)))
samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))]

comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_PROJECTION1+projection-1])

############ 



######### RESULTS ALL (projection1) ######### 
report_outliers_pcas <- NULL

## tic count of all samples togheter plot
#samplenames <- sampleMetadata$SampleName
TIC_per_sample <- apply(sampleMetadata[,COLLUMN_NR_START_VARIABLES:length(sampleMetadata)], 1, function(x){t(sum(x)) }) #iterate over rows =1 (samples), but saves results as column!!! 

png(paste(name_project,'_3-2TICsignal_per_sample.png', sep=""), width=7, height=5, units="in", res=150)
plot(TIC_per_sample)
dev.off()



### PCA of all samples in projection1

#for PCA's: load needed
suppressMessages(library(pcaMethods))
if(PC_AMOUNT == DEFAULT_PC_AMOUNT){
  number_samples <- nrow(samples_matrix_comp_no0)     #Amount of samples
  number_samples <- number_samples-1
}
if(PC_AMOUNT == FIXED_PC_AMOUNT){
  number_samples = FIXED_PC_AMOUNT
  if(number_samples >= nrow(samples_matrix_comp_no0)){
    stop("ERROR: Part II: multivariate analysis stopped because number of PCs exceeds number of samples-1.")
  }
}

#Log transformation and Pareto scaling for pca (wo funtions because matrix here ipv metadata)
lognormalized_samples_matrix_comp_no0 <-   sapply(samples_matrix_comp_no0, log1p) 
x.centered <- lognormalized_samples_matrix_comp_no0
x.centered <-  apply(lognormalized_samples_matrix_comp_no0, 2, function(x) x - mean(x))
scalednormalizedsamples_matrix_comp_no0 <- x.centered
scalednormalizedsamples_matrix_comp_no0 <- data.frame(apply(x.centered, 2, function(x) x/sqrt(sd(x))))

#remove NaN after scaling (again)
#max(scalednormalizedsamples_matrix_comp_no0)
scalednormalizedsamples_matrix_comp_no0 <- data.frame(sapply(scalednormalizedsamples_matrix_comp_no0, function(x) ifelse(is.nan(x), NA, x)))
scalednormalizedsamples_matrix_comp_no0 <- scalednormalizedsamples_matrix_comp_no0[,which(unlist(lapply(scalednormalizedsamples_matrix_comp_no0, function(x)!all(is.na(x)))))]

#pca score plot
nPCsopt <- 2
PCAoptnumber <- pca(scalednormalizedsamples_matrix_comp_no0,nPcs=nPCsopt,scale="none",completeObs=TRUE,subset=NULL,cv="q2",center=FALSE)
PCAscores_comp <- as.data.frame(PCAoptnumber@scores)

HotellingellipsePCA_comp <- Hotellingellipse(number_samples,PCAscores_comp)  #Calculate Hotelling's T2 ellipse for the PCA
dev.off()
PCAplot_comp <- plot_pcascores(PCAscores_comp, comp)

#add all labels
#pca_with_labels <- PCAplot_comp + geom_text(data=PCAscores_comp, aes(x=PC1, y=PC2, label = samples_metadata_comp$SampleName), hjust=-0.05, vjust=-0.05, color="gray6", alpha = 0.5, size = 3)
#not use, see below, add only if label is outlier


#pca score plot weight values
score_weights_pca_comp <- PCAscores_comp[ ,1:2] #PC1 and PC2 are kept
rownames(score_weights_pca_comp) <- samples_metadata_comp$SampleName
#name_df <- paste(name_project, '_score_weights_projection_', projection, '.txt', sep="")
#write_matrix_as_txt_file(score_weights_pca_comp, name_df) 
#plot(score_weights_pca_comp)

#calculate outlier or not
#https://www.mathwarehouse.com/ellipse/equation-of-ellipse.php#equationOfEllipse
h <- mean(HotellingellipsePCA_comp[,1])
k <- mean(HotellingellipsePCA_comp[,2])
a <- abs(max(HotellingellipsePCA_comp[,1])-min(HotellingellipsePCA_comp[,1]))/2
b <- abs(max(HotellingellipsePCA_comp[,2])-min(HotellingellipsePCA_comp[,2]))/2

score_weights_pca_comp$SampleName <- rownames(score_weights_pca_comp)
score_weights_pca_comp$test_outlier <- ( ((score_weights_pca_comp[,1] - h)^2 / a^2) + ((score_weights_pca_comp[,2] - k)^2 / b^2) )

#plot pca with label of outliers
pca_with_labels <- PCAplot_comp + 
  geom_text(data=score_weights_pca_comp, 
            aes(x=PC1, y=PC2, label = ifelse(score_weights_pca_comp$test_outlier>1,as.character(score_weights_pca_comp$SampleName), '')), 
            hjust=-0.05, vjust=-0.05, color="gray6", alpha = 0.5, size = 3)

png(paste(name_project, "_3-3PCA_scoreplot " , projection, "ALL.png", sep=""), width=7, height=5, units="in", res=150)
plot(pca_with_labels)
dev.off()

#table list ouliers
samplenames_outliers <- score_weights_pca_comp[score_weights_pca_comp$test_outlier > 1,]
samplenames_outliers <- as.character(samplenames_outliers$SampleName)
amount_of_outliers <- length(samplenames_outliers)
#nice print character(0):
if(identical(samplenames_outliers, character(0))){
  samplenames_outliers <- ""
}

report_outliers_pca <- NULL
report_outliers_pca$group_name <- 'Projection 1'
report_outliers_pca$amount_of_outliers <- amount_of_outliers
report_outliers_pca$samplenames_outliers <- samplenames_outliers

report_outliers_pcas <- rbind(report_outliers_pcas, report_outliers_pca)
############ 



###### RESULTS PER GROUP DEFINED IN PROJECTION1 ###### 
report_setup_groups <- NULL
report_features_groups <- NULL
report_tic_groups <- NULL
report_cv_raws <- NULL
report_cv_tics <- NULL

for (group in levels(comp)){
  #print(group)
  #group <- "1"
  name_group <- paste0("Group ", group)
  groupnr <- as.integer(group)
  sampleMetadata <- sampleMetadata[!is.na(sampleMetadata$Projection1), ]    #remove samples not in projection1: NA
  sm_group <- sampleMetadata[sampleMetadata$Projection1 == groupnr & sampleMetadata$Type != "Blank",]
  amount_of_samples_group <- nrow(sm_group)
  amount_of_samples_group_per_sample <- rowSums(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)] != 0) #not zero feat/sample
  mean_amount_of_feat_group <- mean(amount_of_samples_group_per_sample)
  sd_amount_of_feat_group <- sd(amount_of_samples_group_per_sample)

  
  #samplenames_group <- sm_group$SampleName
  #since issue newspaces (w >13 samplenames), but no issue in list itself. make df again from report_setup_groups
  i = 1
  samplenames_group = NULL
  for(i in 1:length(sm_group$SampleName)){
    samplenames_group <- paste(samplenames_group, sm_group$SampleName[i], sep=", ")
    i = i + 1
  }
  samplenames_group <- substr(samplenames_group, 3, nchar(samplenames_group))
  
  tic_group <- apply(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)], 1, function(x){t(sum(x)) }) #iterate over rows =1 (samples), but saves results as column!!! 
  tic_mean <- mean(tic_group)
  sd_tic <- sd(tic_group)
  rsd_tic_perc <- (sd_tic / tic_mean *100)
  
  
  #write to tables
  #1st df 1setup
  report_setup_group <- data.frame(name_group, amount_of_samples_group, samplenames_group)
  
  report_setup_groups <- rbind(report_setup_groups, report_setup_group)
  
  
  
  #2nd df 3-1Numberoffeatures
  report_features_group <- NULL
  report_features_group$name_group <- name_group
  report_features_group$mean_amount_of_feat_group <- format(round(mean_amount_of_feat_group,digits=2),nsmall=2)
  report_features_group$sd_amount_of_feat_group <- format(round(sd_amount_of_feat_group,digits=2),nsmall=2)
  
  report_features_groups <- rbind(report_features_groups, report_features_group)
  
  
  #3rd table 3-2TICsignal
  report_tic_group <- NULL
  report_tic_group$name_group <- name_group
  report_tic_group$tic_mean <- format(round(tic_mean,digits=2),nsmall=2)
  report_tic_group$sd_tic <- format(round(sd_tic,digits=2),nsmall=2)
  report_tic_group$rsd_tic_perc <- format(round(rsd_tic_perc,digits=2),nsmall=2)
  
  report_tic_groups <- rbind(report_tic_groups, report_tic_group)
  
  
  
  #for PCA's: load needed
  suppressMessages(library(pcaMethods))
  if(PC_AMOUNT == DEFAULT_PC_AMOUNT){
    number_samples <- amount_of_samples_group    #Amount of samples
    number_samples <- number_samples-1
  }
  if(PC_AMOUNT == FIXED_PC_AMOUNT){
    number_samples = FIXED_PC_AMOUNT
    if(number_samples >= amount_of_samples_group){
      stop("ERROR: Part II: multivariate analysis stopped because number of PCs exceeds number of samples-1.")
    }
  }
  
  #Log transformation and Pareto scaling for pca
  lognormalized_sm_group <- log_transformation(sm_group) 
  scalednormalized_sm_group <- Pareto_scaling(lognormalized_sm_group)
  
  #remove NaN after scaling (again), as matrix
  #max(scalednormalized_sm_group[,COLLUMN_NR_START_VARIABLES:length(scalednormalized_sm_group)])
  scalednormalized_sm_group <- from_df_to_matrix(scalednormalized_sm_group)
  scalednormalized_sm_group <- data.frame(sapply(scalednormalized_sm_group, function(x) ifelse(is.nan(x), NA, x)))
  scalednormalized_sm_group <- scalednormalized_sm_group[,which(unlist(lapply(scalednormalized_sm_group, function(x)!all(is.na(x)))))]
  
  #PCA per group
  nPCsopt <- 2
  PCAoptnumber <- pca(scalednormalized_sm_group,nPcs=nPCsopt,scale="none",completeObs=TRUE,subset=NULL,cv="q2",center=FALSE)
  PCAscores_comp <- as.data.frame(PCAoptnumber@scores)

  #try({ #error if 3 reps, since n-2 needed to calc and already n-1 = 0 ==> not solve via try: previouw ellipse use, so wrong! keep error in code
  HotellingellipsePCA_comp <- Hotellingellipse(number_samples,PCAscores_comp)  #Calculate Hotelling's T2 ellipse for the PCA
  dev.off()
  #})
  PCAplot_comp <- plot_pcascores(PCAscores_comp, group)

  #add all labels
  #pca_with_labels <- PCAplot_comp + geom_text(data=PCAscores_comp, aes(x=PC1, y=PC2, label = samplenames_group), hjust=-0.05, vjust=-0.05, color="gray6", alpha = 0.5, size = 3)
  #see below, only ouliers on plot + png write  

  
  #pca score plot weight values
  score_weights_pca_comp <- PCAscores_comp[ ,1:2] #PC1 and PC2 are kept
  rownames(score_weights_pca_comp) <- sm_group$SampleName
  #name_df <- paste(name_project, '_score_weights_projection_', projection, '.txt', sep="")
  #write_matrix_as_txt_file(score_weights_pca_comp, name_df) 
  #plot(score_weights_pca_comp)
  
  #calculate outlier or not
  #https://www.mathwarehouse.com/ellipse/equation-of-ellipse.php#equationOfEllipse
  h <- mean(HotellingellipsePCA_comp[,1])
  k <- mean(HotellingellipsePCA_comp[,2])
  a <- abs(max(HotellingellipsePCA_comp[,1])-min(HotellingellipsePCA_comp[,1]))/2
  b <- abs(max(HotellingellipsePCA_comp[,2])-min(HotellingellipsePCA_comp[,2]))/2
  
  score_weights_pca_comp$SampleName <- rownames(score_weights_pca_comp)
  score_weights_pca_comp$test_outlier <- ( ((score_weights_pca_comp[,1] - h)^2 / a^2) + ((score_weights_pca_comp[,2] - k)^2 / b^2) )
  
  
  #plot with labels IF outlier
  pca_with_labels <- PCAplot_comp + 
    geom_text(data=score_weights_pca_comp, 
              aes(x=PC1, y=PC2, label = ifelse(score_weights_pca_comp$test_outlier>1,as.character(score_weights_pca_comp$SampleName), '')), 
              hjust=-0.05, vjust=-0.05, color="gray6", alpha = 0.5, size = 3)
  
  png(paste(name_project, "_3-3PCA_scoreplot" , projection, "group", groupnr, ".png", sep=""), width=7, height=5, units="in", res=150)
  plot(pca_with_labels)
  dev.off()
  
  
  #table list ouliers
  samplenames_outliers <- score_weights_pca_comp[score_weights_pca_comp$test_outlier > 1,]
  samplenames_outliers <- as.character(samplenames_outliers$SampleName)
  amount_of_outliers <- length(samplenames_outliers)
  #nice print character(0):
  if(identical(samplenames_outliers, character(0))){
    samplenames_outliers <- ""
  }
  
  report_outliers_pca <- NULL
  report_outliers_pca$group_name <- name_group
  report_outliers_pca$amount_of_outliers <- amount_of_outliers
  report_outliers_pca$samplenames_outliers <- samplenames_outliers
  
  report_outliers_pcas <- rbind(report_outliers_pcas, report_outliers_pca)
  
  ###save sm_group (=SM for part group) so no issue 80% filter, tic,...
  saved_sm_group <- sm_group

  
  ### CV = 3-4Featurerepeatability
  CV_threshold <- CV_THRESHOLD #percent vb 30
  
  ##only retain features if present in 80% of df
  perc_present <- colSums(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)] != 0)/amount_of_samples_group
  matrix_with_perc_present <- data.frame(x=colnames(sm_group)[COLLUMN_NR_START_VARIABLES:length(sm_group)],
                                         y=perc_present)
  names(matrix_with_perc_present) <- c("CompID", "perc_present")
  retained_variables_group <- matrix_with_perc_present[(matrix_with_perc_present[,2]>=0.8),]
  retained_variables_group <- retained_variables_group[,1] #compIDs of good ones
  
  all_matrix <- from_df_to_matrix(sm_group)
  all_matrix <- as.data.frame(all_matrix)
  
  matrix_filter_group <- all_matrix[, names(all_matrix) %in% (retained_variables_group)]
  
  sampleMetadata_until_start_variables <- subset(sm_group, select = c(1:(COLLUMN_NR_START_VARIABLES-1)))
  sampleMetadata_filter_group <- cbind(sampleMetadata_until_start_variables, matrix_filter_group)
  #write_dataframe_as_txt_file(sampleMetadata_filter_QCpool, 'sampleMetadata_filter_QCpool.txt')
  
  sm_group <- sampleMetadata_filter_group
  
  
  ## raw CV
  CV_values <- apply(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)], 2, function(x) sd(x)/mean(x)*100) #in percent
  matrix_with_CV_values <- data.frame(x=colnames(sm_group)[COLLUMN_NR_START_VARIABLES:length(sm_group)], y=CV_values)
  names(matrix_with_CV_values) <- c("CompID", "CV_values")
  #write_dataframe_as_txt_file(matrix_with_CV_values, paste0(name_project,'_CV_values_raw.txt'))
  
  #skip if amount_of_samples_group <= 5
  if(amount_of_samples_group <= 5){
    perc_of_feat_under_CVthesh <- "not calculated"
    amount_of_feat_under_CVthesh <- "not calculated"
    amount_of_feat_retained_80perc <- "not calculated"
  }
  if(amount_of_samples_group > 5){
    #features under CV threshold
    amount_of_feat_under_CVthesh <- nrow(matrix_with_CV_values[matrix_with_CV_values$CV_values <= CV_threshold,])
    #amount_of_feat_group <- ncol(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)])   
    #no, using retained_variables_group from 80% calc:
    amount_of_feat_retained_80perc <- length(retained_variables_group)
    perc_of_feat_under_CVthesh <- format(round((amount_of_feat_under_CVthesh/amount_of_feat_retained_80perc *100),digits=2),nsmall=2)
  }
  
  report_cv_raw <- NULL
  report_cv_raw$name_group <- name_group
  report_cv_raw$amount_of_samples_group <- amount_of_samples_group
  report_cv_raw$amount_of_feat_under_CVthesh <- amount_of_feat_under_CVthesh
  report_cv_raw$amount_of_feat_retained_80perc <- amount_of_feat_retained_80perc
  report_cv_raw$perc_of_feat_under_CVthesh <- perc_of_feat_under_CVthesh

  #report_cv_raw$histogram <- ""
  report_cv_raws <- rbind(report_cv_raws, report_cv_raw)
  
  ## plot per group cumm CV
  #plot(ecdf(matrix_with_CV_values[,2]))
  #h <- hist(matrix_with_CV_values[,2])
  
  #https://stackoverflow.com/questions/29289046/r-ecdf-over-histogram
  max_count <- max(matrix_with_CV_values[,2])
  breaks <- round(max_count/10, 0) #so show intervals 10-20
  
  png(paste(name_project, "_3-4Featurerepeatability_raw" , projection, "group", groupnr, ".png", sep=""), width=7, height=5, units="in", res=150)
  h <- hist(matrix_with_CV_values[,2], col = "dodgerblue3", xlab = "Bin", xlim = c(0,100), breaks=breaks, main=paste0("Group ", group)) #"#4393C3"
  par(new = T)
  ec <- ecdf(matrix_with_CV_values[,2])
  plot(x = h$mids, y=ec(h$mids)*max(h$counts), col = rgb(0,0,0,alpha=0), axes=F, xlab=NA, ylab=NA)
  lines(x = h$mids, y=ec(h$mids)*max(h$counts), col ='red')
  axis(4, at=seq(from = 0, to = max(h$counts), length.out = 11), labels=seq(0, 100, 10), col = 'red', col.axis = 'red')
  mtext(side = 4, line = -2, 'cummulative %', col = 'red')
  dev.off()
  
  #calc cummulative_perc
  total_count <- sum(unlist(h$counts))
  factor <- total_count/100
  count_perc <- unlist(h$counts) / factor
  cummulative_perc <- cumsum(count_perc)

  #histogram cummulat perc table
  hist_table_group <- NULL
  hist_table_group$bin <- unlist(h$breaks[0:(length(h$breaks)-1)]) #bins, only start show, so wo last elem
  hist_table_group$interval <- paste0(unlist(h$breaks[0:(length(h$breaks)-1)]), "-", (unlist(h$breaks[0:(length(h$breaks)-1)])+10))
  hist_table_group$frequency <- unlist(h$counts)
  hist_table_group$cummulative_perc <- format(round(cummulative_perc, 2),nsmall=2)
  hist_table_group <- as.data.frame(hist_table_group)
  
  #cut-off only <100 in tabel hist keep (if present)
  if(length(which(hist_table_group$bin == 90)) > 0)({
    colnr_bin90 <- which(hist_table_group$bin == 90)
    hist_table_group <- hist_table_group[1:colnr_bin90,]
  })

  #write resutls 3-4Featurerepeatability_raw
  write_dataframe_as_txt_file(hist_table_group, paste0(name_project,'_3-4Featurerepeatability_raw', projection, "group", groupnr, '.txt'))
  

  
  
  ##tic-normalised CV (! last part code, sm_group = tic-normal from here)
  sm_group <- saved_sm_group #no 80%, no other normalise
  
  #TIC or total ion count normalisation accord to literature
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3124646/
  #this normalization will be based on the sum of all intensity values in the mass spectrum (i.e., TIC)
  #= intensity compID 1 / sum(intenstity alls compIDs) for 1 sample/spectrum
  averageTIC_before <- mean(as.matrix(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)])) 
  
  #devide by sum per sample
  normalized <-  apply(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)], 1, function(x){t(x/sum(x)) }) #iterate over rows =1 (samples), but saves results as column!!! 
  sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)] <- t(normalized) #needs to be transposed for correct tic-norm
  
  #must be re-scaled to same order of scale (afrondingsfouten bij kleine getallen tijdens berek PCA, daarom all getallen maal zelfde factor)
  averageTIC_after <- mean(as.matrix(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)]))
  factor <- averageTIC_before/averageTIC_after
  sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)] <- factor * sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)]
  

  ##only retain features if present in 80% of df
  perc_present <- colSums(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)] != 0)/amount_of_samples_group
  matrix_with_perc_present <- data.frame(x=colnames(sm_group)[COLLUMN_NR_START_VARIABLES:length(sm_group)],
                                         y=perc_present)
  names(matrix_with_perc_present) <- c("CompID", "perc_present")
  retained_variables_group <- matrix_with_perc_present[(matrix_with_perc_present[,2]>=0.8),]
  retained_variables_group <- retained_variables_group[,1] #compIDs of good ones
  
  all_matrix <- from_df_to_matrix(sm_group)
  all_matrix <- as.data.frame(all_matrix)
  
  matrix_filter_group <- all_matrix[, names(all_matrix) %in% (retained_variables_group)]
  
  sampleMetadata_until_start_variables <- subset(sm_group, select = c(1:(COLLUMN_NR_START_VARIABLES-1)))
  sampleMetadata_filter_group <- cbind(sampleMetadata_until_start_variables, matrix_filter_group)
  #write_dataframe_as_txt_file(sampleMetadata_filter_QCpool, 'sampleMetadata_filter_QCpool.txt')
  
  sm_group <- sampleMetadata_filter_group
  
  
  ## tic CV
  CV_values <- apply(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)], 2, function(x) sd(x)/mean(x)*100) #in percent
  matrix_with_CV_values <- data.frame(x=colnames(sm_group)[COLLUMN_NR_START_VARIABLES:length(sm_group)], y=CV_values)
  names(matrix_with_CV_values) <- c("CompID", "CV_values")
  #write_dataframe_as_txt_file(matrix_with_CV_values, paste0(name_project,'_CV_values_raw.txt'))
  
  #skip if amount_of_samples_group <= 5
  if(amount_of_samples_group <= 5){
    perc_of_feat_under_CVthesh <- "not calculated"
    amount_of_feat_under_CVthesh <- "not calculated"
    amount_of_feat_retained_80perc <- "not calculated"
  }
  if(amount_of_samples_group > 5){
    #features under CV threshold
    amount_of_feat_under_CVthesh <- nrow(matrix_with_CV_values[matrix_with_CV_values$CV_values <= CV_threshold,])
    #amount_of_feat_group <- ncol(sm_group[,COLLUMN_NR_START_VARIABLES:length(sm_group)])   
    #no, using retained_variables_group from 80% calc:
    amount_of_feat_retained_80perc <- length(retained_variables_group)
    perc_of_feat_under_CVthesh <- format(round((amount_of_feat_under_CVthesh/amount_of_feat_retained_80perc *100),2),nsmall=2)
  }
  
  report_cv_tic <- NULL
  report_cv_tic$name_group <- name_group
  report_cv_tic$amount_of_samples_group <- amount_of_samples_group
  report_cv_tic$amount_of_feat_under_CVthesh <- amount_of_feat_under_CVthesh
  report_cv_tic$amount_of_feat_retained_80perc <- amount_of_feat_retained_80perc
  report_cv_tic$perc_of_feat_under_CVthesh <- perc_of_feat_under_CVthesh
  
  #report_cv_raw$histogram <- ""
  report_cv_tics <- rbind(report_cv_tics, report_cv_tic)
  
  
  ## plot per group cumm CV
  #plot(ecdf(matrix_with_CV_values[,2]))
  #h <- hist(matrix_with_CV_values[,2])
  
  #https://stackoverflow.com/questions/29289046/r-ecdf-over-histogram
  max_count <- max(matrix_with_CV_values[,2])
  breaks <- round(max_count/10, 0) #so show intervals 10-20
  
  png(paste(name_project, "_3-4Featurerepeatability_tic" , projection, "group", groupnr, ".png", sep=""), width=7, height=5, units="in", res=150)
  h <- hist(matrix_with_CV_values[,2], col = "dodgerblue3", xlab = "Bin", xlim = c(0,100), breaks = breaks, main=paste0("Group ", group)) #"#4393C3"
  par(new = T)
  ec <- ecdf(matrix_with_CV_values[,2])
  plot(x = h$mids, y=ec(h$mids)*max(h$counts), col = rgb(0,0,0,alpha=0), axes=F, xlab=NA, ylab=NA)
  lines(x = h$mids, y=ec(h$mids)*max(h$counts), col ='red')
  axis(4, at=seq(from = 0, to = max(h$counts), length.out = 11), labels=seq(0, 100, 10), col = 'red', col.axis = 'red')
  mtext(side = 4, line = -2, 'cummulative %', col = 'red')
  dev.off()
  
  #calc cummulative_perc
  total_count <- sum(unlist(h$counts))
  factor <- total_count/100
  count_perc <- unlist(h$counts) / factor
  cummulative_perc <- cumsum(count_perc)
  
  #histogram cummulat perc table
  hist_table_group <- NULL
  hist_table_group$bin <- unlist(h$breaks[0:(length(h$breaks)-1)]) #bins, only start show, so wo last elem
  hist_table_group$interval <- paste0(unlist(h$breaks[0:(length(h$breaks)-1)]), "-", (unlist(h$breaks[0:(length(h$breaks)-1)])+10))
  hist_table_group$frequency <- unlist(h$counts)
  hist_table_group$cummulative_perc <- format(round(cummulative_perc, 2),nsmall=2)
  hist_table_group <- as.data.frame(hist_table_group)
  
  #cut-off only <100 in tabel hist keep (if present)
  if(length(which(hist_table_group$bin == 90)) > 0)({
    colnr_bin90 <- which(hist_table_group$bin == 90)
    hist_table_group <- hist_table_group[1:colnr_bin90,]
  })
  
  #write resutls 3-4Featurerepeatability_raw
  write_dataframe_as_txt_file(hist_table_group, paste0(name_project,'_3-4Featurerepeatability_tic', projection, "group", groupnr, '.txt'))
  
  
}
############ 




#write resutls 1setup
write_dataframe_as_txt_file(report_setup_groups, paste0(name_project,'_1Setup.txt'))

#write resutls 3-1Numberoffeatures
write_dataframe_as_txt_file(report_features_groups, paste0(name_project,'_3-1Numberoffeatures.txt'))

#write resutls 3-2TICsignal
write_dataframe_as_txt_file(report_tic_groups, paste0(name_project,'_3-2TICsignal.txt'))

#write resutls 3-3OutliersPCA
write_dataframe_as_txt_file(report_outliers_pcas, paste0(name_project,'_3-3OutliersPCA.txt'))

#write resutls 3-4Featurerepeatability_raw
write_dataframe_as_txt_file(report_cv_raws, paste0(name_project,'_3-4Featurerepeatability_raw.txt'))

#write resutls 3-4Featurerepeatability_tic
write_dataframe_as_txt_file(report_cv_tics, paste0(name_project,'_3-4Featurerepeatability_tic.txt'))





print("R pipeline - Part II: Validation_statistics - done!")
print(Sys.time())
end_time <- Sys.time()
print(end_time - start_time)
#
#####################
