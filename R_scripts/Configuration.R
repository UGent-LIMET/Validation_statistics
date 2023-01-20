## @knitr INFO
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Configuration


##########Global_settings##########

## @knitr settings
## options
RUN_CODE <- 'run this part of the pipeline'
DONT_RUN_CODE <- 'skip this part of the pipeline, keeps order: pre-processing - targeted analysis - statistical analysis - annotation'

## Adjustments
#Project name:
EXPERIMENT <- 'TEST_validation_statistics' #structured and short, see READ_ME
POLARITY <- "negative" #{"positive", "negative"} #needed this format for annotation to work
# file_conversion and pre-processing need to be performed seperate
# only choose "both" when no pre-processing needed (eg. merge both ionisationmodes)
USER_COMMENT <- "Tutorial comment" #Add info about experiment, eg. explain (Multiple)Comparisons, to include in reports

RUN_PART_VALIDATION <- DONT_RUN_CODE

#
#####################



##########Validation_statistics##########
if(RUN_PART_VALIDATION == RUN_CODE){
  
  ## options
  #Source of file1 variableMetadata:
  VARIABLEMETADATA_FROM_PIPELINE <- 'automatically finds variableMetadata from R pipeline code (based on projectname), present in input folder'
  VARIABLEMETADATA_EXTERN <- 'you choose name of variableMetadata.txt below, present in input folder' #see READ_ME for format file
  
  
  ## Adjustments
  #If you choose 'VARIABLEMETADATA_EXTERN', add additional info here:
  VARIABLEMETADATA_EXTERN <- 'VM.txt'  #'name.txt' of file. Ignore if file created from pipeline 
  COLLUMN_NR_START_SAMPLES <- 20  #always 20 (auto and manual must be same format); unless extra col merged!
  
  #Source of variableMetadata:
  INPUT_VARIABLES <- VARIABLEMETADATA_EXTERN
  
  #input file2 sampleMetadata:
  INPUT_SAMPLES <- 'SM.txt' #don't forget .txt
  COLLUMN_ORDER <- 8
  COLLUMN_NR_TYPE <- 9          #column number of the column 'Type'
  ORDER_NR_OF_FIRST_QC <- 1
  COLLUMN_NR_LAST_BEFORE_COMPARISONS <- 11
  AMOUNT_OF_COMPARISONS <- 1      #amount of pairwise comparisons, =0 if none
  AMOUNT_OF_MULTIPLE_COMPARISONS <- 3         #amount of multiple comparisons, =0 if none present
  AMOUNT_OF_PROJECTIONS <- 4      #amount of PCAs (or alternative projection methods in future), =0 if none present. 
  #! validation statistics are performed using Projection1, other columns are ignored here
  
  #theshold of coeficent of variation (CV) values in percentage (e.g., 30%)
  CV_THRESHOLD <- 30
  
}
#
#####################
