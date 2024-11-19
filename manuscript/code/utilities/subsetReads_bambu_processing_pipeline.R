#!/usr/bin/env Rscript


rm(list = ls())
# set working directory
wkdir <- './chenying'
setwd(wkdir)
library(readxl)
library(data.table)


require(AnnotationDbi)


require(docopt)
'Usage:
subsetReads_bambu_processing_pipeline.R [-g <g> -r <r>]

Options:
-g index
-r subsetType
]' -> doc

opts <- docopt(doc)
print(opts)
nnn <- as.integer(opts$g) # sample id
branch_names <- c("unique","best","noClose","primary")
subsetType <- branch_names[as.integer(opts$r)]
print(nnn)
print(subsetType)


####################
## readclass      ##
####################
sampleData <- data.table(as.data.frame(read_xlsx(paste0('ONT Master Table.xlsx'), sheet = 1))) ## need to convert from tibble to data.frame
sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
                                                                                 runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
#sampleData$runName_combined <- gsub('-pre|-Pre','',sampleData$runName) # there are runs with multiple datasets that should be combined together
sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(SG_NextData_Release=="Yes"),
                                       runName,
                                       `GIS Library ID`)]
sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]
#sampleNames <- unique(sampleData$runName_combined)##
#sampleNames_old <- unique(sampleData[grepl("ON002-RNA-R00177|ON002-RNA-R00178",name)]$runName)
sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:112]

bam.file <- sapply(sampleNames, function(x) paste0("s3://sg-nex-data/data/sequencing_data_ont/bam/genome/"x,"/",x,".bam"))
#


local_path <- "./RunBambu12June_check/"
if(!dir.exists(local_path)) dir.create(local_path)
rcSaveDir <- paste0(local_path,"rc/",subsetType)
if(!dir.exists(rcSaveDir)) dir.create(rcSaveDir)


# bambu_package <- "/mnt/dataSSD/chenying/bambu"
# devtools::load_all(bambu_package)
library(bambu)

# # download bam file
system(paste0("aws s3 cp --no-sign-request", bam.file[nnn], " ",local_path))
local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
# 
anno.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs_v5_reformatted.gtf"
anno_rds_file <- "bambuAnnotations.rds"
if(nnn==1){
    bambuAnnotations <- prepareAnnotations(anno.file)
    saveRDS(bambuAnnotations, file = anno_rds_file)
}else{
     bambuAnnotations <- readRDS(anno_rds_file)
 }

genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"

seNoPut <- bambu(reads = local_bam_file,
                  rcOutDir = rcSaveDir,
                  annotations = bambuAnnotations,
                  genome = genome.file,
                  discovery = FALSE,
                  quant = FALSE,
                  yieldSize = 1000000,
                 trackReads = TRUE,
                  #opt.discovery = list(min.primarySecondaryDistStartEnd2 = 100000),
                  verbose=TRUE,
                 subsetType = subsetType)

system(paste0("rm ",local_bam_file))