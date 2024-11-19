#!/usr/bin/env Rscript
##################################
# re-run bambu for all samples   #
##################################

###########################
## set-up                ##
###########################
rm(list = ls())
# set working directory
library(readxl)
library(data.table)


###########################
## long read all samples ##
###########################
require(AnnotationDbi)


require(docopt)
'Usage:
RunBambu.R [-g <g> ]

Options:
-g index
]' -> doc

opts <- docopt(doc)
print(opts)
nnn <- as.integer(opts$g)


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
sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:112]

bam.file <- sapply(sampleNames, function(x) paste0("s3://sg-nex-data/data/sequencing_data_ont/bam/genome/"x,"/",x,".bam"))


local_path <- "RunBambu5Jul/"
if(!dir.exists(local_path)) dir.create(local_path)
seSaveDir <- paste0(local_path,"se/")
if(!dir.exists(seSaveDir)) dir.create(seSaveDir)
rcSaveDir <- paste0(local_path,"rc")
if(!dir.exists(rcSaveDir)) dir.create(rcSaveDir)
if(!dir.exists(paste0(local_path,"raw_reads"))) dir.create(paste0(local_path,"raw_reads"))
system(paste0("aws s3 cp --no-sign-request ", bam.file[nnn], " ",local_path))
local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
anno_rds_file <- "bambuAnnotations.rds"
bambuAnnotations <- readRDS(anno_rds_file)
genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"


library(bambu)
seNoPut <- bambu(reads = local_bam_file,
                 rcOutDir = rcSaveDir,
                 annotations = bambuAnnotations,
                 genome = genome.file,
                 discovery = FALSE,
                 quant = FALSE,
                 yieldSize = 1000000,
                 trackReads = TRUE,
                 #opt.discovery = list(min.primarySecondaryDistStartEnd2 = 100000),
                 verbose=TRUE)

system(paste0("rm ",local_bam_file))


## run on server per sample
rcSaveDir <- "RunBambu5Jul/rc"
library(BiocFileCache)
rcfiles <- bfcinfo(BiocFileCache(rcSaveDir))$rpath #,
rcfiles <- rcfiles[file.exists(rcfiles)]
anno_rds_file <- "bambuAnnotations.rds"
bambuAnnotations <- readRDS(anno_rds_file)
genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
se <- bambu(reads = rcfiles[nnn],
            annotations = bambuAnnotations,
            genome = genome.file,
            ncore = 1,
            returnDistTable = TRUE,
            trackReads = TRUE,
            opt.em = list(degradationBias = FALSE),
            #opt.discovery = list(min.primarySecondaryDistStartEnd2 = 100000),
            verbose=TRUE)
saveRDS(se, file = paste0(seSaveDir,"/bambuOutput_5Jul_trackReads.rds")



