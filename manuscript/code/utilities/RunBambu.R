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
sampleData <- data.table(as.data.frame(read_xlsx(paste0('.'), sheet = 1))) ## need to convert from tibble to data.frame
sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
                                                                                 runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(SG_NextData_Release=="Yes"),
                                       runName,
                                       `GIS Library ID`)]
sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]
sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:111]

bam.file <- sapply(sampleNames, function(x) paste0("s3://sg-nex-data/data/sequencing_data_ont/bam/genome/"x,"/",x,".bam"))


local_path <- "RunBambu5Jul/"
if(!dir.exists(local_path)) dir.create(local_path)
rcSaveDir <- paste0(local_path,"rc")
if(!dir.exists(rcSaveDir)) dir.create(rcSaveDir)
if(!dir.exists(paste0(rcSaveDir,"/raw_reads")) dir.create(paste0(rcSaveDir,"/raw_reads"))

system(paste0("aws s3 cp --no-sign-request ", bam.file[nnn], " ",local_path))
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
# 
seNoPut <- bambu(reads = local_bam_file,
                  rcOutDir = rcSaveDir,
                  annotations = bambuAnnotations,
                  genome = genome.file,
                  discovery = FALSE,
                  quant = FALSE,
                  yieldSize = 1000000,
                  verbose=TRUE)

system(paste0("rm ",local_bam_file))

library(bambu)
rcSaveDir <- "RunBambu22Apr/rc"
library(BiocFileCache)
rcfiles <- bfcinfo(BiocFileCache(rcSaveDir))$rpath
rcfiles <- rcfiles[file.exists(rcfiles)]
bambuAnnotations <- readRDS("bambuAnnotations.rds")
genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
se <- bambu(reads = rcfiles,
                 annotations = bambuAnnotations,
                 genome = genome.file,
                 ncore = 4,
                 returnDistTable = TRUE,
                 opt.em = list(degradationBias = FALSE),
                 verbose=TRUE)
saveRDS(se, file = "bambuOutput_May25.rds")

## bambu on one sample with read to transcript map##================
library(bambu)
seNoPut <- bambu(reads = local_bam_file,
                  rcOutDir = rcSaveDir,
                  annotations = bambuAnnotations,
                  genome = genome.file,
                  discovery = FALSE,
                  quant = FALSE,
                  yieldSize = 1000000,
                  trackReads = TRUE,
                  verbose=TRUE)


############################
## Pacbio samples         ##
############################
pacbio_data <-  data.table(as.data.frame(read_xlsx(paste0('.'), sheet = 2)))

bam.file <- pacbio_data$`bam-genome.path`)

library(bambu)

local_path <- "RunBambu17Apr_PacBio/"
if(!dir.exists(local_path)) dir.create(local_path)
rcSaveDir <- paste0(local_path,"rc")
if(!dir.exists(rcSaveDir)) dir.create(rcSaveDir)

# download bam file
system(paste0("aws s3 cp ", bam.file[nnn], " ",local_path))
local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)

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
                 verbose=TRUE)

system(paste0("rm ",local_bam_file))


# locally with read class files 
library(bambu)
rcSaveDir <- "RunBambu22May_PacBio/rc"
library(BiocFileCache)
rcfiles <- bfcinfo(BiocFileCache(rcSaveDir))$rpath
bambuAnnotations <- readRDS("bambuAnnotations.rds")
genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
se <- bambu(reads = rcfiles,
            annotations = bambuAnnotations,
            genome = genome.file,
            ncore = 6,
            returnDistTable = TRUE,
            opt.em = list(degradationBias = FALSE),
            verbose=TRUE)
saveRDS(se, file = "bambuOutput_PacBio_May22.rds")

se <- bambu(reads = rcfiles,
            annotations = bambuAnnotations,
            genome = genome.file,
            ncore = 4,
            returnDistTable = TRUE,
            opt.em = list(degradationBias = FALSE),
            NDR = 0.1,
            verbose=TRUE)
saveRDS(se, file = "bambuOutput_PacBioNDR0.1_Aug18.rds")


se <- bambu(reads = rcfiles,
            annotations = bambuAnnotations,
            genome = genome.file,
            ncore = 4,
            returnDistTable = TRUE,
            opt.em = list(degradationBias = FALSE),
            NDR = 1,
            verbose=TRUE)
saveRDS(se, file = "bambuOutput_PacBioNDR1_May22.rds")


############################
## spikein samples        ##
############################
library(bambu)
anno_rds_file <- "bambuAnnotations.rds"
bambuAnnotations <- readRDS(anno_rds_file)
genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
bam.file <- dir("spikein_bam_Apr17/", pattern = ".bam$", full.names = TRUE)
bam.file <- bam.file[!grepl("PacBio",bam.file)]
rcOutDir <- "spikein_bam_Apr17/rc_ont_May22/"
if(!dir.exists(rcOutDir)) dir.create(rcOutDir)
system.time(bambuOutput <- bambu(reads = bam.file,
                                 rcOutDir = rcOutDir,
                                 annotations = bambuAnnotations,
                                 genome = genome.file,
                                 returnDistTable = TRUE,
                                 opt.em = list(degradationBias = FALSE),
                                 stranded = FALSE, ncore = 6,
                                 NDR = 0, 
                                 yieldSize = 1000000, verbose = TRUE))
saveRDS(bambuOutput, file = paste0("bambuOutput_spikein_bam_ont_May22.rds"))

bam.file <- dir("spikein_bam_Apr17/", pattern = ".bam$", full.names = TRUE)
bam.file <- bam.file[grepl("PacBio",bam.file)]
rcOutDir <- "spikein_bam_Apr17/rc_pacbio/"
if(!dir.exists(rcOutDir)) dir.create(rcOutDir)
system.time(bambuOutput <- bambu(reads = bam.file,
                                 rcOutDir = rcOutDir,
                                 annotations = bambuAnnotations,
                                 genome = genome.file,
                                 returnDistTable = TRUE,
                                 opt.em = list(degradationBias = FALSE),
                                 stranded = FALSE, ncore = 6,
                                 NDR = 0, 
                                 yieldSize = 1000000, verbose = TRUE))
saveRDS(bambuOutput, file = paste0("bambuOutput_spikein_bam_pacbio_May22.rds"))


library(BiocFileCache)
rcfiles <- c(bfcinfo(BiocFileCache("spikein_bam_Apr17/rc_ont_May22/"))$rpath,
             bfcinfo(BiocFileCache("spikein_bam_Apr17/rc_pacbio/"))$rpath)
system.time(bambuOutput <- bambu(reads = rcfiles,
                                 annotations = bambuAnnotations,
                                 genome = genome.file,
                                 returnDistTable = TRUE,
                                 opt.em = list(degradationBias = FALSE),
                                 stranded = FALSE, ncore = 6,
                                 NDR = 0,
                                 yieldSize = 1000000, verbose = TRUE))
saveRDS(bambuOutput, file = paste0("bambuOutput_spikein_bam_May22.rds"))
