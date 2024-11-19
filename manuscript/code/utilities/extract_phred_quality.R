#!/usr/bin/env Rscript
#############################################################
# re-run bambu for all samples where filter quality number  #
#############################################################

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
extract_phred_quality.R [-g <g> ]

Options:
-g index
]' -> doc

opts <- docopt(doc)
print(opts)
nnn <- as.integer(opts$g)


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



local_path <- "FilterBamByQuality/"
if(!dir.exists(local_path)) dir.create(local_path)
rcSaveDir <- paste0(local_path,"rc")
if(!dir.exists(rcSaveDir)) dir.create(rcSaveDir)
if(!dir.exists(paste0(rcSaveDir,"/raw_reads"))) dir.create(paste0(rcSaveDir,"/raw_reads"))
   # 
   # # download bam file
#nnn <- 1
system(paste0("aws s3 cp --no-sign-request ", bam.file[nnn], " ",local_path))
local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)

library(bambu)

anno_rds_file <- "bambuAnnotations.rds"
bambuAnnotations <- readRDS(anno_rds_file)

genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"


seNoPut <- bambu(reads = local_bam_file,
                  rcOutDir = rcSaveDir,
                  annotations = bambuAnnotations,
                  genome = genome.file,
                  discovery = FALSE,
                  quant = FALSE,
                  yieldSize = 1000000,
                  #opt.discovery = list(min.primarySecondaryDistStartEnd2 = 100000),
                  verbose=TRUE)
system(paste0("rm ",local_bam_file))
# 

library(bambu)
rcSaveDir <- "FilterBamByQuality/rc"
library(BiocFileCache)
rcfiles <- bfcinfo(BiocFileCache(rcSaveDir))$rpath #,
rcfiles <- rcfiles[file.exists(rcfiles)]
bambuAnnotations <- readRDS("bambuAnnotations.rds")
genome.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
se <- bambu(reads = rcfiles,
            annotations = bambuAnnotations,
            genome = genome.file,
            ncore = 4,
            returnDistTable = TRUE,
            NDR = 0.1, 
            opt.em = list(degradationBias = FALSE),
            #opt.discovery = list(min.primarySecondaryDistStartEnd2 = 100000),
            verbose=TRUE)
saveRDS(se, file = "bambuOutput_21May2024_NDR0.1.rds")


# salmon-lr filtered=======================
tx2gene <- txLengths[,c(2,3)]
x <- 1

salmon_lr.dir <- "salmon_filter/count/"

sampleNames <- dir(salmon_lr.dir)
salmon_lr <- do.call('rbind',lapply(sampleNames,function(k){
    print(k)
    filePath <- sort(dir(paste0(salmon_lr.dir,k),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
    
    
    if(length(filePath)==0){
        return(NULL)
    }
    print(filePath)
    txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
    names(txi)
    
    salmon_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                              abundance = txi$abundance[,1],
                              counts = txi$counts[,1],
                              length = txi$length[,1],
                              countsFromAbundance = txi$countsFromAbundance)
    # short_read <- fread(filePath, header = TRUE)
    salmon_read[, runname:=k]
    return(salmon_read)
}))




salmon_lr[, method:='salmon_lr_q7filter']
salmon_lr[, ntotal:=sum(counts), by = runname]
setnames(salmon_lr, 'abundance','estimates')
salmon_lr[, `:=`(#counts = NULL,
    length = NULL,
    countsFromAbundance = NULL)]

salmon_lr <- geneTxTable[salmon_lr, on = 'tx_name']
salmon_lr[, TPM:=estimates]
salmon_lr[, estimates:=TPM/1000000*ntotal]
saveRDS(salmon_lr,file = "salmon_lr_Q7filtered.rds")
