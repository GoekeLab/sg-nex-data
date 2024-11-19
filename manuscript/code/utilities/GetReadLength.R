#!/usr/bin/env Rscript
##################################
# get read length for each sample #
##################################

###########################
## set-up                ##
###########################
rm(list = ls())
# set working directory
library(readxl)
library(data.table)

require(docopt)
'Usage:
GetReadLength.R [-g <g> ]

Options:
-g index
]' -> doc

opts <- docopt(doc)
print(opts)
nnn <- as.integer(opts$g)

###########################
## long read all samples ##
###########################

# sample information ==================
# sampleData <- data.table(as.data.frame(read_xlsx(paste0('/mnt/projectsInstanceStore1/chenying/ONT Master Table.xlsx'), sheet = 1))) ## need to convert from tibble to data.frame
# sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
#                                                                                  runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
# #sampleData$runName_combined <- gsub('-pre|-Pre','',sampleData$runName) # there are runs with multiple datasets that should be combined together 
# sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(SG_NextData_Release=="Yes"),
#                                        runName, 
#                                        `GIS Library ID`)]
# sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]
# #sampleNames <- unique(sampleData$runName_combined)##
# #sampleNames_old <- unique(sampleData[grepl("ON002-RNA-R00177|ON002-RNA-R00178",name)]$runName)
# sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:111]
# 
# 
# bam.file <- gsub('s3://ontdata.store.genome.sg','/mnt/ontdata',
#                          paste0("s3://ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38/guppy-6.4.2-updated/transcriptome/",sampleNames,".bam"))#
# 
# bam.file.basenames <- gsub("_(R1.sorted|sorted)","",gsub(".genome_alignment","",tools::file_path_sans_ext(BiocGenerics::basename(bam.file))))
# 
# 
# yieldSize <- 1000000
# save.dir <- "/mnt/projectsInstanceStore1/chenying/readLength17Apr/"
# if(!dir.exists(save.dir)){
#     dir.create(save.dir)
# }
# library(GenomicAlignments)
# # rdfile <- dir("/mnt/data/chenying/readLength/", full.names = TRUE)
# # rdfile_name <- gsub("readLength_","",tools::file_path_sans_ext(BiocGenerics::basename(rdfile)))
# 
# local_path <- "/mnt/projectsInstanceStore2/chenying/RunBambu17Apr/"
# noprint <- lapply(seq_along(bam.file.basenames)[nnn], function(x){
#     bam_file_temp <- bam.file[x]
#     system(paste0("aws s3 cp ", gsub('/mnt/ontdata','s3://ontdata.store.genome.sg',bam_file_temp), " ",local_path, " --profile ontdata.store.genome.sg "))
#     local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
#     print(x)
#     print(local_bam_file)
#     bf <- open(Rsamtools::BamFile(local_bam_file, yieldSize = yieldSize))
#     counter <- 1
#     temp_sampleData <- list()
#     while (Rsamtools::isIncomplete(bf)) {
#         rr <- readGAlignments(bf,param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE, 
#                                                                      isDuplicate = FALSE),
#                                                     what=c("qual", "flag","mapq")),use.names=T)
#         rr <- rr[mcols(rr)$flag %in% c(0,16)]
#         temp_sampleData[[counter]] <- data.table(read_len = as.numeric(qwidth(rr)),
#                                                  aligned_len = as.numeric(width(rr)),
#                                                  tx_name = gsub('\\..*','',as.character(seqnames(rr))),
#                                                  runname = bam.file.basenames[x])
#         print(min(length(rr),
#                   counter * yieldSize, na.rm = TRUE))
#         counter <- counter + 1
#     }
#     on.exit(close(bf))
#     temp_sampleData <- do.call("rbind", temp_sampleData)
#     saveRDS(temp_sampleData, file = paste0(save.dir,"readLength_",bam.file.basenames[x],".rds"))
#     system(paste0("rm -v ",local_bam_file))
# })


###########################
## long read pacbio samples ##
###########################
pacbio_data <-  data.table(as.data.frame(read_xlsx(paste0('ONT Master Table.xlsx'), sheet = 2)))

bam.file <- gsub('s3://ontdata.store.genome.sg','/mnt/ontdata',pacbio_data$`bam-tx.path`)#
bam.file.basenames <- pacbio_data$public_name

yieldSize <- 1000000
save.dir <- "readLength17Apr_pacbio/"
if(!dir.exists(save.dir)){
    dir.create(save.dir)
}
library(GenomicAlignments)
# rdfile <- dir("readLength/", full.names = TRUE)
# rdfile_name <- gsub("readLength_","",tools::file_path_sans_ext(BiocGenerics::basename(rdfile)))

local_path <- "RunBambu17Apr/"
noprint <- lapply(seq_along(bam.file.basenames)[nnn], function(x){
    bam_file_temp <- bam.file[x]
    system(paste0("aws s3 cp ", gsub('/mnt/ontdata','s3://ontdata.store.genome.sg',bam_file_temp), " ",local_path, " --profile ontdata.store.genome.sg "))
    local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
    print(x)
    print(local_bam_file)
    bf <- open(Rsamtools::BamFile(local_bam_file, yieldSize = yieldSize))
    counter <- 1
    temp_sampleData <- list()
    while (Rsamtools::isIncomplete(bf)) {
        rr <- readGAlignments(bf,param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE,
                                                                     isDuplicate = FALSE),
                                                    what=c("qual", "flag","mapq")),use.names=T)
        rr <- rr[mcols(rr)$flag %in% c(0,16)]
        temp_sampleData[[counter]] <- data.table(read_len = as.numeric(qwidth(rr)),
                                                 aligned_len = as.numeric(width(rr)),
                                                 tx_name = gsub('\\..*','',as.character(seqnames(rr))),
                                                 runname = bam.file.basenames[x])
        print(min(length(rr),
                  counter * yieldSize, na.rm = TRUE))
        counter <- counter + 1
    }
    on.exit(close(bf))
    temp_sampleData <- do.call("rbind", temp_sampleData)
    saveRDS(temp_sampleData, file = paste0(save.dir,"readLength_",bam.file.basenames[x],".rds"))
    system(paste0("rm -v ",local_bam_file))
})

############################
## short read all samples ##
############################
# sampleData_sr <- data.table(as.data.frame(read_xlsx('ONT Master Table.xlsx', sheet = 3))) ## need to convert from tibble to data.frame
# sampleData_sr <- sampleData_sr[!grepl('#',`ELM library ID`) &(!is.na(runName))&(!grepl("HEYA8.*H9",runName))]
# sr_runNames <- sampleData_sr$runName
# chrm_names <- c(1:22,'X','Y')
# 
# bam.file <- unlist(lapply(sr_runNames, function(r){
#     rnames <- sampleData_sr[runName == r]$runName
#     bam.file <- unlist(lapply(rnames, function(k){
#         if(k == "GIS_HEYA8_Illumina_Rep2-Run1"){
#             bam.file <- "sr_bam/GIS_HEYA8_Illumina_Rep2-Run1.bam"
#         }else if(k == "GIS_MCF7_Illumina_Rep2-Run1"){
#             bam.file <- "GIS_MCF7_Illumina_Rep2-Run1_transcriptome.bam"
#         }else{
#             bam.file <- gsub('GRCh','Grch',gsub('2.1-','2.17-',gsub('(s3://ontdata.store.genome.sg)|(s3://ontdata.store.transcript.sg)','/mnt/ontdata/', 
#                                                                     sampleData_sr[runName == k]$`star_map_txBam.path`)))
#         }
#         
#         if(file_test('-d',bam.file)){
#             bam.file <- dir(bam.file, full.names = TRUE)
#             bam.file <- bam.file[grepl('.bam$',bam.file)]
#         }
#         if(length(bam.file)==0){
#             print(r)
#             print(which(sr_runNames == r))
#         }
#         return(bam.file)
#     }))
#     return(bam.file)
# }))
# bam.file.basenames <- gsub("(_R1.sorted)|(_sorted)|(_transcriptome)","",gsub(".genome_alignment","",tools::file_path_sans_ext(BiocGenerics::basename(bam.file))))
# 
# 
# yieldSize <- 1000000
# save.dir <- "readLength/"
# if(!dir.exists(save.dir)){
#     dir.create(save.dir)
# }
# library(GenomicAlignments)
# # rdfile <- dir("readLength/", full.names = TRUE)
# # rdfile_name <- gsub("readLength_","",tools::file_path_sans_ext(BiocGenerics::basename(rdfile)))
# 
# 
# noprint <- lapply(seq_along(bam.file.basenames), function(x){
#     bam_file_temp <- bam.file[x]
#     if(!grepl(".bam$",bam_file_temp)){
#         return(NULL)
#     }
#     print(x)
#     print(bam_file_temp)
#     bf <- open(Rsamtools::BamFile(bam_file_temp, yieldSize = yieldSize))
#     counter <- 1
#     temp_sampleData <- list()
#     while (Rsamtools::isIncomplete(bf)) {
#         rr <- readGAlignments(bf,param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE, 
#                                                                      isDuplicate = FALSE),
#                                                     what=c("qual", "flag","mapq")),use.names=T)
#         temp_sampleData[[counter]] <- data.table(read_len = as.numeric(qwidth(rr)),
#                                                  aligned_len = as.numeric(width(rr)),
#                                                  tx_name = gsub('\\..*','',as.character(seqnames(rr))),
#                                                  runname = bam.file.basenames[x])
#         print(min(length(rr),
#                   counter * yieldSize, na.rm = TRUE))
#         counter <- counter + 1
#     }
#     on.exit(close(bf))
#     temp_sampleData <- do.call("rbind", temp_sampleData)
#     saveRDS(temp_sampleData, file = paste0(save.dir,"readLength_",bam.file.basenames[x],".rds"))
# })




############################
## spikein samples        ##
############################
bam.file <- dir("spikein_bam_tx_Apr17/", pattern = ".bam$", full.names = TRUE)
bam.file.basenames <- gsub("_(R1.sorted|sorted)","",gsub(".genome_alignment","",tools::file_path_sans_ext(BiocGenerics::basename(bam.file))))


yieldSize <- 1000000
save.dir <- "readLength17Apr_spikein/"
if(!dir.exists(save.dir)){
    dir.create(save.dir)
}
library(GenomicAlignments)
# 
# 
noprint <- lapply(seq_along(bam.file.basenames)[12], function(x){
    bam_file_temp <- bam.file[x]
    if(!grepl(".bam$",bam_file_temp)){
        return(NULL)
    }
    print(x)
    print(bam_file_temp)
    bf <- open(Rsamtools::BamFile(bam_file_temp, yieldSize = yieldSize))
    counter <- 1
    temp_sampleData <- list()
    while (Rsamtools::isIncomplete(bf)) {
        rr <- readGAlignments(bf,param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE,
                                                                     isDuplicate = FALSE),
                                                    what=c("qual", "flag","mapq")),use.names=T)
        if(!grepl("Illumina",bam_file_temp)){
            rr <- rr[mcols(rr)$flag %in% c(0,16)]
        }
        temp_sampleData[[counter]] <- data.table(read_len = as.numeric(qwidth(rr)),
                                                 aligned_len = as.numeric(width(rr)),
                                                 tx_name = gsub('\\..*','',as.character(seqnames(rr))),
                                                 runname = bam.file.basenames[x])
        print(min(length(rr),
                  counter * yieldSize, na.rm = TRUE))
        counter <- counter + 1
    }
    on.exit(close(bf))
    temp_sampleData <- do.call("rbind", temp_sampleData)
    saveRDS(temp_sampleData, file = paste0(save.dir,"readLength_",bam.file.basenames[x],".rds"))
})
