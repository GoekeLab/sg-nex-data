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



###########################
## long read pacbio samples ##
###########################
pacbio_data <-  data.table(as.data.frame(read_xlsx(paste0('data_information.xlsx'), sheet = 2)))

bam.file <- pacbio_data$`bam-tx.path`
bam.file.basenames <- pacbio_data$public_name

yieldSize <- 1000000
save.dir <- "readLength17Apr_pacbio/"
if(!dir.exists(save.dir)){
    dir.create(save.dir)
}
library(GenomicAlignments)

local_path <- "RunBambu17Apr/"
noprint <- lapply(seq_along(bam.file.basenames)[nnn], function(x){
    bam_file_temp <- bam.file[x]
    system(paste0("aws s3 cp ", bam_file_temp, " ",local_path))
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
