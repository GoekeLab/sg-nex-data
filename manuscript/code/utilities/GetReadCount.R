##################################
# get read count for each sample #
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

# sample information ==================
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
sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[109]



bam.file <- sapply(sampleNames, function(x) paste0("s3://sg-nex-data/data/sequencing_data_ont/bam/genome/"x,"/",x,".bam"))


bam.file.basenames <- gsub("_(R1.sorted|sorted)","",gsub(".genome_alignment","",tools::file_path_sans_ext(BiocGenerics::basename(bam.file))))

yieldSize <- 1000000
local_path <- "RunBambu17Apr_readcount/"
if(!dir.exists(local_path)) dir.create(local_path)

readCount <- do.call("rbind",lapply(seq_along(bam.file.basenames), function(x){
    bam_file_temp <- bam.file[x]
    system(paste0("aws s3 cp --no-sign-request ",bam_file_temp, " ",local_path))
    local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
    
    ## number of reads
    statsDt <- data.table(runname = bam.file.basenames[x],
                          total_reads = as.integer(system(paste('samtools view -@ 24 -F 2304 -c ', local_bam_file, sep=''),intern=TRUE)), #number of reads
                          aligned_reads = as.integer(system(paste('samtools view -@ 24 -F 0x904 -c ',  local_bam_file, sep=''),intern=TRUE))) #number of reads aligned
    
    system(paste0("rm -v ",local_bam_file))
    return(statsDt)
}))
saveRDS(readCount, file = "readCount_guppy_6_4_2_updatedsample109.rds") # readCount.rds for combined runs

###########################
## long read pacbio samples ##
###########################
pacbio_data <-  data.table(as.data.frame(read_xlsx(paste0('ONT Master Table.xlsx'), sheet = 2)))

bam.file <- gsub('s3://ontdata.store.genome.sg','/mnt/ontdata',pacbio_data$`bam-genome.path`)#
bam.file.basenames <- pacbio_data$public_name
yieldSize <- 1000000
local_path <- "RunBambu17Apr_readcount/"
if(!dir.exists(local_path)) dir.create(local_path)

readCount <- do.call("rbind",lapply(seq_along(bam.file.basenames), function(x){
    bam_file_temp <- bam.file[x]
    system(paste0("aws s3 cp ", gsub('/mnt/ontdata','s3://ontdata.store.genome.sg', bam_file_temp), " ",local_path, " --profile ontdata.store.genome.sg "))
    local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
    
    ## number of reads
    statsDt <- data.table(runname = bam.file.basenames[x],
                          total_reads = as.integer(system(paste('samtools view -@ 24 -F 2304 -c ', local_bam_file, sep=''),intern=TRUE)), #number of reads
                          aligned_reads = as.integer(system(paste('samtools view -@ 24 -F 0x904 -c ',  local_bam_file, sep=''),intern=TRUE))) #number of reads aligned
    
    system(paste0("rm -v ",local_bam_file))
    return(statsDt)
}))
saveRDS(readCount, file = "readCount_pacbio.rds") # readCount.rds for combined runs


############################
## short read all samples ##
############################
sampleData_sr <- data.table(as.data.frame(read_xlsx('ONT Master Table.xlsx', sheet = 2))) ## need to convert from tibble to data.frame
sampleData_sr <- sampleData_sr[!grepl('#',`ELM library ID`) &(!is.na(runName))&(!grepl("HEYA8.*H9",runName))]
sr_runNames <- sampleData_sr$runName
chrm_names <- c(1:22,'X','Y')

bam.file <- unlist(lapply(sr_runNames, function(r){
    
    bam.file <- gsub('GRCh','Grch',gsub('2.1-','2.17-',gsub('s3://ontdata.store.genome.sg','/mnt/ontdata', 
                                                            sampleData_sr[runName == r]$`star_map_genomeBam.path`)))#paste0("s3://ontdata.store.genome.sg/Nanopore/03_Mapping/Grch38/minimap2-2.17-cDNA/",r,"/")
    if(file_test('-d',bam.file)){
        bam.file <- dir(bam.file, full.names = TRUE)
        bam.file <- bam.file[grepl('.bam$',bam.file)]
    }
    return(bam.file)
}))
bam.file.basenames <- gsub("STAR_alignment","",gsub("_genome","",tools::file_path_sans_ext(BiocGenerics::basename(bam.file))))

yieldSize <- 1000000
readCount <- do.call("rbind",lapply(seq_along(bam.file.basenames), function(x){
    bam_file_temp <- bam.file[x]
    ## number of reads
    statsDt <- data.table(runname = bam.file.basenames[x],
                          total_reads = as.integer(system(paste('samtools view -F 2304 -c ', bam_file_temp, sep=''),intern=TRUE)), #number of reads
                          aligned_reads = as.integer(system(paste('samtools view -F 0x904 -c ',  bam_file_temp, sep=''),intern=TRUE))) #number of reads aligned
    return(statsDt)
}))
saveRDS(readCount, file = "readCountSR.rds")



############################
## spikein samples        ##
############################
bam.file <- dir("spikein_bam_Apr17/", pattern = ".bam$", full.names = TRUE)
bam.file.basenames <- gsub("_(R1.sorted|sorted)","",gsub(".genome_alignment","",tools::file_path_sans_ext(BiocGenerics::basename(bam.file))))

yieldSize <- 1000000
readCount <- do.call("rbind",lapply(seq_along(bam.file.basenames), function(x){
    bam_file_temp <- bam.file[x]
    ## number of reads
    statsDt <- data.table(runname = bam.file.basenames[x],
                          total_reads = as.integer(system(paste('samtools view -@ 24 -F 2304 -c ', bam_file_temp, sep=''),intern=TRUE)), #number of reads
                          aligned_reads = as.integer(system(paste('samtools view -@ 24 -F 0x904 -c ',  bam_file_temp, sep=''),intern=TRUE))) #number of reads aligned
    return(statsDt)
}))
saveRDS(readCount, file = "readCount_spikein_May5.rds") # readCount.rds for combined runs

