#####################################################
# Get number of junctions                           #
#####################################################


.libPaths("/mnt/dataSSD/software/R/site-library")
rm(list = ls())
require(GenomicFeatures)
require(GenomicAlignments) ##readGAlignments
require(AnnotationDbi)#loadDb
require(data.table)#fast large dataset manipulation
require(readxl)
sampleData <- data.table(as.data.frame(read_xlsx(paste0('/mnt/projectsInstanceStore1/chenying/ONT Master Table.xlsx'), sheet = 1))) ## need to convert from tibble to data.frame
sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
                                                                                 runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
#sampleData$runName_combined <- gsub('-pre|-Pre','',sampleData$runName) # there are runs with multiple datasets that should be combined together
sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(SG_NextData_Release=="Yes"),
                                       runName,
                                       `GIS Library ID`)]
sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]
#sampleNames <- unique(sampleData$runName_combined)##
#sampleNames_old <- unique(sampleData[grepl("ON002-RNA-R00177|ON002-RNA-R00178",name)]$runName)
sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[109]#[1:112]
local_path <- "/mnt/projectsInstanceStore2/chenying/RunBambu17Apr_junction/"
if(!dir.exists(local_path)) dir.create(local_path)

bam.file <- sapply(sampleNames, function(x) paste0("s3://sg-nex-data/data/sequencing_data_ont/bam/genome/"x,"/",x,".bam"))

juncData <- lapply(seq_along(sampleNames), function(r){
    system(paste0("aws s3 cp ", gsub('/mnt/ontdata','s3://ontdata.store.genome.sg',bam.file[r]), " ",local_path, " --profile ontdata.store.genome.sg "))
    local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
    bf <- open(BamFile(local_bam_file, yieldSize=500000))#10000000
    tmp_all <- NULL
    while(Rsamtools::isIncomplete(bf)){
        d <- readGAlignments(bf,
                             #index = bamIndicies(fileName),
                             param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)),
                             use.names=TRUE)
        tmp<- list(names = names(d),
                   njunc = njunc(d))#,
        # tx = as.character(gsub('\\..*','',seqnames(d))))
        rm(d)
        gc()
        tmp <- data.table(do.call('cbind',tmp))
        tmp[, runName := sampleNames[r]]
        # number of junction
        tmpJunc <- tmp[,list(nread = .N), by = list(njunc,runName)]
        rm(tmp)
        gc()
        tmp_all <- rbind(tmp_all, tmpJunc)
    }
    close(bf)
    system(paste0("rm -rvf ", local_bam_file))
    return(tmp_all)
})
data_lr <- do.call('rbind',juncData)
saveRDS(data_lr, file = paste0("/mnt/projectInstanceStore2/chenying/juncLR_May5_updatedsample109.rds"))


.libPaths("/mnt/dataSSD/software/R/site-library")
rm(list = ls())
require(GenomicFeatures)
require(GenomicAlignments) ##readGAlignments
require(AnnotationDbi)#loadDb
require(data.table)#fast large dataset manipulation
require(readxl)

pacbio_data <-  data.table(as.data.frame(read_xlsx(paste0('/mnt/projectsInstanceStore1/chenying/ONT Master Table.xlsx'), sheet = 2)))

bam.file <- gsub('s3://ontdata.store.genome.sg','/mnt/ontdata',pacbio_data$`bam-genome.path`)#
local_path <- "RunBambu17Apr_junction_pacbio/"
if(!dir.exists(local_path)) dir.create(local_path)

sampleNames <- pacbio_data$public_name
juncData <- lapply(seq_along(sampleNames), function(r){
    system(paste0("aws s3 cp ", gsub('/mnt/ontdata','s3://ontdata.store.genome.sg',bam.file[r]), " ",local_path, " --profile ontdata.store.genome.sg "))
    local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
    bf <- open(BamFile(local_bam_file, yieldSize=500000))#10000000
    tmp_all <- NULL
    while(Rsamtools::isIncomplete(bf)){
        d <- readGAlignments(bf,
                             #index = bamIndicies(fileName),
                             param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)),
                             use.names=TRUE)
        tmp<- list(names = names(d),
                   njunc = njunc(d))#,
        # tx = as.character(gsub('\\..*','',seqnames(d))))
        rm(d)
        gc()
        tmp <- data.table(do.call('cbind',tmp))
        tmp[, runName := sampleNames[r]]
        # number of junction
        tmpJunc <- tmp[,list(nread = .N), by = list(njunc,runName)]
        rm(tmp)
        gc()
        tmp_all <- rbind(tmp_all, tmpJunc)
    }
    close(bf)
    system(paste0("rm -rvf ", local_bam_file))
    return(tmp_all)
})
data_lr <- do.call('rbind',juncData)
saveRDS(data_lr, file = paste0("juncLR_pacbio_Apr25.rds"))
## short read data
juncData <- lapply(sr_runNames, function(r){
    print(r)
    print(which(sampleNames == r))
    sampleData_runName <- sampleData_sr[sampleData_sr$runName==r,]
    if(grepl("H9|HEYA8",r)){
        bamFile = paste0('/mnt/ontdata/Illumina/02_Mapping/',r) # genome bam file
    }else{
        bamFile = gsub('s3://ontdata.store.genome.sg','/mnt/ontdata',sampleData_runName$star_map_genomeBam.path) # genome bam file
    }
    if(length(bamFile)==0){
        return(NULL)
    }
    if(file_test('-d',bamFile)){
        bamFile <- dir(bamFile, full.names = TRUE)
        bamFile <- bamFile[grepl('.bam$',bamFile)]
    }
    print(bamFile)
    bf <- open(Rsamtools::BamFile(bamFile, yieldSize=1000000))#10000000
    tmp_all <- NULL
    while(Rsamtools::isIncomplete(bf)){
        d <- readGAlignments(bf,
                             #index = bamIndicies(fileName),
                             param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)),
                             use.names=TRUE)
        tmp<- list(names = names(d),
                   njunc = njunc(d))#,
        rm(d)
        gc()
        tmp <- data.table(do.call('cbind',tmp))
        tmp[, runName := r]
        tmp[, cellLine := sampleData_runName$cellLine]
        tmp[, protocol := "short_read"]
        
        # number of junction
        tmpJunc <- tmp[,list(nread = .N), by = list(njunc,runName, cellLine, protocol)]
        rm(tmp)
        gc()
        tmp_all <- rbind(tmp_all, tmpJunc)
    }
    close(bf)
    return(tmp_all)
})
data_sr <- do.call('rbind',juncData)
saveRDS(data_sr, file = paste0(wkdir,"juncSR.rds"))


rm(list = ls())
require(GenomicFeatures)
require(GenomicAlignments) ##readGAlignments
require(AnnotationDbi)#loadDb
require(data.table)#fast large dataset manipulation
require(readxl)
bam.file <- dir("spikein_bam_Apr17/", pattern = ".bam$", full.names = TRUE)
bam.file <- bam.file[!grepl("Illumina",bam.file)]
juncData <- lapply(bam.file, function(r){
    bf <- open(BamFile(r, yieldSize=500000))#10000000
    tmp_all <- NULL
    while(Rsamtools::isIncomplete(bf)){
        d <- readGAlignments(bf,
                             #index = bamIndicies(fileName),
                             param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)),
                             use.names=TRUE)
        tmp<- list(names = names(d),
                   njunc = njunc(d))#,
        # tx = as.character(gsub('\\..*','',seqnames(d))))
        rm(d)
        gc()
        tmp <- data.table(do.call('cbind',tmp))
        tmp[, runName := gsub(".bam","",basename(r))]
        
        # number of junction
        tmpJunc <- tmp[,list(nread = .N), by = list(njunc,runName)]
        rm(tmp)
        gc()
        tmp_all <- rbind(tmp_all, tmpJunc)
    }
    close(bf)
    return(tmp_all)
})
juncData <- do.call("rbind",juncData)
saveRDS(juncData, file = paste0("juncLR_spikein_May5.rds"))

bam.file <- dir("spikein_bam/", pattern = ".bam$", full.names = TRUE)
bam.file <- bam.file[grepl("Illumina",bam.file)]
juncData <- lapply(bam.file, function(r){
    bf <- open(Rsamtools::BamFile(r, yieldSize=1000000))#10000000
    tmp_all <- NULL
    while(Rsamtools::isIncomplete(bf)){
        d <- readGAlignments(bf,
                             #index = bamIndicies(fileName),
                             param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)),
                             use.names=TRUE)
        tmp<- list(names = names(d),
                   njunc = njunc(d))#,
        rm(d)
        gc()
        tmp <- data.table(do.call('cbind',tmp))
        tmp[, runName := gsub(".bam","",basename(r))]
        # number of junction
        tmpJunc <- tmp[,list(nread = .N), by = list(njunc,runName)]
        rm(tmp)
        gc()
        tmp_all <- rbind(tmp_all, tmpJunc)
    }
    close(bf)
    return(tmp_all)
})
data_sr <- do.call('rbind',juncData)
saveRDS(data_sr, file = paste0(wkdir,"juncSR_spikein.rds"))