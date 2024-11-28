#!/usr/bin/env Rscript
.libPaths("/mnt/dataSSD/software/R/site-library")
rm(list = ls())

library(GenomicAlignments) ##readGAlignments
library(AnnotationDbi)#loadDb
library(data.table)


require(docopt)
'Usage:
run_trim_reads_with_error.R [-g <g> -t <t> -s <s>]

Options:
-g index
-t bp
-s times
]' -> doc

opts <- docopt(doc)
print(opts)
g <- as.numeric(opts$g)
trim_bp <- as.numeric(opts$t)
trim_times <- as.numeric(opts$s)


wkdir <- '.'
setwd(wkdir)
## need to manually edit the fasta file
txdbEnsembl91 <- loadDb('hg38_sequins_SIRV_ERCCs_longSIRVs-txdb.sqlite')
txLengths <- transcriptLengths(txdbEnsembl91)
txLengths <- data.table(txLengths)
txRefFile <- "Homo_sapiens.GRCh38.cdna.ncrna_wtChrIs_modified.fa"
txRefFile_matched <- gsub("wtChrIs","wtChrIs_matchedToGTF",txRefFile)


annotationDir <- "transcriptome-index/salmon_index_hg38_sirv_longsirv_ercc_sequin"



setwd(wkdir)
library(readxl)

sampleData <- data.table(as.data.frame(read_xlsx(paste0('.'), sheet = 1))) ## need to convert from tibble to data.frame
sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
                                                                                 runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(SG_NextData_Release=="Yes"),
                                       runName, 
                                       `GIS Library ID`)]
sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]

sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:111]


sampleDataSR <- as.data.frame(read_xlsx(paste0('.'), sheet = 3))## need to convert from tibble to data.frame
sampleDataSR <- sampleDataSR[!grepl('#',sampleDataSR[,1]) &(!is.na(sampleDataSR$runName)),]

sampleNamesSR <- sampleDataSR$runName

# utility_functions
trim_function <- function(trim_times, trim_bp, rnames, sampleData, prefix = "sr",wkdir,annotationDir){
    
    np <- lapply(rnames, function(k){
        outList <- get_fastq(k,sampleData,wkdir, prefix, trim_bp)
       
        outList_trim <- trim_process(prefix,outList[[1]],trim_bp, outList[[2]],wkdir, trim_times)
        
        ## create map directory
        mapDir <- paste0(wkdir,'/02_Mapping_matchedToGTF/',outList[[3]])
        map_function(mapDir, annotationDir,outList_trim[[1]],trim_bp, outList[[2]])
        system(paste0('rm -rvf ',outList[[3]])) 
    })
}

get_fastq <- function(k,sampleData,wkdir, prefix, trim_bp){
    
    if(prefix == "lr"){
        rname <- k
    }else{
        sampleData_runName <- sampleData[sampleData$runName==k,]
        rname <- sampleData_runName$public_name
    }
    fastqDir <- paste0(wkdir,'/01_Fastq/',rname)
    
    if(!dir.exists(fastqDir)){
        dir.create(fastqDir, recursive = TRUE)
    }
    ## download fastq file first
    if(prefix == "lr"){
        s3_fastq_dir <-paste0('s3://sg-nex-data/data/sequencing_data_ont/fastq/',rname)
        system(paste0('aws s3 cp ',s3_fastq_dir,' ', fastqDir))
    }
    
    if(prefix == "sr"){
        s3_fastq_dir <- paste0('s3://sg-nex-data/data/sequencing_data_illumina/fastq/',rname)
        system(paste0('aws s3 sync --no-sign-request ',s3_fastq_dir," ", fastqDir))
    }
    
     
    
    if(prefix == "sr"){
        fileList1 <- paste0(fastqDir,'/',rname,"_R1.fastq.gz")
        fileList2 <- paste0(fastqDir,'/',rname,"_R2.fastq.gz")
        for(i in seq_along(fileList2)){
            system(paste0('rm -rvf ',fileList2[i])) 
        }
    }else{
        fileList1 <- paste0(fastqDir,'/',rname,".fastq.gz")
    }
    print(paste0("finish downloading fastq"))
    return(list(fileList1,fastqDir,rname))
}

trim_process <- function(prefix,fileList1,trim_bp, fastqDir,wkdir, trim_times){
    unzip_fileList1 <- gsub(".gz$","",fileList1)
    for(i in seq_along(fileList1)){
        system(paste0("gunzip < ",fileList1[i]," > ", unzip_fileList1[i]))
        system(paste0('rm -rvf ',fileList1[i])) 
        for(kk in seq_len(trim_times)){
            unzip_trim_fileList1 <- gsub(".fastq$",paste0("_",trim_bp,"bp_",kk,".fastq"),unzip_fileList1)
            renamed_unzip_trim_fileList1 <- gsub(".fastq$","_renamed.fastq",unzip_trim_fileList1)
            gzip_trim_fileList1 <- gsub(".fastq$",".fastq.gz",renamed_unzip_trim_fileList1)
            system(paste0("seqtk trimfq ",c("-L ","-l ")[(prefix=="lr")+1],
                                       trim_bp+1, 
                                        " -q ",
                                       rep(seq(0.01, 0.1, by = 0.01),4)[kk]," ",
                                       unzip_fileList1[i]," > ",unzip_trim_fileList1[i])) #
            
            if(trim_times>1){
                system(paste0("cat ",unzip_trim_fileList1[i]," | sed 's/ runid/",kk," runid/g' > ",
                              renamed_unzip_trim_fileList1[i]))
                system(paste0("gzip < ",renamed_unzip_trim_fileList1[i], " > ",gzip_trim_fileList1[
                    i]))
            }else{
                system(paste0("gzip < ",unzip_trim_fileList1[i], " > ",gzip_trim_fileList1[i]))
            }
            
            system(paste0('rm -rvf ',unzip_trim_fileList1[i]))
            system(paste0('rm -rvf ',renamed_unzip_trim_fileList1[i])) 
        }
        system(paste0('rm -rvf ',unzip_fileList1[i])) 
      
    }
    if(trim_times>1) {
        setwd(fastqDir)
        gzip_trim_fileList1 <- gsub(".fastq$",
                                    paste0("_",trim_bp,"bp.fastq.gz"),unzip_fileList1)
        system(paste0("cat * > ",gzip_trim_fileList1))
        setwd(wkdir)
    }
    print("finishing trimming process")
    return(gzip_trim_fileList1)
}
map_function <- function(mapDir, annotationDir,gzip_trim_fileList1,trim_bp,fastqDir){
    if(!dir.exists(mapDir)){
        dir.create(mapDir, recursive = TRUE)
    }
    salmonPath <- "salmon"
    system(paste0(salmonPath,' quant -i ',annotationDir,
                              ' -l A -r ',
                              paste(gzip_trim_fileList1, collapse = ' '),
                              ' --validateMappings ',
                              ' --fldMean  ', trim_bp+1,' ', # 
                              ' --fldSD 1 ',
                              ' --seqBias ', ' --gcBias ', ' --posBias ',
                              ' -o ', mapDir,'/transcripts_quant_biasCorrected'))
    system(paste0("rm -rvf ",fastqDir))
    print("finish mapping")
}



# run it 
set.seed(1)
prefix <- "sr"
trim_bp <- 75
trim_times <- 1
wkdir <- paste0(wkdir,'trim_reads/',prefix,'_',trim_bp,'bpSingleEnd_',trim_times,'ts')
trim_function(trim_times, trim_bp, sampleNamesSR[g], sampleDataSR, prefix = prefix,wkdir,annotationDir)





