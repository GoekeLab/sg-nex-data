#!/usr/bin/env Rscript
rm(list = ls())

library(GenomicAlignments) ##readGAlignments
library(AnnotationDbi)#loadDb
library(data.table)


require(docopt)
'Usage:
run_trim_reads_pacbio.R [-g <g> ]

Options:
-g index
]' -> doc

opts <- docopt(doc)
print(opts)
g <- as.numeric(opts$g)



wkdir <- 'wkdir'
annoDir <- "anno"
setwd(wkdir)
## need to manually edit the fasta file
txdbEnsembl91 <- loadDb(paste0(annoDir,"hg38_sequins_SIRV_ERCCs_longSIRVs-txdb.sqlite"))
txLengths <- transcriptLengths(txdbEnsembl91)
txLengths <- data.table(txLengths)
txRefFile <- paste0(annoDir,"Homo_sapiens.GRCh38.cdna.ncrna_wtChrIs_modified.fa")
txRefFile_matched <- gsub("wtChrIs","wtChrIs_matchedToGTF",txRefFile)


annotationDir <- paste0("transcriptome-index/salmon_index_hg38_sirv_longsirv_ercc_sequin")



setwd(wkdir)
library(readxl)


sampleData <- data.table(as.data.frame(read_xlsx(paste0(annoDir,"."), sheet = 2))) ## need to convert from tibble to data.frame
sampleNames <- sampleData$public_name
library(BiocParallel)


# utility_functions
trim_function <- function(trim_bp, rnames, sampleNames, sampleData, prefix = "lr",wkdir,annotationDir){
    
    np <- lapply(rnames, function(k){
       
        outList <- get_fastq(k,sampleData,wkdir, prefix, trim_bp)
        outList_trim <- trim_process(prefix,outList[[1]],trim_bp, outList[[2]],wkdir)
        
        ## create map directory
        mapDir <- paste0(wkdir,'trim_reads/',prefix,'_',trim_bp,'bpSingleEnd_1ts/02_Mapping_matchedToGTF/',outList[[3]])
        map_function(mapDir, annotationDir,outList_trim[[1]],trim_bp, outList[[2]])
        system(paste0('rm -rvf ',outList[[3]])) 
    })
}

get_fastq <- function(k,sampleData,wkdir, prefix, trim_bp){
    sampleData_runName <- sampleData[sampleData$public_name==k,]
    rname <- sampleData_runName$name 
    fastqDir <- paste0(wkdir,'trim_reads/',prefix,'_',trim_bp,'bpSingleEnd_1ts/01_Fastq/',k,'/')
    
    if(!dir.exists(fastqDir)){
        dir.create(fastqDir, recursive = TRUE)
    }
    ## download fastq file first
    s3_fastq_dir <- paste0('s3://',rname,'.fastq.gz')
    system(paste0('aws s3 cp ',s3_fastq_dir," ", fastqDir))
    fileList1 <- paste0(fastqDir,k,".fastq.gz")
    fileList_old <- dir(fastqDir, full.names = TRUE)
    system(paste0("mv -v ", fileList_old, " ",fileList1))
    print(paste0("finish downloading fastq"))
    return(list(fileList1,fastqDir,rname))
}




trim_process <- function(prefix,fileList1,trim_bp, fastqDir,wkdir){
    
        unzip_fileList1 <- gsub(".gz$","",fileList1)
        system(paste0("gunzip < ",fileList1," > ", unzip_fileList1))
        system(paste0('rm -rvf ',fileList1)) 
        
        remain_fileList <-  unzip_fileList1
        trim_times <- 1
        while(file.size(remain_fileList)>0){
            ## check if remain_fileList is empty
            
            unzip_trim_fileList1 <- gsub(".fastq$",paste0("_",trim_bp,"bp_",trim_times,".fastq"),unzip_fileList1)
            renamed_unzip_trim_fileList1 <- gsub(".fastq$","_renamed.fastq",unzip_trim_fileList1)
            gzip_trim_fileList1 <- gsub(".fastq$",".fastq.gz",renamed_unzip_trim_fileList1)
            system(paste0("seqtk trimfq ","-L ",
                          trim_bp+1, " ",
                          remain_fileList," > ",unzip_trim_fileList1))
            system(paste0("cat ",unzip_trim_fileList1," | sed 's/ccs/ccs_",trim_times,"/g' > ",renamed_unzip_trim_fileList1))
            system(paste0("gzip < ",renamed_unzip_trim_fileList1, " > ",gzip_trim_fileList1))
            system(paste0('rm -rvf ',unzip_trim_fileList1))
            system(paste0('rm -rvf ',renamed_unzip_trim_fileList1)) 
            
            temp_fileList <- gsub(".fastq$","_remain.fastq",unzip_trim_fileList1)
            system(paste0("seqtk trimfq ","-b ",
                          trim_bp+1, " ",
                          remain_fileList," > ",temp_fileList))
            system(paste0("rm -rvf ", remain_fileList))
            
            ## filter out reads with less than trim_bp+1 so that it can be trimmed again 
            # https://www.biostars.org/p/66996/ 
            remain_fileList <- gsub("_remain.fastq","_filtered.fastq",temp_fileList)
            system(paste0("seqkit seq -m ", trim_bp+1," -g  ", temp_fileList, " > ", remain_fileList))
            system(paste0("rm -rvf ", temp_fileList))
            trim_times <- trim_times + 1
           
        }
        system(paste0("rm -rvf ", remain_fileList))
        
        if(trim_times>1) {
            setwd(fastqDir)
            gzip_trim_fileList1 <- gsub(".fastq$",
                                        paste0("_",trim_bp,"bp.fastq.gz"),unzip_fileList1)#gsub(".fastq$",".fastq.gz",renamed_unzip_trim_fileList1)
            system(paste0("cat * > ",gzip_trim_fileList1))
            setwd(wkdir)
        }
        # this happens for 
        unzip_trim_fileList1 <- gsub(".gz$","",gzip_trim_fileList1)
        filtered_fastq <- gsub(".fastq$","_filtered.fastq",unzip_trim_fileList1)
        prob_fastq <- gsub(".fastq$","_prob.fastq",unzip_trim_fileList1)
        trimmed_prob_fastq <- gsub(".fastq$","_trimmed_prob.fastq",unzip_trim_fileList1)
        system(paste0("filter_fastq.sh -i ", unzip_trim_fileList1,
                      " -o ",filtered_fastq, " -r ", prob_fastq))
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




set.seed(1)
prefix <- 'lr'
trim_bp <- 150
trim_function(trim_bp, sampleNames[g], sampleNames,sampleData, prefix = prefix,wkdir,annotationDir)

