#!/usr/bin/env Rscript

rm(list = ls())

library(GenomicAlignments) ##readGAlignments
library(AnnotationDbi)#loadDb
library(data.table)


require(docopt)
'Usage:
run_trim_reads_with_error_incompatible.R [-g <g> -t <t> -s <s>]

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
sampleData <- data.table(as.data.frame(read_xlsx(paste0('ONT Master Table.xlsx'), sheet = 1))) ## need to convert from tibble to data.frame
sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
                                                                                 runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
#sampleData$runName_combined <- gsub('-pre|-Pre','',sampleData$runName) # there are runs with multiple datasets that should be combined together 
sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(SG_NextData_Release=="Yes"),
                                       runName, 
                                       `GIS Library ID`)]
sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]

#sampleNames <- unique(sampleData$runName_combined)#
sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:111]


gtf.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs_v5_reformatted.gtf"
fasta.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
anno_rds_file <- "bambuAnnotations.rds"
bambuAnnotations <-  readRDS(anno_rds_file)

# utility_functions
trim_function <- function(trim_times, trim_bp, runnames, sampleData, prefix = "sr",wkdir,annotationDir){
    
    np <- lapply(runnames, function(k){
        outList <- get_fastq(k,sampleData,wkdir, prefix, trim_bp)
       
        outList_trim <- trim_process(prefix,outList[[1]],trim_bp, outList[[2]],wkdir, trim_times)
        
        ## create map directory
        mapDir <- paste0(wkdir,'/02_MappingStar_matchedToGTF/',outList[[3]])
        map_function_star(mapDir, annotationDir,outList_trim[[1]],trim_bp, outList[[2]])
        
        ## run bambu to get se with distTable 
        bam.file <- dir(dirname(mapDir), pattern = ".bam$", full.names = TRUE)
        bam.file <- bam.file[grep(k, bam.file)]
        save.dir_rc <- paste0(wkdir,'/03_Bambu/rc/')
        save.dir_se <- paste0(wkdir,'/03_Bambu/se/')
        if(!dir.exists(save.dir_rc)) dir.create(save.dir_rc, recursive = TRUE)
        if(!dir.exists(save.dir_se)) dir.create(save.dir_se, recursive = TRUE)
        run_bambu(bam.file, bambuAnnotations, fasta.file, save.dir_rc, save.dir_se)
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
        s3_fastq_dir <- paste0('s3://sg-nex-data/data/sequencing_data_ont/fastq/',rname)
        system(paste0('aws s3 cp --no-sign-request ',s3_fastq_dir,' ', fastqDir))
    }
    
    if(prefix == "sr"){
        s3_fastq_dir <- paste0('s3://sg-nex-data/data/sequencing_data_illumina/fastq/',rname)
        system(paste0('aws s3 sync --no-sign-request ',s3_fastq_dir," ", fastqDir))
    }
    
    #cmdTMP <- paste0('aws s3 cp ',sampleData_runName$fastq.path,'/ ', fastqDir,'/ --profile ontdata.store.genome.sg --recursive --exclude "*.md5"')
    
    
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
    system(paste0("gunzip < ",fileList1," > ", unzip_fileList1))
    system(paste0('rm -rvf ',fileList1)) 
    for(kk in seq_len(trim_times)){
        unzip_trim_fileList1 <- gsub(".fastq$",paste0("_",trim_bp,"bp_",kk,".fastq"),unzip_fileList1)
        renamed_unzip_trim_fileList1 <- gsub(".fastq$","_renamed.fastq",unzip_trim_fileList1)
        gzip_trim_fileList1 <- gsub(".fastq$",".fastq.gz",renamed_unzip_trim_fileList1)
        system(paste0("seqtk trimfq ",c("-L ","-l ")[(prefix=="lr")+1],
                      trim_bp+1, " -q ",
                      #0.01," ",
                      rep(seq(0.01, 0.1, by = 0.01),4)[kk]," ", # it turns out that the default error rate threshold will give the same sequence all the time, so I want to see if change the error threshold will give different trimmed sequences
                      unzip_fileList1," > ",unzip_trim_fileList1)) #
        # extra step to confirm all reads less than 150bp, if not, extract non-150bp reads and trim again with -L
        # this happens for 
        filtered_fastq <- gsub(".fastq$","_filtered.fastq",unzip_trim_fileList1)
        prob_fastq <- gsub(".fastq$","_prob.fastq",unzip_trim_fileList1)
        trimmed_prob_fastq <- gsub(".fastq$","_trimmed_prob.fastq",unzip_trim_fileList1)
        #empty_fastq <- gsub(".fastq$","_empty.fastq",unzip_trim_fileList1)
        system(paste0("/home/cheny1/filter_fastq.sh -i ", unzip_trim_fileList1,
                      " -o ",filtered_fastq, " -r ", prob_fastq))
       
        if(file.size(prob_fastq)>0){
            system(paste0("rm -rvf ", unzip_trim_fileList1))
            system(paste0("seqtk trimfq -L ",
                          trim_bp+1, " ", # it turns out that the default error rate threshold will give the same sequence all the time, so I want to see if change the error threshold will give different trimmed sequences
                          prob_fastq," > ",trimmed_prob_fastq))
            system(paste0("cat ", filtered_fastq, " ", trimmed_prob_fastq, "  > ", unzip_trim_fileList1))
            system(paste0("rm -rvf ",filtered_fastq))
            system(paste0("rm -rvf ",prob_fastq))
            system(paste0("rm -rvf ",trimmed_prob_fastq))
            
        }else{
            system(paste0("rm -rvf ", filtered_fastq))
        }
        
        if(trim_times>1){
            system(paste0("cat ",unzip_trim_fileList1," | sed 's/ runid/",kk,
                          " runid/g' > ",renamed_unzip_trim_fileList1))
            system(paste0("gzip < ",renamed_unzip_trim_fileList1, " > ",gzip_trim_fileList1))
        }else{
            system(paste0("gzip < ",unzip_trim_fileList1, " > ",gzip_trim_fileList1))
        }
        
        system(paste0('rm -rvf ',unzip_trim_fileList1))
        system(paste0('rm -rvf ',renamed_unzip_trim_fileList1)) 
    }
    system(paste0('rm -rvf ',unzip_fileList1)) 
    if(trim_times>1) {
        setwd(fastqDir)
        gzip_trim_fileList1 <- gsub(".fastq$",
                                    paste0("_",trim_bp,"bp.fastq.gz"),unzip_fileList1)#gsub(".fastq$",".fastq.gz",renamed_unzip_trim_fileList1)
        system(paste0("cat * > ",gzip_trim_fileList1))
        setwd(wkdir)
    }
    print("finishing trimming process")
    return(gzip_trim_fileList1)
    
}

map_function_star <- function(mapDir, annotationDir,gzip_trim_fileList1,trim_bp,fastqDir){
    if(!dir.exists(mapDir)){
        dir.create(mapDir, recursive = TRUE)
    }
    parametersList <- list(indexDir = 'STAR-2.7.10a-Index/',
                           nThread = 24,
                           annotationFaFile = 'hg38_sequins_SIRV_ERCCs_longSIRVs.fa',
                           sjFile = 'hg38_sequins_SIRV_ERCCs_longSIRVs_v5_reformatted.gtf',
                           samTool = 'samtools',
                           star = "STAR"
    )
    
    system(paste0(parametersList[['star']],
                  ' --runThreadN ', parametersList[['nThread']],
                  ' --genomeDir ', parametersList[['indexDir']],
                  ' --genomeLoad LoadAndKeep',
                  ' --readFilesIn ',gzip_trim_fileList1,
                  ' --readFilesCommand gunzip -c ',
                  ' --outFileNamePrefix ',mapDir, 
                  #' --outSAMprimaryFalg AllBestScore',
                  ' --outMultimapperOrder Random ',
                  ' --outSAMattributes NH HI NM MD AS nM jM jI XS ',
                  ' --outSAMtype BAM Unsorted '))
    system(paste0("rm -rvf ",fastqDir))
    print("finish mapping")
}


run_bambu <- function(bam.file, bambuAnnotations, fasta.file, save.dir_rc, save.dir_se){
    library(bambu)
    system.time(bambuOutput <- bambu(reads = bam.file, #rcfile_check,#
                                     rcOutDir = save.dir_rc,
                                     annotations = bambuAnnotations,
                                     genome = fasta.file,
                                     returnDistTable = TRUE,
                                     opt.em = list(degradationBias = FALSE),
                                     stranded = FALSE, ncore = 1,
                                     #NDR = 0.247, # bambu recommended NDR is 0.311
                                     yieldSize = 1000000, verbose = TRUE))
    saveRDS(bambuOutput, file = paste0(save.dir_se,gsub("Aligned.out.bam$","",basename(bam.file)),"_seOutput.rds"))
}

set.seed(1)
prefix <- 'lr'
#sim_LR_from_SR(sampleNames[g],wkdir,txSeqDt)
wkdir <- paste0(wkdir,'trim_reads/',prefix,'_',trim_bp,'bpSingleEnd_',trim_times,'ts_incompatible')
trim_function(trim_times, trim_bp, sampleNames[g], sampleData, prefix = prefix,wkdir,annotationDir)


