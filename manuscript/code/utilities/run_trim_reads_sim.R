#!/usr/bin/env Rscript
rm(list = ls())

library(GenomicAlignments) ##readGAlignments
library(AnnotationDbi)#loadDb
library(data.table)


require(docopt)
'Usage:
run_trim_reads.R [-g <g> ]

Options:
-g index
]' -> doc

opts <- docopt(doc)
print(opts)
g <- as.numeric(opts$g)



wkdir <- '.'
setwd(wkdir)
## need to manually edit the fasta file
#anno.file_wtChrIS <- "/mnt/ont/annotations/Grch38/ensembl-91/Homo_sapiens.GRCh38.91_wtChrIs.gtf"
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

library(BiocParallel)

sim_LR_from_SR <- function(runnames,wkdir,txSeqDt){
    np <- lapply(runnames, function(k){
        bam_path <- get_bam_file(k,wkdir)
        print("finish downloading bam file")
        bam_ranges <- read_in_bam_file(bam_path)
        print("finish reading bam file")
        fqDir <- paste0(wkdir,"trim_reads/simLR/fastq/",k,"/")
        if(!dir.exists(fqDir)) dir.create(fqDir, recursive = TRUE)
        # noprint <- lapply(1:10, function(t){
        #     sim_pos_data <- 
        #     generate_fastq(sim_pos_data,txSeq, t, fqFile)
        # })
        start_end_data <- process_data(bam_ranges)
        print("finish processing data")
        rm(bam_ranges)
        gc()
        myparameters <- bpparam()
        myparameters$workers = 8
        
        np <- bplapply(seq_len(max(start_end_data$seq_times)), function(t){
            fqFile <- paste0(fqDir,t,".fastq")
            
            sim_data <- sim_pos(start_end_data[seq_times>=t], t, txSeqDt,trim_bp)
            print(paste0("finish simulating sequence positions for ",t))
            #!#$"%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
            # quality_string <- paste(sample(c(letters, LETTERS, 0:9,
            # "*","&","^","[","]","%","$","#","!","+","~",":","j","{","}",">","<","|","?","="),
            # trim_bp+1,replace = TRUE), collapse = "")
            quality_string <- paste(rep("~",trim_bp+1),collapse = "")
            sim_data[, quality := quality_string]
            final_data <- sim_data[,c("qname","seqChar","quality"),with = FALSE]
            rm(sim_data)
            gc()
            setnames(final_data, 1:3, c("Header","Sequence","Quality"))
            microseq::writeFastq(setDF(final_data), fqFile) # library(microseq)
        }, BPPARAM = myparameters)
        fqFiles <- dir(fqDir)
        setwd(fqDir)
        fqFile_final <- paste0(wkdir,"trim_reads/simLR/fastq/",k,".fastq")
        system(paste0("cat * > ", fqFile_final))
        system(paste0("rm -rvf ",fqDir))
        setwd(wkdir)
        system(paste0("gzip ", fqFile_final))
        print("finish generating simulated fastq file")
        mapDir <- paste0(wkdir,"trim_reads/simLR/map/",k,"/")
        if(!dir.exists(mapDir)) dir.create(mapDir, recursive = TRUE)
        salmonPath <- "/mnt/dataSSD/software/salmon-1.9.0_linux_x86_64/bin/salmon"
        system(paste0(salmonPath," quant -p 48 -i ",annotationDir,
                      " -l A -r ",
                      fqFile_final,".gz ",
                      " --validateMappings ",
                      " --fldMean  ", trim_bp+1," ", # 
                      " --fldSD 1 ",
                      #" --seqBias ", " --gcBias ", " --posBias ",
                      " -o ", mapDir,"/transcripts_quant_biasCorrected"))
        print("finish salmon mapping")
	system(paste0("rm -vf ",fqFile_final,".gz"))
    })
    
}

get_bam_file <- function(k,wkdir){
    # sampleData_runName <- sampleData[sampleData$runName==k,]
    # rname <- sampleData_runName$`publicName (SGNex_CellLine_protocol_replicate1_run1)`
    bamDir <- paste0(wkdir,'trim_reads/bam/',k,'/')
    if(!dir.exists(bamDir)){
        dir.create(bamDir, recursive = TRUE)
    }
    s3_bam_dir <- paste0('s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/',rname)
    
    system(paste0('aws s3 cp --no-sign-request ',s3_bam_dir," ", bamDir))
    bam_path <- dir(bamDir, pattern = ".bam$", full.names = TRUE)
    return(bam_path)
}


read_in_bam_file <- function(bam_path){
    bamFile <- BamFile(bam_path, yieldSize = 1000000)
    bf <- open(bamFile)
    readGrgList <- NULL
    counter <- 1
    while (isIncomplete(bf)) {
        print(counter)
        rgl <-
            readGAlignments(bf,
                            param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment=FALSE,
                                                                    isDuplicate = FALSE),
                                                 what=c("flag","mapq")), #"qual", 
                            use.names = TRUE)
        rgl.dt <- data.table(qname =as.numeric(as.factor(names(rgl))),
                                     tx_name = gsub("\\..*","",as.character(seqnames(rgl))),
                                     start = start(rgl),
                                     strand = as.character(strand(rgl)),
                                     mapq = rgl@elementMetadata@listData$mapq,
                                     flag = rgl@elementMetadata@listData$flag,
                                     end = end(rgl))
        readGrgList <- do.call("rbind", list(readGrgList,rgl.dt))
        rm(rgl)
        rm(rgl.dt)
        gc()
        counter <- counter + 1
    }
    on.exit(close(bf))
    system(paste0("rm -vf ",bam_path))
    system(paste0("rm -vf ",bam_path,".bai"))
    return(readGrgList)
}

# the final fastq should have trim_bp+1
process_data <- function(start_end_data){
    start_end_data <- start_end_data[(end-start)>150]
    start_end_data[, last_possible_start_position := end-trim_bp]
    start_end_data[, seq_times :=floor((end-start+1)/150)]
    # only use primary alignments
    start_end_data <- start_end_data[flag %in% c(0,16)]
    
    return(start_end_data)
}

sim_pos <- function(data, t,txSeqDt,trim_bp){
    
    data[,sim_start_pos := round(runif(1)*(last_possible_start_position-start)+start)]
    start_end_data_sim <- unique(data[,.(qname, sim_start_pos, end, tx_name, strand)])
    start_end_data_sim[, qname := paste0(qname,".",t)]
    start_end_data_sim[, sim_end_pos := pmin(end, sim_start_pos + trim_bp)]
   
    start_end_data_sim <- txSeqDt[start_end_data_sim, on = "tx_name"]
    start_end_data_sim[, seqChar := substring(seq, sim_start_pos, sim_end_pos)]
    start_end_data_sim[strand == "=", seqChar := reverseComplement(seqChar)] 
    # if negative strand, is this typo? or there is something special about this?
    start_end_data_sim <- start_end_data_sim[!is.na(seq)]
    start_end_data_sim <- start_end_data_sim[,.(qname, seqChar)]
    return(start_end_data_sim)
}

# long read vs short read: two sources of difference, read length vs error rate, we want show that read length
# the thing is we don't know which one is more correct, but we know that if the error rate of long read is not a big issue, 
# the trimmed reads from long read should be more similar to short read as compared to original long read to short read 

# seq_pos <- min(start(tmp_range)):max(end(tmp_range))
# seqChar <- geneSeq[[match(as.character(unique(seqnames(tmp_range))), listNames)]][seq_pos]
generate_fastq <- function(sim_data, txSeq,t,fqFile){
    if(t == 1){
        file.create(fqFile)
    }
    np <- lapply(seq_len(nrow(sim_data)), function(x){
        tmp_data <- sim_data[x]
        seqChar <- txSeq[[match(tmp_data$tx_name, listNames)]][tmp_data$sim_start_pos:tmp_data$sim_end_pos]
        cat(paste0("@",tmp_data$qname,"-",t, " \n"), file=fqFile,append=TRUE)
        cat(paste0(seqChar, " \n"), file=fqFile,append=TRUE)
        cat("+ \n", file=fqFile,append=TRUE)
        cat("&^***  \n", file=fqFile,append=TRUE)
    })
}


# run it 
cat('Load transcript sequence information')
txSeq <- readDNAStringSet(file='hg38_sequins_SIRV_ERCCs_longSIRVs_cdna.fa')
listNames <- unlist(lapply(strsplit(names(txSeq)," "),'[[',1))
txSeqDt <- data.table(tx_name = listNames, 
                      seq = as.character(txSeq))
set.seed(1)
prefix <- 'lr'
trim_bp <- 150
sim_LR_from_SR(sampleNames[g],wkdir,txSeqDt)
#trim_function(trim_bp, sampleNames[g], sampleData, prefix = prefix,wkdir,annotationDir)


