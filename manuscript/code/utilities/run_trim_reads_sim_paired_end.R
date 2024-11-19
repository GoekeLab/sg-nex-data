#!/usr/bin/env Rscript
rm(list = ls())

library(GenomicAlignments) ##readGAlignments
library(AnnotationDbi)#loadDb
library(data.table)


require(docopt)
'Usage:
run_trim_reads_sim_paired_end.R [-g <g> ]

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


sampleData <- data.table(as.data.frame(read_xlsx(paste0(annoDir,"ONT Master Table.xlsx"), sheet = 1))) ## need to convert from tibble to data.frame
sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
                                                                                 runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
#sampleData$runName_combined <- gsub('-pre|-Pre','',sampleData$runName) # there are runs with multiple datasets that should be combined together 
sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(SG_NextData_Release=="Yes"),
                                       runName, 
                                       `GIS Library ID`)]
sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]

#sampleNames <- unique(sampleData$runName_combined)#
sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:112]

library(BiocParallel)

insertsizes_distribution <- fread("sr_stats.txt", skip = 49, fill = TRUE)
is <- insertsizes_distribution[V1 == "IS",1:5, with =FALSE]
isVec <- rep(as.numeric(is$V2), is$V3)
sim_LR_from_SR <- function(runnames,wkdir,txSeqDt,seq_size, insert_size,paired){
    np <- lapply(runnames, function(k){
        bam_path <- get_bam_file(k,wkdir)
        print("finish downloading bam file")
        bam_ranges <- read_in_bam_file(bam_path)
        print("finish reading bam file")
        
        # noprint <- lapply(1:10, function(t){
        #     sim_pos_data <- 
        #     generate_fastq(sim_pos_data,txSeq, t, fqFile)
        # })
        start_end_data <- process_data(bam_ranges, insert_size, seq_size)
        print("finish processing data")
        rm(bam_ranges)
        gc()
        myparameters <- bpparam()
        myparameters$workers = 8
        fqDir1 <- paste0(wkdir,"trim_reads/simLR/fastq/",k,"_1/")
        if(!dir.exists(fqDir1)) dir.create(fqDir1, recursive = TRUE)
        if(paired){
            fqDir2 <- paste0(wkdir,"trim_reads/simLR/fastq/",k,"_2/")
            if(!dir.exists(fqDir2)) dir.create(fqDir2, recursive = TRUE)
        }
        np <- bplapply(seq_len(max(start_end_data$seq_times)), function(t){
            sim_data <- sim_pos(start_end_data[seq_times>=t], t, txSeqDt,seq_size,paired)
            print(paste0("finish simulating sequence positions for ",t))
            #!#$"%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
            # quality_string <- paste(sample(c(letters, LETTERS, 0:9,
            # "*","&","^","[","]","%","$","#","!","+","~",":","j","{","}",">","<","|","?","="),
            # trim_bp+1,replace = TRUE), collapse = "")
            quality_string <- paste(rep("~",seq_size),collapse = "")
            sim_data[, quality := quality_string]
            final_data <- sim_data[,c("qname","seqChar","quality"),with = FALSE]
            rm(sim_data)
            gc()
            setnames(final_data, 1:3, c("Header","Sequence","Quality"))
            if(paired){
                fqFile1 <- paste0(fqDir1,t,"_1.fastq")
                microseq::writeFastq(setDF(final_data[grepl("_1",Header)]), fqFile1) # library(microseq)
                fqFile2 <- paste0(fqDir2,t,"_2.fastq")
                microseq::writeFastq(setDF(final_data[grepl("_2",Header)]), fqFile2) # library(microseq)
            }else{
                fqFile <- paste0(fqDir1,t,".fastq")
                microseq::writeFastq(setDF(final_data), fqFile) # library(microseq)
            }
            
        }, BPPARAM = myparameters)
        setwd(fqDir1)
        fqFile_final <- paste0(wkdir,"trim_reads/simLR/fastq/",k,"_1.fastq")
        system(paste0("cat * > ", fqFile_final))
        system(paste0("rm -rvf ",fqDir1))
        setwd(wkdir)
        system(paste0("gzip ", fqFile_final))
        fqFile_final <- paste0(fqFile_final,".gz")
        if(paired){
            setwd(fqDir2)
            fqFile_final2 <- paste0(wkdir,"trim_reads/simLR/fastq/",k,"_2.fastq")
            system(paste0("cat * > ", fqFile_final2))
            system(paste0("rm -rvf ",fqDir2))
            setwd(wkdir)
            system(paste0("gzip ", fqFile_final2))
            fqFile_final2 <- paste0(fqFile_final2,".gz")
        }
        print("finish generating simulated fastq file")
        
        mapDir <- paste0(wkdir,"trim_reads/simLR/map/",k,"/")
        if(!dir.exists(mapDir)) dir.create(mapDir, recursive = TRUE)
        salmonPath <- "salmon" # nolint: infix_spaces_linter.
        if(paired){
            system(paste0(salmonPath," quant -p 24 -i ",annotationDir,
                          " -l A -1 ",
                          fqFile_final,
                          " -2 ", fqFile_final2,
                          " --validateMappings ",
                          # " --fldMean  ", trim_bp+1," ", # 
                          # " --fldSD 1 ",
                          #" --seqBias ", " --gcBias ", " --posBias ",
                          " -o ", mapDir,"/transcripts_quant_biasCorrected"))
        }else{
            system(paste0(salmonPath," quant -p 24 -i ",annotationDir,
                          " -l A -r ",
                          fqFile_final,
                          " --validateMappings ",
                          " --fldMean  ", seq_size, # 
                          " --fldSD 1 ",
                          #" --seqBias ", " --gcBias ", " --posBias ",
                          " -o ", mapDir,"/transcripts_quant_biasCorrected"))
        }
        
        print("finish salmon mapping")
	system(paste0("rm -vf ",fqFile_final))
	if(paired) system(paste0("rm -vf ",fqFile_final2))
    })
    
}

get_bam_file <- function(k,wkdir){
      bamDir <- paste0(wkdir,'trim_reads/bam/',k,'/')
    if(!dir.exists(bamDir)){
        dir.create(bamDir, recursive = TRUE)
    }
    s3_bam_dir <- paste0('s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/',k,"/",k,".bam")
    system(paste0("aws s3 cp --no-sign-request ",s3_bam_dir," ", bamDir))
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
        rgl.dt <- data.table(qname =paste0(as.numeric(as.factor(names(rgl))),".",counter),
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
    readGrgList[, qname := as.numeric(as.factor(qname))]
    return(readGrgList)
}

# the final fastq should have trim_bp+1
process_data <- function(start_end_data,insert_size, seq_size,isVec){
    
    start_end_data <- start_end_data[(end-start)>(seq_size-1)]
    #start_end_data[, last_possible_start_position := pmax(start,end-(insert_size-1))]
    start_end_data[, seq_times := ceiling((end-start+1)/(seq_size-1))] # maximum number of reads can be sequenced 
    # only use primary alignments
    start_end_data <- start_end_data[flag %in% c(0,16)]
    return(start_end_data)
}

sim_pos <- function(data, t,txSeqDt,seq_size,paired){
    data[, width := end - start+1]
    # order simulated insert sizes 
    sim_insert_size <- sort(sample(isVec[(isVec>seq_size)],nrow(data), replace = TRUE))

    # assign according to width rank 
    data[order(width), insert_size := sim_insert_size]

    # for entries with insert size greater than width, resample according to a reduced set 
    if(nrow(data[width < insert_size])>0){
        data[which(width < insert_size), insert_size := sample(isVec[(isVec>seq_size)&(isVec<width)],1), by = qname]
    }
    data[, last_possible_start_position := pmax(start,end-(insert_size-1))] # update last_possible_start_position
    data[,sim_start_pos := round(runif(1)*(last_possible_start_position-start)+start)] # simulate start position
    start_end_data_sim <- unique(data[,.(qname, sim_start_pos, end, tx_name, strand, insert_size, width)])
    
    if(paired){# if paired_end is generated, take the begin and end sequences and name it as two 
        start_end_data_sim2 <- copy(start_end_data_sim)
        start_end_data_sim1 <- copy(start_end_data_sim)
        start_end_data_sim1[, qname := paste0(qname,"_1.",t)] # forward
        start_end_data_sim1[, sim_end_pos := pmin(end, sim_start_pos + seq_size-1)]
        
        start_end_data_sim2[, qname := paste0(qname,"_2.",t)] # backward
        start_end_data_sim2[, sim_end_pos := pmin(end, sim_start_pos + insert_size - 1)]
        start_end_data_sim2[, sim_start_pos := pmax(0,sim_end_pos - seq_size + 1)]
        start_end_data_sim2[, strand := ifelse(strand == "-", "+","-")]
        
        start_end_data_sim <- do.call("rbind", list(start_end_data_sim1, start_end_data_sim2))
    }else{
        start_end_data_sim[, qname := paste0(qname,".",t)]
        start_end_data_sim[, sim_end_pos := pmin(end, sim_start_pos + seq_size-1)]
    }
    start_end_data_sim <- txSeqDt[start_end_data_sim, on = "tx_name"]
    start_end_data_sim[, seqChar := substring(seq, sim_start_pos, sim_end_pos)]
    start_end_data_sim[strand == "-", seqChar := microseq::reverseComplement(seqChar)] 
    # if negative strand, is this typo? or there is something special about this?
    start_end_data_sim <- start_end_data_sim[!is.na(seq)]
    start_end_data_sim <- start_end_data_sim[,.(qname, seqChar)]
    return(start_end_data_sim)
}

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
txSeq <- readDNAStringSet(file=paste0(annoDir,"hg38_sequins_SIRV_ERCCs_longSIRVs_cdna.fa"))
listNames <- unlist(lapply(strsplit(names(txSeq)," "),'[[',1))
txSeqDt <- data.table(tx_name = listNames, 
                      seq = as.character(txSeq))
set.seed(1)
prefix <- 'lr'
seq_size <- 151
paired <- TRUE
sim_LR_from_SR(sampleNames[g],wkdir,txSeqDt, seq_size, insert_size, paired)

