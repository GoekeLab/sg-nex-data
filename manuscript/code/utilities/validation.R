# prepare blat sequences
library(readxl)
library(data.table)
dt <- read_xlsx("Copy of candidate_for_novel_transcripts_16May2024.xlsx")
dt <- data.table(as.data.frame(dt))


sink("validation.txt")
for(ii in  which(!is.na(dt$Correct))){
    seq1 <- ifelse(!is.na(dt[ii]$extended_seq_before),
                   dt[ii]$extended_seq_before,
                   dt[ii]$before_junc_seq)
    seq2 <- ifelse(!is.na(dt[ii]$extended_seq_after),
                   dt[ii]$extended_seq_after,
                   dt[ii]$after_junc_seq)
    cat("#### ",dt[ii]$tx_name," \n")
    cat(">seq1 \n")
    cat(seq1, " \n")
    cat(">seq2  \n")
    cat(seq2, " \n")
    cat(">primerF \n")
    cat(" \n")
    cat(">primerR \n")
    cat(" \n")
    cat(" \n")
}
sink()



seq_dir <- "SeqResult-0-212555/"
seqFiles <- list.files(seq_dir, pattern = "seq", full.names = TRUE)
txnames <- unique(unlist(lapply(seqFiles, function(x) paste0("BambuTx",gsub("TX","",unlist(strsplit(basename(x),"_"))[4])))))
txnames <- c(txnames[c(1:6)], paste0(txnames[7],"_1"),paste0(txnames[7],"_2"),txnames[8])
sink("validation_newset.txt")
for(ii in  seq_along(txnames)){
    temp_txname <- txnames[ii]
    print(temp_txname)
    seq1 <- ifelse(!is.na(dt[tx_name == temp_txname]$extended_seq_before),
                   dt[tx_name == temp_txname]$extended_seq_before,
                   dt[tx_name == temp_txname]$before_junc_seq)
    seq2 <- ifelse(!is.na(dt[tx_name == temp_txname]$extended_seq_after),
                   dt[tx_name == temp_txname]$extended_seq_after,
                   dt[tx_name == temp_txname]$after_junc_seq)
    cat("#### ",temp_txname," \n")
    cat(">seq1 \n")
    cat(seq1, " \n")
    cat(">seq2  \n")
    cat(seq2, " \n")
    cat(">primerF \n")
    cat(fread(seqFiles[grepl(paste0("_",gsub("BambuTx","", temp_txname),"_"),seqFiles)&grepl("F.seq",seqFiles)], header = FALSE)$V1," \n")
    cat(" \n")
    cat(">primerR \n")
    cat(fread(seqFiles[grepl(paste0("_",gsub("BambuTx","", temp_txname),"_"),seqFiles)&grepl("R.seq",seqFiles)], header = FALSE)$V1," \n")
    cat(" \n")
}
sink()



## lr vs sr check
dt <- read_xlsx("Copy of candidateTable_YF_28052024.xlsx")
dt <- data.table(as.data.frame(dt))


sink("validation_dpcr.txt")
for(ii in  which(!is.na(dt$Comment))){
    seq1 <- dt[ii]$uniqueSeq
    cat("#### ",dt[ii]$tx_name," \n")
    cat(">seq \n")
    cat(seq1, " \n")
    cat(">primerF \n")
    cat(" \n")
    cat(">primerR \n")
    cat(" \n")
    cat(" \n")
}
sink()
