########################################
# get quantification expression values #
########################################

# library set-up =============================
require(GenomicAlignments) ##readGAlignments
require(AnnotationDbi)#loadDb
require(data.table)#fast large dataset manipulation
require(readxl)


library(tximport) #
library(bambu)
general_list <- readRDS("general_list2023-04-27.rds")

samples <- general_list$samples
txLengths <- general_list$txLengths


sampleData_sr <- data.table(as.data.frame(read_xlsx('.', sheet = 3))) ## need to convert from tibble to data.frame
sampleData_sr <- sampleData_sr[!grepl('#',`ELM library ID`) &(!is.na(runName))&(!grepl("HEYA8.*H9",runName))]
sr_runNames <- sampleData_sr$runName

# bambu-lr =====================
seOutput <- readRDS("bambuOutput_May25.rds")
bambu_lr <- as.data.table(assays(seOutput)$counts, keep.rownames = TRUE)
geneTxTable <- as.data.table(rowData(seOutput))
setnames(bambu_lr, "rn","tx_name")
setnames(geneTxTable, c("GENEID","TXNAME"),c("gene_name","tx_name"))

bambu_lr <- melt(bambu_lr, id.vars = "tx_name", measure.vars = colnames(bambu_lr)[-1])
setnames(bambu_lr, c("variable","value"),c("runname","estimates"))

seGene <- transcriptToGeneExpression(seOutput)
sumCount <- apply(assays(seGene)$counts, 2,sum)
ntotalDt <- data.table(runname = names(sumCount),
                       ntotal = as.numeric(sumCount))


bambu_lr <- ntotalDt[bambu_lr, on = "runname"]
bambu_lr <- geneTxTable[bambu_lr, on = "tx_name"]
bambu_lr[, method := "bambu_lr"]
saveRDS(bambu_lr,file = "bambu_lr.rds")

# bambu-lr gene ===============================
seGene <- transcriptToGeneExpression(seOutput)
bambu_lr <- as.data.table(assays(seGene)$counts, keep.rownames = TRUE)
setnames(bambu_lr, "rn","gene_name")
bambu_lr <- melt(bambu_lr, id.vars = "gene_name", measure.vars = colnames(bambu_lr)[-1])
setnames(bambu_lr, c("variable","value"),c("runname","estimates"))
bambu_lr[, ntotal:=sum(estimates), by = runname]
bambu_lr[, method := "bambu_lr"]
saveRDS(bambu_lr,file = paste0("bambu_lr_gene.rds"))


# bambu-lr gene excluding non-full-length reads ===============================
counts <- as.data.table(assays(seOutput)$fullLengthCounts, keep.rownames = TRUE)
runnames <- colnames(counts)[-1]
colnames(counts)[-1] <- rename_duplicatedNames(runnames)
colData(seOutput)@rownames <- rename_duplicatedNames(colData(seOutput)@rownames)
counts <- melt(counts, id.vars = "rn", measure.vars = colnames(counts)[-1])
setnames(counts, "rn", "TXNAME")
rowDataSe <- as.data.table(rowData(seOutput))
counts <- rowDataSe[, .(TXNAME, GENEID)][counts, on = "TXNAME"]

incompatibleCounts <- metadata(seOutput)$incompatibleCounts
incompatibleCounts[, TXNAME := "incompatible"]
counts_incompatible <- melt(incompatibleCounts, id.vars = c("GENEID","TXNAME"), 
                            measure.vars = setdiff(colnames(incompatibleCounts), c("GENEID","TXNAME")))
# GENEID, TXNAME, variable, value
counts <- rbind(counts, counts_incompatible[variable %in% unique(counts$variable)])
counts[, valueGene := sum(value), by = list(variable, GENEID)]
seGene <- transcriptToGeneExpression(seOutput)
sumCount <- apply(assays(seGene)$counts, 2,sum)
ntotalDt <- data.table(runname = names(sumCount),
                       ntotal = as.numeric(sumCount))
setnames(counts, c("GENEID","variable"),c("gene_name","runname"))
counts <- ntotalDt[counts, on = "runname"]
counts[, valueGeneCPM := valueGene / ntotal * 10^6]

## counts
bambu_lr <- unique(counts[, .(gene_name, runname,valueGeneCPM)])
setnames(bambu_lr, c("valueGeneCPM"),c("normEst"))
bambu_lr[, method := "bambu_lr"]
saveRDS(bambu_lr,file = paste0("bambu_lr_gene_excluding_non_full_length_reads.rds"))


# salmon-lr=======================
tx2gene <- txLengths[,c(2,3)]
x <- 1

salmon_lr.dir <- "salmon_fastq6.4.2/count/"

sampleNames <- dir(salmon_lr.dir)
salmon_lr <- do.call('rbind',lapply(sampleNames,function(k){
        print(k)
        filePath <- sort(dir(paste0(salmon_lr.dir,k),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
        
        
        if(length(filePath)==0){
            return(NULL)
        }
        print(filePath)
        txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
        names(txi)
        
        salmon_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                                  abundance = txi$abundance[,1],
                                  counts = txi$counts[,1],
                                  length = txi$length[,1],
                                  countsFromAbundance = txi$countsFromAbundance)
        # short_read <- fread(filePath, header = TRUE)
        salmon_read[, runname:=k]
        return(salmon_read)
    }))




salmon_lr[, method:='salmon_lr']
salmon_lr[, ntotal:=sum(counts), by = runname]
setnames(salmon_lr, 'abundance','estimates')
salmon_lr[, `:=`(#counts = NULL,
                 length = NULL,
                 countsFromAbundance = NULL)]

salmon_lr <- geneTxTable[salmon_lr, on = 'tx_name']
salmon_lr[, TPM:=estimates]
salmon_lr[, estimates:=TPM/1000000*ntotal]
saveRDS(salmon_lr,file = "salmon_lr.rds")



# nanocount-lr====================================
nanocount.file <- dir(paste0("nanocount_fastq6.4.2/count/"), full.names = TRUE)
sampleNames <- tools::file_path_sans_ext(BiocGenerics::basename(nanocount.file))
nanoCountData <- do.call("rbind",lapply(seq_along(nanocount.file), function(i){
    nanoCountData <- fread(nanocount.file[i])
    
    nanoCountData[, runname:= tools::file_path_sans_ext(BiocGenerics::basename(nanocount.file[i]))]
    return(nanoCountData)
}))
setnames(nanoCountData, "transcript_name","tx_name")
nanoCountData <- geneTxTable[nanoCountData, on = "tx_name"]
nanoCountData[, tx_name := gsub("\\..*","", tx_name)]
nanoCountData[, method := "NanoCount_lr"]
nanoCountData[, ntotal := sum(est_count), by = list(runname)]
nanoCountData[, estimates:=tpm/1000000*ntotal]
saveRDS(nanoCountData, file = "nanocount_lr.rds")



# salmon-sr========================================
tx2gene <- txLengths[,c(2,3)]
x <- 1

salmon_sr <- do.call('rbind',lapply(sr_runNames,function(k){
    print(k)
    filePath <- sort(dir(paste0('.',k,'/transcripts_quant'),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
    if(length(filePath)==0){
        return(NULL)
    }
    print(filePath)
    txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
    names(txi)
    
    short_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                             abundance = txi$abundance[,1],
                             counts = txi$counts[,1],
                             length = txi$length[,1],
                             countsFromAbundance = txi$countsFromAbundance)
    short_read[, runname:=k]
    return(short_read)
}))

salmon_sr[, method:='salmon_sr']
salmon_sr[, ntotal:=sum(counts), by = runname]
setnames(salmon_sr, 'abundance','estimates')
salmon_sr[, `:=`(length = NULL,
                 countsFromAbundance = NULL)]

salmon_sr <- geneTxTable[salmon_sr, on = 'tx_name']
salmon_sr[, TPM:=estimates]
salmon_sr[, estimates:=TPM/1000000*ntotal]

saveRDS(salmon_sr,file = "salmon_sr.rds")


# salmon-sr bambuAnnotation ==========================
tx2gene <- txLengths[,c(2,3)]
x <- 1
sampleDir <- '.'
sr_runNames <- dir(sampleDir)
salmon_sr <- do.call('rbind',lapply(sr_runNames,function(k){
    print(k)
    filePath <- sort(dir(paste0(sampleDir,k,'/transcripts_quant'),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
    
    print(filePath)
    if(length(filePath)==0){
        return(NULL)
    }
    print(filePath)
    txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
    names(txi)
    
    short_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                             abundance = txi$abundance[,1],
                             counts = txi$counts[,1],
                             length = txi$length[,1],
                             countsFromAbundance = txi$countsFromAbundance)
    short_read[, runname:=k]
    return(short_read)
}))

salmon_sr[, method:='salmon_sr']
salmon_sr[, ntotal:=sum(counts), by = runname]
setnames(salmon_sr, 'abundance','estimates')
salmon_sr[, `:=`(counts = NULL,
                 length = NULL,
                 countsFromAbundance = NULL)]

salmon_sr <- geneTxTable[salmon_sr, on = 'tx_name']
salmon_sr[, TPM:=estimates]
salmon_sr[, estimates:=TPM/1000000*ntotal]

saveRDS(salmon_sr,file = "salmon_sr_bambuAnnotations.rds")




# rsem-sr====================
rsem.dir <- "."
rnames <- dir(rsem.dir)
rsem_sr_tx <- do.call("rbind",lapply(rnames, function(r){
    rsem_sr_tx <- fread(dir(paste0(rsem.dir,r), full.names = TRUE, pattern = "isoforms.results"))
    rsem_sr_tx <- rsem_sr_tx[,.(transcript_id, expected_count, TPM)]
    rsem_sr_tx[, runname := r]
    rsem_sr_tx[, ntotal := sum(expected_count)]
    setnames(rsem_sr_tx, c("transcript_id","TPM"), c("tx_name","estimates"))
    rsem_sr_tx[, expected_count := NULL]
    return(rsem_sr_tx)
}))

rsem_sr_tx <- geneTxTable[rsem_sr_tx, on = 'tx_name']
rsem_sr_tx[, method := "rsem_sr"]
rsem_sr_tx[, TPM:=estimates]
rsem_sr_tx[, estimates:=TPM/1000000*ntotal]
saveRDS(rsem_sr_tx,"rsem_sr_tx.rds")

rsem_sr_gene <- do.call("rbind",lapply(rnames, function(r){
    rsem_sr_gene <- fread(dir(paste0(rsem.dir,r), full.names = TRUE, pattern = "genes.results"))
    rsem_sr_gene <- rsem_sr_gene[,.(gene_id, expected_count, TPM)]
    rsem_sr_gene[, runname := r]
    rsem_sr_gene[, ntotal := sum(expected_count)]
    setnames(rsem_sr_gene, c("gene_id","TPM"), c("gene_name","estimates"))
    rsem_sr_gene[, expected_count := NULL]
    return(rsem_sr_gene)
}))


rsem_sr_gene[, method := "rsem_sr"]
rsem_sr_gene[, TPM:=estimates]
rsem_sr_gene[, estimates:=TPM/1000000*ntotal]
saveRDS(rsem_sr_gene,"rsem_sr_gene.rds")


# rsem-sr bambuAnnotation ==========================
rsem.dir <- "."
rnames <- dir(rsem.dir)
rsem_sr_tx <- do.call("rbind",lapply(rnames, function(r){
    rsem_sr_tx <- fread(dir(paste0(rsem.dir,r), full.names = TRUE, pattern = "isoforms.results"))
    rsem_sr_tx <- rsem_sr_tx[,.(transcript_id, expected_count, TPM)]
    rsem_sr_tx[, runname := r]
    rsem_sr_tx[, ntotal := sum(expected_count)]
    setnames(rsem_sr_tx, c("transcript_id","TPM"), c("tx_name","estimates"))
    rsem_sr_tx[, expected_count := NULL]
    return(rsem_sr_tx)
}))

rsem_sr_tx <- geneTxTable[rsem_sr_tx, on = 'tx_name']
rsem_sr_tx[, method := "rsem_sr"]
rsem_sr_tx[, TPM:=estimates]
rsem_sr_tx[, estimates:=TPM/1000000*ntotal]
saveRDS(rsem_sr_tx,"rsem_sr_bambufasta_tx.rds")

rsem_sr_gene <- do.call("rbind",lapply(rnames, function(r){
    rsem_sr_gene <- fread(dir(paste0(rsem.dir,r), full.names = TRUE, pattern = "genes.results"))
    rsem_sr_gene <- rsem_sr_gene[,.(gene_id, expected_count, TPM)]
    rsem_sr_gene[, runname := r]
    rsem_sr_gene[, ntotal := sum(expected_count)]
    setnames(rsem_sr_gene, c("gene_id","TPM"), c("gene_name","estimates"))
    rsem_sr_gene[, expected_count := NULL]
    return(rsem_sr_gene)
}))


rsem_sr_gene[, method := "rsem_sr"]
rsem_sr_gene[, TPM:=estimates]
rsem_sr_gene[, estimates:=TPM/1000000*ntotal]
saveRDS(rsem_sr_gene,"rsem_sr_bambufasta_gene.rds")


# trim-reads================================
dir_path <- c(
              "lr_150bpSingleEnd_1ts",
              "lr_300bpSingleEnd_1ts",
              "simLR",
              "simLR_PE",
              "simLR_pb")
np <- lapply(seq_along(dir_path)[5], function(path){
    temp_path <- dir_path[path] 
    if(temp_path == "simLR"){
        trim_path <- paste0("trim_reads/",temp_path,"/map/")
    }else if(temp_path == "simLR_PE"){
        trim_path <- paste0("trim_reads/simLR/map/")
    }else if(temp_path == "simLR_pb"){
        trim_path <- paste0("trim_reads/lr_150bpSingleEnd_1ts/")
    }else{
        trim_path <- paste0("trim_reads/",temp_path,"/")
    }
    
    filePathList <- dir(trim_path, full.names = TRUE, recursive = TRUE, pattern = "quant.sf")
    print(length(filePathList))
    tx2gene <- txLengths[,c(2,3)]
    x <- 1
    if(temp_path %in% c("simLR","simLR_PE")){
        trim_names <- dir(paste0(trim_path))
    }else{
        trim_names <- dir(paste0(trim_path,"/02_Mapping_matchedToGTF"))
    }
    
    salmon_sr <- do.call("rbind",lapply(trim_names,function(k){
            print(k)
        filePath <- filePathList[grep(k, filePathList)]
        if(length(filePath)==0){
            return(NULL)
        }
        print(filePath)
        txi <- try(tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)), TRUE)
        names(txi)
        if(class(txi)=="try-error") return(NULL)
        short_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                                 abundance = txi$abundance[,1],
                                 counts = txi$counts[,1],
                                 length = txi$length[,1],
                                 countsFromAbundance = txi$countsFromAbundance)
       
        short_read[, runname:=k]
        return(short_read)
    }))
    
    salmon_sr[, method:=paste0('trim_',temp_path)]
    salmon_sr[, ntotal:=sum(counts), by = runname]
    setnames(salmon_sr, 'abundance','estimates')
    salmon_sr[, `:=`(counts = NULL,
                     length = NULL,
                     countsFromAbundance = NULL)]
    
    salmon_sr <- geneTxTable[salmon_sr, on = 'tx_name']
    salmon_sr[, TPM:=estimates]
    salmon_sr[, estimates:=TPM/1000000*ntotal]
   
    saveRDS(salmon_sr,file = paste0("trim_",temp_path,".rds"))
})



trim_path <- "trim_reads/lr_150bp_sim_goodquality/"
filePathList <- dir(trim_path, full.names = TRUE, recursive = TRUE, pattern = "quant.sf")

tx2gene <- txLengths[,c(2,3)]
x <- 1
trim_names <- dir(trim_path)
salmon_sr <- do.call('rbind',lapply(trim_names,function(k){+
    print(k)
    filePath <- filePathList[grep(k, filePathList)]
    if(length(filePath)==0){
        return(NULL)
    }
    print(filePath)
    txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
    names(txi)
    
    short_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                             abundance = txi$abundance[,1],
                             counts = txi$counts[,1],
                             length = txi$length[,1],
                             countsFromAbundance = txi$countsFromAbundance)
    short_read[, runname:=k]
    return(short_read)
}))

salmon_sr[, method:='trim_lr_sim_150bp_salmon']
salmon_sr[, ntotal:=sum(counts), by = runname]
setnames(salmon_sr, 'abundance','estimates')
salmon_sr[, `:=`(counts = NULL,
                 length = NULL,
                 countsFromAbundance = NULL)]

salmon_sr <- geneTxTable[salmon_sr, on = 'tx_name']
salmon_sr[, TPM:=estimates]
salmon_sr[, estimates:=TPM/1000000*ntotal]

saveRDS(salmon_sr,file = "trim_lr_150bp_sim_goodquality.rds")

bp <- c(100,75,50,150)
np <- lapply(bp[2], function(b){
    trim_path <- paste0("trim_reads/sr_",b,"bp/")
    if(b != 100){
        trim_path <- paste0(trim_path,"02_Mapping_matchedToGTF/")
    }
    filePathList <- dir(trim_path, full.names = TRUE, recursive = TRUE, pattern = "quant.sf")
    
    tx2gene <- txLengths[,c(2,3)]
    x <- 1
    trim_names <- dir(trim_path)
    salmon_sr <- do.call('rbind',lapply(trim_names,function(k){
        print(k)
        filePath <- filePathList[grep(k, filePathList)]
        if(length(filePath)==0){
            return(NULL)
        }
        print(filePath)
        txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
        names(txi)
        
        short_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                                 abundance = txi$abundance[,1],
                                 counts = txi$counts[,1],
                                 length = txi$length[,1],
                                 countsFromAbundance = txi$countsFromAbundance)
        short_read[, runname:=k]
        return(short_read)
    }))
    
    salmon_sr[, method:=paste0('trim_sr_',b,'bp_salmon')]
    salmon_sr[, ntotal:=sum(counts), by = runname]
    setnames(salmon_sr, 'abundance','estimates')
    salmon_sr[, `:=`(counts = NULL,
                     length = NULL,
                     countsFromAbundance = NULL)]
    
    salmon_sr <- geneTxTable[salmon_sr, on = 'tx_name']
    salmon_sr[, TPM:=estimates]
    salmon_sr[, estimates:=TPM/1000000*ntotal]
    
    saveRDS(salmon_sr,file = paste0("trim_sr_",b,"bp_updated.rds"))
    
})


# spikein bambu-lr =================
seSpikein <- readRDS("bambuOutput_spikein_bam_ont_May22.rds")
bambu_lr_spikein <- as.data.table(assays(seSpikein)$counts, keep.rownames = TRUE)
geneTxTable <- as.data.table(rowData(seSpikein))
setnames(bambu_lr_spikein, "rn","tx_name")
setnames(geneTxTable, c("GENEID","TXNAME"),c("gene_name","tx_name"))

bambu_lr_spikein <- melt(bambu_lr_spikein, id.vars = "tx_name", measure.vars = colnames(bambu_lr_spikein)[-1])
setnames(bambu_lr_spikein, c("variable","value"),c("runname","estimates"))


seSpikeinGene <- transcriptToGeneExpression(seSpikein)


sumCount <- apply(assays(seSpikeinGene)$counts, 2,sum)
ntotalDt <- data.table(runname = names(sumCount),
                       ntotal = as.numeric(sumCount))
bambu_lr_spikein <- ntotalDt[bambu_lr_spikein, on = "runname"]



#bambu_lr_spikein[, ntotal:=sum(estimates), by = runname]
bambu_lr_spikein <- geneTxTable[bambu_lr_spikein, on = "tx_name"]
bambu_lr_spikein[, method := "bambu_lr"]
saveRDS(bambu_lr_spikein,file = "bambu_lr_spikein.rds")



# spikein bambu-lr gene ===============================
seSpikeinGene <- transcriptToGeneExpression(seSpikein)
bambu_lr_spikein <- as.data.table(assays(seSpikeinGene)$counts, keep.rownames = TRUE)
setnames(bambu_lr_spikein, "rn","gene_name")
bambu_lr_spikein <- melt(bambu_lr_spikein, id.vars = "gene_name", measure.vars = colnames(bambu_lr_spikein)[-1])
setnames(bambu_lr_spikein, c("variable","value"),c("runname","estimates"))
bambu_lr_spikein[, ntotal:=sum(estimates), by = runname]
bambu_lr_spikein[, method := "bambu_lr"]
saveRDS(bambu_lr_spikein,file = paste0("bambu_lr_spikein_gene.rds"))




# spikein salmon-lr =================
tx2gene <- txLengths[,c(2,3)]
x <- 1
spikein_salmon.dir <- 'spikein_fastq_Apr17/map_salmon_lr/'
spikein_samples <- dir(spikein_salmon.dir)
salmon_lr <- do.call('rbind',lapply(spikein_samples,function(k){
    print(k)
    filePath <- sort(dir(paste0(spikein_salmon.dir,k),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
    if(length(filePath)==0){
        return(NULL)
    }
    print(filePath)
    txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
    names(txi)
    
    salmon_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                              abundance = txi$abundance[,1],
                              counts = txi$counts[,1],
                              length = txi$length[,1],
                              countsFromAbundance = txi$countsFromAbundance)
    salmon_read[, runname:=k]
    return(salmon_read)
}))




salmon_lr[, method:='salmon_lr']
salmon_lr[, ntotal:=sum(counts), by = runname]
setnames(salmon_lr, 'abundance','estimates')
salmon_lr[, `:=`(counts = NULL,
                 length = NULL,
                 countsFromAbundance = NULL)]

salmon_lr <- geneTxTable[salmon_lr, on = 'tx_name']
salmon_lr[, TPM:=estimates]
salmon_lr[, estimates:=TPM/1000000*ntotal]
saveRDS(salmon_lr,file = "salmon_lr_spikein.rds")

# spikein nanocount lr ======================
nanocount.file <- dir(paste0("spikein_fastq_Apr17/map_nanocount_lr/"), full.names = TRUE, recursive = TRUE, pattern = "tsv")
sampleNames <- tools::file_path_sans_ext(BiocGenerics::basename(nanocount.file))
nanoCountData <- do.call("rbind",lapply(seq_along(nanocount.file), function(i){
    nanoCountData <- fread(nanocount.file[i])
    
    nanoCountData[, runname:= tools::file_path_sans_ext(BiocGenerics::basename(nanocount.file[i]))]
    return(nanoCountData)
}))
setnames(nanoCountData, "transcript_name","tx_name")
nanoCountData <- geneTxTable[nanoCountData, on = "tx_name"]
nanoCountData[, tx_name := gsub("\\..*","", tx_name)]
nanoCountData[, method := "NanoCount_lr"]
nanoCountData[, ntotal := sum(est_count)]
nanoCountData[, estimates:=tpm/1000000*ntotal]
saveRDS(nanoCountData,file = "nanocount_lr_spikein.rds")


# spikein salmon-sr ==================
tx2gene <- txLengths[,c(2,3)]
x <- 1
spikein_salmon.dir <- 'spikein/map/'
spikein_samples <- dir(spikein_salmon.dir)
salmon_sr <- do.call('rbind',lapply(spikein_samples,function(k){
    print(k)
    filePath <- sort(dir(paste0(spikein_salmon.dir,k,'/transcripts_quant'),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]

    if(length(filePath)==0){
        return(NULL)
    }
    print(filePath)
    txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
    names(txi)
    
    short_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                             abundance = txi$abundance[,1],
                             counts = txi$counts[,1],
                             length = txi$length[,1],
                             countsFromAbundance = txi$countsFromAbundance)
    short_read[, runname:=k]
    return(short_read)
}))

salmon_sr[, method:='salmon_sr']
salmon_sr[, ntotal:=sum(counts), by = runname]
setnames(salmon_sr, 'abundance','estimates')
salmon_sr[, `:=`(counts = NULL,
                 length = NULL,
                 countsFromAbundance = NULL)]

salmon_sr <- geneTxTable[salmon_sr, on = 'tx_name']
salmon_sr[, TPM:=estimates]
salmon_sr[, estimates:=TPM/1000000*ntotal]

saveRDS(salmon_sr,file = "salmon_sr_spikein.rds")

# spikein rsem-sr ===================
rsem.dir <- "spikein/map_rsem/"
rnames <- dir(rsem.dir)
rsem_sr_tx <- do.call("rbind",lapply(rnames, function(r){
    rsem_sr_tx <- fread(dir(paste0(rsem.dir,r), full.names = TRUE, pattern = "isoforms.results"))
    rsem_sr_tx <- rsem_sr_tx[,.(transcript_id, expected_count, TPM)]
    rsem_sr_tx[, runname := r]
    rsem_sr_tx[, ntotal := sum(expected_count)]
    setnames(rsem_sr_tx, c("transcript_id","TPM"), c("tx_name","estimates"))
    rsem_sr_tx[, expected_count := NULL]
    return(rsem_sr_tx)
}))

rsem_sr_tx <- geneTxTable[rsem_sr_tx, on = 'tx_name']
rsem_sr_tx[, method := "rsem_sr"]
rsem_sr_tx[, TPM:=estimates]
rsem_sr_tx[, estimates:=TPM/1000000*ntotal]
saveRDS(rsem_sr_tx,"rsem_sr_tx_spikein.rds")

rsem_sr_gene <- do.call("rbind",lapply(rnames, function(r){
    rsem_sr_gene <- fread(dir(paste0(rsem.dir,r), full.names = TRUE, pattern = "genes.results"))
    rsem_sr_gene <- rsem_sr_gene[,.(gene_id, expected_count, TPM)]
    rsem_sr_gene[, runname := r]
    rsem_sr_gene[, ntotal := sum(expected_count)]
    setnames(rsem_sr_gene, c("gene_id","TPM"), c("gene_name","estimates"))
    rsem_sr_gene[, expected_count := NULL]
    return(rsem_sr_gene)
}))


rsem_sr_gene[, method := "rsem_sr"]
rsem_sr_gene[, TPM:=estimates]
rsem_sr_gene[, estimates:=TPM/1000000*ntotal]
saveRDS(rsem_sr_gene,"rsem_sr_gene_spikein.rds")



# bambu-lr subset bam files =====================
type <- c("unique","best","seconBestAlignDiff10percent","onePrimary")
np <- lapply(type, function(k){
    seOutput <- readRDS(paste0("bambuOutput_LR_",k,"_30Jan2023.rds"))
    bambu_lr <- as.data.table(assays(seOutput)$counts, keep.rownames = TRUE)
    geneTxTable <- as.data.table(rowData(seOutput))
    setnames(bambu_lr, "rn","tx_name")
    setnames(geneTxTable, c("GENEID","TXNAME"),c("gene_name","tx_name"))
    
    bambu_lr <- melt(bambu_lr, id.vars = "tx_name", measure.vars = colnames(bambu_lr)[-1])
    setnames(bambu_lr, c("variable","value"),c("runname","estimates"))
    
    
    bambu_lr[, runname:=gsub("HCT116","Hct116",gsub("(.genome_alignment.sorted)|(_R1.sorted)|(_sorted)","",gsub("-pre|-Pre","",runname))), by = runname]
    bambu_lr[runname == "GIS_Hct116_cDNA_Rep2_Run4", runname:="GIS_Hct116_cDNA_Rep2_Run5"]
    bambu_lr[, estimates:=sum(estimates), by = list(runname, tx_name)]
    bambu_lr <- unique(bambu_lr[,.(tx_name, estimates, runname)])
    
    seGene <- transcriptToGeneExpression(seOutput)
    sumCount <- apply(assays(seGene)$counts, 2,sum)
    ntotalDt <- data.table(runname = names(sumCount),
                           ntotal = as.numeric(sumCount))
    ntotalDt[, runname:=gsub("HCT116","Hct116",gsub("(.genome_alignment.sorted)|(_R1.sorted)|(_sorted)","",gsub("-pre|-Pre","",runname))), by = runname]
    ntotalDt[runname == "GIS_Hct116_cDNA_Rep2_Run4", runname:="GIS_Hct116_cDNA_Rep2_Run5"]
    ntotalDt <- ntotalDt[, ntotal := sum(ntotal), by = list(runname)]
    ntotalDt <- unique(ntotalDt, by = NULL)
    bambu_lr <- ntotalDt[bambu_lr, on = "runname"]
    
    bambu_lr <- geneTxTable[bambu_lr, on = "tx_name"]
    bambu_lr[, method := "bambu_lr"]
    saveRDS(bambu_lr,file = paste0("bambu_lr",k,".rds"))
    
})


type <- c("unique","best","seconBestAlignDiff10percent","onePrimary")
np <- lapply(type, function(k){
    seOutput <- readRDS(paste0("bambuOutput_LR_",k,"_30Jan2023.rds"))
    seGene <- transcriptToGeneExpression(seOutput)
    bambu_lr <- as.data.table(assays(seGene)$counts, keep.rownames = TRUE)
    setnames(bambu_lr, "rn","gene_name")
    bambu_lr <- melt(bambu_lr, id.vars = "gene_name", measure.vars = colnames(bambu_lr)[-1])
    setnames(bambu_lr, c("variable","value"),c("runname","estimates"))
    bambu_lr[, runname := gsub(paste0("_",k,"_"), "", runname), by = runname]
    
    bambu_lr[, runname:=gsub("HCT116","Hct116",gsub("(.genome_alignment.sorted)|(_R1.sorted)|(_sorted)","",gsub("-pre|-Pre","",runname))), by = runname]
    bambu_lr[runname == "GIS_Hct116_cDNA_Rep2_Run4", runname:="GIS_Hct116_cDNA_Rep2_Run5"]
    bambu_lr[, estimates:=sum(estimates), by = list(runname, gene_name)]
    bambu_lr <- unique(bambu_lr[,.(gene_name, estimates, runname)])
    bambu_lr[, ntotal:=sum(estimates), by = runname]
    bambu_lr[, method := "bambu_lr"]
    saveRDS(bambu_lr,file = paste0("bambu_lr",k,"_gene.rds"))
    
})


# pacbio quantification results: by itself ===================
seOutput <- readRDS("bambuOutput_PacBio_May22.rds")
bambu_lr_pacbio <- as.data.table(assays(seOutput)$counts, keep.rownames = TRUE)
geneTxTable <- as.data.table(rowData(seOutput))
setnames(bambu_lr_pacbio, "rn","tx_name")
setnames(geneTxTable, c("GENEID","TXNAME"),c("gene_name","tx_name"))

bambu_lr_pacbio <- melt(bambu_lr_pacbio, id.vars = "tx_name", measure.vars = colnames(bambu_lr_pacbio)[-1])
setnames(bambu_lr_pacbio, c("variable","value"),c("runname","estimates"))

seGene <- transcriptToGeneExpression(seOutput)
sumCount <- apply(assays(seGene)$counts, 2,sum)
ntotalDt <- data.table(runname = names(sumCount),
                       ntotal = as.numeric(sumCount))


bambu_lr_pacbio <- ntotalDt[bambu_lr_pacbio, on = "runname"]
bambu_lr_pacbio <- geneTxTable[bambu_lr_pacbio, on = "tx_name"]
bambu_lr_pacbio[, method := "bambu_lr"]
saveRDS(bambu_lr_pacbio,file = "bambu_lr_pacbio.rds")

# pacbio gene ===============================
seGene <- transcriptToGeneExpression(seOutput)
bambu_lr_pacbio <- as.data.table(assays(seGene)$counts, keep.rownames = TRUE)
setnames(bambu_lr_pacbio, "rn","gene_name")
bambu_lr_pacbio <- melt(bambu_lr_pacbio, id.vars = "gene_name", measure.vars = colnames(bambu_lr_pacbio)[-1])
setnames(bambu_lr_pacbio, c("variable","value"),c("runname","estimates"))
bambu_lr_pacbio[, ntotal:=sum(estimates), by = runname]
bambu_lr_pacbio[, method := "bambu_lr"]
saveRDS(bambu_lr_pacbio,file = paste0("bambu_lr_pacbio_gene.rds"))


# salmon-lr-pacbio=======================
tx2gene <- txLengths[,c(2,3)]
x <- 1

salmon_lr.dir <- "salmon_fastq6.4.2/count/"

sampleNames <- dir(salmon_lr.dir, pattern = "PacBio")
salmon_lr <- do.call('rbind',lapply(sampleNames,function(k){
    print(k)
    filePath <- sort(dir(paste0(salmon_lr.dir,k),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
    
    
    if(length(filePath)==0){
        return(NULL)
    }
    print(filePath)
    txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
    names(txi)
    
    salmon_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                              abundance = txi$abundance[,1],
                              counts = txi$counts[,1],
                              length = txi$length[,1],
                              countsFromAbundance = txi$countsFromAbundance)
    salmon_read[, runname:=k]
    return(salmon_read)
}))




salmon_lr[, method:='salmon_lr']
salmon_lr[, ntotal:=sum(counts), by = runname]
setnames(salmon_lr, 'abundance','estimates')
salmon_lr[, `:=`(
    length = NULL,
    countsFromAbundance = NULL)]

salmon_lr <- geneTxTable[salmon_lr, on = 'tx_name']
salmon_lr[, TPM:=estimates]
salmon_lr[, estimates:=TPM/1000000*ntotal]
saveRDS(salmon_lr,file = "salmon_lr_pacbio.rds")



# nanocount-lr====================================
nanocount.file <- dir(paste0("nanocount_fastq6.4.2/count/"), full.names = TRUE, pattern = "PacBio")
sampleNames <- tools::file_path_sans_ext(BiocGenerics::basename(nanocount.file))
nanoCountData <- do.call("rbind",lapply(seq_along(nanocount.file), function(i){
    nanoCountData <- fread(nanocount.file[i])
    
    nanoCountData[, runname:= tools::file_path_sans_ext(BiocGenerics::basename(nanocount.file[i]))]
    return(nanoCountData)
}))
setnames(nanoCountData, "transcript_name","tx_name")
nanoCountData <- geneTxTable[nanoCountData, on = "tx_name"]
nanoCountData[, tx_name := gsub("\\..*","", tx_name)]
nanoCountData[, method := "NanoCount_lr"]
nanoCountData[, ntotal := sum(est_count), by = list(runname)]
nanoCountData[, estimates:=tpm/1000000*ntotal]
saveRDS(nanoCountData, file = "nanocount_lr_pacbio.rds")

# spikein bambu-lr-pacbio =================
seSpikein <- readRDS("bambuOutput_spikein_bam_May22.rds")
bambu_lr_spikein <- as.data.table(assays(seSpikein)$counts, keep.rownames = TRUE)
geneTxTable <- as.data.table(rowData(seSpikein))
setnames(bambu_lr_spikein, "rn","tx_name")
setnames(geneTxTable, c("GENEID","TXNAME"),c("gene_name","tx_name"))

bambu_lr_spikein <- melt(bambu_lr_spikein, id.vars = "tx_name", measure.vars = colnames(bambu_lr_spikein)[-1])
setnames(bambu_lr_spikein, c("variable","value"),c("runname","estimates"))

seSpikeinGene <- transcriptToGeneExpression(seSpikein)

sumCount <- apply(assays(seSpikeinGene)$counts, 2,sum)
ntotalDt <- data.table(runname = names(sumCount),
                       ntotal = as.numeric(sumCount))
bambu_lr_spikein <- ntotalDt[bambu_lr_spikein, on = "runname"]

bambu_lr_spikein <- geneTxTable[bambu_lr_spikein, on = "tx_name"]
bambu_lr_spikein[, method := "bambu_lr"]
saveRDS(bambu_lr_spikein,file = "bambu_lr_spikein_wtpacbio.rds")



# spikein bambu-lr gene ===============================
seSpikeinGene <- transcriptToGeneExpression(seSpikein)
bambu_lr_spikein <- as.data.table(assays(seSpikeinGene)$counts, keep.rownames = TRUE)
setnames(bambu_lr_spikein, "rn","gene_name")
bambu_lr_spikein <- melt(bambu_lr_spikein, id.vars = "gene_name", measure.vars = colnames(bambu_lr_spikein)[-1])
setnames(bambu_lr_spikein, c("variable","value"),c("runname","estimates"))
bambu_lr_spikein[, ntotal:=sum(estimates), by = runname]
bambu_lr_spikein[, method := "bambu_lr"]
saveRDS(bambu_lr_spikein,file = paste0("bambu_lr_spikein_gene_wtpacbio.rds"))



# public-data-set: short read ====================


#### processing in linux as local processing using links with fread fails very often but when process on server, very fast and smooth 
library(readxl)
library(data.table)
public_filepath <- "public_read_download_links.xlsx"
public_dt <- data.table(as.data.frame(read_xlsx(public_filepath)))

rsem_sr_tx <- do.call("rbind",lapply(seq_len(10), function(r){
    tmp_dt <-  public_dt[2*r-1]
    rsem_sr_tx <- fread(input = tmp_dt$transcript_tsv_rsem)
    rsem_sr_tx <- rsem_sr_tx[,.(transcript_id, gene_id, expected_count, TPM)]
    rsem_sr_tx[, runname := paste0("ENCODE_",tmp_dt$cellLine,"_",tmp_dt$rep)]
    rsem_sr_tx[, ntotal := sum(expected_count)]
    setnames(rsem_sr_tx, c("transcript_id","gene_id","TPM"), c("tx_name","gene_name","estimates"))
    rsem_sr_tx[, expected_count := NULL]
    return(rsem_sr_tx)
}))
saveRDS(rsem_sr_tx,"rsem_sr_tx_encode.rds")
rsem_sr_gene <- do.call("rbind",lapply(seq_len(10), function(r){
    tmp_dt <-  public_dt[2*r-1]
    rsem_sr_gene <- fread(input = tmp_dt$gene_tsv_rsem)
    rsem_sr_gene <- rsem_sr_gene[,.(gene_id, expected_count, TPM)]
    rsem_sr_gene[, runname := paste0("ENCODE_",tmp_dt$cellLine,"_",tmp_dt$rep)]
    rsem_sr_gene[, ntotal := sum(expected_count)]
    setnames(rsem_sr_gene, c("gene_id","TPM"), c("gene_name","estimates"))
    rsem_sr_gene[, expected_count := NULL]
    return(rsem_sr_gene)
}))
saveRDS(rsem_sr_gene,"rsem_sr_gene_encode.rds")


kallisto_sr_tx <- do.call("rbind",lapply(seq_len(10), function(r){
    tmp_dt <-  public_dt[2*r-1]
    kallisto_sr_tx <- fread(input = tmp_dt$transcript_tsv_kallisto)
    kallisto_sr_tx <- kallisto_sr_tx[,.(target_id, est_counts, tpm)]
    kallisto_sr_tx[, target_id := gsub("\\..*","",target_id)]
    kallisto_sr_tx[, runname := paste0("ENCODE_",tmp_dt$cellLine,"_",tmp_dt$rep)]
    kallisto_sr_tx[, ntotal := sum(est_counts)]
    setnames(kallisto_sr_tx, c("target_id","tpm"), c("tx_name","estimates"))
    kallisto_sr_tx[, est_counts := NULL]
    return(kallisto_sr_tx)
}))
saveRDS(kallisto_sr_tx,"kallisto_sr_tx_encode.rds")



## further processing locally 

rsem_sr_tx <- readRDS("rsem_sr_tx_encode.rds")
rsem_sr_tx[, tx_name :=gsub("tSpikein_|\\..*","",tx_name), by = tx_name]
rsem_sr_tx[, gene_name :=gsub("gSpikein_|\\..*","",gene_name), by = gene_name]
# there are 1600 transcripts with multiple version after ., checking one showing it is _PAR_Y version
# thus one additional step is done here to take the maximum for each transcript, it can also be done later 
rsem_sr_tx[, estimates := sum(estimates), by = list(tx_name, gene_name, runname)]
rsem_sr_tx <- unique(rsem_sr_tx, by = NULL)
rsem_sr_tx <- rsem_sr_tx[grep("ENST|ERCC|phiX",tx_name)] # ntotal is based on total sequencing depth including tRNAs 
rsem_sr_tx <- geneTxTable[rsem_sr_tx, on = c("tx_name","gene_name")]
rsem_sr_tx[, method := "rsem_sr_public"]
rsem_sr_tx[, TPM:=estimates]
rsem_sr_tx[, estimates:=TPM/1000000*ntotal] # strange thing though, using estimate sum to approximate the total sequencing depth 
saveRDS(rsem_sr_tx,"rsem_sr_tx_encode_final.rds")


rsem_sr_gene <- readRDS("rsem_sr_gene_encode.rds")
rsem_sr_gene[, gene_name :=gsub("gSpikein_|\\..*","",gene_name)]
rsem_sr_gene[, estimates := sum(estimates), by = list(gene_name, runname)]
rsem_sr_gene <- unique(rsem_sr_gene, by = NULL)
rsem_sr_gene <- rsem_sr_gene[grep("ENSG|ERCC|phiX",gene_name)] # ntotal is based on total sequencing depth including tRNAs 
rsem_sr_gene[, method := "rsem_sr_public"]
rsem_sr_gene[, TPM:=estimates]
rsem_sr_gene[, estimates:=TPM/1000000*ntotal]
saveRDS(rsem_sr_gene,"rsem_sr_gene_encode_final.rds")


kallisto_sr_tx <- readRDS("kallisto_sr_tx_encode.rds")
kallisto_sr_tx <- unique(kallisto_sr_tx, by = NULL)
kallisto_sr_tx <- geneTxTable[kallisto_sr_tx, on = 'tx_name']
kallisto_sr_tx[, method := "kallisto_sr_public"]
kallisto_sr_tx[, TPM:=estimates]
kallisto_sr_tx[, estimates:=TPM/1000000*ntotal]
saveRDS(kallisto_sr_tx,"kallisto_sr_tx_encode_final.rds")

# salmon lr post-q7 filtering ==============
tx2gene <- txLengths[,c(2,3)]
x <- 1

salmon_lr.dir <- "salmon_filter/count/"

sampleNames <- dir(salmon_lr.dir)
salmon_lr <- do.call('rbind',lapply(sampleNames,function(k){
    print(k)
    filePath <- sort(dir(paste0(salmon_lr.dir,k),pattern = 'quant.sf', recursive = TRUE, full.names = TRUE),decreasing = TRUE)[1]
    if(length(filePath)==0){
        return(NULL)
    }
    print(filePath)
    txi <- tximport(filePath, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = as.logical(x)) # requires 'rjson'
    names(txi)
    salmon_read <- data.table(tx_name = gsub('\\..*','',rownames(txi$abundance)),
                              abundance = txi$abundance[,1],
                              counts = txi$counts[,1],
                              length = txi$length[,1],
                              countsFromAbundance = txi$countsFromAbundance)
    salmon_read[, runname:=k]
    return(salmon_read)
}))


salmon_lr[, method:='salmon_lr_q7']
salmon_lr[, ntotal:=sum(counts), by = runname]
setnames(salmon_lr, 'abundance','estimates')
salmon_lr[, `:=`(
    length = NULL,
    countsFromAbundance = NULL)]

salmon_lr <- geneTxTable[salmon_lr, on = 'tx_name']
salmon_lr[, TPM:=estimates]
salmon_lr[, estimates:=TPM/1000000*ntotal]
saveRDS(salmon_lr,file = "salmon_lr_q7_filtering.rds")

