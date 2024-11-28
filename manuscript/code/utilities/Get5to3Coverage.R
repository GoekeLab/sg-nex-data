#!/usr/bin/env Rscript
##################################
# get 5to3 coverage for each sample #
##################################

###########################
## set-up                ##
###########################

require(docopt)
'Usage:
coverage5To3.R [-g <g>  -l <l>]

Options:
-g index
-l sampleIndex
]' -> doc

opts <- docopt(doc)

g <- as.integer(opts$g)
l <- as.integer(opts$l)

print(paste(g, l))
#=============================================================
suppressMessages(require(GenomicAlignments)) ##readGAlignments
suppressMessages(require(AnnotationDbi))#loadDb
suppressMessages(require(data.table))#fast large dataset manipulation
suppressMessages(require(readxl))

txdbEnsembl91 <- loadDb('hg38_sequins_SIRV_ERCCs_longSIRVs-txdb.sqlite')
tx <- transcripts(txdbEnsembl91, use.names = TRUE)
txLengths <- transcriptLengths(txdbEnsembl91)
gr <- data.frame(chr = txLengths$tx_name, start = 1, end = txLengths$tx_len, strand = "+")
gr <- makeGRangesFromDataFrame(gr)
n_gr <- 1
gr_seqnames <- as.character(seqnames(gr))
if(g %in% c(110, 111)) { # inly sequins
    chris_gr <- grep("^R|ENST",gr_seqnames)
    gr <- gr[chris_gr]
    seqlevels(gr) <- as.character(seqnames(gr)) # reset seqlevels 
    n_gr <- length(gr)
    n_chunk <- ceiling(n_gr/2000)
    splitIndex_gr <- split(1:n_gr, rep(1:n_chunk, ceiling(n_gr/n_chunk))[1:n_gr])
}
if(g %in% c(98,99,100,104,105,106)){ # Heya8 big samples, only sequin and spliced sirvs 
     chris_gr <- c(grep("ENST|^R",gr_seqnames),
                  which(grepl("SIRV",gr_seqnames)&(nchar(gr_seqnames)<8)))
    gr <- gr[chris_gr]
    seqlevels(gr) <- as.character(seqnames(gr)) # reset seqlevels 
    n_gr <- length(gr)
    n_chunk <- ceiling(n_gr/2000)
    splitIndex_gr <- split(1:n_gr, rep(1:n_chunk, ceiling(n_gr/n_chunk))[1:n_gr])
}

if(g %in% c(101,102,103,107,108,109)){ #H9 big samples: no sequins
   chris_gr <- grep("^R",gr_seqnames)
    gr <- gr[-chris_gr]
    seqlevels(gr) <- as.character(seqnames(gr)) # reset seqlevels 
    n_gr <- length(gr)
    n_chunk <- ceiling(n_gr/2000)
    splitIndex_gr <- split(1:n_gr, rep(1:n_chunk, ceiling(n_gr/n_chunk))[1:n_gr])
}

save.dir <- "coverageSaveRDSApr17"
if(!dir.exists(save.dir)){
    dir.create(save.dir)
}

txInfo <- data.table(tx_name = tx$tx_name, 
                     strand = as.character(strand(tx)))


###########################
# functions               #
###########################
read_in_bam_file <- function(bam_path){
    bamFile <- BamFile(bam_path, yieldSize = 1000000)
    bf <- open(bamFile)
    readGrgList <- list()
    counter <- 1
    while (isIncomplete(bf)) {
        readGrgList[[counter]] <-
            readGAlignments(bf,
                            param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment=FALSE,
                                                                    isDuplicate = FALSE),
                                                 what=c("qual", "flag","mapq")),
                            use.names = TRUE)
        counter <- counter + 1
    }
    on.exit(close(bf))
    if (length(readGrgList) > 1) {
        readGrgList <- do.call(c, readGrgList)
    } else if(length(readGrgList) == 1){
        readGrgList <- readGrgList[[1]]
    }
    return(readGrgList)
}

###########################
## long read all samples ##
###########################


###########################
## long read pacbio samples ##
###########################

pacbio_data <-  data.table(as.data.frame(read_xlsx(paste0('.'), sheet = 2)))

bam.file <- pacbio_data$`bam-tx.path`
sampleNames <- pacbio_data$public_name
names(bam.file) <- sampleNames

local_path <- "RunBambu17Apr_coverage/"
if(!dir.exists(local_path)) dir.create(local_path)

noprint <- lapply(sampleNames[g], function(r){
    print(r)
    bam_file_temp <- bam.file[r]
    system(paste0("aws s3 cp ", bam_file_temp, " ",local_path))
    system(paste0("aws s3 cp ", bam_file_temp), ".bai ",local_path))
    local_bam_file <- dir(local_path, pattern = ".bam$", full.names = TRUE)
    print(local_bam_file)
        bam_ranges <- read_in_bam_file(local_bam_file)
        coverageV <- coverage(bam_ranges, drop.D.ranges=FALSE)
        names(coverageV) <- gsub("\\..*","",names(coverageV))
        nread_each <- length(bam_ranges)
        txvec <- gsub("\\..*","",unique(as.character(seqnames(bam_ranges))))
        n <- length(txvec)
        n_chunk <- ceiling(n/2000)
        splitIndex <- split(1:n, rep(1:n_chunk, ceiling(n/n_chunk))[1:n])
        pb <- progress::progress_bar$new(
            format = "  Genes [:bar] :percent in :elapsed",
            total = n_chunk, clear = FALSE, width= 60)
        pb$tick(0)
        
        covdt <- do.call("rbind",lapply(as.list(1:n_chunk), function(i){
            sid_chunk <- txvec[splitIndex[[i]]]
            covVec <- do.call("rbind",parallel::mclapply(sid_chunk,function(x){
                covNum <- as.double(coverageV[[x]])
                dt <- data.table(position = 1:length(covNum),
                                 coverage = covNum)
                dt[, rel_pos:=position/max(position)]
                dt[, pos_bin:=ceiling(rel_pos*100)]
                dt[, bin_count:=mean(coverage), by = pos_bin]
                dt[,tx_len:=length(covNum)]
                dt <- unique(dt[,.(pos_bin, bin_count,tx_len)])
                dt[, tx_name:=x]
                dt[, strand:=txInfo[tx_name == x]$strand]
                dt[, rel_bin_count:=bin_count/max(bin_count)]
                return(dt)
            },mc.set.seed = TRUE,mc.preschedule = TRUE,
            mc.silent = FALSE, mc.cores = 24))#parallel::detectCores()
            pb$tick()
            return(covVec)
        }))
        covdt[, nread:=nread_each]
        coverageV <-  NULL
        bam_ranges <- NULL
        covdt[, run_name:=r]
        if(!dir.exists(paste0(save.dir,"/",r))) dir.create(paste0(save.dir,"/",r))
        saveRDS(covdt, file = paste0(save.dir,"/",r,"/covDT",r,".rds"))
        print(file.exists(paste0(save.dir,"/",r,"/covDT",r,".rds")))
    system(paste0("rm -v ", local_bam_file,"*"))
})




############################
## spikein samples        ##
############################
bam.file <- dir("spikein_bam_tx_Apr17/", pattern = ".bam$", full.names = TRUE)
bam.file.basenames <- gsub("_(R1.sorted|sorted)","",gsub(".genome_alignment","",tools::file_path_sans_ext(BiocGenerics::basename(bam.file))))
names(bam.file) <- bam.file.basenames
noprint <- lapply(bam.file.basenames[g], function(r){
    print(r)
    bamFiles <- bam.file[[r]]
    if(length(bamFiles)==0|(all(!grepl(".bam$",bamFiles)))){
        return(NULL)
    }else{
        bamFiles <- bamFiles[grepl(".bam$",bamFiles)]
    }

    bam_ranges <- read_in_bam_file(bamFiles)
    coverageV <- coverage(bam_ranges, drop.D.ranges=FALSE)
    names(coverageV) <- gsub("\\..*","",names(coverageV))
    nread_each <- length(bam_ranges)
    txvec <- gsub("\\..*","",unique(as.character(seqnames(bam_ranges))))
    n <- length(txvec)
    n_chunk <- ceiling(n/2000)
    splitIndex <- split(1:n, rep(1:n_chunk, ceiling(n/n_chunk))[1:n])
    pb <- progress::progress_bar$new(
        format = "  Genes [:bar] :percent in :elapsed",
        total = n_chunk, clear = FALSE, width= 60)
    pb$tick(0)

    covdt <- do.call("rbind",lapply(as.list(1:n_chunk), function(i){
        sid_chunk <- txvec[splitIndex[[i]]]
        covVec <- do.call("rbind",parallel::mclapply(sid_chunk,function(x){
            covNum <- as.double(coverageV[[x]])
            dt <- data.table(position = 1:length(covNum),
                             coverage = covNum)
            dt[, rel_pos:=position/max(position)]
            dt[, pos_bin:=ceiling(rel_pos*100)]
            dt[, bin_count:=mean(coverage), by = pos_bin]
            dt[,tx_len:=length(covNum)]
            dt <- unique(dt[,.(pos_bin, bin_count,tx_len)])
            dt[, tx_name:=x]
            dt[, strand:=txInfo[tx_name == x]$strand]
            dt[, rel_bin_count:=bin_count/max(bin_count)]
            return(dt)
        },mc.set.seed = TRUE,mc.preschedule = TRUE,
        mc.silent = FALSE, mc.cores = 24))#parallel::detectCores()
        pb$tick()
        return(covVec)
    }))
    covdt[, nread:=nread_each]
    coverageV <-  NULL
    bam_ranges <- NULL
    covdt[, run_name:=r]
    saveRDS(covdt, file = paste0(save.dir,"/covDT",r,".rds"))
    print(file.exists(paste0(save.dir,"/covDT",r,".rds")))

})
