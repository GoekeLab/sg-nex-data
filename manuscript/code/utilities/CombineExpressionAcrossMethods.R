##################################
# combine quantification across methods #
##################################



require(GenomicFeatures)
require(GenomicAlignments) ##readGAlignments
require(AnnotationDbi)#loadDb
require(data.table)#fast large dataset manipulation
require(readxl)


library(tximport) # 
source("utility_function.R")

cat('Setting working directory')
wkdir <- 'wkdir'
general_list <- readRDS("general_list2023-04-20.rds")

samples <- general_list$samples
sr_runNames <- general_list$sr_runNames
txLengths <- general_list$txLengths


## overall count data ==============================
bambu_lr <- readRDS("bambu_lr.rds")
bambu_lr_gene <- readRDS("bambu_lr_gene.rds")
bambu_lr_pacbio <- readRDS("bambu_lr_pacbio.rds")
bambu_lr_gene_pacbio <- readRDS("bambu_lr_pacbio_gene.rds")

nanocount_lr <- readRDS("nanocount_lr.rds")
salmon_lr <- readRDS("salmon_lr.rds") 
nanocount_lr_pacbio <- readRDS("nanocount_lr_pacbio.rds")
salmon_lr_pacbio <- readRDS("salmon_lr_pacbio.rds")
salmon_sr <- readRDS("salmon_sr.rds")
rsem_sr_gene <- readRDS("rsem_sr_gene.rds")
rsem_sr_tx <- readRDS("rsem_sr_tx.rds")
merge.colnames <- c('tx_name','gene_name','estimates','runname','method','ntotal')
com_data <- rbind(bambu_lr[, merge.colnames, with =FALSE],
                  bambu_lr_pacbio[, merge.colnames, with =FALSE],
                  nanocount_lr[, merge.colnames, with =FALSE],
                  salmon_lr_pacbio[, merge.colnames, with =FALSE],
                  nanocount_lr_pacbio[, merge.colnames, with =FALSE],
                  salmon_lr[, merge.colnames, with =FALSE],
                  salmon_sr[, merge.colnames, with =FALSE],
                  rsem_sr_tx[, merge.colnames, with =FALSE])
com_data[, runname := as.character(runname)]
com_data[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]

com_data[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
com_data[, protocol_method:=paste0(protocol_general,'.',method)]
com_data[, short_read:=as.numeric(protocol_general == 'Illumina')]
com_data[, runname_method:=paste0(runname,'.',method)]
com_data[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]
com_data[, normEst:=estimates/ntotal*1000000]



rsem_sr <- rsem_sr_gene[, c("gene_name","estimates","runname","ntotal","method"), with =FALSE]
rsem_sr[, runname := as.character(runname)]
rsem_sr[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]

rsem_sr[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
rsem_sr[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]

bambu_lr <- bambu_lr_gene[, c("gene_name","estimates","runname","ntotal","method"), with =FALSE]
bambu_lr[, runname := as.character(runname)]
bambu_lr[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]

bambu_lr[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
bambu_lr[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]

bambu_lr_pacbio <- bambu_lr_gene_pacbio[, c("gene_name","estimates","runname","ntotal","method"), with =FALSE]
bambu_lr_pacbio[, runname := as.character(runname)]
bambu_lr_pacbio[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]

bambu_lr_pacbio[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
bambu_lr_pacbio[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]



com_data_gene <- unique(com_data[!(method %in% c("rsem_sr","bambu_lr")), 
    list(estimates=sum(estimates), ntotal = ntotal), by = list(gene_name, runname, protocol_general, cellLine, method)])
com_data_gene <- rbind(com_data_gene, 
                       rsem_sr[,.(gene_name, runname, protocol_general, cellLine, estimates, ntotal, method)],
                       bambu_lr[,.(gene_name, runname, protocol_general, cellLine, estimates, ntotal, method)],
                       bambu_lr_pacbio[,.(gene_name, runname, protocol_general, cellLine, estimates, ntotal, method)])
com_data_gene[, normEst:=estimates/ntotal*1000000, by = runname]
com_data_gene[, runname := gsub("k562","K562",runname)]

rm(list = c("bambu_lr","bambu_lr_gene","bambu_lr_pacbio","bambu_lr_gene_pacbio",
            "salmon_lr","rsem_sr_tx","rsem_sr_gene","salmon_sr","rsem_sr","nanocount_lr",
            "salmon_lr_pacbio","nanocount_lr_pacbio"))
gc()

saveRDS(com_data, file = paste0(wkdir, "combinedExpressionDataTranscript_19June2023.rds"))
saveRDS(com_data_gene, file = paste0(wkdir, "combinedExpressionDataGene_19June2023.rds"))

### trim reads ===========================
trim_lr_150bp_1ts <- readRDS("trim_lr_150bpSingleEnd_1ts.rds")
trim_lr_150bp_1ts[, method := "trim_lr_150bp_1ts"]
trim_lr_300bp_1ts <- readRDS("trim_lr_300bpSingleEnd_1ts.rds")
trim_lr_300bp_1ts[, method := "trim_lr_300bp_1ts"]
trim_lr_sim <- readRDS("trim_simLR.rds")
trim_lr_sim[, method := "trim_lr_sim"]
trim_lr_sim_pe <- readRDS("trim_simLR_PE.rds")
trim_lr_sim_pe[, method := "trim_lr_sim_pe"]
trim_lr_sim_pb <- readRDS("trim_simLR_pb.rds")
trim_lr_sim_pb[, method := "trim_lr_sim_pb"]
trim_sr_125bp <- readRDS("trim_sr_125bp.rds")
trim_sr_150bp <- readRDS("trim_sr_150bp.rds")
trim_sr_100bp <- readRDS("trim_sr_100bp.rds")
trim_sr_75bp <- readRDS("trim_sr_75bp_updated.rds")
trim_sr_50bp <- readRDS("trim_sr_50bp.rds")
rsem_sr_tx_public <- readRDS("rsem_sr_tx_encode_final.rds")
rsem_sr_gene_public <- readRDS("rsem_sr_gene_encode_final.rds")
kallisto_sr_tx_public <- readRDS("kallisto_sr_tx_encode_final.rds")
merge.colnames <- c('tx_name','gene_name','estimates','runname','method','ntotal')
com_data <- rbind(salmon_lr_pacbio[, merge.colnames, with =FALSE],
                  salmon_lr[, merge.colnames, with =FALSE],
                  salmon_sr[, merge.colnames, with =FALSE],
                  rsem_sr_tx_public[, merge.colnames, with = FALSE],
                  kallisto_sr_tx_public[, merge.colnames, with = FALSE],
                  trim_lr_150bp_1ts[, merge.colnames, with =FALSE],
                  trim_lr_300bp_1ts[, merge.colnames, with =FALSE],
                  trim_lr_sim[, merge.colnames, with =FALSE],
                  trim_lr_sim_pe[, merge.colnames, with =FALSE],
                  trim_lr_sim_pb[, merge.colnames, with =FALSE],
                  trim_sr_125bp[, merge.colnames, with =FALSE],
                  trim_sr_150bp[, merge.colnames, with =FALSE],
                  trim_sr_100bp[, merge.colnames, with =FALSE],
                  trim_sr_75bp[, merge.colnames, with =FALSE],
                  trim_sr_50bp[, merge.colnames, with =FALSE])
com_data[, runname := as.character(runname)]
com_data[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]
com_data[grep("Rep",protocol), protocol := "encode"]
com_data[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
com_data[, protocol_method:=paste0(protocol_general,'.',method)]
com_data[, short_read:=as.numeric(protocol_general %in% c('Illumina','encode'))]
com_data[, runname_method:=paste0(runname,'.',method)]
com_data[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]
com_data[, normEst:=estimates/ntotal*1000000]




rsem_sr <- rsem_sr_gene_public[, c("gene_name","estimates","runname","ntotal","method"), with =FALSE]
rsem_sr[, runname := as.character(runname)]
rsem_sr[, protocol:="encode"]
rsem_sr[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
rsem_sr[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]


com_data_gene <- unique(com_data[, 
                                 list(estimates=sum(estimates), ntotal = ntotal), by = list(gene_name, runname, protocol_general, cellLine, method)])
com_data_gene <- rbind(com_data_gene, 
                       rsem_sr[,.(gene_name, runname, protocol_general, cellLine, estimates, ntotal, method)])
com_data_gene[, normEst:=estimates/ntotal*1000000, by = runname]
com_data_gene[, runname := gsub("k562","K562",runname)]

rm(list = c("salmon_lr","rsem_sr_tx_public","rsem_sr_gene_public","salmon_sr","rsem_sr",
            "trim_lr_150bp_1ts","trim_lr_300bp_1ts","trim_lr_150bp","trim_lr_300bp",
            "trim_lr_sim","trim_lr_sim_pe","trim_lr_sim_pb",
            "trim_sr_125bp","trim_sr_150bp","trim_sr_100bp",
            "trim_sr_75bp","trim_sr_50bp","nanocount_lr"))
gc()

saveRDS(com_data, file = paste0(wkdir, "combinedExpressionDataTranscript_trimreads_20May2024.rds"))
saveRDS(com_data_gene, file = paste0(wkdir, "combinedExpressionDataGene_trimreads_20May2024.rds"))


## counts data for deseq2 only =====================
salmon_lr <- readRDS("salmon_lr.rds") 
salmon_sr <- readRDS("salmon_sr.rds")
merge.colnames <- c('tx_name','gene_name','estimates','runname','method','ntotal','counts')
com_data <- rbind(salmon_lr[, merge.colnames, with =FALSE],
                  salmon_sr[, merge.colnames, with =FALSE])
com_data[, runname := as.character(runname)]
com_data[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]

com_data[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
com_data[, protocol_method:=paste0(protocol_general,'.',method)]
com_data[, short_read:=as.numeric(protocol_general == 'Illumina')]
com_data[, runname_method:=paste0(runname,'.',method)]
com_data[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]
com_data[, normEst:=estimates/ntotal*1000000]





com_data_gene <- unique(com_data[, 
                                 list(estimates=sum(estimates), 
                                      counts = sum(counts),
                                      ntotal = ntotal), by = list(gene_name, runname, protocol_general, cellLine, method)])

com_data_gene[, normEst:=estimates/ntotal*1000000, by = runname]
com_data_gene[, runname := gsub("k562","K562",runname)]

rm(list = c("salmon_lr","salmon_sr"))
gc()

saveRDS(list(com_data,com_data_gene), file = paste0(wkdir, "combinedExpressionDataCountsList_5May2023.rds"))

## spikein-count data==================
bambu_lr <- readRDS("bambu_lr_spikein_wtpacbio.rds")
bambu_lr_gene <- readRDS("bambu_lr_spikein_gene_wtpacbio.rds")
salmon_lr <- readRDS("salmon_lr_spikein.rds") 
nanocount_lr <- readRDS("nanocount_lr_spikein.rds")
salmon_sr <- readRDS("salmon_sr_spikein.rds")
rsem_sr_gene <- readRDS("rsem_sr_gene_spikein.rds")
rsem_sr_tx <- readRDS("rsem_sr_tx_spikein.rds")
merge.colnames <- c('tx_name','gene_name','estimates','runname','method','ntotal')
com_data <- rbind(bambu_lr[, merge.colnames, with =FALSE],
                  salmon_lr[, merge.colnames, with =FALSE],
                  nanocount_lr[, merge.colnames, with =FALSE],
                  salmon_sr[, merge.colnames, with =FALSE],
                  rsem_sr_tx[, merge.colnames, with =FALSE])
com_data[, runname := gsub("__","_",gsub("_sorted.bam","",gsub(",","",runname)))]
com_data[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]
com_data[, normEst:=estimates/ntotal*1000000, by = runname]




rsem_sr <- rsem_sr_gene[, c("gene_name","estimates","runname","ntotal","method"), with =FALSE]
rsem_sr[, runname := gsub("__","_",gsub("_sorted.bam","",gsub(",","",runname)))]
rsem_sr[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]

bambu_lr <- bambu_lr_gene[, c("gene_name","estimates","runname","ntotal","method"), with =FALSE]
bambu_lr[, runname := gsub("__","_",gsub("_sorted.bam","",gsub(",","",runname)))]
bambu_lr[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]

com_data_gene <- unique(com_data[!(method %in% c("rsem_sr","bambu_lr")), 
                                 list(estimates=sum(estimates), ntotal = ntotal), 
                                 by = list(gene_name, runname, protocol, method)])
com_data_gene <- rbind(com_data_gene, 
                       rsem_sr[,.(gene_name, runname, protocol, estimates, ntotal, method)],
                       bambu_lr[,.(gene_name, runname, protocol, estimates, ntotal, method)])
com_data_gene[, normEst:=estimates/ntotal*1000000, by = runname]
saveRDS(list(com_data, com_data_gene), file = paste0(wkdir, "combinedExpressionDataList_spikein_25May.rds"))


## subset-count data ==================================
salmon_lr <- readRDS("salmon_lr_ONT_flag.rds") 
salmon_sr <- readRDS("salmon_sr.rds")
merge.colnames <- c('tx_name','gene_name','estimates','runname','method','ntotal')
type <- c("unique","best","seconBestAlignDiff10percent","onePrimary")
bambuList <- lapply(type, function(k){
    bambu_lr <- readRDS(paste0("bambu_lr",k,".rds"))
    bambu_lr[, runname := gsub(paste0("_",k,"_"),"", runname), by = runname]
    bambu_lr[, method := paste0("bambu_lr",k)]
    return(bambu_lr[, merge.colnames, with =FALSE])
})

bambuGeneList <- lapply(type, function(k){
    bambu_lr_gene <- readRDS(paste0("bambu_lr",k,"_gene.rds"))
    bambu_lr_gene[, runname := gsub(paste0("_",k,"_"),"", runname), by = runname]
    bambu_lr_gene[, method := paste0("bambu_lr",k)]
    bambu_lr <- bambu_lr_gene[, c("gene_name","estimates","runname","ntotal","method"), with =FALSE]
    bambu_lr[, runname := as.character(runname)]
    bambu_lr[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]
    
    bambu_lr[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
    bambu_lr[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]
    
    return(bambu_lr[,.(gene_name, runname, protocol_general, cellLine, estimates, ntotal, method)])
})


com_data <- list(salmon_lr[, merge.colnames, with =FALSE],
                  salmon_sr[, merge.colnames, with =FALSE])
com_data <- do.call("rbind", c(bambuList,
                               com_data))
com_data[, runname := as.character(runname)]
com_data[, protocol:=unlist(strsplit(runname, '_'))[3],by= runname]

com_data[, protocol_general:=gsub('PromethionDirect','direct',gsub('cDNAStranded','cDNA',protocol)), by = protocol]
com_data[, protocol_method:=paste0(protocol_general,'.',method)]
com_data[, short_read:=as.numeric(protocol_general == 'Illumina')]
com_data[, runname_method:=paste0(runname,'.',method)]
com_data[, cellLine := gsub("-EV","",gsub('k562','K562',unlist(strsplit(runname,'_'))[2])), by = runname]
com_data[, normEst:=estimates/ntotal*1000000]




com_data_gene <- unique(com_data[!grepl("bambu_lr",method), 
    list(estimates=sum(estimates), ntotal = ntotal), by = list(gene_name, runname, protocol_general, cellLine, method)])
com_data_gene <- do.call("rbind", c(list(com_data_gene), bambuGeneList))
com_data_gene[, normEst:=estimates/ntotal*1000000, by = runname]
com_data_gene[, runname := gsub("k562","K562",runname)]

rm(list = c("bambuList","bambuGeneList","salmon_lr","salmon_sr"))
gc()

saveRDS(list(com_data, com_data_gene), file = paste0(wkdir, "combinedExpressionDataList_subset.rds"))
