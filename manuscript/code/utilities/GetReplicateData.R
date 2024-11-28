#!/usr/bin/env Rscript
###########################################
# get replicate data for each metric type #
###########################################

###########################
## set-up                ##
###########################
library(readxl)
library(data.table)

require(docopt)
'Usage:
GetReplicateData.R [-g <g> -t <t>]

Options:
-g index
-t threshold
]' -> doc

opts <- docopt(doc)
print(opts)
nnn <- as.integer(opts$g)
ttt <- as.integer(opts$t)
###########################
## set-up libraries      ##
###########################
require(GenomicFeatures)
require(GenomicAlignments) ##readGAlignments
require(AnnotationDbi)#loadDb
require(data.table)#fast large dataset manipulation
require(readxl)
require(dplyr)


require(ggplot2)
require(RColorBrewer)
require(gridExtra)

library(gplots)
library(RColorBrewer)
library(limma)
library(ggpubr)

library(ComplexHeatmap)
library(circlize)
library(viridis)
gctorture(FALSE)
library(tximport) # 
source("utility_function.R")


###########################
## load general data     ##
###########################
cat('Setting working directory')
wkdir <- '.'
general_list <- readRDS("general_list2023-04-27.rds")
samples_wSpikein <- general_list$samples_wSpikein
cellLines <- general_list$cellLines
protocolCol <-  general_list$protocolCol
protocolVec <-  general_list$protocolVec
protocolLabel <- general_list$protocolLabel

txvec <- fread(paste0(".txList_matchingToGTF_wtChrIs.txt"), header = FALSE)

txvec <- gsub("\\..*","",txvec$V1)
ensemblAnnotations.transcripts <- copy(general_list$ensemblAnnotations.transcripts)
setnames(ensemblAnnotations.transcripts, "ensembl_gene_id","gene_name")
ensemblAnnotations.transcripts <- data.table(tx_name = txvec, status = TRUE)[ensemblAnnotations.transcripts, on = "tx_name"]
ensemblAnnotations.transcripts[is.na(status), status := FALSE]
ensemblAnnotations.transcripts[, all_in := all(status), by = gene_name]
genevec <- unique(ensemblAnnotations.transcripts[which(all_in)]$gene_name)
samples <- general_list$samples

################################
## get metrics for replicates ##
################################

source("utility_function.R")
methodNamesList <- CJ(lr = c("bambu_lr","salmon_lr"),
                      sr = c("rsem_sr","salmon_sr"),
                      gene = c(TRUE, FALSE))

library(BiocParallel)
bpParameters <- bpparam()
bpParameters$workers <- 1

dominant_typeData <- readRDS("dominant_typeData_25May2023.rds")
metric_type_id <- nnn
methodNameId <- ttt
runnamevec <- readRDS("runnamevec.rds")
methodNames <- c(methodNamesList[methodNameId]$lr, methodNamesList[methodNameId]$sr)
gene <- methodNamesList[methodNameId]$gene
print(paste0(c(methodNames,gene)))


#############################
## load total count data   ##
#############################
temp <- samples_wSpikein[,.(old_runname, runname)]
setnames(temp, "runname", "new_name")
if(gene){
    comDataGene <- readRDS(paste0(wkdir, "combinedExpressionDataGene_19June2023.rds"))
    comDataGene[, old_runname := runname]
    comDataGene <- temp[comDataGene, on = "old_runname"]
    comDataGene[!is.na(new_name), runname := new_name]
    comDataGene[, `:=`(old_runname = NULL, new_name = NULL)]
    replicateData <- process_replicate(dt = comDataGene[!grepl("PacBio",runname)], methodNames, gene,
                        samples, ensemblAnnotations.transcripts, genevec, runnamevec, majorMinor = FALSE, 
                        scatterPlot = FALSE, complexity = FALSE, expressionLevel = FALSE, metric_type_id = metric_type_id, bpParameters) # focus on protein coding genes
    
}else{
    comDataTranscript <- readRDS(paste0(wkdir, "combinedExpressionDataTranscript_19June2023.rds"))
    comDataTranscript[, old_runname := runname]
    comDataTranscript <- temp[comDataTranscript, on = "old_runname"]
    comDataTranscript[!is.na(new_name), runname := new_name]
    comDataTranscript[, `:=`(old_runname = NULL, new_name = NULL)]
    
    complexityInfo <- unique(comDataTranscript[,.(gene_name, tx_name)])
    complexityInfo[!is.na(gene_name), ntx := length(unique(tx_name)), by = gene_name]
    
    comDataTranscript <- complexityInfo[comDataTranscript, on = c("gene_name","tx_name")]
    replicateData <- process_replicate(dt = comDataTranscript[!grepl("PacBio",runname)], methodNames, gene, 
                                        samples, ensemblAnnotations.transcripts, genevec,
                                       runnamevec, majorMinor = TRUE, scatterPlot = FALSE, 
                                       complexity = FALSE, expressionLevel = FALSE, metric_type_id = metric_type_id,bpParameters) # focus on protein coding genes
    
}

saveRDS(replicateData, file = paste0("replicateResults/replicate_",c("transcript","gene")[as.numeric(gene)+1],"_expression_comparison_",methodNames[1],"_",methodNames[2],"_",
                                     c("cor","mae","mard","mard_mod","rmse")[metric_type_id],"_1Aug2023.rds"))









