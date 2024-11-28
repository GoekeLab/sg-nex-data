# novel transcripts 
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


seList <- dir(".",
              pattern = "_20Jul2023_NDR", full.names = TRUE)

novelTxCount <- do.call("rbind",lapply(seList, function(r){
    se <- readRDS(r)
    txData <- data.table(as.data.frame(rowData(se)))
    txData[, tx_type := novelGene+novelTranscript]
    txData[, ndr_value := as.numeric(gsub(".*NDR|\\.rds","",r))]
    return(unique(txData[,list(n = .N), by = list(tx_type,ndr_value)]))
}))

novelTxCount$ndr_value <- c(rep(as.numeric(gsub(".*NDR|\\.rds",
      "",seList))[1:9],each = 3),
      gsub(".*NDR|\\.rds","",seList)[10],
      rep(as.numeric(gsub(".*NDR|\\.rds","",
        seList))[11],each = 3))

se <-  readRDS("bambuOutput_May25.rds")
txData <- data.table(as.data.frame(rowData(se)))
txData[, tx_type := novelGene+novelTranscript]
txData[, ndr_value := 0.316]
recommendedData <- unique(txData[,list(n = .N), by = list(tx_type,ndr_value)])
saveRDS(novelTxCount, file = "novelTxCount.rds")

# novelTxCount unique alignments
seList <- dir(".",
              pattern = "_20Jul2023_NDR", full.names = TRUE)

novelTxCount <- do.call("rbind",lapply(seList, function(r){
    se <- readRDS(r)
    txData <- data.table(as.data.frame(rowData(se)))
    txData[, tx_type := novelGene+novelTranscript]
    txData[, ndr_value := as.numeric(gsub(".*NDR|\\.rds","",r))]
    return(unique(txData[,list(n = .N), by = list(tx_type,ndr_value)]))
}))

novelTxCount$ndr_value <- c(rep(as.numeric(gsub(".*NDR|\\.rds",
                                                "",seList))[1:9],each = 3),
                            gsub(".*NDR|\\.rds","",seList)[10],
                            rep(as.numeric(gsub(".*NDR|\\.rds","",
                                                seList))[11],each = 3))

se <-  readRDS("bambuOutput_May25.rds")
txData <- data.table(as.data.frame(rowData(se)))
txData[, tx_type := novelGene+novelTranscript]
txData[, ndr_value := 0.316]
recommendedData <- unique(txData[,list(n = .N), by = list(tx_type,ndr_value)])
saveRDS(novelTxCount, file = "novelTxCount.rds")

lineData <- unique(novelTxCount[tx_type >0,list(n = sum(n)), 
                                by = list(ndr_value)])
lineData[, tx_type := "all"]
lineData <- do.call("rbind",list(lineData, novelTxCount[tx_type != 0]))
p_novelTx <- ggplot(novelTxCount[tx_type!=0], aes(x = ndr_value, y = n))+
   
    geom_bar(aes(fill = factor(tx_type, 
                                 levels = c(2,1),
                               labels = c("novel gene",
                                          "novel isoform"))),stat = "identity", 
    position = position_stack())+
    geom_line( data = lineData[tx_type == "all"], aes(group = tx_type),
               linetype = 2, col = "grey")+
    geom_point( data = lineData[tx_type == "all"], 
                aes(x = ndr_value, y = n), col = "grey", shape = 16)+
    geom_text( data = lineData[tx_type == "all"], 
               aes(label = n))+
    ylab("Number of novel transcripts")+
    xlab("NDR")+
    scale_fill_brewer(type = "qual", palette = 2, name = "")+
    theme_classic()+
    theme(legend.position = "top")

pdf(paste0("novelTxCount_varyNDR.pdf"),
    width = 6, height = 4)
p_novelTx
dev.off()


