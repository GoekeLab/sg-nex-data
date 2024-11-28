##only m6A sites with nreads >= 50 are considered unless specified
##consistency barplots -> done (nreads >= 50)
#$python create_matrix_ensg_enst.py new_SGNEx_bedgraphs_noreadfilter 50 
tab <- read.csv("Downloads/50_filtered_m6anet_ensg_enst_prob_withm6ACE_matrix.csv")
names <- paste(tab$gene_id,tab$chromosome,tab$genomic_position,sep="_")
tab <- tab[,c(8:30)]
ntab<-data.frame(names,tab)
write.table(ntab,"Downloads/50_filtered.csv",quote=F,row.names=F,sep=",")
#$python find_consistency/mean.py 50_filtered.csv 0.9
tab <- read.csv("Downloads/0.9_avg_bycellline_matrix.csv",header=T)
pdf("pt9_barplot_nread50_new_SGNEx_bedgraphs_noreadfilter.pdf", height=10, width=20)
par(mfrow=c(1,2))
test <- table(rowSums(tab >= 0.9,na.rm=F))
barplot(test[c(2,3,4,5,6,7,8)],main=">0.9 in replicates from the same cell line (NA eliminated)",xlab="number of cell lines",ylab="total number of positions")
test <- table(rowSums(tab >= 0.9,na.rm=T))
barplot(test[c(2,3,4,5,6,7,8)],main=">0.9 in replicates from the same cell line (NA kept)",xlab="number of cell lines",ylab="total number of positions")
dev.off()
## visualize cell line-specific 
celllinespecific <- na.omit(tab[rowSums(tab >= 0.9,na.rm=F) == 2,])
celllinespecific_counts <- colSums(celllinespecific >= 0.9)
labels <- c("A549","H9","HEYA8","Hct116","HepG2","K562","MCF7")
pdf("pt9_piechart_nread50_new_SGNEx_bedgraphs_noreadfilter.pdf", height=10, width=20)
pie(celllinespecific_counts[c(2,3,4,5,6,7,8)],labels=labels)
dev.off()
## remove H9 and HEYA8 test
# tab <- tab[,c(1,2,5,6,7,8)]
# pdf("pt9_barplot_nread50_new_SGNEx_bedgraphs_noreadfilter_noH9HEYA8.pdf", height=10, width=20)
# par(mfrow=c(1,2))
# test <- table(rowSums(tab >= 0.9,na.rm=F)) #there are more m6a sites than with H9 and HEYA8 b/c there are less NA being eliminated
# barplot(test[c(2,3,4,5,6)],main=">0.9 in replicates from the same cell line (NA eliminated)",xlab="number of cell lines",ylab="total number of positions")
# test <- table(rowSums(tab >= 0.9,na.rm=T))
# barplot(test[c(2,3,4,5,6)],main=">0.9 in replicates from the same cell line (NA kept)",xlab="number of cell lines",ylab="total number of positions")
# dev.off()
## use nreads >= 100 for testing
#$python create_matrix_ensg_enst.py new_SGNEx_bedgraphs_noreadfilter 100
# tab <- read.csv("Downloads/100_filtered_m6anet_ensg_enst_prob_withm6ACE_matrix.csv")
# names <- paste(tab$gene_id,tab$chromosome,tab$genomic_position,sep="_")
# tab <- tab[,c(8:30)]
# ntab<-data.frame(names,tab)
# write.table(ntab,"Downloads/100_filtered_for_testing.csv",quote=F,row.names=F,sep=",")
# #$python find_consistency/mean.py 100_filtered_for_testing.csv 0.9
# tab <- read.csv("Downloads/0.9_avg_bycellline_matrix.csv",header=T)
# pdf("pt9_barplot_nread100_new_SGNEx_bedgraphs_noreadfilter.pdf", height=10, width=20)
# par(mfrow=c(1,2))
# tab <- tab[,c(1,2,5,6,7,8)] ##exclude H9 and HEYA8
# test <- table(rowSums(tab >= 0.9,na.rm=F))
# barplot(test[c(2,3,4,5,6)],main=">0.9 in replicates from the same cell line (NA eliminated)",xlab="number of cell lines",ylab="total number of positions")
# test <- table(rowSums(tab >= 0.9,na.rm=T))
# barplot(test[c(2,3,4,5,6)],main=">0.9 in replicates from the same cell line (NA kept)",xlab="number of cell lines",ylab="total number of positions")
# dev.off()
##recreate the original barplot in the preprint
#$python create_matrix_ensg_enst.py bedgraphs_for_supplementary 100 
# tab <- read.csv("Downloads/100_filtered_m6anet_ensg_enst_prob_withm6ACE_matrix.csv")
# names <- paste(tab$gene_id,tab$chromosome,tab$genomic_position,sep="_")
# tab <- tab[,c(8:19)]
# ntab<-data.frame(names,tab)
# write.table(ntab,"Downloads/100_filtered_for_testing.csv",quote=F,row.names=F,sep=",")
# #$python find_consistency/mean.py 100_filtered_for_testing.csv 0.9
# tab <- read.csv("Downloads/0.9_avg_bycellline_matrix.csv",header=T)
# pdf("pt9_original_bedgraphs_for_supplementary.pdf", height=10, width=20)
# par(mfrow=c(1,2))
# test <- table(rowSums(tab >= 0.9,na.rm=F))
# barplot(test[c(2,3,4,5,6)],main=">0.9 in replicates from the same cell line (NA eliminated)",xlab="number of cell lines",ylab="total number of positions")
# test <- table(rowSums(tab >= 0.9,na.rm=T))
# barplot(test[c(2,3,4,5,6)],main=">0.9 in replicates from the same cell line (NA kept)",xlab="number of cell lines",ylab="total number of positions")
# dev.off()

##heatmap positons vs samples -> done (no nreads filter b/c nreads >=50 leaves too many positions in certain cell lines empty)
#$python heatmap_top50_from_rank.py new_SGNEx_bedgraphs_noreadfilter rank_m6apositions.csv 0.9
library(gplots)
library(colorspace)
dat <- na.omit(read.csv("Downloads/top50_rank_heatmap.csv",header=T,row.names=1))
y<- as.matrix(dat)
#y<- as.matrix(read.table("ngsfcHeatmap.csv",header=T,row.names = 1, sep=","))
# hc <- hclust(as.dist(1-cor(y, method="pearson")), method="median")
hc <- hclust(dist(cor(y, method="pearson")), method="median")
mycol <- sequential_hcl(50000, "Blues 3",rev=T)
pdf("heatmap.pdf", height=500, width=20)
heatmap.2(y,dendrogram="row",col=mycol,density.info="none", trace="none",key=T, cexRow = 0.45, cexCol=0.75,ylab = "combinations", xlab = "GI scores")
heatmap.2(y, density.info="none",key=F,
          trace="none", cexRow = 0.45, cexCol=0.45,  key.xlab = "m6A probabilty",
          key.title = "Color Scale",margins=c(6,6),breaks=seq(0,1,0.00002), col=mycol)
dev.off()
# heatmap.2(y, dendrogram="row",  Colv=FALSE, col=mycol,scale="row", density.info="none", trace="none",key=T, cexRow = 0.45, cexCol=0.75,  key.xlab = "log2 fold change",key.title = "Color Scale",ylab = "combinations", xlab = "GI scores")

##rank plot -> done
#$python make_bedgraph_ensg.py m6anet_outputs (nreads >= 50)
#$python rank_numModPostions.py new_SGNEx_bedgraphs 0.9
rank_tab <- read.csv("Downloads/rank_m6apositions.csv",header=T)
library(ggplot2)
pdf("rank_plot.pdf", height=10, width=20)
rank_tab[rank_tab$gene == "ENSG00000136997",]
L <- ifelse(rank_tab$rank==40, "MYC", NA)
ggplot(aes(x=rank, y=number_of_sites_in_all_celllines, label=L), data=rank_tab) +
    geom_point() +
    geom_label() +
    xlab("rank") +
    ylab("number of m6A sites with a proability over 0.9") +
    ggtitle("Total number of m6A sites with probabily over 0.9 found in all cell lines") +
    theme_classic()
dev.off()

##raw current signal/dwell time boxplot -> done (nread >= 50)
# setwd("Downloads/separate_kmer/separated/")
# fns <- list.files(".")
library(ggplot2)
library(ggpubr)
# for (i in 1:length(fns)){
#   ofn<-gsub("csv", "png", fns[i])
#   tab <- read.table(fns[i],sep=",",header=T)
#   tab$m6a <- as.factor(tab$m6a)
#   p <- ggplot(tab, aes(x=sample, y=current, fill=m6a)) + 
#     geom_boxplot()
#   ggsave(ofn,width=45,height=20,units="in")
# }
tab <- read.table("Downloads/separate_kmer/separated/GGACT.csv",sep=",",header=T)
colnames(tab)<-c("cell_line","m6A","current")
tab$m6a <- as.factor(tab$m6A)

GGACT <- ggplot(tab, aes(x=cell_line, y=current, fill=m6a)) + geom_boxplot()+ggtitle("GGACT")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")+
    scale_fill_manual(values=c("lightgreen", "steelblue1"))
ggsave("GGACT_boxplot.pdf",width=5,height=5,units="in")

tab <- read.table("Downloads/separate_kmer/separated/GAACT.csv",sep=",",header=T)
colnames(tab)<-c("cell_line","m6A","current")
tab$m6a <- as.factor(tab$m6A)
GAACT <- ggplot(tab, aes(x=cell_line, y=current, fill=m6a)) + geom_boxplot()+ggtitle("GAACT")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")+
    scale_fill_manual(values=c("lightgreen", "steelblue1"))
tab <- read.table("Downloads/separate_kmer/separated/AGACT.csv",sep=",",header=T)
colnames(tab)<-c("cell_line","m6A","current")
tab$m6a <- as.factor(tab$m6A)
AGACT <- ggplot(tab, aes(x=cell_line, y=current, fill=m6a)) + geom_boxplot() +ggtitle("AGACT")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")+
    scale_fill_manual(values=c("lightgreen", "steelblue1"))
figure <- ggarrange(GAACT, AGACT,
                    ncol = 2, nrow = 1)
ggsave("GAACT_AGACT_boxplot.pdf",width=10,height=5,units="in")

##sequence logo -> done
#$ python seqlogo_all_samples.py m6anet_outputs 0.9
library(ggseqlogo)
library(ggplot2)
mod_tab <- read.table("Downloads/seqlogos/modified.csv",header=T)
unmod_tab <- read.table("Downloads/seqlogos/unmodified.csv",header=T)
kmer_list <- list(modified_DRACH = as.vector(mod_tab$motif),
                  unmodified_DRACH= as.vector(unmod_tab$motif))
pdf("sequence_logo.pdf")
ggseqlogo(kmer_list,method="prob")+ggtitle("modified (probability >= 0.9) vs unmodified DRACH")+xlab("position")
dev.off()

##correlation heatmap -> done ##excluding the Run2 replicates (low throughput)
#$python create_matrix.py m6anet_outputs 50
# setwd("Downloads/")
dat <- na.omit(read.table("Downloads/m6anet_prob_matrix.csv",sep=",",header=TRUE,row.names = 1))
# colnames(dat) <- c("A549_Rep5.Run1","A549_Rep6.Run1","Hct116_Rep2.Run1","Hct116_Rep2.Run4","Hct116_Rep3.Run3",
#                    "HepG2_Rep5.Run2","HepG2_Rep6.Run1","K562_Rep4.Run1","K562_Rep5.Run1","K562_Rep6.Run1",
#                    "MCF7_Rep3.Run1","MCF7_Rep4.Run1")
y <- cor(dat)
pdf("sample_correlation_heatmap.pdf", height=20,width=20, title = "sample correlation heatmap")
mycol <- sequential_hcl(50000, "Blues 3",rev=T)
cl <- c("A549","A549","H9","H9","H9","HEYA8","HEYA8","HEYA8","HEYA8","HEYA8",
        "Hct116","Hct116","Hct116","HepG2","HepG2","K562","K562","K562","MCF7","MCF7")
ha <- HeatmapAnnotation(cell_line=cl,col=list(cell_line=c("A549"="pink","K562"="palegreen","Hct116"="skyblue","MCF7"="limegreen","HepG2"="steelblue","H9"="tan1","HEYA8"="khaki")))
ra <- rowAnnotation(cell_line=cl,col=list(cell_line=c("A549"="pink","K562"="palegreen","Hct116"="skyblue","MCF7"="limegreen","HepG2"="steelblue","H9"="tan1","HEYA8"="khaki")))
Heatmap(y,col=mycol,bottom_annotation = ha,right_annotation = ra)
dev.off()

##odds matrix -> done ##excluding the Run2 replicates (low throughput)
#$python create_matrix.py m6anet_outputs 50
library(colorspace)
library(ComplexHeatmap)
data <- read.table("Downloads/m6anet_prob_matrix.csv",sep=",",header=TRUE)
# data <- read.table("Downloads/20_filtered_m6anet_prob_matrix.csv",sep=",",header=TRUE)
oddsMatrix <- matrix(NA, ncol=20, nrow=20, dimnames=list(colnames(data)[-1], colnames(data)[-1]))
for(i in 1:20){for(j in 1:20){ oddsMatrix[i,j] <- fisher.test(table(data[, i+1] > 0.9, data[, j+1] > 0.9))$estimate}}
oddsMatrix[oddsMatrix > 1000] <- 1000
diag(oddsMatrix) <- max(oddsMatrix[upper.tri(oddsMatrix)])
col <- sequential_hcl(50000, "Blues 3",rev=T)
pdf("samples_heatmap.pdf", height=20,width=20, title = "odd ratios heatmap")
cl <- c("A549","A549","H9","H9","H9","HEYA8","HEYA8","HEYA8","HEYA8","HEYA8",
        "Hct116","Hct116","Hct116","HepG2","HepG2","K562","K562","K562","MCF7","MCF7")
ha <- HeatmapAnnotation(cell_line=cl,col=list(cell_line=c("A549"="pink","K562"="palegreen","Hct116"="skyblue","MCF7"="limegreen","HepG2"="steelblue","H9"="tan1","HEYA8"="khaki")))
ra <- rowAnnotation(cell_line=cl,col=list(cell_line=c("A549"="pink","K562"="palegreen","Hct116"="skyblue","MCF7"="limegreen","HepG2"="steelblue","H9"="tan1","HEYA8"="khaki")))
Heatmap(oddsMatrix,col=col,bottom_annotation = ha,right_annotation = ra)
dev.off()

#plot tracks for MYC -> done (no read filter b/c nreads >=50 filters out all positions in region1 and region2)
#$python summarized_bedgraph.py m6anet_outputs [cell line]
library(Sushi)
library(gridExtra)
full_bed <- read.table("Downloads/Homo_sapiens.GRCh38.91.exon.bed",header=T)
# Import bedgraphs
AC1 <- read.table("Downloads/m6ace/collapsed_Hct116_m6ACE.bedGraph", sep = " ", skip=1)
colnames(AC1)<- c("chrom","chromstart","chromend","prob")
AC1<-AC1[order(AC1$chromstart),]
JM1 <- read.table("Downloads/new_SGNEx_bedgraphs_noreadfilter/A549.summarized.bedGraph", sep = " ", skip=1)
colnames(JM1)<- c("chrom","chromstart","chromend","prob")
JM1<-JM1[order(JM1$chromstart),]
JM2 <- read.table("Downloads/new_SGNEx_bedgraphs_noreadfilter/Hct116.summarized.bedGraph", sep = " ", skip=1)
colnames(JM2)<- c("chrom","chromstart","chromend","prob")
JM2<-JM2[order(JM2$chromstart),]
JM3 <- read.table("Downloads/new_SGNEx_bedgraphs_noreadfilter/HepG2.summarized.bedGraph", sep = " ", skip=1)
colnames(JM3)<- c("chrom","chromstart","chromend","prob")
JM3<-JM3[order(JM3$chromstart),]
JM4 <- read.table("Downloads/new_SGNEx_bedgraphs_noreadfilter/K562.summarized.bedGraph", sep = " ", skip=1)
colnames(JM4)<- c("chrom","chromstart","chromend","prob")
JM4<-JM4[order(JM4$chromstart),]
JM5 <- read.table("Downloads/new_SGNEx_bedgraphs_noreadfilter/MCF7.summarized.bedGraph", sep = " ", skip=1)
colnames(JM5)<- c("chrom","chromstart","chromend","prob")
JM5<-JM5[order(JM5$chromstart),]
JM6 <- read.table("Downloads/new_SGNEx_bedgraphs_noreadfilter/H9.summarized.bedGraph", sep = " ", skip=1)
colnames(JM6)<- c("chrom","chromstart","chromend","prob")
JM6<-JM6[order(JM6$chromstart),]
JM7 <- read.table("Downloads/new_SGNEx_bedgraphs_noreadfilter/HEYA8.summarized.bedGraph", sep = " ", skip=1)
colnames(JM7)<- c("chrom","chromstart","chromend","prob")
JM7<-JM7[order(JM7$chromstart),]
###############################################################################
chrom <- "8"
chromstart <- 127735400
chromend <- 127743000
##remove potential overlaps in collapsed m6ACE bedGraph
check<-AC1[AC1$chrom == chrom,]
check<-check[check$chromstart >= chromstart,]
check<-check[check$chromend <= chromend,]
check<- check[check$chromstart != 127736229,]
subset_bed <- full_bed[full_bed$chrom==paste0("chr",chrom,""),]
subset_bed <- subset_bed[subset_bed$start>=chromstart,]
subset_bed <- subset_bed[subset_bed$end<=chromend,]
pdf("zoomtracks_MYC.pdf", width = 100, height = 300)
## three regions
layout.matrix <- matrix(c(1,1,1,
                          2,10,18,
                          3,11,19,
                          4,12,20,
                          5,13,21,
                          6,14,22,
                          7,15,23,
                          8,16,24,
                          9,17,25)
                        ,9, 3, byrow = TRUE)
layout(mat = layout.matrix,
       heights = matrix(c(1,1,1,
                          0.75,0.75,0.75,
                          0.75,0.75,0.75,
                          0.75,0.75,0.75,
                          0.75,0.75,0.75,
                          0.75,0.75,0.75,
                          0.75,0.75,0.75,
                          0.75,0.75,0.75,
                          0.75,0.75,0.75)
                        ,9, 3, byrow = TRUE))
#par(mar=c(1,1,1,1))
plotGenes(subset_bed,chrom,chromstart,chromend,labeltext=T,bheight=0.1,
          plotgenetype="box",bentline=F, col = "black",fontsize = 1.2)
#MYC
zoomregion1 = c(127736050,127736700)
zoomregion2 = c(127738200,127739100)
zoomregion3 = c(127740300,127741400)
# zoomsregion(zoomregion1,wideextend=0.01)
## three zoomregions
zoomsregion(zoomregion1,wideextend=0.01,offsets=c(0,0.504))
zoomsregion(zoomregion2,wideextend=0.01,offsets=c(0.503,0))
zoomsregion(zoomregion3,wideextend=0.16,offsets=c(0,0))
plotBedgraph(signal=check, chrom=chrom, chromstart=zoomregion1[1], chromend=zoomregion1[2], lwd = 1, main="Hct116_m6ACE",color = "darkgreen",range=c(0,1))
zoombox(zoomregion1)
plotBedgraph(signal=JM1, chrom=chrom, chromstart=zoomregion1[1], chromend=zoomregion1[2], lwd = 1, main="A549_m6Anet predictions",range=c(0,1))
zoombox(zoomregion1)
plotBedgraph(signal=JM2, chrom=chrom, chromstart=zoomregion1[1], chromend=zoomregion1[2], lwd = 1, main="Hct116 m6Anet predictions",range=c(0,1))
zoombox(zoomregion1)
plotBedgraph(signal=JM3, chrom=chrom, chromstart=zoomregion1[1], chromend=zoomregion1[2], lwd = 1, main="HepG2 m6Anet predictions",range=c(0,1))
zoombox(zoomregion1)
plotBedgraph(signal=JM4, chrom=chrom, chromstart=zoomregion1[1], chromend=zoomregion1[2], lwd = 1, main="K562 m6Anet predictions",range=c(0,1))
zoombox(zoomregion1)
plotBedgraph(signal=JM5, chrom=chrom, chromstart=zoomregion1[1], chromend=zoomregion1[2], lwd = 1, main="MCF7 m6Anet predictions",range=c(0,1))
zoombox(zoomregion1)
plotBedgraph(signal=JM6, chrom=chrom, chromstart=zoomregion1[1], chromend=zoomregion1[2], lwd = 1, main="H9 m6Anet predictions",range=c(0,1))
zoombox(zoomregion1)
plotBedgraph(signal=JM7, chrom=chrom, chromstart=zoomregion1[1], chromend=zoomregion1[2], lwd = 1, main="HEYA8 m6Anet predictions",range=c(0,1))
zoombox(zoomregion1)

plotBedgraph(signal=check, chrom=chrom, chromstart=zoomregion2[1], chromend=zoomregion2[2], lwd = 1, main="Hct116_m6ACE",color = "darkgreen")
zoombox(zoomregion2)
plotBedgraph(signal=JM1, chrom=chrom, chromstart=zoomregion2[1], chromend=zoomregion2[2], lwd = 1, main="A549_m6Anet predictions")
zoombox(zoomregion2)
plotBedgraph(signal=JM2, chrom=chrom, chromstart=zoomregion2[1], chromend=zoomregion2[2], lwd = 1, main="Hct116 m6Anet predictions")
zoombox(zoomregion2)
plotBedgraph(signal=JM3, chrom=chrom, chromstart=zoomregion2[1], chromend=zoomregion2[2], lwd = 1, main="HepG2 m6Anet predictions")
zoombox(zoomregion2)
plotBedgraph(signal=JM4, chrom=chrom, chromstart=zoomregion2[1], chromend=zoomregion2[2], lwd = 1, main="K562 m6Anet predictions")
zoombox(zoomregion2)
plotBedgraph(signal=JM5, chrom=chrom, chromstart=zoomregion2[1], chromend=zoomregion2[2], lwd = 1, main="MCF7 m6Anet predictions")
zoombox(zoomregion2)
plotBedgraph(signal=JM6, chrom=chrom, chromstart=zoomregion2[1], chromend=zoomregion2[2], lwd = 1, main="H9 m6Anet predictions")
zoombox(zoomregion2)
plotBedgraph(signal=JM7, chrom=chrom, chromstart=zoomregion2[1], chromend=zoomregion2[2], lwd = 1, main="HEYA8 m6Anet predictions")
zoombox(zoomregion2)

plotBedgraph(signal=check, chrom=chrom, chromstart=zoomregion3[1], chromend=zoomregion3[2], lwd = 1, main="Hct116_m6ACE",color = "darkgreen")
zoombox(zoomregion3)
plotBedgraph(signal=JM1, chrom=chrom, chromstart=zoomregion3[1], chromend=zoomregion3[2], lwd = 1, main="A549_m6Anet predictions")
zoombox(zoomregion3)
plotBedgraph(signal=JM2, chrom=chrom, chromstart=zoomregion3[1], chromend=zoomregion3[2], lwd = 1, main="Hct116 m6Anet predictions")
zoombox(zoomregion3)
plotBedgraph(signal=JM3, chrom=chrom, chromstart=zoomregion3[1], chromend=zoomregion3[2], lwd = 1, main="HepG2 m6Anet predictions")
zoombox(zoomregion3)
plotBedgraph(signal=JM4, chrom=chrom, chromstart=zoomregion3[1], chromend=zoomregion3[2], lwd = 1, main="K562 m6Anet predictions")
zoombox(zoomregion3)
plotBedgraph(signal=JM5, chrom=chrom, chromstart=zoomregion3[1], chromend=zoomregion3[2], lwd = 1, main="MCF7 m6Anet predictions")
zoombox(zoomregion3)
plotBedgraph(signal=JM6, chrom=chrom, chromstart=zoomregion3[1], chromend=zoomregion3[2], lwd = 1, main="H9 m6Anet predictions")
zoombox(zoomregion3)
plotBedgraph(signal=JM7, chrom=chrom, chromstart=zoomregion3[1], chromend=zoomregion3[2], lwd = 1, main="HEYA8 m6Anet predictions")
zoombox(zoomregion3)
dev.off()

# plot coverage plots -> done
library("reshape2")
library("raster")
library("Rsamtools")
library("GenomicAlignments")
library("ggplot2")
library("ggpubr")
CO1 <- c("Downloads/genomecov/bam/SGNex_A549_directRNA_replicate5_run1_genome.bam",
         "Downloads/genomecov/bam/SGNex_A549_directRNA_replicate6_run1_genome.bam")
CO2 <- c("Downloads/genomecov/bam/SGNex_Hct116_directRNA_replicate2_run1_genome.bam",
         "Downloads/genomecov/bam/SGNex_Hct116_directRNA_replicate2_run4_genome.bam")
CO3 <- c("Downloads/genomecov/bam/SGNex_HepG2_directRNA_replicate5_run2_genome.bam",
         "Downloads/genomecov/bam/SGNex_HepG2_directRNA_replicate6_run1_genome.bam")
CO4 <- c("Downloads/genomecov/bam/SGNex_K562_directRNA_replicate4_run1_genome.bam",
         "Downloads/genomecov/bam/SGNex_K562_directRNA_replicate5_run1_genome.bam",
         "Downloads/genomecov/bam/SGNex_K562_directRNA_replicate6_run1_genome.bam")
CO5 <- c("Downloads/genomecov/bam/SGNex_MCF7_directRNA_replicate3_run1_genome.bam",
         "Downloads/genomecov/bam/SGNex_MCF7_directRNA_replicate4_run1_genome.bam")
CO6 <- c("Downloads/genomecov/newSGNEx_bam/GIS_H9_directRNA_Rep1_Run1_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_H9_directRNA_Rep1_Run2_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_H9_directRNA_Rep2_Run1_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_H9_directRNA_Rep2_Run2_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_H9_directRNA_Rep3_Run1_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_H9_directRNA_Rep3_Run2_R1.sorted.bam")
CO7 <- c("Downloads/genomecov/newSGNEx_bam/GIS_HEYA8_directRNA_Rep1_Run1_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_HEYA8_directRNA_Rep1_Run2_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_HEYA8_directRNA_Rep2_Run1_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_HEYA8_directRNA_Rep2_Run2_R1.sorted.bam",
         "Downloads/genomecov/newSGNEx_bam/GIS_HEYA8_directRNA_Rep3_Run1_R1.sorted.bam")
###
plotBamCoverages <- function(bam.files,chr,start_pos,end_pos,ymax){
    geneGRangesObject <- GRanges(chr, IRanges(start_pos,end_pos))
    v<-BamViews(bam.files,bamRanges=geneGRangesObject)
    gappedAlign=readGAlignments(v)
    gappedAlign<- lapply(gappedAlign, keepSeqlevels ,value = seqnames(geneGRangesObject))
    geneSeqnames = seqnames(geneGRangesObject)
    geneStartIndex = start(geneGRangesObject)
    geneEndIndex = end(geneGRangesObject)
    geneRanges = geneStartIndex:geneEndIndex
    normalizedCoverageVector <- vector('list', length(gappedAlign))
    for (i in 1:length(gappedAlign)){
        normalizedCoverageVector[[i]] = as.double(coverage(gappedAlign[[i]], drop.D.ranges=F)[[1]][geneRanges])
    }
    xlimRange = range(geneRanges)
    conditionIdentifiers <- rep("none",times=length(bam.files))
    uniqueConditionIdentifiers = unique(conditionIdentifiers)
    ggplotObjectVector <- vector('list', length(uniqueConditionIdentifiers))
    for(i in 1:length(uniqueConditionIdentifiers)){
        sameConditionIndex = which(conditionIdentifiers %in% uniqueConditionIdentifiers[i])
        ylimRange = c(0, max(unlist(lapply(normalizedCoverageVector[sameConditionIndex],max))))
        tempCoverageVector = normalizedCoverageVector[sameConditionIndex]
        tempDataFrame = data.frame(geneRanges, tempCoverageVector)
        tempColumnNames = c('geneRangesInfo', paste0(as.character(conditionIdentifiers[sameConditionIndex]),seq(from=1,to=length(sameConditionIndex))))
        colnames(tempDataFrame) = tempColumnNames
        if(length(sameConditionIndex)==1){
            covMean = tempDataFrame[,tempColumnNames[2:length(tempColumnNames)]]
        }
        else{
            covMean =  rowMeans(tempDataFrame[,tempColumnNames[2:length(tempColumnNames)]])
        }
        averageTempDataFrame = data.frame(tempDataFrame$geneRangesInfo, covMean)
        colnames(averageTempDataFrame) = c('geneRangesInfo', 'averageValues')
        tempDataFrameMelt <- melt(tempDataFrame, id.vars="geneRangesInfo")
        ggplotObjectVector[[i]] <- ggplot()
        ggplotObjectVector[[i]] <- ggplot(tempDataFrameMelt, aes(x=geneRangesInfo,  y = value, col=variable)) + theme_bw() + coord_cartesian(ylim = c(0,ymax)) + geom_ribbon(aes(ymin=0,ymax=value, linetype=NA), fill="lightblue3", alpha=0.2) 
        ggplotObjectVector[[i]] <- ggplotObjectVector[[i]] + geom_line(data=averageTempDataFrame, aes(x=geneRangesInfo, y=averageValues), color="black", size=0.1)
        ggplotObjectVector[[i]] <- ggplotObjectVector[[i]] + labs(y=paste0('Normalized coverages'),x="")
        ggplotObjectVector[[i]] <- ggplotObjectVector[[i]] + theme(legend.position='none',panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
        return(ggplotObjectVector[[i]])
    }
}

##coverage plots
p1_1 <- plotBamCoverages(CO1,"8",127736050,127736700,150)
p1_2 <- plotBamCoverages(CO1,"8",127738200,127739100,150)
p1_3 <- plotBamCoverages(CO1,"8",127740300,127741400,150)

p2_1 <- plotBamCoverages(CO2,"8",127736050,127736700,25)
p2_2 <- plotBamCoverages(CO2,"8",127738200,127739100,25)
p2_3 <- plotBamCoverages(CO2,"8",127740300,127741400,25)

p3_1 <- plotBamCoverages(CO3,"8",127736050,127736700,400)
p3_2 <- plotBamCoverages(CO3,"8",127738200,127739100,400)
p3_3 <- plotBamCoverages(CO3,"8",127740300,127741400,400)

p4_1 <- plotBamCoverages(CO4,"8",127736050,127736700,250)
p4_2 <- plotBamCoverages(CO4,"8",127738200,127739100,250)
p4_3 <- plotBamCoverages(CO4,"8",127740300,127741400,250)

p5_1 <- plotBamCoverages(CO5,"8",127736050,127736700,150)
p5_2 <- plotBamCoverages(CO5,"8",127738200,127739100,150)
p5_3 <- plotBamCoverages(CO5,"8",127740300,127741400,150)

p6_1 <- plotBamCoverages(CO6,"8",127736050,127736700,50)
p6_2 <- plotBamCoverages(CO6,"8",127738200,127739100,50)
p6_3 <- plotBamCoverages(CO6,"8",127740300,127741400,50)

p7_1 <- plotBamCoverages(CO7,"8",127736050,127736700,350)
p7_2 <- plotBamCoverages(CO7,"8",127738200,127739100,350)
p7_3 <- plotBamCoverages(CO7,"8",127740300,127741400,350)
gl <- c(p1_1,p1_2,p1_3,p2_1,p2_2,p2_3,p3_1,p3_2,p3_3,p4_1,p4_2,p4_3,p5_1,p5_2,p5_3,
        p6_1,p6_2,p6_3,p7_1,p7_2,p7_3)
pdf("Downloads/coverage_plots_MYC.pdf", width = 100, height = 125)
ggarrange(p1_1,p1_2,p1_3,p2_1,p2_2,p2_3,p3_1,p3_2,p3_3,p4_1,p4_2,p4_3,p5_1,p5_2,p5_3,
          p6_1,p6_2,p6_3,p7_1,p7_2,p7_3,
          ncol = 3, nrow = 7)
dev.off()

