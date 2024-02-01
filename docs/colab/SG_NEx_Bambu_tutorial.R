#==================#
# R script version #
#==================#
# Notes:
# In the case if google colab is not working, to be run locally to demonstrate how it works



### Installation =============================
# library installation 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("bambu", update = FALSE)

## later we will need to do some plotting. which require some plotting packages
BiocManager::install("ggbio")
BiocManager::install("ggplot2")
BiocManager::install("gridExtra")
BiocManager::install("circlize")
BiocManager::install("ComplexHeatmap")
BiocManager::install("grid")

# For post-bambu downstream analysis, we will need the following packages.
BiocManager::install("DESeq2", update = FALSE)
BiocManager::install("apeglm", update = FALSE)
BiocManager::install("DEXSeq", update = FALSE)
### Data access and preparation ==============
# install AWS for R
install.packages("aws.s3")
# set region for AWS
Sys.setenv("AWS_DEFAULT_REGION" = 'ap-southeast-1')



# create a directory to store the data
dir.create("bambu_tutorial")
setwd("./bambu_tutorial")


# download genome fasta file
aws.s3::save_object(
  object="data/data_tutorial/annotations/hg38_chr22.fa",
  bucket="sg-nex-data",
  region="ap-southeast-1",
  file="hg38_chr22.fa")
# download genome index fastai file
aws.s3::save_object(
  object="data/data_tutorial/annotations/hg38_chr22.fa.fai",
  bucket="sg-nex-data",
  region="ap-southeast-1",
  file="hg38_chr22.fa.fai")
# download gtf file
aws.s3::save_object(
  object="data/data_tutorial/annotations/hg38_chr22.gtf",
  bucket="sg-nex-data",
  region="ap-southeast-1",
  file="hg38_chr22.gtf")



### Download raw fastq files and perform alignment ================
fastq_list=c("A549_directRNA_sample1.fastq.gz",
           "A549_directRNA_sample2.fastq.gz",
           "A549_directRNA_sample3.fastq.gz",
           "HepG2_directRNA_sample1.fastq.gz",
           "HepG2_directRNA_sample2.fastq.gz",
           "HepG2_directRNA_sample3.fastq.gz")
for (fastq in fastq_list[1]){ # here we try with one sample
  aws.s3::save_object(
    object=paste0("data/data_tutorial/fastq/",fastq),
    bucket = "sg-nex-data",
    region="ap-southeast-1",
    file=basename(fastq))
}


# visualizing fastq file
system("zcat A549_directRNA_sample1.fastq.gz | head -8", intern = TRUE)


# install minimap2
cat(system("curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf - ", intern = TRUE), sep = "\n")
cat(system("./minimap2-2.26_x64-linux/minimap2 --version", intern = TRUE), sep = "\n")


# install samtools
cat(system("sudo apt install samtools", intern = TRUE), sep = "\n")
cat(system("samtools --version-only", intern = TRUE), sep = "\n")

# perform alignment
cat(system("./minimap2-2.26_x64-linux/minimap2 -ax splice -uf -k14 hg38_chr22.fa A549_directRNA_sample1.fastq.gz > A549_directRNA_sample1_test.sam ", intern = TRUE), sep = "\n")
cat(system("samtools view -Sb A549_directRNA_sample1_test.sam | samtools  sort -o A549_directRNA_sample1_test.bam ", intern = TRUE), sep = "\n")
cat(system("samtools index A549_directRNA_sample1_test.bam", intern = TRUE), sep = "\n")




### Download bam files from S3 bucket directly ===========
# download aligned bam files for A549 samples and HepG2 samples
bam_list=c("A549_directRNA_sample1.bam",
           "A549_directRNA_sample2.bam",
           "A549_directRNA_sample3.bam",
           "HepG2_directRNA_sample1.bam",
           "HepG2_directRNA_sample2.bam",
           "HepG2_directRNA_sample3.bam")
for (bam in bam_list){
  aws.s3::save_object(
    object=paste0("data/data_tutorial/bam/",bam),
    bucket = "sg-nex-data",
    region="ap-southeast-1",
    file=basename(bam))
}


# just had a quick check at the generated files, sligh change in bam file size might be due to minimap2 version
file.size("A549_directRNA_sample1_test.bam")
file.size("A549_directRNA_sample1.bam")

# Do not run unless you have specified the sample_alias with the sample of interest in the SG-NEx bukcet as recommended.
#aws.s3::save_object(
#  object="data/sequencing_data_ont/bam/genome/<sample_alias>",
#  bucket="sg-nex-data",
#  region="ap-southeast-1",
#  file="<sample_alias>")

#### Quality control ===================
cat(system("samtools stats A549_directRNA_sample1.bam | head -46 ", intern = TRUE), sep = "\n")



#### Preparing data for bambu ===============
# set work directory if you are in a different directory
setwd("bambu_tutorial")

# data preparation
library(bambu)
fa.file <- "hg38_chr22.fa"
gtf.file <- "hg38_chr22.gtf"
annotations <- prepareAnnotations(gtf.file) # This function creates a reference annotation object which is used for transcript discovery and quantification in Bambu.
samples.bam <- list.files("./", pattern = "[0-9].bam$", full.names = TRUE)

# before running commands, good to check if file all exists 
file.exists(fa.file)
file.exists(gtf.file)
print(samples.bam)


# running Bambu
se <- bambu(reads = samples.bam, annotations = annotations, genome = fa.file, ncore = 2)

# check the output object
assays(se) #returns the transcript abundance estimates as counts or CPM.

# visualize the quantification results for first 5 transcripts 
head(assays(se)$counts)


rowRanges(se) #returns a GRangesList (with genomic coordinates) with all annotated and newly discovered transcripts.
rowData(se) #returns additional information about each transcript such as the gene name and the class of the newly discovered transcript.

writeBambuOutput(se, path = "./output")



# Identify sample-specific novel transcripts============
se.filtered <- se[mcols(rowRanges(se))$novelTranscript == TRUE,]
#investigate high confidence novel isoforms
head(mcols(rowRanges(se.filtered))[order(mcols(rowRanges(se.filtered[]))$NDR),])


plotBambu(se, type = "annotation", gene_id = "BambuGene1")


se.filtered <- se[mcols(rowRanges(se))$novelTranscript == TRUE,]
#investigate high confidence novel isoforms in novel genes
se.filtered.NovelGene = se.filtered[rowData(se.filtered)$novelGene == TRUE,]
head(mcols(rowRanges(se.filtered.NovelGene))[order(mcols(rowRanges(se.filtered.NovelGene))$NDR),])


plotBambu(se, type = "annotation", gene_id = "BambuGene1")

#investigate high confidence novel isoforms in annotated genes
se.filtered.noNovelGene <- se.filtered[rowData(se.filtered)$novelGene == FALSE,]
head(mcols(rowRanges(se.filtered.noNovelGene))[order(mcols(rowRanges(se.filtered.noNovelGene))$NDR),])

plotBambu(se, type = "annotation", gene_id = "ENSG00000272779")


#find which ones are unique based on counts
expression.A549 <- apply(assays(se.filtered)$counts[,grep("A549",colnames(se.filtered))],1,mean)
expression.HepG2 <- apply(assays(se.filtered)$counts[,grep("HepG2",colnames(se.filtered))],1,mean)


se.filtered[expression.A549>=1 &(expression.HepG2==0)] # unique in A549
se.filtered[expression.A549==0 &(expression.HepG2>=1)] # unique in HepG2


annotations.UCSC = rowRanges(se.filtered)
seqlevelsStyle(annotations.UCSC) <- "UCSC" #this reformats the chromosome names
annotations.UCSC = keepStandardChromosomes(annotations.UCSC, pruning.mode="coarse") #removes chromosomes UCSC doesn't have
writeToGTF(annotations.UCSC, "novel_annotations.UCSC_10July.gtf")


http://sg-nex-data.s3.amazonaws.com/data/data_tutorial/annotations/novel_annotations.UCSC_10July.gtf


### Downstream analysis ==============================
seGene <- transcriptToGeneExpression(se)
colData(seGene)$condition <- as.factor(rep(c("A549","HepG2"),each = 3))


library(DESeq2)
dds <- DESeqDataSetFromMatrix(round(assays(seGene)$counts),
                                    colData = colData(seGene),
                                    design = ~ condition)
dds.deseq <- DESeq(dds)
deGeneRes <- DESeq2::results(dds.deseq, independentFiltering = FALSE)
head(deGeneRes[order(deGeneRes$padj),])


# summary of differentially expressed genes
summary(deGeneRes)

library(apeglm)
resLFC <- lfcShrink(dds.deseq, coef = "condition_HepG2_vs_A549", type = "apeglm")
plotMA(resLFC, ylim = c(-3,3))


# summary of differentially used transcripts
library(DEXSeq)
colData(se)$condition <- as.factor(rep(c("A549","HepG2"),each = 3))
dxd <- DEXSeqDataSet(countData = round(assays(se)$counts),
sampleData = as.data.frame(colData(se)),
design = ~sample + exon + condition:exon,
featureID = rowData(se)$TXNAME,
groupID = rowData(se)$GENEID)
dxr <- DEXSeq(dxd)
head(dxr)

plotMA(dxr, cex = 0.8)

