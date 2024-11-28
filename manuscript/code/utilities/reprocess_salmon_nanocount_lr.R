#!/usr/bin/env Rscript

rm(list = ls())
# set working directory
library(readxl)
library(data.table)

require(docopt)
'Usage:
reprocess_salmon_nanocount_lr.R [-g <g> ]

Options:
-g sample fastq file
]' -> doc

opts <- docopt(doc)
print(opts)
fastq_file <- as.character(opts$g)

save.dir <- paste0("salmon_filter/")
if(!dir.exists(save.dir)){
    dir.create(save.dir)
    dir.create(paste0(save.dir,"bam"))
    dir.create(paste0(save.dir,"count"))
    dir.create(paste0(save.dir,"index"))
}

save.dir_nanocount <- paste0("nanocount_fastq6.4.2/")
if(!dir.exists(save.dir_nanocount)){
    dir.create(save.dir_nanocount)
    dir.create(paste0(save.dir_nanocount,"bam"))
    dir.create(paste0(save.dir_nanocount,"count"))
}


system(paste0("aws s3 cp --no-sign-request ", fastq_file, " ", save.dir))

fastqFileFinal <- paste0(save.dir, fastq_file)

minimap2Path <- "minimap2"  ## 
salmonPath <- "salmon"
chopperPath <- "chopper"
nanoqPath <- "nanoq"
tx_ref <- "hg38_sequins_SIRV_ERCCs_longSIRVs_cdna.fa"
nanocountPath <- "NanoCount" # NanoCount v1.0.0.post3



nthreads <- 24
indexFile <- paste0(save.dir,"index/index.mmi")
system(paste0(minimap2Path," -t ",nthreads," -I 1000G -d ", indexFile, " ", tx_ref))#tx_ref_spikein

r <- gsub(".fastq.gz","",basename(fastq_file))

bam.file <- paste0(save.dir,"bam/",r,".bam")
output.file <- paste0(save.dir,"count/",r)

fastqFileFinal_filtered <- gsub(".fastq.gz","_filtered.fastq.gz", fastqFileFinal)
filter_report <- paste0(save.dir,"bam/",r,"_nanoq_filter_report.txt")
system(paste0("cat ",fastqFileFinal," | ",nanoqPath," -q 7 -O g -r ",filter_report," > ",fastqFileFinal_filtered))
system(paste0("rm -v ", fastqFileFinal))
system(paste0(minimap2Path," -t ",nthreads," -ax map-ont -p 1.0 -N 100 ",
             indexFile," ", fastqFileFinal_filtered," | samtools view -@ ",nthreads," -Sb > ",bam.file))
system(paste0(salmonPath," quant --ont -p ",nthreads," -t ",tx_ref," -q -l U -a ",bam.file," -o ", output.file))
system(paste0("rm -rvf ",bam.file,"*"))



bam.file <- paste0(save.dir_nanocount,"bam/",r,".bam")
output.file <- paste0(save.dir_nanocount,"count/",r,".tsv")
system(paste0(minimap2Path, " -t ",nthreads," -ax map-ont -p 0 -N 10 ", tx_ref, " ", fastqFileFinal,
              " | samtools view -@ ",nthreads," -bh > ",bam.file))
system(paste0(nanocountPath, " -i ",bam.file," -o ", output.file))
system(paste0("rm -rvf ",bam.file,"*"))


system(paste0("rm -v ", fastqFileFinal_filtered))





