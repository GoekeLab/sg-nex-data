##################################
# get spikein fastq from bam and perform salmon/rsem/nanocount #
##################################




rm(list = ls())
bam.file <- dir("spikein_bam_Apr17/", pattern = ".bam$", full.names = TRUE)

noprint <- lapply(bam.file[12], function(k){
    fastq.dir <- paste0("spikein_fastq_Apr17/",basename(k),"/")
    if(!dir.exists(fastq.dir)) dir.create(fastq.dir, recursive = TRUE)
    fastq.file <- paste0(fastq.dir,gsub(".bam",".fastq.gz",basename(k)))
    
    if(!grepl("Illumina", k)){# for long read 
        system(paste0("samtools fastq -@ 24 ",k, " > ",fastq.file))
    }else{# for short read 
        fastq.file1 <- gsub(".fastq.gz","_R1.fastq.gz",fastq.file)
        fastq.file2 <- gsub(".fastq.gz","_R2.fastq.gz",fastq.file)
        system(paste0("samtools fastq -@ 24 -1 ",fastq.file1, " -2 ",fastq.file2," ",k))
    }
    
    
})


# after converting run salmon
# short read cannot directly convert as they are paired reads 

# short read salmon===============
fastq.dir <- "spikein_fastq_Apr17/"
short_read_samples <- dir(fastq.dir, pattern = "Illumina")
salmonPath <- "salmon"

noprint <- lapply(short_read_samples[-1], function(k){
    fastq.dir_k <- paste0(fastq.dir,k)
    if(!dir.exists(fastq.dir_k)){
        dir.create(fastq.dir_k, recursive = TRUE)
    }
    mapDir <- paste0(fastq.dir,"map/",k)
    if(!dir.exists(mapDir)){
        dir.create(mapDir, recursive = TRUE)
    }
    
    fileList1 <- list.files(fastq.dir_k, pattern = glob2rx("*R1*.fastq.gz$*"), full.names = TRUE)
    fileList2 <- list.files(fastq.dir_k, pattern = glob2rx("*R2*.fastq.gz$*"), full.names = TRUE)
    
    system(paste0(salmonPath, ' quant -i ',annotationDir,
                  ' -p ', nthreads,' ',
                  ' -l A -1 ',
                  paste(fileList1, collapse = ' '), ' -2 ',paste(fileList2, collapse = ' '),
                  ' --validateMappings ',
                  ' --seqBias ', ' --gcBias ', ' --posBias ',
                  ' -o ', mapDir,'/transcripts_quant'))
    
})


# run rsem===================
rsemPath <- "RSEM-1.3.3/"
starPath <- "STAR-2.7.10a/bin/Linux_x86_64/"
anno.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs_v5_reformatted.gtf"
annotationDir <- 'rsem-ref'
gene_ref <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"

if(!dir.exists(annotationDir)){
    dir.create(annotationDir, recursive = TRUE)
}
system(paste0(rsemPath,"rsem-prepare-reference --gtf ",anno.file, 
              " --star ",
              " --star-path ", starPath,
              " -p 8 ",gene_ref, " ",annotationDir,"/rsem-ref"))


nthreads <- 24
noprint <- lapply(short_read_samples[1], function(k){
    fastq.dir_k <- paste0(fastq.dir,k)
    if(!dir.exists(fastq.dir_k)){
        dir.create(fastq.dir_k, recursive = TRUE)
    }
    mapDir <- paste0(fastq.dir,"map_rsem/",k)
    if(!dir.exists(mapDir)){
        dir.create(mapDir, recursive = TRUE)
    }
    
    fileList1 <- list.files(fastq.dir_k, pattern = glob2rx("*R1*.fastq.gz$*"), full.names = TRUE)
    fileList2 <- list.files(fastq.dir_k, pattern = glob2rx("*R2*.fastq.gz$*"), full.names = TRUE)
    
    tempDir <- paste0(fastq.dir,'map_rsem/temp/',k)
    if(!dir.exists(tempDir)){
        dir.create(tempDir, recursive = TRUE)
    }
    setwd(mapDir)
    system(paste0(rsemPath, 'rsem-calculate-expression --paired-end  ',
                  ' --star ', ' --star-path ', starPath,
                  ' ',paste(fileList1, collapse = ','), 
                  ' ',paste(fileList2, collapse = ','),
                  ' ',annotationDir,"/rsem-ref",
                  ' -p ', nthreads,
                  ' --star-gzipped-read-file ',
                  ' --no-bam-output ',
                  ' --calc-pme ',
                  ' --calc-ci ',
                  ' --temporary-folder ',tempDir, #temp dir will be removed automatically 
                  ' --time ',
                  k))
    
})


# run salmon for long read======================
# same bam2fastq first, then runsalmon
fastq.dir <- "spikein_fastq_Apr17/"
short_read_samples <- dir(fastq.dir, pattern = "Illumina")
long_read_samples <- setdiff(dir(fastq.dir), c(short_read_samples,"map","map_rsem"))
salmonPath <- "salmon"
nthreads <- 24

minimap2Path <- "minimap2"  ## 
tx_ref <- "hg38_sequins_SIRV_ERCCs_longSIRVs_cdna.fa"

indexFile <- paste0(fastq.dir,"index.mmi")
system(paste0(minimap2Path," -t ",nthreads," -I 1000G -d ", indexFile, " ", tx_ref))#tx_ref_spikein



noprint <- lapply(long_read_samples[12], function(k){
    fastq.dir_k <- paste0(fastq.dir,k)
    if(!dir.exists(fastq.dir_k)){
        dir.create(fastq.dir_k, recursive = TRUE)
    }
    mapDir <- paste0(fastq.dir,"map_salmon_lr/",k)
    if(!dir.exists(mapDir)){
        dir.create(mapDir, recursive = TRUE)
    }
    fastqFileFinal <- dir(fastq.dir_k, full.names = TRUE)
    print(fastqFileFinal)
    bam.file <- paste0(mapDir,"/",k,".bam")
    #sbam.file <- paste0(save.dir,"bam/",r,".sorted.bam")
    output.file <- paste0(mapDir,"/","count/",k)
    
    system(paste0(minimap2Path," -t ",nthreads," -ax map-ont -p 1.0 -N 100 ",
                  indexFile," ", fastqFileFinal," | samtools view -@ ",nthreads,
                  " -Sb > ",bam.file))
    system(paste0(salmonPath," quant --ont -p ",nthreads," -t ",
                  tx_ref," -q -l U -a ",bam.file," -o ", output.file))
    system(paste0("rm -v ",bam.file))
    #system(paste0("rm ",sbam.file,".bai"))
})



## nanocount =======================
fastq.dir <- "spikein_fastq_Apr17/"
short_read_samples <- dir(fastq.dir, pattern = "Illumina")
long_read_samples <- setdiff(dir(fastq.dir), c(short_read_samples,"map","map_rsem","index.mmi"))
salmonPath <- "salmon"
nthreads <- 24

minimap2Path <- "minimap2"  ## 
tx_ref <- "hg38_sequins_SIRV_ERCCs_longSIRVs_cdna.fa"
nanocountPath <- "NanoCount" # NanoCount v1.0.0.post3

noprint <- lapply(long_read_samples[12], function(k){
    fastq.dir_k <- paste0(fastq.dir,k)
    if(!dir.exists(fastq.dir_k)){
        dir.create(fastq.dir_k, recursive = TRUE)
    }
    mapDir <- paste0(fastq.dir,"map_nanocount_lr/",k)
    if(!dir.exists(mapDir)){
        dir.create(mapDir, recursive = TRUE)
        dir.create(paste0(mapDir,"/count"))
    }
    fastqFileFinal <- dir(fastq.dir_k, full.names = TRUE)
    print(fastqFileFinal)
    bam.file <- paste0(mapDir,"/",k,".bam")
    output.file <- paste0(mapDir,"/","count/",k,".tsv")
    system(paste0(minimap2Path, " -t ",nthreads," -ax map-ont -p 0 -N 10 ", tx_ref, " ", fastqFileFinal,
                  " | samtools view -@ ",nthreads," -bh > ",bam.file))
    system(paste0(nanocountPath, " -i ",bam.file," -o ", output.file))
    system(paste0("rm -rvf ",bam.file))
})
