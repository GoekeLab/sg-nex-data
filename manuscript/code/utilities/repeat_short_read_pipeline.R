rm(list = ls())
# set working directory
wkdir <- '.'
setwd(wkdir)
library(readxl)
library(data.table)

sampleDataSR <- as.data.frame(read_xlsx(paste0('ONT Master Table.xlsx'), sheet = 3))## need to convert from tibble to data.frame
sampleDataSR <- sampleDataSR[!grepl('#',sampleDataSR[,1]) &(!is.na(sampleDataSR$public_name)),]

sampleNamesSR <- sampleDataSR$public_name


# short read pipeline
anno.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs_v5_reformatted.gtf"
# stringtie2-step1, discovery by each sample
stringtie2_path <- "stringtie"
ncpu <- 12
bamDir <- "sr_stringtie2/bam"
dir.create(bamDir, recursive = TRUE)
saveDir <- "stringtie2_results"
dir.create(saveDir, recursive = TRUE)
np <- lapply(sampleNamesSR, function(rname){
    bamDir_temp <- paste0(bamDir, "/", rname)
    dir.create(bamDir_temp)
    s3_bam_dir <- paste0('s3://sg-nex-data/data/sequencing_data_illumina/bam/genome/',rname)
    system(paste0('aws s3 sync --no-sign-request ',s3_bam_dir," ", bamDir_temp))
    bam.file <- dir(bamDir_temp, full.names = TRUE, pattern = ".bam$")
    system(paste0(stringtie2_path, " -p ",ncpu," -G ", # for short read 
                  anno.file, " -o ", saveDir,"/", #anno.file_halfspikein
                  gsub(".bam","",basename(bam.file)),".out.gtf ", bam.file))
    system(paste0("rm -rvf ",bamDir_temp))
})

# stringtie2-step2, merge annotation across samples 
gtf_list_file <- dir(saveDir, full.names = TRUE, pattern  = ".gtf")
merged.gtf <- "sr_stringtie2/stringtie2_merged.gtf"
system(paste0(stringtie2_path," --merge -o ",
              merged.gtf, " -G ", anno.file, "  ", gtf_list_file))


tx_ref <- "sr_stringtie2/stringtie2_merge.fasta"
gene_ref <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
# to obtain tx_ref 
library(Biostrings)
cat('Load transcript sequence information')
geneSeq <- readDNAStringSet(file=gene_ref)
listNames <- unlist(lapply(strsplit(names(geneSeq)," "),'[[',1))
stringtie2Annotations <- bambu::readFromGTF(merged.gtf)
txNames <- names(stringtie2Annotations)
library(BiocParallel)
bpParameters <- bpparam()
#===# set parallel options: otherwise use parallel to distribute samples
bpParameters$workers <- 4
txSeq <- unlist(bplapply(seq_along(txNames), function(s){
    print(paste(s, txNames[s]))
    exonrange <- stringtie2Annotations[[txNames[s]]]
    
    seqCharList <- lapply(seq_along(exonrange), function(g){
        
        tmp_range <- exonrange[g]
        seq_pos <- start(tmp_range):end(tmp_range)
        seqChar <- geneSeq[[match(as.character(unique(seqnames(tmp_range))), listNames)]][seq_pos]
        if(as.character(unique(strand(tmp_range))) == "-"){
            seqChar <- reverseComplement(seqChar)
        }
        as.character(seqChar)
    })
    
    if(as.character(unique(strand(exonrange))) == "-"){
        seqCharList <- rev(seqCharList)
    }
    seqChar <- paste(unlist(seqCharList), collapse = "")
}, BPPARAM = bpParameters))


sink(tx_ref)
noprint <- lapply(seq_along(txNames), function(s){
    cat(paste0(">",txNames[s]), " \n")
    cat(txSeq[s]," \n")
})
sink()

## create salmon annotations
genetx_ref <- "stringtie2_merge_genetx.fa"
gene_ref <- "hg38_sequins_SIRV_ERCCs_longSIRVs.fa"
decoyFile <- "decoys.txt"
salmonPath <- "salmon"
annotationDir <- 'transcriptome-index/salmon_index_stringtie2_merge'
system(paste0('grep "^>" < ',gene_ref,' | cut -d " " -f 1 > ', decoyFile))
system(paste0('sed -i.bak -e "s/>//g" ',decoyFile))
system(paste0("cat ",tx_ref," ", gene_ref, " > ",genetx_ref))
system(paste0(salmonPath," index -t ",genetx_ref," -d ",decoyFile," -p 12 -i ",annotationDir))

# salmon
fastqDir <- "sr_stringtie2/fastq"
dir.create(fastqDir, recursive = TRUE)
nthreads <- 12

np <- lapply(sampleNamesSR[-1], function(rname){
    fastqDir_temp <- paste0(fastqDir, "/", rname)
    dir.create(fastqDir_temp)
    s3_fastq_dir <- paste0('s3://sg-nex-data/data/sequencing_data_illumina/fastq/',rname)
    system(paste0('aws s3 sync --no-sign-request ',s3_fastq_dir," ", fastqDir_temp))
    fileList1 <- paste0(fastqDir_temp,'/',rname,"_R1.fastq.gz")
    fileList2 <- paste0(fastqDir_temp,'/',rname,"_R2.fastq.gz")
    mapDir <- paste0("sr_stringtie2/salmon/",rname)
    dir.create(mapDir, recursive = TRUE)
    system( paste0(salmonPath, ' quant -i ',annotationDir,
                   ' -p ', nthreads, ' ',
                   ' -l A -1 ',
                   paste(fileList1, collapse = ' '), ' -2 ',paste(fileList2, collapse = ' '),
                   ' --validateMappings ',
                   ' --seqBias ', ' --gcBias ', ' --posBias ',
                   ' -o ', mapDir,'/transcripts_quant'))
    system(paste0("rm -rvf ",fastqDir_temp))
})



