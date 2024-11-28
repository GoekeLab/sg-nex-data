rm(list = ls())
# set working directory
wkdir <- '.'
setwd(wkdir)
library(readxl)
library(data.table)

sampleData <- data.table(as.data.frame(read_xlsx(paste0('.'), sheet = 1))) ## need to convert from tibble to data.frame
sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
                                                                                 runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(SG_NextData_Release=="Yes"),
                                       runName,
                                       `GIS Library ID`)]
sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]
sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:112]

bam.file <- sapply(sampleNames, function(x) paste0("s3://sg-nex-data/data/sequencing_data_ont/bam/genome/"x,"/",x,".bam"))

names(bam.file) <- sampleNames
# short read pipeline
anno.file <- "hg38_sequins_SIRV_ERCCs_longSIRVs_corrected.gtf"
# stringtie2-step1, discovery by each sample
stringtie2_path <- "stringtie"
ncpu <- 12
bamDir <- "lr_stringtie2/bam"
dir.create(bamDir, recursive = TRUE)
saveDir <- "lr_stringtie2/stringtie2_results"
dir.create(saveDir, recursive = TRUE)
np <- lapply(sampleNames, function(rname){
    bamDir_temp <- paste0(bamDir, "/", rname)
    dir.create(bamDir_temp)
    s3_bam_dir <- bam.file[rname]
    system(paste0('aws s3 sync --no-sign-request ',s3_bam_dir," ", bamDir_temp))
    bam.file_temp <- dir(bamDir_temp, full.names = TRUE, pattern = ".bam$")
    system(paste0(stringtie2_path, " -p ",ncpu," -G ", 
                  anno.file, " -o ", saveDir,"/", #anno.file_halfspikein
                  rname,".out.gtf ", bam.file_temp))
    system(paste0("rm -rvf ",bamDir_temp))
})

# stringtie2-step2, merge annotation across samples 
gtf_list_file <- dir(saveDir, full.names = TRUE, pattern  = ".gtf")
merged.gtf <- "lr_stringtie2/stringtie2_merged.gtf"
system(paste0(stringtie2_path," --merge -o ",
              merged.gtf, " -G ", anno.file, "  ", gtf_list_file))
