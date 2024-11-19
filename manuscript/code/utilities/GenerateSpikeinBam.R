#!/usr/bin/env Rscript
##################################
# generate spikein bam file      #
##################################

###########################
## set-up                ##
###########################
rm(list = ls())
wkdir <- "."
setwd(wkdir)
# set working directory
library(readxl)
library(data.table)

#samtools idxstats NA18152.bam | cut -f 
samtools_path <- "samtools"
n_threads <- 48

###########################
# functions               #
###########################
get_spikein_bam <- function(sampleNames, bam.file,sampleData,save.dir,n_threads, samtools_path){
    np <- lapply(sampleNames, function(r){
        r_path <- bam.file[[r]]
       
        system(paste0("aws s3 cp --no-sign-request ",  r_path, " ", save.dir))
         
        in.bam <- paste0(save.dir, basename(r_path))
        out.bam <- gsub(".bam$","_chrIS_only.bam", in.bam)
        system(paste0(samtools_path, " index -@  ",n_threads," ", in.bam))
        system(paste0(samtools_path, " view -@  ",n_threads," -b -h -o  ", out.bam," ",in.bam," chrIS "))
        system(paste0(samtools_path," index -@  ",n_threads," ",out.bam))
        system(paste0(samtools_path," idxstats -@  ",n_threads," ", out.bam, " > ",gsub(".bam$",".idxstats",out.bam)))
        out.bam <- gsub(".bam$","_SIRVomeERCCome_only.bam", in.bam)
        system(paste0(samtools_path, " index -@  ",n_threads," ", in.bam))
        system(paste0(samtools_path, " view -@  ",n_threads," -b -h -o  ", out.bam," ",in.bam," SIRVomeERCCome "))
        system(paste0(samtools_path," index -@  ",n_threads,"  ",out.bam))
        system(paste0(samtools_path," idxstats -@  ",n_threads," ", out.bam, " > ",gsub(".bam$",".idxstats",out.bam)))
        system(paste0("rm -vf ", in.bam,"*"))
    })
}

get_spikein_bam_tx <- function(sampleNames, bam.file,sampleData,save.dir,n_threads, samtools_path){
    np <- lapply(sampleNames, function(r){
        r_path <- bam.file[[r]]
        system(paste0("aws s3 cp --no-sign-request ",  r_path, " ", save.dir))
        
        in.bam <- paste0(save.dir, basename(r_path))
        out.bam <- gsub(".bam$","_chrIS_only.bam", in.bam)
        system(paste0(samtools_path, " index -@  ",n_threads," ", in.bam))
        system(paste0(samtools_path, " view -@  ",n_threads," -N spikein_bam_Apr17/read.nameschrIS.txt -b -h -o  ", out.bam," ",in.bam))
        system(paste0(samtools_path," index -@  ",n_threads," ",out.bam))
        out.bam <- gsub(".bam$","_SIRVomeERCCome_only.bam", in.bam)
        system(paste0(samtools_path, " index -@  ",n_threads," ", in.bam))
        system(paste0(samtools_path, " view -@  ",n_threads," -N spikein_bam_Apr17/read.namesSIRVomeERCCome.txt  -b -h -o  ", out.bam," ",in.bam))
        system(paste0(samtools_path," index -@  ",n_threads,"  ",out.bam))
        system(paste0(samtools_path," idxstats -@  ",n_threads," ", out.bam, " > ",gsub(".bam$",".idxstats",out.bam)))
        system(paste0("rm -vf ", in.bam,"*"))
    })
}
merge_bam_file <- function(bam.path, save.dir_all, protocol_spikein){
    np <- lapply(seq_len(nrow(protocol_spikein)), function(k){
        si_type <- protocol_spikein[k]$spikein_type
        chromosome_name <- ifelse(grepl("sequin",si_type),"chrIS","SIRVomeERCCome")
        protocol_type <- protocol_spikein[k]$protocol
        p_bam_path <- bam.path[grepl(protocol_type, bam.path)&grepl(chromosome_name, 
                                                                    bam.path)&(gsub("_chrIS_only|_SIRVomeERCCome_only|_sorted|_R1_sorted|\\.bam|\\.tx","",
                                                                                    basename(bam.path)) %in% c(samples_wSpikein[grepl(si_type,RNAcontent)]$runname,
                                                                                                               samples_wSpikein[grepl(si_type,RNAcontent)]$old_runname))]
        if(length(p_bam_path)==0){
            return(NULL)
        }
        out.bam <- paste0(save.dir_all,"allSpikinReadsCombined_",gsub(",","",gsub(" ","",si_type)),"_",protocol_type,".bam")
        sout.bam <- gsub(".bam$","_sorted.bam",out.bam)
        system(paste0(samtools_path," merge -@  ",n_threads," -o ",out.bam, " ", paste(p_bam_path, collapse = " ")))
        system(paste0(samtools_path," sort -@  ",n_threads,"  -o ",sout.bam, " ", out.bam))
        system(paste0(samtools_path," index -@  ",n_threads," ",sout.bam))
        system(paste0("rm -vf ",out.bam))
    })
}


###########################
## sample information    ##
###########################

sampleData <- data.table(as.data.frame(read_xlsx(paste0('ONT Master Table.xlsx'), sheet = 1))) ## need to convert from tibble to data.frame
sampleData <- sampleData[(grepl("H9|HEYA8",runName)&(grepl("ON00",name))&(!grepl("HEYA8.*H9",
                                                                                 runName)))|(SG_NextData_Release=="Yes"&(!is.na(SG_NextData_Release))&(!grepl("CRC",runName)))|(grepl("HCT116",runName))]
#sampleData$runName_combined <- gsub('-pre|-Pre','',sampleData$runName) # there are runs with multiple datasets that should be combined together 
sampleData[,runName_combined := ifelse(grepl("directRNA",runName)|(!grepl("H9|HEYA8",runName))|(grepl("WINSTON",name)),
                                       runName, 
                                       `GIS Library ID`)]
sampleData[, runName_combined := gsub("HCT","Hct", runName_combined)]

sampleData[, runName := gsub("HCT","Hct", runName)]
sampleData[runName_combined != runName, runName_combined := paste0(gsub("_Run.*","",runName),"_",runName_combined)]
sampleData[, runName_combined := gsub('-pre|-Pre','',runName_combined)]
sampleData$demultiplexed <- grepl("NB", sampleData$name)|(grepl("barcode",sampleData$name))|(!is.na(sampleData$barcoding.kit))


sampleNames <- unique(sampleData$`publicName (SGNex_CellLine_protocol_replicate1_run1)`)[1:111]

sampleData_sr <- data.table(as.data.frame(read_xlsx("ONT Master Table.xlsx", sheet = 3))) ## need to convert from tibble to data.frame
sampleData_sr <- sampleData_sr[!grepl('#',`ELM library ID`) &(!is.na(runName))&(!grepl("HEYA8.*H9",runName))]
sr_runNames <- sampleData_sr$runName
chrm_names <- c(1:22,'X','Y')


pacbio_data <-  data.table(as.data.frame(read_xlsx(paste0('ONT Master Table.xlsx'), sheet = 2)))



samples <- unique(as.data.table(sampleData)[,.(runName_combined,`publicName (SGNex_CellLine_protocol_replicate1_run1)`,
                                               `replicate-id`,Platform,demultiplexed)],by=NULL)
setnames(samples, 1:3, c("runname","publicName","replicate_id"))
samples <- rbindlist(list(samples,
                          data.table(runname = sr_runNames,
                                     publicName = sampleData_sr$public_name,
                                     replicate_id = sampleData_sr$`replicate-id`,
                                     Platform = rep("Illumina",21)),
                          data.table(runname = pacbio_data$name,
                                     publicName = pacbio_data$public_name,
                                     replicate_id = pacbio_data$Replicate,
                                     Platform = "PacBio")), fill= TRUE)


samples[, cellLine:=gsub("-EV","",gsub('k562','K562',strsplit(runname, '\\_')[[1]][2])),by = runname]
samples[, protocol:=strsplit(runname, '\\_')[[1]][3], by = runname]
samples[, cDNAstranded:=ifelse(protocol %in% c('cDNA','cDNAStranded'), protocol=='cDNAStranded',NA)]
samples[, randomPrimer:=grepl('RandomPrimer',protocol)]
samples[, protocol_type:=gsub('Stranded|RandomPrimer','',gsub('PromethionD','d', protocol))]

## add replicate ids for new cell lines
samples[cellLine %in% c("H9","HEYA8") &(is.na(replicate_id)), replicate_id := gsub("Rep", " rep ",gsub("GIS_|(_RHH.*)|(_Run.*)|(_cDNA_)|(_directcDNA_)|(_directRNA_)","",runname)), by = runname]

samples[protocol != "Illumina", bioRep:=strsplit(publicName, '\\_')[[1]][4], by = runname]
samples[, old_runname := runname]
samples[, runname := publicName]
samples[, publicName := NULL]




cellLines <- c('Hct116','HepG2','K562','A549','MCF7',"H9","HEYA8")
spike_in_info <- unique(sampleData[,.(runName_combined, RNAcontent)])
setnames(spike_in_info, "runName_combined","old_runname")
samples_wSpikein <- spike_in_info[samples, on = "old_runname"]
# RNA content information for short read missing
samples_wSpikein[grepl("Illumina", runname)&(cellLine %in% cellLines[1:5]), RNAcontent := "1% RNA sequin Mix A v1.0"]
samples_wSpikein[grepl("Illumina", runname)&(cellLine == "H9"), RNAcontent := "1% spike-in of 6ng SIRV-4"]
samples_wSpikein[grepl("Illumina", runname)&(cellLine == "HEYA8"), RNAcontent := "1% spike-in of 6ng RNA sequin Mix A, V2, 6ng SIRV-1 E2"]
samples_wSpikein[grepl("PacBio", runname), RNAcontent := "sequin Mix A v1.0 sequin MixA V2 SIRV-1 E2 SIRV-4"]
samples_wSpikein[, RNAcontent := gsub("SIRV-1 \\(E2","SIRV-1 E2",gsub("A\\,","A",RNAcontent))]



#sampleNames <- samples_wSpikein[grepl("sequin|SIRV",RNAcontent)&(!grepl("Illumina",runname))]$runname
sampleNames <- samples_wSpikein[grepl("SIRV-4",RNAcontent)&(grepl("_cDNA",runname))]$runname
#sampleNames <- samples_wSpikein[grepl("sequin|SIRV",RNAcontent)&(grepl("PacBio",runname))]$runname

##########################
# genome bam file        #
##########################
bam.file <- unlist(lapply(sampleNames, function(r){ # sampleNames for all
    if(grepl("PacBio",r)){
        bam.file <- gsub('s3://ontdata.store.genome.sg','/mnt/ontdata',pacbio_data[public_name == r]$`bam-genome.path`)
    }else if(!grepl("Illumina",r)){
        bam.file <- paste0("s3://sg-nex-data/data/sequencing_data_ont/bam/genome/",r,"/",r,".bam")#
        
    }else{
        rname <- unique(sampleData_sr[runName == r]$public_name)
        bam.file <- paste0("s3://sg-nex-data/data/sequencing_data_illumina/bam/genome/",rname,"/",rname,".bam")
    }
    
    return(bam.file)
}))
names(bam.file) <- sampleNames


save.dir <- "spikein_bam_Apr17/bysample/"
if(!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)
get_spikein_bam(sampleNames, bam.file,sampleData,save.dir,n_threads, samtools_path)


bam.path <- dir(save.dir, pattern = ".bam$", full.names = TRUE)
save.dir_all <- paste0(wkdir,"/spikein_bam_Apr17/")
protocolVec <- c("directcDNA","_cDNA","directRNA","Illumina","PacBio")[2]
spikein_type <- c("sequin Mix A v1.0","sequin Mix A V2","SIRV-1","SIRV-4")[4]
#spikein_type <- c("sequin Mix A v1.0","sequin Mix A V2","SIRV-1","SIRV-4")

protocol_spikein <- CJ(protocol = protocolVec, spikein_type = spikein_type)
merge_bam_file(bam.path, save.dir_all, protocol_spikein)

##########################
# trancriptome bam file        #
##########################
# can actually use the seqnames in gene bam files and then susbet transcriptome bam files 
chrnames <- c("chrIS","SIRVomeERCCome")
path.dir <- "spikein_bam_Apr17/"
setwd(path.dir)
for(i in seq_along(chrnames)){
    bam.file <- list.files(".", pattern = paste0(c("sequin","SIRV")[i],".*bam$"), full.names = TRUE)
    noprint <- lapply(seq_along(bam.file), function(r){
        system(paste0("samtools view -@  ",n_threads," ",bam.file[r],
                      " | cut -f 1 | awk '!x[$0]++' > x",r,"read.names_May5.txt"))
    })
    qname_file <- paste0("read.names",chrnames[i],"_May5.txt")
    qname_file_list <- list.files(".", pattern = "read.names_May5.txt",full.names = TRUE)
    system(paste0("cat ",paste(qname_file_list, collapse = " "), " > ",qname_file))
    system(paste0("rm -rvf *read.names_May5.txt"))
}


bam.file <- unlist(lapply(sampleNames, function(r){ # sampleNames for all
    if(grepl("PacBio",r)){
        bam.file <- gsub('s3://ontdata.store.genome.sg','/mnt/ontdata',pacbio_data[public_name == r]$`bam-tx.path`)
    }else if(!grepl("Illumina",r)){
        bam.file <- paste0("s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/",r,"/",r,".bam")#
    }else{
        rname <- unique(sampleData_sr[runName == r]$public_name)
        bam.file <- paste0("s3://sg-nex-data/data/sequencing_data_illumina/bam/transcriptome/",rname,"/",rname,".bam")
    }
        
    return(bam.file)
}))
names(bam.file) <- sampleNames
    
    
save.dir <- "spikein_bam_tx_Apr17/bysample/"
if(!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)
get_spikein_bam_tx(sampleNames, bam.file,sampleData,save.dir,n_threads, samtools_path) #[37:42]

bam.path <- dir(save.dir, pattern = ".bam$", full.names = TRUE)
save.dir_all <- paste0(wkdir,"/spikein_bam_tx_Apr17/")
protocolVec <- c("directcDNA","_cDNA","directRNA","Illumina","PacBio")[2]
spikein_type <- c("sequin Mix A v1.0","sequin Mix A V2","SIRV-1","SIRV-4")[4]

protocol_spikein <- CJ(protocol = protocolVec, spikein_type = spikein_type)

merge_bam_file(bam.path, save.dir_all, protocol_spikein)
    
    
    