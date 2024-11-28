# pychopper analysis on 2 samples: 1 latest cDNA sample and 1 direct cDNA sample 
# can work for both cDNA and direct cDNA
# cDNA sample
aws s3 cp --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fastq/SGNex_Hct116_cDNA_replicate3_run5/SGNex_Hct116_cDNA_replicate3_run5.fastq.gz  ./adpater_analysis/cDNA/ 
gunzip SGNex_Hct116_cDNA_replicate3_run5.fastq.gz
pychopper -k PCS110  -r report.pdf -A aln_hits.bed -S statistics.tsv -u unclassified.fq -w rescued.fq -t 8  SGNex_Hct116_cDNA_replicate3_run5.fastq full_length_output.fq
# direct cDNA sample:SGNex_Hct116_directcDNA_replicate5_run1.fastq.gz
aws s3 cp --no-sign-request  s3://sg-nex-data/data/sequencing_data_ont/fastq/SGNex_Hct116_directcDNA_replicate5_run1/SGNex_Hct116_directcDNA_replicate5_run1.fastq.gz   ./adpater_analysis/dcDNA/  
cd ./adpater_analysis/dcDNA
gunzip SGNex_Hct116_directcDNA_replicate5_run1.fastq.gz
{ /usr/bin/time -v  pychopper -r report.pdf -b dcs109_primers.fas -m edlib -A aln_hits.bed -S statistics.tsv -u unclassified.fq -w rescued.fq  -t 16 SGNex_Hct116_directcDNA_replicate5_run1.fastq full_length_output.fq  ; } 2> log.txt

## for aligned reads 
# for cdna sample only 
aws s3 cp --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/genome/SGNex_Hct116_cDNA_replicate3_run5/SGNex_Hct116_cDNA_replicate3_run5.bam  ./adpater_analysis/cDNA/ 
# samtools view -@ 48 -b -F 4  > mapped.bam
rm SGNex_Hct116_cDNA_replicate3_run5.bam

samtools fastq -@ 48 -F 0x904  SGNex_Hct116_cDNA_replicate3_run5.bam > mapped_primary_alignment.fastq
{ /usr/bin/time -v pychopper -k PCS110  -r report.pdf -A aln_hits.bed -S statistics.tsv -u unclassified.fq -w rescued.fq -t 48  mapped_primary_alignment.fastq full_length_output.fq   ; } 2> log.txt

#0x904

## processing output
library(data.table)

alnHits <- fread("adapter_analysis/pychopper_results/cdna/aln_hits.bed", header = FALSE)
read_vec <- unique(alnHits$V1)
alnHits[, read_id := match(V1, read_vec)]
alnHits[, V1 := NULL]
alnHit[, V6_mod := ifelse(V6=="+",1,ifelse(V6=="-",-1,NA))]
alnHits[, `:=`(count = .N), by = list(read_id, V4)]
alnHits_summary <- unique(alnHits[,list(total_count = .N,
                                        count_diff = max(diff(count)),
                                        strand_prod = prod(V6_mod)), by = list(read_id)])
cdna_results <- copy(alnHits_summary)

alnHits <- fread("adapter_analysis/pychopper_results/dcdna/aln_hits.bed", header = FALSE)
read_vec <- unique(alnHits$V1)
alnHits[, read_id := match(V1, read_vec)]
alnHits[, V1 := NULL]
alnHits[, V6_mod := ifelse(V6=="+",1,ifelse(V6=="-",-1,NA))]
alnHits <- unique(alnHits)

alnHits[V4=="VNP"&(V6_mod == 1), status := 3] 
alnHits[V4=="VNP"&(V6_mod == -1), status := -3] 
alnHits[V4=="SSP"&(V6_mod == -1), status := -5] 
alnHits[V4=="SSP"&(V6_mod == 1), status := 5] 
alnHits[,primer_combination := paste(paste0(V6,V4), collapse = ","), by = read_id]
alnHits_summary <- unique(alnHits[, list(fwd = all(c(3,-5) %in% unique(status)),
                                  rvs = all(c(-3,5) %in% unique(status)),
                                  vnp = 3 %in% unique(status)), by = read_id])

length(unique(cdna_raw[primer_combination %in% c("+VNP,-SSP","-SSP,+VNP","+SSP,-VNP","-VNP,+SSP")]$read_id))
length(unique(alnHits[primer_combination %in% c("+VNP","+SSP,-VNP","-VNP,+SSP")]$read_id))
# 1474224/6598848 reverse
# 1161132/6598848  forward


## installation =========================
# use mamba
mamba install -c nanoporetech -c anaconda -c bioconda "nanoporetech::pychopper" 
# this would cause error: what() could not unlink 
# fixed by using just one channel
mamba install -c bioconda "nanoporetech::pychopper"  ## successful 





