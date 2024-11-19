#!/bin/bash

cat sqanti3_annotation_test.txt | while read line
do  
     
      echo "$line"
      dirname=$(echo $line | cut -d"_" -f3 | sed -e 's/.gtf//g')
      echo "$dirname"
     ./sqanti3_qc.py  $line   hg38_sequins_SIRV_ERCCs_longSIRVs_v5_reformatted.gtf hg38_sequins_SIRV_ERCCs_longSIRVs.fa    \
--CAGE_peak ./SQANTI3-5.2/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed    \
--polyA_motif_list ./SQANTI3-5.2/data/polyA_motifs/mouse_and_human.polyA_motif.txt    \
-o sqanti3_filter -d  $dirname \
 -t 2 --report both --force_id_ignore  --isoAnnotLite
      echo "$line"
done 
