#!/usr/bin/env nextflow
// Author: Jonathan Göke
// Contact: gokej@gis.a-star.edu.sg
// This script is based on https://github.com/GoekeLab/bioinformatics-workflows/tree/master/nextflow

// Define pipeline input parameters 
// note: input files require the use of absolute paths
params.refFa = '/path/to/ref.fa'
params.refGtf = '/path/to/ref.gtf'
params.reads = '/path/to/reads.fq'
params.outdir = 'results'

// Align reads to the genome using Minimap2
process MINIMAP2_ALIGN {

  input: 
    path refFa
    path reads
  output:
    path "aligned_reads.sam"

  """
    minimap2 -ax splice -uf -k14 $refFa $reads > aligned_reads.sam
  """
}
// Option to reduce memory usage in Minimap2 -I 2G 
// Option to use multiple threads: -t 8

// Convert sam to bam file using Samtools
process SAM_TO_BAM {
  publishDir params.outdir, mode: 'copy'

  input:
    path reads_sam
  output:
    path "aligned_reads.bam"
  
  """
    samtools view -b $reads_sam > aligned_reads.bam
  """
}
// The bam index/sorting is not required for this workshop
// For visualisation and specific access operations it will be required
// samtools sort reads_unsorted.bam -o aligned_reads.bam
// samtools index aligned_reads.bam


// Transcript discovery and quantification with Bambu
process BAMBU {
  publishDir params.outdir, mode: 'copy'

  input:
    path refFa
    path refGtf
    path reads_bam
  output:
    path "counts_transcript.txt"
    path "counts_gene.txt"
    path "extended_annotations.gtf"

    """
    #!/usr/bin/env Rscript --vanilla
    library(bambu)
    annotations <- prepareAnnotations("$refGtf")
    se     <- bambu(reads = "$reads_bam", 
                    annotations = annotations,
                    genome = "$refFa",
                    NDR=1,
                    ncore = 1)
    writeBambuOutput(se, path = "./")
    """
}

workflow {
  MINIMAP2_ALIGN(params.refFa, params.reads) 
  SAM_TO_BAM(MINIMAP2_ALIGN.out)
  BAMBU(params.refFa, params.refGtf, SAM_TO_BAM.out)
}

/*
 nextflow run nextflow/workflow_longReadRNASeq.nf -with-report -resume \
      --reads $PWD/fastq/A549_directRNA_sample2.fastq.gz \
      --refFa $PWD/reference/hg38_chr22.fa \
      --refGtf $PWD/reference/hg38_chr22.gtf \
      --outdir $PWD/results/

*/
