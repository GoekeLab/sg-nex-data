# **Isoform quantification from direct RNA data with the SG-NEx samples** 

In this tutorial, we will perform isoform quantification on a direct RNA sequencing SG-NEx sample. We will use the Oxford Nanopore direct RNA sequencing data from the A549 cell line from the SG-NEx project for this analysis. The A549 cell line was extracted from the lung tissues of a patient with lung cancer. The package **NanoCount** will be used to obtain isoform counts in the cell line.

## **Content**

- [Installation](#installation)
- [Data Access](#data-access)
- [Running NanoCount](#running-nanocount)
- [Reference](#reference)

## **Installation**

Install NanoCount with pip from GitHub.

```
pip install git+https://github.com/a-slide/NanoCount.git
```

NanoCount repository: https://github.com/a-slide/NanoCount

Detailed NanoCount documentation: https://a-slide.github.io/NanoCount/


## **Data Access**
Next, we will need to download the required data to run NanoCount. The required data for this tutorial:

-   aligned reads to the transcriptome from the A549 cell line direct RNA run (bam file and bam index file)

### **Download Data**
The data files can be downloaded from the [SG-NEx AWS S3 bucket](http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/). Please refer to this [page](AWS_data_access_tutorial.md) for a comprehensive tutorial to access the SG-NEx dataset.

To begin we will create a directory to store the required files.
```bash
# create a directory to store the files for this sample
mkdir A549_directRNA_replicate1_run1 
```

Then we will access the S3 bucket to download the data into the created directory. There are two methods to download the data. 
<br>

#### **Method One (Recommended):** 
This command requires you to have [AWS CLI](https://aws.amazon.com/cli/) installed.
```bash
# To download the processed data
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/SGNex_A549_directRNA_replicate1_run1/SGNex_A549_directRNA_replicate1_run1.bam ./A549_directRNA_replicate1_run1

aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/SGNex_A549_directRNA_replicate1_run1/SGNex_A549_directRNA_replicate1_run1.bam.bai ./A549_directRNA_replicate1_run1

# The directory should now contain two files, the BAM file and BAM index file. 
```

#### **Method Two:**
```bash
cd A549_directRNA_replicate1_run1 

# To download the BAM file
wget https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/sequencing_data_ont/bam/transcriptome/SGNex_A549_directRNA_replicate1_run1/SGNex_A549_directRNA_replicate1_run1.bam

# To download the BAM index file
wget https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/sequencing_data_ont/bam/transcriptome/SGNex_A549_directRNA_replicate1_run1/SGNex_A549_directRNA_replicate1_run1.bam.bai

# The directory should now contain two files, the BAM file and BAM index file. 
```

## **Running NanoCount**

For this tutorial we will start with the BAM file of aligned reads, however if you are starting from the FASTQ file the reads need to be mapped to the transcriptome with minimap2 with the following command:

```
minimap2 -t 4 -ax map-ont -N 10 transcriptome.fasta reads.fastq | samtools view -bh > aligned_reads.bam
```

We recommend using the -N 10 option to retain at least 10 secondary mappings.

### Run NanoCount

In this tutorial we will run NanoCount on the A549_directRNA_replicate1_run1.bam file. 

```
# within the A549_directRNA_replicate1_run1 directory
NanoCount -i SGNex_A549_directRNA_replicate1_run1.bam -b filtered_SGNex_A549_directRNA_replicate1_run1.bam --extra_tx_info -o SGNex_A549_directRNA_replicate1_run1_counts.tsv
```

#### Expected output while running NanoCount

```
## Checking options and input files ##
## Initialise Nanocount ##
	Parse Bam file and filter low quality alignments
[W::hts_idx_load3] The index file is older than the data file: SGNex_A549_directRNA_replicate1_run1.bam.bai
	Summary of alignments parsed in input bam file
		Valid alignments: 978,468
		Discarded alignment with invalid 3 prime end: 447,613
		Discarded unmapped alignments: 115,299
		Discarded negative strand alignments: 7,638
		Discarded supplementary alignments: 3,899
		Discarded short alignments: 10
	Summary of reads filtered
		Valid secondary alignments: 360,755
		Reads with valid best alignment: 335,179
		Invalid secondary alignments: 275,778
		Reads with low query fraction aligned: 3,613
	Write selected alignments to BAM file
	Summary of alignments written to bam
		Alignments skipped: 856,993
		Alignments to select: 695,934
		Alignments written: 695,934
	Generate initial read/transcript compatibility index
## Start EM abundance estimate ##
	Progress: 14.0 rounds [00:09, 1.51 rounds/s]
	Exit EM loop after 14 rounds
	Convergence value: 0.004933839788023745
## Summarize data ##
	Convert results to dataframe
	Compute estimated counts and TPM
	Write file
```

#### Output TSV file

NanoCount returns a file containing count data per transcript isoform. 

- **transcript_name**: Transcript name as in source bam file.
- **raw**: Raw abundance estimates. The sum of all abundance values is 1.
- **est_count**: Estimated counts obtained by multiplying the raw abundance by the number of primary alignments.
- **tpm**: Estimated counts obtained by multiplying the raw abundance by 1M.
- **transcript_length**: Optional column included with the option extra_tx_info.

```
# the first 10 lines of the TSV file produced from the sample SGNex_A549_directRNA_replicate1_run1

head SGNex_A549_directRNA_replicate1_run1_counts.tsv

transcript_name	raw	est_count	tpm	transcript_length
ENST00000331825.10	0.019231515100894832	6446.0000000028285	19231.51510089483	878
ENST00000327741.9	0.009842502066064544	3299.000000001448	9842.502066064544	1929
ENST00000359579.4	0.00772651054799532	2589.7640789665234	7726.51054799532	1590
ENST00000297785.7	0.0075303047028635675	2524.0000000011078	7530.304702863567	2107
ENST00000534180.1	0.006889745554571151	2309.2980252356037	6889.745554571151	983
ENST00000251453.7	0.006321448488095429	2118.8167827913376	6321.448488095429	560
ENST00000380554.4	0.005240445849168986	1756.4873992786115	5240.445849168986	1807
ENST00000451311.6	0.0052348178074999255	1754.6009979000175	5234.817807499926	628
ENST00000387347.2	0.00505998287482433	1696.0000000007442	5059.982874824331	1559
```

#### Output BAM file

In the above command we have included the -b flag which outputs the alignments that pass NanoCount's filtering thresholds into a new BAM for downstream analysis such ad QC or visualisation. The alignments are written in the same order as the source BAM file.

## **Reference**

If you use NanoCount please cite as follows:

Josie Gleeson, Adrien Leger, Yair D J Prawer, Tracy A Lane, Paul J Harrison, Wilfried Haerty, Michael B Clark, Accurate expression quantification from nanopore direct RNA sequencing with NanoCount, Nucleic Acids Research, Volume 50, Issue 4, 28 February 2022, Page e19, https://doi.org/10.1093/nar/gkab1129

If you use the dataset from SG-NEx in your work, please cite the following paper:

Chen, Ying, et al. “A systematic benchmark of Nanopore long read RNA
sequencing for transcript level analysis in human cell lines.” bioRxiv
(2021). doi: <https://doi.org/10.1101/2021.04.21.440736>
