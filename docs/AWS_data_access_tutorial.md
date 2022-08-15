# Accessing the SG-NEx dataset

SG-NEx data source contains long read (Oxford Nanopore) RNA sequencing data for commonly used cell lines. The data is hosted by AWS on S3 and can be accessed using direct links or the aws CLI.

The SG-NEx S3 bucket contains the following types of data:

   - [Raw sequencing signal (fast5)](#raw-sequencing-signal)            
   - [Basecalled sequences (fastq)](#basecalled-sequences)            
   - [Aligned sequences (bam)](#aligned-sequences)     
   - [Data visualisation tracks (bigwig/bigbed)](#data-visualisation-tracks)        
   - [Annotations](#annotations)            
   - [Processed data for RNA modification detection](#processed-data)     
   - [Sample and experiment information](#sample-and-experimental-data)               

 Below is the folder index for the open data bucket:

![folder indexing\!](/images/folder_index.png)

# Raw sequencing signal
To access raw sequencing (fast5) files:

```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fast5/ # list samples 
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fast5/sample_name .    # download fast5 files to your local directory
```

# Basecalled sequences
To access basecalled sequencing (fastq) files:

```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fastq/  # list samples 
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fastq/sample_name .   # download fastq files to your local directory
```
# Aligned sequences

We provide both genome and transcriptome aligned files:

```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/genome  # list samples inside this folder
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/genome/sample_name .   # download bam files that are aligned to genome 

aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome  # list samples inside this folder
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/sample_name .   # download bam files that are aligned to transcriptome
```
# Data visualisation tracks

We provide bigbed and bigwig files which can be directly visualised any genome browser. These files follow the UCSC chromosome naming convention and they can be directly visualised using the UCSC Genome Browser:

- Visualise the SG-NEx data in the UCSC Genome Browser(coming soon)

The files can be accessed and downloaded through S3 as well:
```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/genome_browser_data/bigbed/  # list all bigbed files
aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/genome_browser_data/bigwig/  # list all bigwig files
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/genome_browser_data/bigbed/sample_name.bigbed .   # download bigbed file for the a specific sample
```
# Annotations

The genome and transcriptome fasta files and the gtf file describing the genome annotations and which were used to process thedata can also be accessed. The latest SG-NEx data release used Ensembl version 91 (see [here](/docs/ANNOTATIONS.md) for links to original data). Two sets of annotations are provided in the bucket:

- Grch38 Ensembl annotations (without spike in RNAs) 
- Grch38 Ensembl + Sequin + SIRV and ERCC annotations

```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/annotations/genome_fasta/  # list included genome fasta files used for processing the sequencing data 
aws s3 sync --no-sign-request s3://sg-nex-data/data/annotations/genome_fasta .   # download genome fasta files used for processing the sequencing data 
```


![genome_fasta\!](/images/genome_fasta.png)


```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/annotations/transcriptome_fasta/  # list included transcriptome fasta files used for processing the sequencing data 
aws s3 sync --no-sign-request s3://sg-nex-data/data/annotations/transcriptome_fasta .   # download transcriptome fasta files used for processing the sequencing data 
```

![transcriptome_fasta\!](/images/transcriptome_fasta.png)


```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/annotations/gtf_file/  # list included annotation gtf files used in processing the sequencing data 
aws s3 sync --no-sign-request s3://sg-nex-data/data/annotations/gtf_file .  # download nnotation gtf files used for processing the sequencing data 
```
![gtf_file\!](/images/gtf_file.png)



# Processed data 

## RNA modification detection
 Long read direct RNA sequencing has allows the detection of RNA modification with RNA modification tools, such as [xPore](https://github.com/GoekeLab/xpore) and [m6Anet](https://github.com/GoekeLab/m6anet). To simplify the analysis of RNA modifications using the SG-Nex datasets, you can download the processed files to use with xPore and m6Anet. 
 
 To download the processed data for differential RNA modification analysis with xPore:
 ```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/processed_data/xpore/  # list all samples that have processed data for RNA modification detection using xPore
aws s3 sync --no-sign-request s3://sg-nex-data/data/processed_data/xpore/sample_name .  # download the json and index file needed for running xPore
```
To download the processed data for detection of m6A using m6Anet:
 ```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/processed_data/m6Anet/  # list all samples that have processed data for RNA modification detection using m6Anet
aws s3 sync --no-sign-request s3://sg-nex-data/data/processed_data/m6Anet/sample_name .  # download the json and index file needed for running m6Anet
```

These files are provided for a subset of samples, please see [here](/docs/samples_with_RNAmod_data.tsv) for the sample list with matched processed data for xPore and m6Anet.

# Sample and experimental data 

Detailed information for each sequencing sample is provided [here](/docs/samples.tsv). The data also includes multiplexed samples which share the same fast5 files. The information about the multiplexed samples can be found [here](/docs/multiplexed_samples.tsv). The files can also be accessed directly on S3:


 ```bash
aws s3 ls --no-sign-request s3://sg-nex-data/metadata/  # list metadata files
aws s3 sync --no-sign-request s3://sg-nex-data/metadata/ .  # download the metadata files
```
