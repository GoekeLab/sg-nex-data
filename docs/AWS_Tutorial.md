# On accessing SG-NEx datasets from AWS S3 bucket

### Bucket indexing

SG-NEx data source contains long read (Oxford Nanopore) RNA sequencing data for commonly used cancer cell lines. The annotation files used for processing these sequencing data are also stored in the bucket. Below is the folder index for the open data bucket
![folder indexing\!](/images/folder_index.png)

The bucket contains the following types of data

   - [Raw sequencing signals](#raw-sequencing-signals)            
   - [Basecalled sequences](#basecalled-sequences)            
   - [Aligned sequences](#aligned-sequences)             
   - [Annotations](#annotations)            
   - [Processed data](#processed-data)                    

The sample information is provided [here](/docs/samples.tsv). The data also include multiplexed samples, for those samples, they will share the same fast5 files, to find the sample mapping to mux samples, please refer the multiplexed sample info [here](/docs/multiplexed_samples.tsv)

# Raw sequencing signals
To access raw sequencing fast5 file

```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data/fast5/ # list samples 
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data/fast5/sample_name .    # download fast5 files to your local directory
```

# Basecalled sequences
To access basecalled sequencing fastq file

```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data/fastq/  # list samples 
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data/fastq/sample_name .   # download fastq files to your local directory
```
# Aligned sequences

We provide both genome and transcriptome aligned files

```bash
aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data/bam/genome  # list samples inside this folder
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data/bam/genome/sample_name .   # download bam files that are aligned to genome 

aws s3 ls --no-sign-request s3://sg-nex-data/data/sequencing_data/bam/transcriptome  # list samples inside this folder
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data/bam/transcriptome/sample_name .   # download bam files that are aligned to transcriptome
```

# Annotations

We provide genome fasta, gtf file and transcriptome fasta files to cater for all needs.

Two sets of annotations are provided in the bucket: 

- Grch38 Ensembl annotations: 
- Grch38 Ensembl + Sequin + SIRVERCCome annotations

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
 Long read RNA sequencing has allowed for detection of RNA modification with RNA modification tools, such as xPore and m6Anet. In the SG-Nex datasets, you can also find the processed data for xPore and m6Anet. 
 
 To download xpore processed data
 ```bash

aws s3 ls --no-sign-request s3://sg-nex-data/data/processed_data/xpore/  # list all samples that have processed data for RNA modification detection using xPore
aws s3 sync --no-sign-request s3://sg-nex-data/data/processed_data/xpore/sample_name .  # download the json and index file needed for running xPore
```
To download m6Anet processed data
 ```bash

aws s3 ls --no-sign-request s3://sg-nex-data/data/processed_data/m6Anet/  # list all samples that have processed data for RNA modification detection using m6Anet
aws s3 sync --no-sign-request s3://sg-nex-data/data/processed_data/m6Anet/sample_name .  # download the json and index file needed for running m6Anet
```

[Here](/docs/samples_with_RNAmod_data.tsv) you can find all samples with matched processed data for xPore and m6Anet
