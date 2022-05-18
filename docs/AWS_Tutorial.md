# On accessing SG-NEx datasets from AWS open data.

### Bucket indexing

SG-NEx data source contains long read (Oxford Nanopore) RNA sequencing data for commonly used cancer cell lines. The annotation files used for processing these sequencing data are also stored in the bucket. Below is a indexing figure for the open data bucket
![folder indexing\!](/docs/folder_indexing.png)

The bukcet contains the following types of data

   - [Raw sequencing signals](#fast5-data)            
   - [Basecalled sequences](#fastq-data)            
   - [Aligned sequences](#bam-data)             
   - [Annotations](#annotations)            
   - [Processed data](#processed-data)                    

The sample information are provided [here](/docs/sample_information.tsv). The data also include multiplexed samples, for those samples, they will share the same fast5 files, to find the sample mapping to mux samples, please refer the multiplexed sample info [here](/docs/multiplexed_samples.tsv)

### Raw sequencing signals
To access raw sequencing fast5 file

```bash
aws s3 ls --no-sign-request sg-nex-data/data/sequencing_data/fast5/ # list samples 
aws s3 cp --no-sign-request sg-nex-data/data/sequencing_data/fast5/sample_name . --recursive  # download fast5 files to your local directory
```

### Basecalled sequences
To access basecalled sequencing fastq file

```bash
aws s3 ls --no-sign-request sg-nex-data/data/sequencing_data/fastq/  # list samples 
aws s3 cp --no-sign-request sg-nex-data/data/sequencing_data/fastq/sample_name . --recursive  # download fastq files to your local directory
```
### Aligned sequences

We provide both genome and transcriptome aligned files

```bash
aws s3 ls --no-sign-request sg-nex-data/data/sequencing_data/bam/genome  # list samples inside this folder
aws s3 cp --no-sign-request sg-nex-data/data/sequencing_data/bam/genome/sample_name . --recursive  # download bam files that are aligned to genome 

aws s3 ls --no-sign-request sg-nex-data/data/sequencing_data/bam/transcriptome  # list samples inside this folder
aws s3 cp --no-sign-request sg-nex-data/data/sequencing_data/bam/transcriptome/sample_name . --recursive  # download bam files that are aligned to transcriptome
```


