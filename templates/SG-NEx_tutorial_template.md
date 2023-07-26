# **TITLE** 

One paragraph on the aim of the tutorial

**Optional: can have a note message to describe the expected running time for this tutorial** 

## **Content**

- [Installation](#installation)
- [Data Access and Preparation](#data-access-and-preparation)
- [Running software](#running-software)
- [Reference](#reference)

## **Installation**

this session should describe how the software can be installed 




## **Data Access and Preparation**
Next, we will need to download the required data to run XXXX. The required data include:

-   a set of aligned reads to the genome from the A549 and HepG2 cell lines (bam files),
-   reference human genome annotations (gtf file, TxDb object, or Bambu
    Annotation object),
-   reference human genome sequence (fasta file or BSgenome).

Generally, you may want to learn how to get access to these data using the [data
access
tutorial](https://github.com/GoekeLab/sg-nex-data/blob/updated-documentation/docs/AWS_data_access_tutorial.md). Below we only show the necessary steps to download the required data. The following command requires you to have [AWS CLI](https://aws.amazon.com/cli/) installed.

``` bash
# create a directory to store the data
mkdir tutorial

# download genome fasta file 
aws s3 cp --no-sign-request s3://sg-nex-data/data/data_tutorial/annotations/hg38_chr22.fa ./tutorial

# download genome index fastai file 
aws s3 cp --no-sign-request s3://sg-nex-data/data/data_tutorial/annotations/hg38_chr22.fa.fai ./tutorial

# download gtf file
aws s3 cp --no-sign-request s3://sg-nex-data/data/data_tutorial/annotations/hg38_chr22.gtf ./tutorial

# download aligned bam files for A549 samples and HepG2 samples
aws s3 sync --no-sign-request s3://sg-nex-data/data/data_tutorial/bam ./tutorial --include *.bam 
```

You may also download the required data directly from the [SG-NEx AWS S3
bucket](http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/) if you are unfamiliar with AWS CLI command. They are stored in the `data/data_tutorial/bam` folder.


``` bash
# Note: Please make sure to replace the "sample_alias" with your sample name 
# To download genome bam files
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/genome/<sample_alias> ./tutorial
```

## **Running software**

This section describes how to run the software after downloading the datasets



### **Access to the required files for the complete workflow** 

If you wish to run the complete workflow of m6Anet, you can access all the required files (fast5, fastq and bam files) from the [SG-NEx S3 bucket](https://github.com/GoekeLab/sg-nex-data/blob/update-docs-aws/docs/samples.tsv). 


```bash
# Note: Please make sure to replace the "sample_alias" with your sample name

# To download the fast5 file
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fast5/<sample_alias> ./

# To download the fastq file
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fastq/<sample_alias> ./

# To download the bam file
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/<sample_alias> ./
```

Meanwhile, you can also access the transcriptome fasta files and the gtf file (annotation file) using the following command.

```bash 
# To download the transcriptome fasta file
aws s3 sync --no-sign-request s3://sg-nex-data/data/annotations/transcriptome_fasta ./ --exclude hg38*

# To download the gtf file
aws s3 cp --no-sign-request s3://sg-nex-data/data/annotations/gtf_file/Homo_sapiens.GRCh38.91.gtf ./
```
Alternatively, you can download the files from the [AWS S3 browser](http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/#data/annotations/) directly. 


## **Reference**

References to the software related publications and SGNEx can be added here.

If you use the dataset from SG-NEx in your work, please cite the following paper.

Chen, Ying, et al. “A systematic benchmark of Nanopore long read RNA
sequencing for transcript level analysis in human cell lines.” bioRxiv
(2021). doi: <https://doi.org/10.1101/2021.04.21.440736>
 
