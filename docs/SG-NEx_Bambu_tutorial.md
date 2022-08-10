# **Transcript discovery and quantification of SG-NEx samples**

In this tutorial, we will perform novel transcript discovery and
quantification on the SG-NEx samples. We will be using six Nanopore direct RNA-Sequencing
samples, three replicates each from the A549 and
HepG2 cell lines. The A549 cell line was extracted from lung tissues
from a patient with lung cancer whereas HepG2 was extracted from
hepatocellular carcinoma from a patient with liver cancer. We will use
Bambu, a R package hosted on the Bioconductor platform to identify and
quantify novel isoforms in these cell lines. 

**Note: This tutorial may take 10 minutes to complete.**

## **Content**

- [Installation](#installation)
- [Data Access and Preparation](#data-access-and-preparation) 
- [Running Bambu](#running-bambu)
- [Reference](#reference)

## **Installation**

First, we have to install Bambu. Before that, make sure you have R
(version \>= 4.1) installed on your machine. We can install Bambu using the following command:

``` bash
R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("bambu")

q()
```
If you want a more recent version of Bambu, you may refer to the Bambu Github repository [here](https://github.com/GoekeLab/bambu). 

## **Data Access and Preparation**
### **Download Data for Bambu**
Next, we will need to download the required data to run Bambu. The required data include:

-   a set of aligned reads to the genome from the A549 and HepG2 cell lines (bam files),
-   reference human genome annotations (gtf file, TxDb object, or Bambu
    Annotation object),
-   reference human genome sequence (fasta file or BSgenome).

Generally, you may want to learn how to get access to these data using the [data
access
tutorial](https://github.com/GoekeLab/sg-nex-data/blob/updated-documentation/docs/AWS_data_access_tutorial.md). Below we only show the necessary steps to download the required data. The following command requires you to have [AWS CLI](https://aws.amazon.com/cli/) installed.

``` bash
# create a directory to store the data
mkdir bambu_tutorial

# download genome fasta file 
aws s3 cp --no-sign-request s3://sg-nex-data/data/data_tutorial/annotations/hg38_chr22.fa ./bambu_tutorial

# download genome index fastai file 
aws s3 cp --no-sign-request s3://sg-nex-data/data/data_tutorial/annotations/hg38_chr22.fa.fai ./bambu_tutorial

# download gtf file
aws s3 cp --no-sign-request s3://sg-nex-data/data/data_tutorial/annotations/hg38_chr22.gtf ./bambu_tutorial

# download aligned bam files for A549 samples and HepG2 samples
aws s3 sync --no-sign-request s3://sg-nex-data/data/data_tutorial/bam ./bambu_tutorial --include *.bam 
```

You may also download the required data directly from the [SG-NEx AWS S3
bucket](http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/) if you are unfamiliar with AWS CLI command. They are stored in the `data/data_tutorial/bam` folder.

**NOTE: We have downsampled the Hg38 genome, A549 and HepG2 samples to ensure this tutorial can be completed in 10 minutes. If you want to run Bambu on the original samples, you can find the sample name [here](https://github.com/GoekeLab/sg-nex-data/blob/updated-documentation/docs/samples.tsv) and amend it into the following code chunk:**

```bash
# Note: Please make sure to replace the "sample_alias" with your sample name 
# To download genome bam files
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/genome/<sample_alias> ./bambu_tutorial
```
### **Prepare Data for Bambu**

All required data are now stored in the `bambu_tutorial` folder of the
current working directory. Next, we prepare the data to run Bambu.

``` bash
R  

# set work directory if you are in a different directory
setwd("bambu_tutorial")

# data preparation
library(bambu)
fa.file <- "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
gtf.file <- "Homo_sapiens.GRCh38.91.gtf"
annotations <- prepareAnnotations(gtf.file) # This function creates a reference annotation object which is used for transcript discovery and quantification in Bambu.
samples.bam <- list.files(".", pattern = ".bam$", full.names = TRUE)
```

## **Running Bambu**

Now we can run Bambu with these data. For a
faster running speed, you can increase the `ncore` parameter up
to the total number of samples at your availability. 

``` bash
# running Bambu 
se <- bambu(reads = samples.bam, annotations = annotations, genome = fa.file, ncore = 2)  
```

Bambu returns a `SummarizedExperiment` object with the genomic
coordinates of the annotated & novel transcripts and their expression
estimates. They can be assessed using the following code:

``` bash
assays(se) #returns the transcript abundance estimates as counts or CPM.
rowRanges(se) #returns a GRangesList (with genomic coordinates) with all annotated and newly discovered transcripts.
rowData(se) #returns additional information about each transcript such as the gene name and the class of the newly discovered transcript.
```
This `SummarizedExperiment` object can also be used for further downstream analysis (eg. [DESeq](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)). If you want to save the transcript &  genomic annotations and their expression
estimates, you can then write them into an `output` folder using the `writeBambuOutput` function.

``` bash
writeBambuOutput(se, path = "./output")

q()
```

The files in the `output` folder is described below:

| Output file name                | Description                                                             |
|:----------------------------|:------------------------------------------|
| extended_annotations.gtf        | Extended transcript & gene annotations for the genome using long reads data.        |
| counts_transcript.txt           | Total read counts estimates for each transcript in each sample.        |
| CPM_transcript.txt              | Counts per million (CPM) estimates for each transcript in each sample. |
| fullLengthCounts_transcript.txt | Full length read counts estimates for each transcript in each sample.  |
| uniqueCounts_transcript.txt                | Unique read counts estimates for each transcript in each sample.       |
| counts_gene.txt                 | Gene read counts estimates for each transcript in each sample.         |

**NOTE: This is a short tutorial to demonstrate the usage of Bambu on the SG-NEx data. Please refer to the [Bambu documentation](https://github.com/GoekeLab/bambu) for a more complete workflow in novel transcript discovery and quantification.**

## **Reference**

In this tutorial, we extended the existing transcript & gene annotations
on the [SGNEx](https://github.com/GoekeLab/sg-nex-data) dataset using
[Bambu](https://github.com/GoekeLab/bambu). If you use Bambu and the
dataset from SG-NEx in your work, please cite the following paper.

Chen, Ying, et al. “A systematic benchmark of Nanopore long read RNA
sequencing for transcript level analysis in human cell lines.” bioRxiv
(2021). doi: <https://doi.org/10.1101/2021.04.21.440736>
