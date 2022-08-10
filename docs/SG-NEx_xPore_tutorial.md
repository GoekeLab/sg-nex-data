# **Analysing differential RNA modifications of SG-NEx samples** 

In this tutorial, we will perform differential RNA modification analysis on the SG-NEx samples. We will use the Nanopore direct RNA-sequencing data of A549 and HepG2 cell lines from the SG-NEx project for this analysis. The A549 cell line was extracted from the lung tissues of a patient with lung cancer. Meanwhile, the HepG2 cell line was extracted from hepatocellular carcinoma from a patient with liver cancer. A Python package known as **xPore**, will be used to identify differential RNA modification sites in these two cell lines. <br>
**Note: This tutorial may take 2 hours to complete(due to the use of large datasets).**

## **Content**

- [Installation](#installation)
- [Data Access and Preparation](#data-access-and-preparation)
- [Running xPore](#running-xpore)
- [Reference](#reference)

## **Installation**

xPore requires Python 3 to run. It can be installed through PyPI using the following command. Note: This command requires pip3 to be installed. 


```bash
pip install xpore
```

Alternatively, it can be installed from our GitHub repository.


```bash
git clone https://github.com/GoekeLab/xpore.git
cd xpore
python setup.py install
cd ..  # return to the working directory 
```

## **Data Access and Preparation**

xPore requires a control group and a treatment group to work. For this tutorial, we will use the data from A549_directRNA_replicate6_run1 as the control and HepG2_directRNA_replicate6_run1 as the treatment group. The required data include:

- feature values from both the A549 and HepG2 samples (two json files, one for each sample),
- indexes of the feature values file. (two index files, one from each sample)

### **Download data for xPore**
The data files can be downloaded from the [SG-NEx AWS S3 bucket](http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/). Please refer to this [page](https://github.com/GoekeLab/sg-nex-data/blob/update-docs-aws/docs/AWS_data_access_tutorial.md) for a comprehensive tutorial to access the SG-NEx dataset.

To begin, we will create directories to store the required files for the respective samples. 
```bash
# create directories for A549 and HepG2 samples
mkdir A549_directRNA_replicate6_run1
mkdir HepG2_directRNA_replicate6_run1
```

Then, we will access the S3 bucket to download the data into the created directories. There are two methods to download the data. 
<br>

#### **Method One(Recommended):** 
Note: The command requires you to have [AWS CLI](https://aws.amazon.com/cli/) installed.

**For A549 sample:**
```bash
# To download the processed data for A549 sample
aws s3 sync --no-sign-request s3://sg-nex-data/data/processed_data/xpore/SGNex_A549_directRNA_replicate6_run1 ./A549_directRNA_replicate6_run1
# The directory should now contain two files. They are the json and index files.
```

**For HepG2 sample:**
```bash
# To download the processed data for HepG2 sample
aws s3 sync --no-sign-request s3://sg-nex-data/data/processed_data/xpore/SGNex_HepG2_directRNA_replicate6_run1 ./HepG2_directRNA_replicate6_run1
# The directory should now contain two files. They are the json and index files.
```

Alternatively, we can use the URL to the S3 bucket to download the processed data. 
<br>

#### **Method Two:**
**For A549 sample:**
```bash
cd A549_directRNA_replicate6_run1 

# To download the json file
wget https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/processed_data/xpore/SGNex_A549_directRNA_replicate6_run1/data.json

# To download the index file
wget https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/processed_data/xpore/SGNex_A549_directRNA_replicate6_run1/data.index

# Exit the A549 directory to return to our working directory
cd ..
```

**For HepG2 sample:**
```bash
cd HepG2_directRNA_replicate6_run1 

# To download the json file
wget https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/processed_data/xpore/SGNex_HepG2_directRNA_replicate6_run1/data.json	

# To download the index file
wget https://sg-nex-data.s3.ap-southeast-1.amazonaws.com/data/processed_data/xpore/SGNex_HepG2_directRNA_replicate6_run1/data.index

# Exit the HepG2 directory to return to our working directory 
cd ..
```
<br>

Then, we need to prepare a `YAML` file to specify the experimental design, the data directories, the output directory, and the method options. We can create a `.yml` file using the command below: 

```bash
# To create .yml file
cat > config.yml
```

Next, **copy** and **paste** the command below to the **terminal**!
```yaml
notes: Pairwise comparison without replicates with default parameter setting.

data:
    Control:  # Condition Name 1
        rep1: ./A549_directRNA_replicate6_run1   # Input directory containing the json and index file
    Treatment:  # Condition Name 2
        rep1: ./HepG2_directRNA_replicate6_run1  # Input directory containing the json and index file
    
out: ./out    # Output directory 

method:
    prefiltering:  # To remove positions that are unlikely to be differentially modified
        method: t-test
        threshold: 0.1
```
#### **Please press the `Enter` key followed by the `Ctrl` and `D` keys to save the file. This file is now saved as `config.yml`**
<br>

**Note: If you want xPore to test every genomic/transcriptomic position, you may remove the prefiltering section.** <br>
Please refer to [xPore](https://xpore.readthedocs.io/en/latest/configuration.html) documentation for more information on the .yml configuration file. 
<br>

You may refer to this [page](https://github.com/GoekeLab/sg-nex-data/blob/update-docs-aws/docs/samples_with_RNAmod_data.tsv) for the list of URLs to access the processed data for xPore from different samples. 


## **Running xPore** 

Now, we are all set to run xPore to identify differential RNA modification. 


```bash
# Run xPore
# Before running the command, please ensure that you are at the directory where the .yml configuration file is located.
xpore diffmod --config config.yml --n_processes 8
```
**Note: You can reduce the "n_processes" if you have lesser processes available in your machine. However, lowering the "n_processes" may increase the run time. Alternatively, you may visit the [xPore documentation](https://xpore.readthedocs.io/en/latest/data.html) for downsampled dataset (demo.tar.gz file), which can reduce the total run time to less than ten minutes.**


xPore returns a table `(diffmod.table)` consisting of differential RNA modification rate across all tested positions in the output directory that we specify in the `.yml` file. Please refer to the [xPore documentation](https://xpore.readthedocs.io/en/latest/outputtable.html) for more information regarding the output table. 

<br>


#### **NOTE: This is a short tutorial to demonstrate the usage of xPore on the SG-NEx data. Please refer to the [xPore documentation](https://xpore.readthedocs.io/en/latest/quickstart.html#) for the complete workflow to perform differential RNA modification analysis.** 
<br>

### **Access to the required files for complete workflow** 

If you wish to run the complete workflow of xPore, you can access all the required files (fast5, fastq and bam files) from the [SG-NEx S3 bucket](https://github.com/GoekeLab/sg-nex-data/blob/update-docs-aws/docs/samples.tsv). 


```bash
# Note: Please make sure to replace the "sample_alias" with your sample name
# To download the fast5 file
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fast5/<sample_alias> ./

# To download the fastq file
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fastq/<sample_alias> ./

# To download the bam file
aws s3 sync --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/<sample_alias> ./
```

Meanwhile, you can also access the genome and transcriptome fasta files and the gtf file (annotation file) using the following command.

```bash 
# To download the transcriptome fasta file
aws s3 sync --no-sign-request s3://sg-nex-data/data/annotations/transcriptome_fasta ./ --exclude hg38*

# To download the gtf file
aws s3 cp --no-sign-request s3://sg-nex-data/data/annotations/gtf_file/Homo_sapiens.GRCh38.91.gtf ./
```
Alternatively, you can download the files with their respective [URL](https://github.com/GoekeLab/sg-nex-data/blob/update-docs-aws/docs/samples.tsv). 


## **Reference**
In this tutorial, we performed differential RNA modification analysis on the [SG-NEx](https://github.com/GoekeLab/sg-nex-data) dataset using [xPore](https://github.com/GoekeLab/xpore). If you use xPore and the dataset from SG-NEx in your work, please cite the following paper. 

Pratanwanich, Ploy Naruemon, et al. "Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore." Nature Biotechnology (2021). doi: https://doi.org/10.1038/s41587-021-00949-w

Chen, Ying, et al. "A systematic benchmark of Nanopore long read RNA sequencing for transcript level analysis in human cell lines." bioRxiv (2021). doi: https://doi.org/10.1101/2021.04.21.440736