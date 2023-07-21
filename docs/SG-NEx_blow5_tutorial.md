# **Basecalling and analysing SG-NEx samples in S/BLOW5 format**

In this tutorial, we will learn how to use the BLOW5 format for fast and efficient data manipulation and analysis.
Tutorial contains four parts: 1, basecalling a full SG-NEx sample directly from S/BLOW5 format; 2, extracting a subset of reads and basecall that subset. 3, running nanopolish/f5c eventalign; and, 4, an optional section on converting existing FAST5 data to BLOW5.

In this tutorial, we will be using a Nanopore direct RNA-Sequencing sample, one replicate from the K562 cell line.

## **Content**

- [Installation](#installation)
- [Data Access and Preparation](#data-access-and-preparation)
- [Running Buttery-eel on a full sample](#running-buttery-eel-on-a-full-sample)
- [Running Buttery-eel on a subset](#running-buttery-eel-on-a-subset-of-reads)
- [Running f5c eventalign](#running-f5c-eventalign)
- [Converting FAST5 to BLOW5](#converting-fast5-to-blow5)
- [Reference](#reference)

## **Installation**

### Installing buttery-eel

To directly basecall a S/BLOW5 file, we have to install [buttery-eel](https://github.com/Psy-Fer/buttery-eel), the S/BLOW5 basecaller wrapper for ONT Guppy. As of 18/02/2023, the [multiproc branch](https://github.com/Psy-Fer/buttery-eel/tree/multiproc) is recommended for efficient basecalling. Steps are briefly given below.

1. Download and setup ONT Guppy from https://community.nanoporetech.com/downloads and note the Guppy version. We recommend downloading the Linux 64-bit GPU version and then simply extracting the tarball.
2. Now clone buttery-eel multiproc branch:
```bash
git clone https://github.com/Psy-Fer/buttery-eel.git -b multiproc
cd  buttery-eel
```
3. modify `requirements.txt` in the repository base to have `ont-pyguppy-client-lib` version match the ONT Guppy version you downloaded. For instance, if your ONT Guppy version is 6.4.2, modify `requirements.txt` to have `ont-pyguppy-client-lib==6.4.2`.
4. Setup a Python virtual environment and setup buttery-eel.
```bash
python3 -m venv venv3
source ./venv3/bin/activate
pip install --upgrade pip
pip install --upgrade setuptools wheel
python setup.py install
```
5. Check if buttery-eel works:
```bash
buttery-eel --help
```
Please refer to the buttery-eel Github repository [here](https://github.com/Psy-Fer/buttery-eel) for detailed setup instructions. Open a [issue](https://github.com/Psy-Fer/buttery-eel/issues) if you encounter problems setting up buttery-eel.


### Installing slow5tools

For extracting a subset of reads from a BLOW5 file for the second part of the tutorial, we need slow5tools. We can install slow5tools on Linux x86_64 architecture using the following commands:

```bash
VERSION=v0.8.0
wget "https://github.com/hasindu2008/slow5tools/releases/download/$VERSION/slow5tools-$VERSION-x86_64-linux-binaries.tar.gz" && tar xvf slow5tools-$VERSION-x86_64-linux-binaries.tar.gz && cd slow5tools-$VERSION/
./slow5tools
```
For different systems and architecture, you may refer to the slow5tools GitHub repository [here](https://github.com/hasindu2008/slow5tools).


### Installing f5c

For the third part of the tutorial, we need [f5c](https://github.com/hasindu2008/f5c/), an optimised re-implementation of the index, call-methylation and eventalign modules in Nanopolish. We can install f5c on Linux x86_64 architecture using the following commands:

```
VERSION=v1.1
wget "https://github.com/hasindu2008/f5c/releases/download/$VERSION/f5c-$VERSION-binaries.tar.gz" && tar xvf f5c-$VERSION-binaries.tar.gz && cd f5c-$VERSION/
./f5c_x86_64_linux        # CPU version
./f5c_x86_64_linux_cuda   # cuda supported version
```

Nanopolish from version 0.14.0 onwards supports BLOW5, so you can use that instead of f5c if you wish. However, f5c is efficient and faster and produces the same results as nanopolish.

## **Data Access and Preparation**

Next, we will need to download the required data in BLOW5 format to run buttery-eel.

Generally, you may want to learn how to get access to these data using the [data
access
tutorial](AWS_data_access_tutorial.md). Below we only show the necessary steps to download the required data. The following command requires you to have [AWS CLI](https://aws.amazon.com/cli/) installed.

```bash
# download BLOW5 file for the K562 replicate 4 run 1 to the current directory
aws s3 cp --no-sign-request s3://sg-nex-data-blow5/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.blow5 ./
```

For the second part of the tutorial on [extract a subset of reads](#Running-buttery-eel-subset), you also need the BLOW5 index:

```bash
# download BLOW5 index file for the K562 replicate 4 run 1 to the current directory
aws s3 cp --no-sign-request s3://sg-nex-data-blow5/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.blow5.idx ./
```

You may also download the required data directly from the [SG-NEx AWS S3
bucket](http://sg-nex-data-blow5.s3-website-ap-southeast-1.amazonaws.com/) if you are unfamiliar with AWS CLI command.


## **Running Buttery-eel on a full sample**

Now we can run buttery-eel with these data.

``` bash
# source the buttery-eel virtual environment if you haven't already done so
source ./venv3/bin/activate

# running buttery-eel
buttery-eel  -g /path/to/ont-guppy/bin/  --config rna_r9.4.1_70bps_hac.cfg --device 'cuda:all' -i ./SGNex_K562_directRNA_replicate4_run1.blow5 -o  SGNex_K562_directRNA_replicate4_run1.fastq --port 5555  --use_tcp
```

The above command assumes that you have NVIDIA CUDA enabled GPUs and GPU version of Guppy is installed. Make sure to change `/path/to/ont-guppy/bin/` to where your installed ONT Guppy binary lives. The model `rna_r9.4.1_70bps_hac.cfg` is for direct RNA which is the case for this example. If you sample is cDNA, make sure to change the model to `dna_r9.4.1_450bps_hac.cfg` or `dna_r9.4.1_450bps_sup.cfg`. You can change the port `5555` to whatever port which is free on your system. Optionally you may use `-q NUM` option with buttery-eel to split basecalls to pass and fail based on mean quality score `NUM`.

On a Server with 4 NVIDIA Tesla V100 GPUs using buttery-eel commit [cdbe2eda940b2a4](https://github.com/Psy-Fer/buttery-eel/commit/cdbe2eda940b2a42b6e9f51c809683ba609d9aa4), it took ~10 minutes to basecall this example and consumed ~4GB of RAM.

### **Advanced tricks for the tech-savvy**

If you are doing your compute on AWS or if you have a very-high bandwidth and very-low latency Internet connection, given below is a trick to simply mounting the S3 bucket and then directly basecalling.  This eliminates the need to download the BLOW5 file.

First install s3fs and then mount the s3 public bucket:

```bash
# install s3fs first
sudo apt-get install s3fs
# directory to mount the bucket
mkdir ./s3
# command to mount a public bucket called sg-nex-data onto ./s3 DIRECTORY
s3fs sg-nex-data-blow5 ./s3/ -o public_bucket=1  -o url=http://s3.amazonaws.com/ -o dbglevel=info -o curldbg -o umask=0005 -o  uid=$(id -u)
```

Verify if the bucket is properly mounted by listing contents:
```bash
$ls ./s3/
index.html                               SGNex_HepG2_cDNAStranded_replicate4_run3
SGNex_A549_cDNA_replicate1_run2          SGNex_HepG2_cDNAStranded_replicate5_run4
SGNex_A549_cDNAStranded_replicate3_run3  SGNex_HepG2_directcDNA_replicate1_run1
SGNex_A549_cDNAStranded_replicate5_run2  SGNex_HepG2_directcDNA_replicate1_run2
...
```

Now simply call buttey-eel on the BLOW5 file inside the mounted S3 bucket:

``` bash
buttery-eel  -g /path/to/ont-guppy/bin/  --config rna_r9.4.1_70bps_hac.cfg --device 'cuda:all' -i ./s3/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.blow5 -o  SGNex_K562_directRNA_replicate4_run1.fastq --port 5555  --use_tcp
```

On a system with a single Tesla V100 GPU connected to Internet with a 1 Gbps connection, it took ~30 minutes. Parameters could be further optimised for performance, but not discussed in this tutorial.


## **Running Buttery-eel on a subset of reads**

Now let us see how we can quickly extract a subset of reads using SLOW5 and basecall those extracted reads. In this example, let us extract the reads that map to the transcript 'ENST00000564818.5'.

```bash
# copy the BAM file and the index for the sample
aws s3 cp --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.bam ./
aws s3 cp --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/bam/transcriptome/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.bam.bai ./

# extract the list of read IDs that map to ENST00000564818.5 using samtools and bash commands
samtools view  SGNex_K562_directRNA_replicate4_run1.bam ENST00000564818.5 | cut -f 1 | sort -u > read_ids_ENST00000564818.5.txt

# extract the raw reads from the BLOW5 file and index downloaded at the beginning of this tutorial
slow5tools get --list read_ids_ENST00000564818.5.txt SGNex_K562_directRNA_replicate4_run1.blow5 -o read_ids_ENST00000564818.5.blow5

# basecall those selected reads
buttery-eel  -g /path/to/ont-guppy/bin/  --config rna_r9.4.1_70bps_hac.cfg --device 'cuda:all' -i read_ids_ENST00000564818.5.blow5 -o  selected_reads.fastq --port 5555  --use_tcp
```

## Running f5c eventalign

Now let us execute the event alignment on the subset of data that we extracted and basecalled in the [previous section](#running-buttery-eel-on-a-subset-of-reads). 

```bash
# run minimap2 to align the reads to the reference (assuming you already have gencode.v40.transcripts.fa file as the reference)
minimap2 -ax map-ont -uf --secondary=no gencode.v40.transcripts.fa selected_reads.fastq > selected_reads.sam
samtools sort selected_reads.sam -o selected_reads.bam
samtools index selected_reads.bam

# execute f5c eventalign
f5c index selected_reads.fastq --slow5 read_ids_ENST00000564818.5.blow5
f5c eventalign --min-mapq 0 --rna --slow5 read_ids_ENST00000564818.5.blow5 \
    -r selected_reads.fastq -b selected_reads.bam -g gencode.v40.transcripts.fa \
    --signal-index  --scale-events --summary summary.txt -t 32 > eventalign.txt
```

If you are using nanopolish, the commands are slightly different:
```bash
nanopolish index selected_reads.fastq --slow5 read_ids_ENST00000564818.5.blow5
nanopolish eventalign -r selected_reads.fastq -b selected_reads.bam -g gencode.v40.transcripts.fa \
    --signal-index  --scale-events --summary summary.txt -t 32 > eventalign.txt
```

That is how we run f5c/nanopolish on an extracted subset. Of course, nothing stops us from running f5c on the full dataset if necessary. We already downloaded the full BLOW5 file (and index) and have the full basecalls from the previous sections. Simply do the following: 

```bash
# align
minimap2 -ax map-ont -uf --secondary=no -t32 gencode.v40.transcripts.fa SGNex_K562_directRNA_replicate4_run1.fastq > SGNex_K562_directRNA_replicate4_run1.sam
samtools sort SGNex_K562_directRNA_replicate4_run1.sam -o SGNex_K562_directRNA_replicate4_run1.bam
samtools index SGNex_K562_directRNA_replicate4_run1.bam

# execute f5c eventalign
f5c index SGNex_K562_directRNA_replicate4_run1.fastq --slow5 SGNex_K562_directRNA_replicate4_run1.blow5 \
    -t32 --skip-slow5-idx #we already have a downloaded blow5.idx, so we ask to skip making it
f5c eventalign --min-mapq 0 --rna --slow5 SGNex_K562_directRNA_replicate4_run1.blow5 \
    -r SGNex_K562_directRNA_replicate4_run1.fastq -b SGNex_K562_directRNA_replicate4_run1.bam  -g gencode.v40.transcripts.fa \
    --signal-index  --scale-events --summary summary.txt -t 32 > eventalign.txt
```

Note that f5c is efficient and faster, yet, produces the same results as nanopolish. For best performance with f5c, you can select a profile that suits your computer system using -x option as explained [here](https://hasindu2008.github.io/f5c/docs/profiles).


## Converting FAST5 to BLOW5

It is recommended that you directly download the already converted BLOW5 files from the s3 bucket as explained [before](#data-access-and-preparation). BLOW5 files are much smaller than FAST5 and are much faster to basecall and process. For the sake of completeness, here is how original FAST5 data can be converted to BLOW5.

Download the FAST5 raw signal data for conversion:

``` bash
# create a directory to store the data
mkdir blow5_convert_tutorial

# download the FAST5 archive
aws s3 cp --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fast5/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.tar.gz ./blow5_convert_tutorial

# extract the archive into a directory called fast5
mkdir blow5_convert_tutorial/fast5/
tar xf SGNex_K562_directRNA_replicate4_run1.tar.gz -C blow5_convert_tutorial/fast5/

```

Now we can run slow5tools on these data:

```bash
# convert each individual multi-FAST5 file to BLOW5
slow5tools f2s blow5_convert_tutorial/fast5/ -d blow5_convert_tutorial/slow5_tmp --retain  -p 16 #increase/decrease -p based on your processor count

# Merge the BLOW5 files to a single BLOW5 file
slow5tools merge blow5_convert_tutorial/slow5_tmp -o blow5_convert_tutorial/SGNex_K562_directRNA_replicate4_run1.blow5 -t 16 #increase/decrease -t based on your processor count

# Now index the BLOW5 file
slow5tools index blow5_convert_tutorial/SGNex_K562_directRNA_replicate4_run1.blow5

#clean up
rm -r blow5_convert_tutorial/SGNex_K562_directRNA_replicate4_run1.tar.gz blow5_convert_tutorial/fast5 blow5_convert_tutorial/slow5_tmp
```

You can convert other SG-NEx samples by following the steps above. However, note that those steps only will work for samples with multi-FAST5 files. Some of the samples in SG-NEx are in the early single-FAST5 format. For those samples, slow5tools f2s will fail with an error mentioning of "multiple runIDs in same file". For such samples, you may use the helper script that comes with slow5tools which is available [here](https://github.com/hasindu2008/slow5tools/blob/master/scripts/mixed-single-fast5-to-blow5.sh). Commands for using this script are:

```bash
# down the script
wget https://raw.githubusercontent.com/hasindu2008/slow5tools/master/scripts/mixed-single-fast5-to-blow5.sh

# set executable permission to script
chmod +x mixed-single-fast5-to-blow5.sh

#launch the script
./mixed-single-fast5-to-blow5.sh blow5_convert_tutorial/fast5/
```

If successful, a merged BLOW5 file called reads.blow5 will be created along with its index reads.blow5.idx. You can rename these files to what you want. There are a few samples where some of the original FAST5 files are corrupted where the above mentioned method will fail. Converting such samples need some manual FAST5 file curation which is too advanced to discussed here.

In summary, the easiest way is simply downloading the already converted BLOW5 files from the s3 bucket as explained [before](#data-access-and-preparation).


## **Reference**

If you use the dataset from SG-NEx in your work, please cite:
Chen, Ying, et al. “A systematic benchmark of Nanopore long read RNA sequencing for transcript level analysis in human cell lines.” bioRxiv (2021). doi: <https://doi.org/10.1101/2021.04.21.440736>

If you used S/BLOW5 in your work, please cite:
Gamaarachchi, H., Samarakoon, H., Jenner, S.P. et al. “Fast nanopore sequencing data analysis with SLOW5.” Nat Biotechnol 40, 1026–1029 (2022). https://doi.org/10.1038/s41587-021-01147-4

If you used buttery-eel in your work, please cite:
Samarakoon, Hiruna, et al. "Accelerated nanopore basecalling with SLOW5 data format." bioRxiv (2023). doi: https://doi.org/10.1101/2023.02.06.527365

If you used f5c in your work, please cite:
Gamaarachchi, H., Lam, C.W., Jayatilaka, G. et al. GPU accelerated adaptive banded event alignment for rapid comparative nanopore signal analysis. BMC Bioinformatics 21, 343 (2020). https://doi.org/10.1186/s12859-020-03697-x
