# **Converting G-NEx samples to BLOW5**
=======

In this tutorial, we will convert a SG-NEx sample to S/BLOW5 format.
We will be using a Nanopore direct RNA-Sequencing sample, one replicate from the K562 cell line.

## **Content**

- [Installation](#installation)
- [Data Access and Preparation](#data-access-and-preparation)
- [Running slow5tools](#running-slow5tools)
- [Reference](#reference)

## **Installation**

First, we have to install slow5tools. We can install slow5tools on Linux x86_64 architecture using the following commands:

``` bash
VERSION=v0.8.0
wget "https://github.com/hasindu2008/slow5tools/releases/download/$VERSION/slow5tools-$VERSION-x86_64-linux-binaries.tar.gz" && tar xvf slow5tools-$VERSION-x86_64-linux-binaries.tar.gz && cd slow5tools-$VERSION/
./slow5tools
```
For different systems and architecture, you may refer to the slow5tools Github repository [here](https://github.com/hasindu2008/slow5tools).


## **Data Access and Preparation**

### **Download FAST5 data**

Next, we will need to download the FAST5 raw signal data for conversion. Generally, you may want to learn how to get access to these data using the [data
access
tutorial](https://github.com/GoekeLab/sg-nex-data/blob/updated-documentation/docs/AWS_data_access_tutorial.md). Below we only show the necessary steps to download the required data. The following command requires you to have [AWS CLI](https://aws.amazon.com/cli/) installed.

``` bash
# create a directory to store the data
mkdir blow5_convert_tutorial

# download the FAST5 archive
aws s3 cp --no-sign-request s3://sg-nex-data/data/sequencing_data_ont/fast5/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.tar.gz ./blow5_convert_tutorial

# extract the archive into a directory called fast5
mkdir blow5_convert_tutorial/fast5/
tar xf SGNex_K562_directRNA_replicate4_run1.tar.gz -C blow5_convert_tutorial/fast5/

```

You may also download the required data directly from the [SG-NEx AWS S3
bucket](http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/) if you are unfamiliar with AWS CLI command. They are stored in the `data/sequencing_data_ont/fast5/` folder.

## **Running Slow5tools**

Now we can run slow5tools on these data.

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

Already converted BLOW5 files and an indexes is also available to be directly downloaded.

```
aws s3 cp --no-sign-request s3://sg-nex-data/data/processed_data/blow5/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.blow5
aws s3 cp --no-sign-request s3://sg-nex-data/data/processed_data/blow5/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.blow5.idx
```

You can convert other SG-NEx samples by following the steps in this tutorial above. However, note that those steps only will work for samples with multi-FAST5 files. Some of the samples in SG-NEx are in the early single-FAST5 format. For those samples, slow5tools f2s will fail with an error mentioning of "multiple runIDs in same file". For such samples, you may use the helper script that comes with slow5tools which is available [here](https://github.com/hasindu2008/slow5tools/blob/master/scripts/mixed-single-fast5-to-blow5.sh). Commands for using this script are:

```


# down the script
wget https://raw.githubusercontent.com/hasindu2008/slow5tools/master/scripts/mixed-single-fast5-to-blow5.sh

#set executable permission to script
chmod +x mixed-single-fast5-to-blow5.sh

#launch the script
./mixed-single-fast5-to-blow5.sh blow5_convert_tutorial/fast5/
```

If successful, a merged BLOW5 file called reads.blow5 will be created along with its index reads.blow5.idx. You can rename these files to what you want.

There are a few samples where some of the original FAST5 files are corrupted where the above mentioned method will fail. Converting such samples need some manual FAST5 file curation which is too advanced to discussed here.


## **Reference**

If you use the
dataset from SG-NEx in your work, please cite the following paper.

Chen, Ying, et al. “A systematic benchmark of Nanopore long read RNA
sequencing for transcript level analysis in human cell lines.” bioRxiv
(2021). doi: <https://doi.org/10.1101/2021.04.21.440736>

If you used S/BLOW5 in your work, please cite the following paper.

Gamaarachchi, H., Samarakoon, H., Jenner, S.P. et al. “Fast nanopore sequencing data analysis with SLOW5.” Nat Biotechnol 40, 1026–1029 (2022). https://doi.org/10.1038/s41587-021-01147-4
