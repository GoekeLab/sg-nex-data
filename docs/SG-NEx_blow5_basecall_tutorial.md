# **Directly Basecalling a SG-NEx samples in S/BLOW5 format**

In this tutorial, we will perform basecalling of a SG-NEx sample using S/BLOW5 format.
We will be using a Nanopore direct RNA-Sequencing sample, one replicate from the K562 cell line.

## **Content**

- [Installation](#installation)
- [Data Access and Preparation](#data-access-and-preparation)
- [Running Buttery-eel](#Running-buttery-eel)
- [Reference](#reference)

## **Installation**

To directly basecall a S/BLOW5 file, we have to install [buttery-eel](https://github.com/Psy-Fer/buttery-eel), the S/BLOW5 basecaller wrapper for ONT Guppy. As of 18/02/2023, the [multiproc branch](https://github.com/Psy-Fer/buttery-eel/tree/multiproc) is recommended for efficient basecalling. Steps are briefly given below.

1. Download and setup ONT Guppy from https://community.nanoporetech.com/downloads and note the Guppy version. We recommend downloading the Linux 64-bit GPU version and then simply extracting the tarball.
2. Now clone buttery-eel multiproc branch:
```
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

## **Data Access and Preparation**

Next, we will need to download the required data in BLOW5 format to run buttery-eel.

Generally, you may want to learn how to get access to these data using the [data
access
tutorial](https://github.com/GoekeLab/sg-nex-data/blob/updated-documentation/docs/AWS_data_access_tutorial.md). Below we only show the necessary steps to download the required data. The following command requires you to have [AWS CLI](https://aws.amazon.com/cli/) installed.

```bash
# download BLOW5 file for the K562 replicate 4 run 1 to the current directory
aws s3 cp --no-sign-request s3://sg-nex-data/data/processed_data/blow5/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.blow5 ./
```

You may also download the required data directly from the [SG-NEx AWS S3
bucket](http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/) if you are unfamiliar with AWS CLI command. They are stored in the `processed_data/blow5/` folder.

If you want to convert the FAST5 data yourself, you can refer the conversion tutorial [here](SG-NEx_blow5_conversion_tutorial.md).


## **Running Buttery-eel**

Now we can run buttery-eel with these data.

``` bash
# source the butter-eel virtual environment if you haven't already done so
source ./venv3/bin/activate

# running buttery-eel
buttery-eel  -g /path/to/ont-guppy/bin/  --config rna_r9.4.1_70bps_hac.cfg --device 'cuda:all' -i ./SGNex_K562_directRNA_replicate4_run1.blow5 -o  SGNex_K562_directRNA_replicate4_run1.fastq --port 5555  --use_tcp
```

The above command assumes that you have NVIDIA CUDA enabled GPUs and GPU version of Guppy is installed. Make sure to change `/path/to/ont-guppy/bin/` to where your installed ONT Guppy binary lives. The model `rna_r9.4.1_70bps_hac.cfg` is for direct RNA which is the case for this example. If you sample is cDNA, make sure to change the model to `dna_r9.4.1_450bps_hac.cfg` or `dna_r9.4.1_450bps_sup.cfg`. You can change the port `5555` to whatever port which is free on your system. Optionally you may use `-q NUM` option with buttery-eel to split basecalls to pass and fail based on mean quality score `NUM`.

On a Server with 4 NVIDIA Tesla V100 GPUs using buttery-eel commit [cdbe2eda940b2a4](https://github.com/Psy-Fer/buttery-eel/commit/cdbe2eda940b2a42b6e9f51c809683ba609d9aa4), it took ~10 minutes to basecall this example and consumed ~4GB of RAM.

## **Advanced tricks for the tech-savvy**

If you are doing your compute on AWS or if you have a very-high bandwidth and very-low latency Internet connection, given below is a trick to simply mounting the S3 bucket and then directly basecalling.  This eliminates the need to download the BLOW5 file.

First install s3fs and then mount the s3 public bucket:

```bash
# install s3fs first
sudo apt-get install s3fs
# directory to mount the bucket
mkdir ./s3
# command to mount a public bucket called sg-nex-data onto ./s3 DIRECTORY
s3fs sg-nex-data ./s3/ -o public_bucket=1  -o url=http://s3.amazonaws.com/ -o dbglevel=info -o curldbg -o umask=0005 -o  uid=$(id -u)
```

Verify if the bucket is properly mounted by listing contents:

```bash
ls ./s3/
data  index.html  metadata  README  RELEASE_NOTE
```

Now simply call buttey-eel on the BLOW5 file inside the mounted S3 bucket:

``` bash
buttery-eel  -g /path/to/ont-guppy/bin/  --config rna_r9.4.1_70bps_hac.cfg --device 'cuda:all' -i ./s3/data/processed_data/blow5/SGNex_K562_directRNA_replicate4_run1/SGNex_K562_directRNA_replicate4_run1.blow5 -o  SGNex_K562_directRNA_replicate4_run1.fastq --port 5555  --use_tcp
```

On a system with a single Tesla V100 GPU connected to Internet with a 1 Gbps connection, it took ~30 minutes. Parameters could be further optimised for performance, but not discussed in this tutorial.

## **Reference**

If you use the dataset from SG-NEx in your work, please cite:
Chen, Ying, et al. “A systematic benchmark of Nanopore long read RNA
sequencing for transcript level analysis in human cell lines.” bioRxiv
(2021). doi: <https://doi.org/10.1101/2021.04.21.440736>

If you used S/BLOW5 in your work, please cite:
Gamaarachchi, H., Samarakoon, H., Jenner, S.P. et al. “Fast nanopore sequencing data analysis with SLOW5.” Nat Biotechnol 40, 1026–1029 (2022). https://doi.org/10.1038/s41587-021-01147-4

If you used butter-eel in your work, please cite:
Samarakoon, Hiruna, et al. "Accelerated nanopore basecalling with SLOW5 data format." bioRxiv (2023). doi: https://doi.org/10.1101/2023.02.06.527365
