![The Singapore Nanopore-Expression Project\!](
https://jglaborg.files.wordpress.com/2021/10/sg_nex_textlogo.png)

[![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/GoekeLab/sg-nex-data?color=blue&include_prereleases)](#data-release-and-access) [![cell lines](https://img.shields.io/badge/cell_lines-13-green)](#data-release-and-access) [![Sequencing Experiments](https://img.shields.io/badge/sequencing_runs-112-green)](docs/samples.tsv) 

The SG-NEx project is an international collaboration initiated at the [Genome Institute of Singapore](https://www.a-star.edu.sg/gis/) to provide reference transcriptomes for 5 of the most commonly used cancer cell lines using Nanopore long read RNA-Seq data:

![The Singapore Nanopore-Expression Project - Design\!](
https://jglaborg.files.wordpress.com/2020/10/sg_nex_design-1.png)

Transcriptome profiling is done using PCR-cDNA sequencing ("PCR-cDNA"), amplification-free cDNA sequencing ("direct cDNA"), direct sequencing of native RNA (“direct RNA”), and short read RNA-Seq. All samples are sequenced with at least 3 high quality replicates. For a subset of samples spike-in RNAs are included and matched m6A profiling is available.

The raw, aligned, and processed data is hosted on the AWS open data registry (see below for data access and analysis tutorial).


## Content

- [Email list](#sign-up-for-data-release-notifications-and-updates)
- [Latest Data Release and Access](#data-release-and-access)
- [Browse the data](#browse-the-data)
- [Data Processing](#data-processing)
- [Use Cases and Applications](#use-cases-and-applications)
- [Data Analysis Tutorials](#data-analysis-tutorials-and-workflows)
- [Contributing](#contributing)
- [Acknowledgements](#acknowledgements)
- [Citing the SG-NEx project](#citing-the-sg-nex-project)
- [Contact](#contact)

## Sign up for data release notifications and updates
You can sign up for the sg-nex-updates email list to receive notifications about upcoming data releases: 

https://groups.google.com/forum/#!forum/sg-nex-updates/join

## Data Release and Access

**Latest Release (v0.5)**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10795697.svg)](https://doi.org/10.5281/zenodo.10795697)

This release includes 113 samples from 13 different cell lines. 

**Data Access**

You can access the following data through the [AWS Open Data Registry](https://registry.opendata.aws/sgnex/):

- raw files (fast5)
- raw files (blow5)
- basecalled files (fastq)
- aligned reads (genome and transcriptome) (bam)
- tracks for visualisation (bigwig and bigbed)
- processed data for differential RNA modification analysis (json, for use with xPore)
- processed data for identification of m6A (json, for use with m6Anet)
- annotation files
- detailed sample and experiment information

You can browse the S3 data here: 1) [fast5, fastq, and bam](http://sg-nex-data.s3-website-ap-southeast-1.amazonaws.com/) and 2) [blow5](http://sg-nex-data-blow5.s3-website-ap-southeast-1.amazonaws.com/). 

Please refer to the [data access tutorial](docs/AWS_data_access_tutorial.md) which describes the S3 data structure and how to access files with [AWS CLI](https://docs.aws.amazon.com/cli/latest/reference/s3/). The direct links to the data are listed in the [sample spreadsheet](docs/samples.tsv).

Here are the locations for the spike-in concentrations used in SG-NEx samples:
- [Sequin concentration](docs/RNAsequins_MixA.xlsx)
- [SIRV-1 concentration](https://www.lexogen.com/wp-content/uploads/2021/10/025UI368V0200_SIRV-Calculator_Set-1_extended.xlsx)
- [SIRV-4 concentration](https://www.lexogen.com/wp-content/uploads/2021/10/141UI371V0200_SIRV-Calculator_Set-4_extended.xlsx)

_**Citation**_: Please cite the pre-print describing the SG-NEx data resource when using these data, and add the following details: "The SG-NEx data was accessed on [DATE] at registry.opendata.aws/sg-nex-data".

Chen, Y. _et al._ "A systematic benchmark of Nanopore long read RNA sequencing for transcript level analysis in human cell lines." _bioRxiv_ (2021). doi: https://doi.org/10.1101/2021.04.21.440736

**Release Note & Updates**

Version Number: V0.5.1       
Date: 2024-04-15
Release of new sample     
- new RNA004 sample of Hek293T (SGNex_Hek293T_directRNA_replicate5_run1)
- pod5, fastq, genome and transcriptome aligned bam files are included in this release 

Version Number: V0.5.0       
Date: 2024-03-08         
Release of new samples      
- direct RNA data for H9 and HEYA8 samples       
- cDNA and direct cDNA samples for H9 and HEYA8     
- cDNA promethion samples of Hct116 samples using SQK-PCS110 (100 million reads on average)     
- cDNA sample of Hct116 sampe using the SQK-PCS111     

Update of existing sample files      
- SGNex_MCF7_cDNAStranded_replicate2_run1.fastq.gz additional info characters removed before @ for the first read      
- SGNex_K562_cDNAStranded_replicate3_run3.fastq.gz  line48000 added 1 character of “ for quality to match sequence length        
- SGNex_A549_directRNA_replicate5_run1.tar.gz updated as previous version is incomplete        
- SGNex_MCF7-EV_directRNA_replicate1_run1.fastq.gz updated on ENA as it is a duplicated file        
- SGNex_MCF7_directRNA_replicate2_run2  fixed with this command “zcat SGNex_MCF7_directRNA_replicate2_run2.fastq.gz | sed 's/.*@/@/g' | sed '$d' | gzip > SGNex_MCF7_directRNA_replicate2_run2_fixed.fastq.gz” thanks to Alex 


Version Number: V0.4.0                
Date: 2023-03-06                          
Update of the SG-NEx data on AWS. Includes raw signal data in blow5 format. 

Version Number: V0.3.0               
Date: 2022-07-28                 
Initial release of the SG-NEx data on AWS. Includes Nanopore direct RNA, cDNA, direct cDNA-Seq, short read RNA-Seq and m6ACE-Seq.

**Release History**

You can find previous releases here in the [release history](https://github.com/GoekeLab/sg-nex-data/releases)

## Browse the data

You can now browse the data using the UCSC genome browser:

[View the SG-NEx data in the UCSC Genome Browser](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=https://sg-nex-data.s3.amazonaws.com/data/sequencing_data_ont/genome_browser_data/hub_track/hub-ONT-Grch38-complete-2024-01-19.txt)

By default only selected tracks are shown, but you can visualise all reads (bigbed tracks) and their coverage tracks (bigwig) from each individual sample.

## Data Processing

All data was aligned against the human genome version Grch38 (please refer to the [data access tutorial](docs/AWS_data_access_tutorial.md) for reference files). We collaborated with [nf-core](https://github.com/nf-core) to develop [nanoseq](https://github.com/nf-core/nanoseq), a standardardized pipeline for Nanopore RNA-Seq data processing. 

## Use Cases and Applications

You can browse a list of articles that review or use the SG-NEx data [here](/docs/SGNEx_usecases.md). If you have used the data for your own research, feel free to add a publication entry.

## Data Analysis Tutorials and Workflows

The following short tutorials are available that demonstrate how to analyse the SG-NEx data:

- [Transcript discovery and quantification of SG-NEx samples (using Bambu)](./docs/SG-NEx_Bambu_tutorial_notebook.ipynb) <a target="_blank" href="https://colab.research.google.com/github/GoekeLab/sg-nex-data/blob/master/docs/SG-NEx_Bambu_tutorial_notebook.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

- [Analysing differential RNA modifications of SG-NEx samples (using xPore)](./docs/SG-NEx_xPore_tutorial.md)

- [Identification of m6A with the SG-NEx samples (using m6Anet)](./docs/SG-NEx_m6Anet_tutorial.md)

- [Basecalling and analysing SG-NEx samples in S/BLOW5 format](./docs/SG-NEx_blow5_tutorial.md)
  
- [Isoform discovery and Quantification with FLAIR](./docs/SG_NEx_FLAIR_tutorial_notebook.ipynb) <a target="_blank" href="https://colab.research.google.com/github/GoekeLab/sg-nex-data/blob/master/docs/SG_NEx_FLAIR_tutorial_notebook.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

- [Transcript reconstruction and quantification of SG-NEx samples with IsoTools](./docs/SG_NEx_IsoTools_tutorial.ipynb)<a target="_blank" href="https://colab.research.google.com/github/GoekeLab/sg-nex-data/blob/master/docs/SG_NEx_IsoTools_tutorial.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>

- [Transcript quantification of SG-NEx direct RNA samples with NanoCount](./docs/SG-NEx_NanoCount_tutorial.md)

Additional, more detailed workflows can be found here:

- [Transcript discovery, quantification, and differential transcript expression from long read RNA-Seq data (using Bambu)](https://github.com/GoekeLab/bambu)

- [Identification of differential RNA modifications using a METTL3 knockout cell line (using xPore)](./docs/xPore_ONT_tutorial.ipynb)

- [Analysing transcriptome-wide m6A modifications (using m6Anet)](https://m6anet.readthedocs.io/en/latest/)


## Contributing
We welcome contributions from all long read RNA-seq tool developers! You may follow the steps below to contribute:

- Fork this repository
- Add your tutorial document to the docs folder 
- Adding your tutorial workflow link in the Data Analysis Tutorials and Workflows section in README.md in this format: [tutorial title](path_to_the_tutorial)
- Submit a pull request.


## Acknowledgements
**GIS Sequencing Platform and Data Generation**            
Hwee Meng Low, Yao Fei, Sarah Ng, Wendy Soon, CC Khor   

**Cancer Genomics and RNA Modifications**            
Viktoriia Iakovleva, Puay Leng Lee, Lixia Xin, Hui En Vanessa Ng, Jia Min Loo, Xuewen Ong, Hui Qi Amanda Ng, Suk Yeah Polly Poon, Hoang-Dai Tran, Kok Hao Edwin Lim, Huck Hui Ng, Boon Ooi Patrick Tan, Huck-Hui Ng, N.Gopalakrishna Iyer, Wai Leong Tam, Wee Joo Chng, Leilei Chen, Ramanuj DasGupta, Yun Shen Winston Chan, Qiang Yu, Torsten Wüstefeld, Wee Siong Sho Goh

**Statistical Modeling and Data Analytics**                     
Ying Chen, Nadia M. Davidson, Yuk Kei Wan, Hasindu Gamaarachchi, Andre Sim, Harshil Patel,  Min Hao Ling, Yu Song Chuah, Naruemon Pratanwanich, Christopher Hendra, Laura Watten, Chelsea Sawyer, Dominik Stanojevic, Philip Andrew Ewels, Andreas Wilm, Mile Sikic, Alexandre Thiery, Michael I. Love, Alicia Oshlak, Jonathan Göke
## Citing the SG-NEx project

The SG-NEx resource is described in:

Chen, Ying, et al. "A systematic benchmark of Nanopore long read RNA sequencing for transcript level analysis in human cell lines." _bioRxiv_ (2021). doi: https://doi.org/10.1101/2021.04.21.440736

Please cite this pre-print when using these data, and add the following details: "The SG-NEx data was accessed on [DATE] at registry.opendata.aws/sg-nex-data".

## Contact

Questions about SG-NEx? Please add an entry in the [Discussions Forum](https://github.com/GoekeLab/sg-nex-data/discussions). You can also contact [Jonathan Göke](https://www.a-star.edu.sg/gis/our-people/faculty-staff)

![The Singapore Nanopore-Expression Project\!](https://jglaborg.files.wordpress.com/2020/10/sg_nex_logos-1.png)
