# SG-NEx - The Singapore Nanopore-Expression Project
## The Singapore Nanopore-Expression Project

The SG-NEx project was initiated by the Genome Institute of Singapore with the aim to generate reference transcriptomes for 5 of the most commonly used cancer cell lines using Nanopore RNA-Seq data. The following cell lines are used for sequencing:

* A549 (Lung Cancer)    
* MCF7 (Breast Cancer)   
* Hct116 (Colon Cancer)
* K562 (Leukemia)  
* HepG2 (Liver Cancer)     
   

Transcriptome profiling is done using PCR-cDNA sequencing ("PCR-cDNA"), amplification-free cDNA sequencing ("direct cDNA"), and direct sequencing of native RNA (“direct RNA”). 

## Data download
Data can be downloaded [here](DATA.md)     
Notes on data usage: This site provides early access to the SG-NEx data for research. Please note that the data is under publication embargo until the SG-NEx project is published.

## Basecalling and Alignment

Basecalling was done using albacore-2.3.1. Alignment was done with minimap2 using the following parameters:    
 
* Alignment to genome:  -ax splice  -I 32G   -t 3 
* Alignment to transcriptome: -ax  map-ont  -t 3


## Reference files
Details on reference files can be found [here](ANNOTATIONS.md).


## Contributors
**GIS Sequencing Platform**            
Hwee Meng Low, Wendy Soon, CC Khor     
**GIS Genome Innovation Lab**        
Yao Fei, Sarah Ng    
**Cancer Genomics**            
Torsten Wüstefeld, Viktoriia Iakovleva, Ramanuj DasGupta, Shumei Chia, Lixia Xin, Shyam Prabhakar, Puay Leng Lee, Yu Qiang, Wai Leong Tam, Patrick Tan, Sho Goh     
**Statistical Modeling and Data Analytics**                     
Chen Ying, Naruemon Pratanwanich, Andreas Wilm, Alexandre Thiery, Jonathan Göke


## Contact
[Jonathan Göke](https://www.a-star.edu.sg/gis/Our-People/Investigator-details/source/faculty_member/user_id/160)

[https://www.jglab.org](https://www.jglab.org)


