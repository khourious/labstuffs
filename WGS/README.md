## Bioinformatic pipeline for whole genome sequencing using Illumina and MinION

This repo contains scripts and files to run the bioinformatic analysis of whole genome sequencing using MinION and Illlumina platforms, and was built based on the [CADDE](https://www.caddecentre.org/) scripts and [ARTIC](https://artic.network/) bioinformatics workflow

---

### Setting up the pipeline

Download and install the pipeline from the github repo:
```sh
sudo apt-get install -y npm
sudo npm install -g github-files-fetcher
fetcher --url="https://github.com/lpmor22/scripts/tree/master/WGS"
cd WGS
chmod 777 -R INSTALL SCRIPTS
bash INSTALL
```
---

### Illumina pipeline

#### ILLUMINA script for genome assembly -- it requires 2 parameters:

- The prime scheme information (example: nCoV-2019/FIOCRUZ_2kb_v1 or nCoV-2019/ARTIC_V3)
- The PATH to the directory containing the raw sequencing data downloaded from Illumina's BaseSpace Sequence Hub (fastq.gz files).

Then run the script using:
```sh
ILLUMINA nCoV-2019/FIOCRUZ_2kb_v1 $HOME/WGS/RAW/LIBRARY01_20210123 
```

---

### MinION pipeline

#### Create the sample sheet

Create a csv file in ``CSV_FILES`` directory -- the csv file name **corresponds to the library name**.
	
The csv file contains this format: sample,barcode,virus_reference/version -- NO HEADER!!__
	
You can combine pool A and B if they are on 2 different barcodes, by adding an extra line at the end of the csv file:
```sh
sample01A,BC01,nCoV-2019/ARTIC_V3
sample01B,BC02,nCoV-2019/ARTIC_V3
sample01,BC01-BC02,nCoV-2019/ARTIC_V3
```
	
#### RAMPART script for MinION real-time analysis -- it requires 2 parameters:

- The PATH to the csv file.
- The PATH to the directory containing the fastq pass files obtained from MinKNOW fast basecalling during the sequencing.

Then run the script using:

```sh
RAMPART $HOME/WGS/CSV_FILES/LIBRARY01_20210123.csv $HOME/WGS/RAW/LIBRARY01_20210123/../fastq_pass
```

#### MINION script for genome assembly -- it requires 2 parameters:

- The PATH to the csv file.
- The PATH to the directory containing the raw sequencing data (fast5 files).
Then run the script using:
```sh
MINION $HOME/WGS/CSV_FILES/LIBRARY01_20210123.csv $HOME/WGS/RAW/LIBRARY01_20210123 
```