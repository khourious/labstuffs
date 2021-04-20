# WGS script
Based on the [CADDE](https://www.caddecentre.org/) scripts and [ARTIC](https://artic.network/) bioinformatics workflow

---

## Setting up the pipeline

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

## MinION pipeline

### Create the sample sheet

1. Create a csv file in ``CSV_FILES`` directory -- the csv file name corresponds to the library name.
The csv file contains this format: sample,barcode,virus_reference/version -- NO HEADER!!

```sh
LIBRARY01_nCOV-19_20210123.csv
```

### RAMPART script for MinION real-time analysis -- it requires 2 parameters:

- The PATH to the csv file.

- The PATH to the directory containing the fastq pass files obtained from MinKNOW fast basecalling during the sequencing.

Then run the script using:

```sh
RAMPART $HOME/WGS/CSV_FILES/LIBRARY01_nCOV-19_20210123.csv $HOME/WGS/RAW/LIBRARY01_nCOV-19_20210123/../fastq_pass
```

### MINION script for genome assembly -- it requires 2 parameters:

- The PATH to the csv file

- The PATH to the directory containing the raw files (fast5 files)

Then run the script using:

```sh
ONT $HOME/WGS/CSV_FILES/LIBRARY01_nCOV-19_20210123.csv $HOME/WGS/RAW/LIBRARY01_nCOV-19_20210123 
```

The consensus and stats results are in the ``CONSENSUS`` directory of the library.

#### This pipeline can:

- Perform real-time analysis using RAMPART.
- Perform high accuracy basecalling using guppy_basecaller.
- Do demultiplexing using guppy_barcoder.
- Estimate the min and max read filtering lengths automatically.

- Combine pool A and B if they are on 2 different barcodes, by adding an extra line at the end of the csv file:
```sh
	sample01A,BC01,nCoV-2019/ARTIC_V3
	sample01B,BC02,nCoV-2019/ARTIC_V3
	sample01,BC01-BC02,nCoV-2019/ARTIC_V3
```

- Generate consensus sequences using medaka.
- Do assembly statistics.