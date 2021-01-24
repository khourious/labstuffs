# [CADDE](https://www.caddecentre.org/) script
Modified by lpmor22

---

## Setting up and running the pipeline

1. Download and install the pipeline from the github repo:

```sh
sudo apt-get install -y npm
sudo npm install -g github-files-fetcher
fetcher --url="https://github.com/lpmor22/scripts/tree/master/nanopore/CADDE"
cd CADDE
chmod 700 -R INSTALL SCRIPTS
bash INSTALL
```

---

2. Create a csv file in ``CSV_FILES`` directory -- the csv file name corresponds to the library name.

```sh
SSA_FIOCRUZ_LIBRARY01_20210112.csv
```

The csv file contains this format: sample,barcode,virus_reference/version -- NO HEADER!!

---

3. Then run the CADDE script -- it requires 2 parameters:

- The PATH to the directory containing the raw files (fast5 files)

- The PATH to the csv file

```sh
CADDE $HOME/ONT/RAW/SSA_FIOCRUZ_LIBRARY01_20210112 $HOME/ONT/CSV_FILES/SSA_FIOCRUZ_LIBRARY01_20210112.csv
```

The consensus and stats results are in the ``CONSENSUS`` directory of the library.

---

#### This pipeline can:

- Basecall using guppy_basecaller.

- Demultiplex using guppy_barcoder.

- Estimate the min and max read filtering lengths automatically.

- Combine pool A and B if they are on 2 different barcodes, by adding an extra line at the end of the csv file:
```sh
	sample01A,BC01,nCoV-2019/V3
	sample01B,BC02,nCoV-2019/V3
	sample01,BC01-BC02,nCoV-2019/V3
```

- Generate consensus sequences using medaka.

- Do assembly statistics.
