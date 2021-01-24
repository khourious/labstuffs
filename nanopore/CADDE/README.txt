
###################
###### SETUP ######
###################

## I assume the directories LIBRARIES, CSV_FILES, SCRIPTS and PRIMER_SCHEMES 
## are already created on the Desktop from the first CADDE pipleine setup
## Copy all the ".sh" files to the SCRIPTS directory
## then type in a terminal

cd ~/Desktop/SCRIPTS

git clone --recursive https://github.com/artic-network/artic-ncov2019.git

conda env create -f artic-ncov2019/environment.yml
conda install -c bioconda exonerate -n artic-ncov2019
conda env create -f artic-ncov2019/environment-medaka.yml
conda install -c bioconda exonerate -n artic-ncov2019-medaka

cp -r ~/Desktop/SCRIPTS/artic-ncov2019/primer_schemes/nCoV-2019 ~/Desktop/PRIMER_SCHEMES

cd ~/Desktop

###################
###################



###################
###### USAGE ######
###################

## Like for the first CADDE pipieline, create a csv file in CSV_FILES directory 
## csv file name corresponds to the library name with this format NameVirus(es)_LocationSequencing_libraryNumber_Date   
## example: nCov19_SP_library4_20200305.csv
## csv file contains this format: sample,barcode,virus_reference/version (NO HEADER !!)	 
## example: sample01,BC01,ncov2019/V2
##			sample02,BC02,ncov2019/V2
##			...

## then run the script CADDE_pipeline.NEW.sh
## It requires 3 parameters
## 1. the csv file
## 2. the path to the directory containing the fastq files
## 3. the path to the directory containing the fast5 files

## example: ~/Desktop/SCRIPTS/CADDE_pipeline.NEW.sh ~/Desktop/CSV_FILES/nCov19_SP_library4_20200305.csv ~/Desktop/LIBRARIES/nCov19_SP_library4_20200305/XXX/fastq_pass ~/Desktop/LIBRARIES/nCov19_SP_library4_20200305/XXX/fast5_pass

## The consensus and stats results are in the CONSENSUS directory of the library

 
## The pipeline can:
## 	- demultiplex using guppy_barcoder, give good results like porechop but is much faster
## 	- estimate the min and max read filtering lengths automatically
## 	- combine pool A and B if they are on 2 different barcodes, by adding an extra line at the end of the csv file, for example:
## 		sample01A,BC01,ncov2019/V2
## 		sample01B,BC02,ncov2019/V2
## 		...
## 		sample01,BC01-BC02,ncov2019/V2
## 	- generate consensus sequences using nanopolish (It can also generate consensus using medaka, you just have to remove the ## at the beginning of the lines starting with:  "source activate" and "artic minion" and add ## to the lines just above starting with the same "source activate" and "artic minion")
## 	- do assembly statistics





