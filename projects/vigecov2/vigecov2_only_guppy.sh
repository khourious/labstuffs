#!/bin/bash

FAST5=$1

LIBNAME=$2

mkdir $LIBNAME $LIBNAME/reads $LIBNAME/barcodes; cd $LIBNAME

guppy_basecaller -r -i $FAST5 -s reads/ -c dna_r9.4.1_450bps_fast.cfg -x auto; guppy_barcoder -r -i reads/pass/ -s barcodes/ --require_barcodes_both_ends -c configuration.cfg --barcode_kits "EXP-NBD196"; cd barcodes/; list=$(ls -l | grep ^d | awk '$9 ~/^barcode/ {print $9}'); for i in $list; do cat $i"/"*".fastq" > "../"$LIBNAME"_"$i".fastq"; gzip "../"$LIBNAME"_"$i".fastq"; done 
