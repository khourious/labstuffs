#!/bin/bash

LIBNAME=$1

cd $LIBNAME/barcodes

gdp.py -c /home/fiocruz-luiz/Documents/GenomeDetective/watcher3/config/config.yaml upload-files -u /home/fiocruz-luiz/Desktop/Resultado_sequenciamento/$LIBNAME/*.fastq.gz
