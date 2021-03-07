#!/bin/bash

#basecalling data from .fast5 files with guppy

guppy_basecaller --input_path /path/to/fast5 --save_path /path/for/output --flowcell FLO-MIN106 --kit SQK-LSK109 --qscore_filtering --min_qscore 8

#demultiplexing, trimming, and filtering reads with guppy_barcoder

guppy_barcoder --input_path /path/to/fastq --save_path /path/for/output --require_barcode_both_ends --trim_barcodes --num_extra_bases_trim 18 --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"