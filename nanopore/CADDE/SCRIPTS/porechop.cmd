#porechop --verbosity 1 -t 2 -i "$1" -b "$1"_DIR --barcode_threshold 75 --check_reads `grep ">" "$1" | wc -l` --barcode_diff 5 --require_two_barcodes > "$1".demultiplexing.txt

porechop --verbosity 1 -t 2 -i "$1" -b "$1"_DIR --barcode_threshold 80 --native_barcodes --discard_middle --discard_unassigned --barcode_diff 5 --require_two_barcodes > "$1".demultiplexing.txt
