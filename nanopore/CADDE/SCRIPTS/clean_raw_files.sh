find RAW_FILES -name "*.fast5" | parallel mv {} RAW_FILES

find RAW_FILES -type d | parallel rmdir {}
