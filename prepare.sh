#!/bin/bash

# Check for correct usage
if [ "$#" -ne 1 ] || [ ! -d "$1" ]; then
    printf "Usage: %s data_dir\n" "$0"
    exit 1
fi

# Create necessary directories
directories=("logs" "results" "rawdata" "benchmarks" "tmp")
for dir in "${directories[@]}"; do
    if [ ! -d "$dir" ]; then
        mkdir "$dir" || { printf "Error: Could not create directory '%s'.\n" "$dir"; exit 1; }
    fi
done

# Change to rawdata directory
cd rawdata || { printf "Error: Could not change to 'rawdata' directory.\n"; exit 1; }

# Use find to create symbolic links to .fq.gz files in the input directory and its subdirectories
while IFS= read -r fq_file; do
    ln -s "$fq_file" . || { printf "Error: Could not create symbolic link for '%s' in 'rawdata' directory.\n" "$fq_file"; exit 1; }
done < <(find "$1" -type f -name "*.fq.gz")

# Check if at least one .fq.gz file was found and linked
fq_files=(./*.fq.gz)
if [ ${#fq_files[@]} -eq 0 ]; then
    printf "Error: No .fq.gz files found in directory '%s'.\n" "$1"
    exit 1
fi

ls *.fq.gz >> ../samples
