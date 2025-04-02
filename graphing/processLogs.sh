#!/bin/bash

for file in /home/ket581/skyline/SUMMER/skyline-SRFBD-validation/logs/*; do
    # Extract the part after the last underscore
    new_name=$(echo "$file" | sed -E 's/slurm-[0-9]+_([0-9]+)\.(err|out)/\1.\2/')
    
    # Rename the file if the new name is different
    if [[ "$file" != "$new_name" ]]; then
        mv "$file" "$new_name"
    fi
done

for file in /home/ket581/skyline/SUMMER/skyline-SRFBD-validation/logs/*.out; do
    # Extract the index (filename without extension)
    index=$(basename "$file" .out)
    
    # Define the source and destination paths
    src="/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/inf/$index/sRanges.42.log"
    dest="/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/logs/$index.log"
    
    # Check if the source file exists before copying
    if [[ -f "$src" ]]; then
        cp "$src" "$dest"
        echo "Copied $src to $dest"
    else
        echo "Warning: $src not found!"
    fi
done
