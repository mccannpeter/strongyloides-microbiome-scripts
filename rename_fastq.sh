#!/bin/bash

# Path to SampleSheet.csv
SAMPLESHEET="SampleSheet.csv"

# Declare an associative array to store mapping
declare -A sample_map

# Convert the SampleSheet.csv to Unix-style line endings to avoid '\r' issues
sed -i 's/\r//' "$SAMPLESHEET"

# Read SampleSheet.csv to populate the sample_map
while IFS=',' read -r sample_id i7_index_id index i5_index_id index2 sample_project description; do
    # Skip header or empty lines
    if [[ "$sample_id" == "Sample_ID" || -z "$sample_id" ]]; then
        continue
    fi
    # Populate associative array
    sample_map["$sample_id"]="$description"
done < "$SAMPLESHEET"

# Loop through all fastq.gz files in the current directory
for file in *.fastq.gz; do
    # Extract the Sample_ID from the filename
    sample_id=$(echo "$file" | grep -oE 'PN0545_[0-9]{4}')
    
    # If the Sample_ID is found in the mapping, rename the file
    if [[ -n "${sample_map[$sample_id]}" ]]; then
        # Extract other parts of the filename
        rest=$(echo "$file" | sed -E "s/^${sample_id}_//")
        
        # New filename
        new_name="${sample_map[$sample_id]}_${rest}"
        
        # Rename the file
        mv "$file" "$new_name"
        echo "Renamed $file to $new_name"
    else
        echo "No mapping found for $file"
    fi
done
