#!/bin/bash
#
DATASET=$1
#  
printf "" > $DATASET.fasta
# 
if [ -d "${DATASET}_complete/ncbi_dataset/data" ]; then
  echo "Directory exists: ${DATASET}_complete/ncbi_dataset/data"
  
  # Loop over each subdirectory in the current directory
  for dir in "${DATASET}_complete/ncbi_dataset/data/"*/; do
    echo "Directory: $dir"
 
    # Loop over each file in that subdirectory
    for file in "$dir"*; do
      if [[ -f "$file" ]]; then
        echo "  File: $file"
        CODE=$(head -1 $file | cut -d' ' -f1 | tr -d ">" )
       
        echo "$CODE"
        HEADER=$(python3 Label_datasets.py $CODE)
        
        printf "$HEADER\n" >> $DATASET.fasta
        tail -n +2 $file >> $DATASET.fasta
         
        fi
      done
    done
  else
  echo "Directory does NOT exist: ${DATASET}_complete/ncbi_dataset/data"
fi
#
