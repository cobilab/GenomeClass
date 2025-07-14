#!/bin/bash
#
declare -a DATASETS_AVAILABLE=( "Archaea" "Viruses" "Unclassified" "Eukaryota" "2"); #"Other"
NUMBER_DS=5000;
#
wget https://www.genome.jp/ftp/db/virushostdb/virushostdb.formatted.genomic.fna.gz
gunzip virushostdb.formatted.genomic.fna.gz
#
for DATASET in "${DATASETS_AVAILABLE[@]}"
  do
  
  datasets summary genome taxon "${DATASET}" --as-json-lines  > ${DATASET}_summary.jsonl

  shuf -n $NUMBER_DS ${DATASET}_summary.jsonl > new_${DATASET}_summary.jsonl
  mv new_${DATASET}_summary.jsonl ${DATASET}_summary.jsonl
  
  
  jq -r '.accession' ${DATASET}_summary.jsonl > ${DATASET}_accessions.txt
  datasets download genome accession --inputfile ${DATASET}_accessions.txt --dehydrated --filename ${DATASET}_complete.zip

  unzip ${DATASET}_complete.zip -d ${DATASET}_complete
  datasets rehydrate --directory ${DATASET}_complete
    
  ./Label_datasets.sh $DATASET
done
