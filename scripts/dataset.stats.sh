#!/bin/bash

in_file=$1
source_data=$2
assembly_info=$3
out_file=$4

echo "Total number of proteins blasted" > ${out_file}
cat ${source_data} | cut -f 2 -d "," | sort -u | wc -l >> ${out_file}

echo "Total number of genomes found" > ${out_file}
cat ${assembly_info} | cut -f 1 -d "," |  wc -l >> ${out_file}

echo "Total number of genomes downloaded" >> ${out_file}
cat ${in_file} | wc -l >> ${out_file}

echo "Number of unique phage genomes" >> ${out_file}
cat ${in_file} | cut -f 1 | sort -u | wc -l >> ${out_file}