#!/bin/bash

in_file=some_ids/protein_ids.txt
out_file=metadata/assembly_ids.txt
for protein in $(cat ${in_file})
do
  efetch -db ipg -id ${protein} -format ipg  -email kotovskaya.oa@outlook.com | grep ${protein} | \
  rev |  cut -f 1,3 | rev | \
  grep -v "Assembly" >> ${out_file}
done