#!/bin/env bash

regions_file="$1"
haplotig_fasta_file="$2"
out_file="$3"

if [ -e ${regions_file} ]; then
  gunzip -c "${haplotig_fasta_file}" > "temp.fa"
  samtools faidx "temp.fa" -r "${regions_file}" -o ${out_file}
  rm "temp.fa"
fi
