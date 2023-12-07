#!/bin/env bash

regions_file="$1"
haplotig_fasta_file="$2"
out_file="$3"

if [ -e ${regions_file} ]; then
  gunzip -c "${haplotig_fasta_file}" > "temp.fa"
  sed 's/\/rc//g' "${regions_file}" > "regions_file_no_rc.txt"
  samtools faidx "temp.fa" -r "regions_file_no_rc.txt" -o ${out_file}
  for rc_ctg in $( grep -o '\S*/rc\S*' "${regions_file}" ); do
      non_rc_ctg=$( echo "$rc_ctg" | sed 's|/rc||g' )
      sed -i "s|$non_rc_ctg|$rc_ctg|g" "${out_file}"
  done
  rm temp.fa
fi
