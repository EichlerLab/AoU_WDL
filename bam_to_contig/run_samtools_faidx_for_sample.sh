#!/bin/env bash

regions_file="$1"
haplotig_fasta_file="$2"
out_file="$3"

temp_fasta="temp.fasta"

if [ -e ${regions_file} ]; then
  if file "${haplotig_fasta_file}" | grep -q gzip; then
    gunzip -c "${haplotig_fasta_file}" > "${temp_fasta}"
  else
    cp "${haplotig_fasta_file}" "${temp_fasta}"
  fi
  sed 's/\/rc//g' "${regions_file}" | sed 's/ .*$//' > "regions_file_no_rc.txt"
  sed -i 's/#.*//' ${temp_fasta}
  samtools faidx "${temp_fasta}" -r "regions_file_no_rc.txt" -o ${out_file}

  for rc_ctg in $( grep -o '\S*/rc\S*' "${regions_file}" ); do
      non_rc_ctg=$( echo "$rc_ctg" | sed 's|/rc||g' )
      sed -i "s|$non_rc_ctg|$rc_ctg|g" "${out_file}"
  done

  while IFS=' ' read -r ctg locus; do
     if [[ -n "$locus" ]]; then
       sed -i "s|${ctg}|${locus}_${ctg}|g" "${out_file}" 
     fi
  done < ${regions_file}
  rm "${temp_fasta}"
  rm "${temp_fasta}.fai"
fi
