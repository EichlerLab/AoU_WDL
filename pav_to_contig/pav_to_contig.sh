#!/bin/env bash
gunzip -c $1 | \
bedtools intersect -header -a stdin -b $2 -wa | \
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/TIG_REGION\t[%GT]\t[%SAMPLE]\n' | \
python $3 | \
sort -k1,1 -k2,2n | \
bedtools merge -c 1,2,3,4,5,6 -o distinct -d 1000000 | \
awk -v flank="$4" '{ $2=$2-flank; $3=$3+flank } 1'
