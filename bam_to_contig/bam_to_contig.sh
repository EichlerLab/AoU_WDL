#!/bin/env bash
python $1 $2 $3 | \
sort -k1,1 -k2,2n | \
bedtools merge -d 20 | \
awk -v flank="$4" '{ $2=($2-flank<0)?0:$2-flank; $3=$3+flank } 1' | \
sed 's/ /:/; s/ /-/'
