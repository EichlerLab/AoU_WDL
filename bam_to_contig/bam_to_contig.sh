#!/bin/env bash
python $1 $2 $3 $4 $5 | \
sort -k1,1 -k2,2n | \
bedtools merge -d 1000000 | \
awk -v flank="0" '{ $2=($2-flank<0)?0:$2-flank; $3=$3+flank } 1' | \
sed 's/ /:/; s/ /-/'
