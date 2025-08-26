#!/bin/env bash
python $1 $2 $3 $4 | \
sort -k1,1 -k2,2n | \
awk 'NR==1{cols=NF} {if(cols>=4) print | "bedtools merge -d 100 -c 4 -o distinct"; else print | "bedtools merge -d 100"}' | \
awk -v flank="0" '{ $2=($2-flank<0)?0:$2-flank; $3=$3+flank } 1' | \
sed 's/ /:/; s/ /-/'