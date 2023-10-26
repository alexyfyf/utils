#!/bin/bash

## 1.01
## for paired end, and only 2 or 4 lanes
## when used in single end, it should still work, with R2 command reporting error
## 20211108
## To be used following get_files_basespace.sh

mkdir merge

for i in *L001_R1_001.fastq.gz
do
BASE=$(basename $i _L001_R1_001.fastq.gz)
echo $BASE
cat ${BASE}_L00[1-4]_R1_001.fastq.gz > merge/${BASE}_R1.fastq.gz 
cat ${BASE}_L00[1-4]_R2_001.fastq.gz > merge/${BASE}_R2.fastq.gz
done

