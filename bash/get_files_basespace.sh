#!/bin/bash

## version 1.0 
## 20211108

## from Nectar server to access basespace Project folders (containing fastq files)
## To be used with merge_pe.sh

DIR=$1
OLDIFS=$IFS
IFS="
"

for i in ${DIR}/Samples/*/Files/*.fastq.gz; do 
	ln -s $i; 
done

IFS=$OLDIFS

