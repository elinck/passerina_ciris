#!/bin/bash

cd ~/Dropbox/ruhu/Passerina_ciris/stacks/

files=~/Dropbox/Passerina_ciris/pyrad/fastq/*
for f in $files
do
  echo "Processing $f"
  
  bowtie2 -p 8 -x /Users/cj/Desktop/Setophaga_corona_genome/Scoronata -U $f \
  -S ${f}.sam
  
  samtools view ${f}.sam \
  -b -q 10 -o ${f%_R1.fq.gz}.bam
  
  rm ${f}.sam
  
done

