#!/bin/bash
cd ~/Desktop/

refSeqs=tgut2_split/*

for chr in $refSeqs

	do
	echo "aligning to $chr"
	
	nucmer $chr ./Scoronata.fa
	
	show-coords ./out.delta > ~/Desktop/mummer_out/coords/${chr##*/}.coords #Check regex if changing paths.
	
	mv out.delta ~/Desktop/mummer_out/alignments/${chr##*/}.delta
	
done
