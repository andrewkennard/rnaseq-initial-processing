#!/bin/bash 

module load qualimap/2.2.1

prefix=$( basename "$1" .bam)

unset DISPLAY

qualimap bamqc \
	-bam $1 \
	-gff "../../genomes/Naegr1_best_models_gff3.gff" \
	-outdir "./4_qualimap/" \
	-outfile "$prefix" \
	-outformat "HTML" \
	--java-mem-size=4G >& "${prefix}-qualimap.log"
