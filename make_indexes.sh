module load star/2.7.6a

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ~/nl/genomes/Naegr1 --genomeFastaFiles ~/nl/genomes/Naegr1_scaffolds.fasta --sjdbGTFfile ~/nl/genomes/Naegr1_best_models_gff3.gff --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 11
