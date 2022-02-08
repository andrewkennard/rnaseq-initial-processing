module load star/2.7.6a
module load samtools/1.9

STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S1_S1_L001_R1_001.fastq,fastq/G141-S1_S1_L002_R1_001.fastq,fastq/G141-S1_L003_R1_001.fastq,fastq/G141-S1_S1_L004_R1_001.fastq fastq/G141-S1_S1_L001_R2_001.fastq,fastq/G141-S1_S1_L002_R2_001.fastq,fastq/G141-S1_S1_L003_R2_001.fastq,fastq/G141-S1_S1_L004_R2_001.fastq --outFileNamePrefix star/G141-S1.star

samtools view -S -b star/G141-S1.starAligned.out.sam > star/G141-S1.star.bam

mkdir -p sorted
samtools sort -o sorted/G141-S1.star.sorted.bam star/G141-S1.star.bam
samtools index sorted/G141-S1.star.sorted.bam
