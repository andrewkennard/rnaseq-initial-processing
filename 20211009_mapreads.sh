module load star/2.7.6a
module load samtools/1.9
module load qualimap/2.2.1

#STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S6_S6_L001_R1_001.fastq,fastq/G141-S6_S6_L002_R1_001.fastq,fastq/G141-S6_S6_L003_R1_001.fastq,fastq/G141-S6_S6_L004_R1_001.fastq fastq/G141-S6_S6_L001_R2_001.fastq,fastq/G141-S6_S6_L002_R2_001.fastq,fastq/G141-S6_S6_L003_R2_001.fastq,fastq/G141-S6_S6_L004_R2_001.fastq  --outFileNamePrefix star/G141-S6.star
#samtools view -S -b star/G141-S6.starAligned.out.sam > star/G141-S6.star.bam
#samtools sort -o sorted/G141-S6.star.sorted.bam star/G141-S6.star.bam
#samtools index sorted/G141-S6.star.sorted.bam
qualimap bamqc -bam sorted/G141-S6.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S6/" -outfile "G141-S6.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S6.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S11_S11_L001_R1_001.fastq,fastq/G141-S11_S11_L002_R1_001.fastq,fastq/G141-S11_S11_L003_R1_001.fastq,fastq/G141-S11_S11_L004_R1_001.fastq fastq/G141-S11_S11_L001_R2_001.fastq,fastq/G141-S11_S11_L002_R2_001.fastq,fastq/G141-S11_S11_L003_R2_001.fastq,fastq/G141-S11_S11_L004_R2_001.fastq  --outFileNamePrefix star/G141-S11.star
samtools view -S -b star/G141-S11.starAligned.out.sam > star/G141-S11.star.bam
samtools sort -o sorted/G141-S11.star.sorted.bam star/G141-S11.star.bam
samtools index sorted/G141-S11.star.sorted.bam
qualimap bamqc -bam sorted/G141-S11.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S11/" -outfile "G141-S11.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S11.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S2_S2_L001_R1_001.fastq,fastq/G141-S2_S2_L002_R1_001.fastq,fastq/G141-S2_S2_L003_R1_001.fastq,fastq/G141-S2_S2_L004_R1_001.fastq fastq/G141-S2_S2_L001_R2_001.fastq,fastq/G141-S2_S2_L002_R2_001.fastq,fastq/G141-S2_S2_L003_R2_001.fastq,fastq/G141-S2_S2_L004_R2_001.fastq  --outFileNamePrefix star/G141-S2.star
samtools view -S -b star/G141-S2.starAligned.out.sam > star/G141-S2.star.bam
samtools sort -o sorted/G141-S2.star.sorted.bam star/G141-S2.star.bam
samtools index sorted/G141-S2.star.sorted.bam
qualimap bamqc -bam sorted/G141-S2.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S2/" -outfile "G141-S2.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S2.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S10_S10_L001_R1_001.fastq,fastq/G141-S10_S10_L002_R1_001.fastq,fastq/G141-S10_S10_L003_R1_001.fastq,fastq/G141-S10_S10_L004_R1_001.fastq fastq/G141-S10_S10_L001_R2_001.fastq,fastq/G141-S10_S10_L002_R2_001.fastq,fastq/G141-S10_S10_L003_R2_001.fastq,fastq/G141-S10_S10_L004_R2_001.fastq  --outFileNamePrefix star/G141-S10.star
samtools view -S -b star/G141-S10.starAligned.out.sam > star/G141-S10.star.bam
samtools sort -o sorted/G141-S10.star.sorted.bam star/G141-S10.star.bam
samtools index sorted/G141-S10.star.sorted.bam
qualimap bamqc -bam sorted/G141-S10.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S10/" -outfile "G141-S10.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S10.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S5_S5_L001_R1_001.fastq,fastq/G141-S5_S5_L002_R1_001.fastq,fastq/G141-S5_S5_L003_R1_001.fastq,fastq/G141-S5_S5_L004_R1_001.fastq fastq/G141-S5_S5_L001_R2_001.fastq,fastq/G141-S5_S5_L002_R2_001.fastq,fastq/G141-S5_S5_L003_R2_001.fastq,fastq/G141-S5_S5_L004_R2_001.fastq  --outFileNamePrefix star/G141-S5.star
samtools view -S -b star/G141-S5.starAligned.out.sam > star/G141-S5.star.bam
samtools sort -o sorted/G141-S5.star.sorted.bam star/G141-S5.star.bam
samtools index sorted/G141-S5.star.sorted.bam
qualimap bamqc -bam sorted/G141-S5.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S5/" -outfile "G141-S5.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S5.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S3_S3_L001_R1_001.fastq,fastq/G141-S3_S3_L002_R1_001.fastq,fastq/G141-S3_S3_L003_R1_001.fastq,fastq/G141-S3_S3_L004_R1_001.fastq fastq/G141-S3_S3_L001_R2_001.fastq,fastq/G141-S3_S3_L002_R2_001.fastq,fastq/G141-S3_S3_L003_R2_001.fastq,fastq/G141-S3_S3_L004_R2_001.fastq  --outFileNamePrefix star/G141-S3.star
samtools view -S -b star/G141-S3.starAligned.out.sam > star/G141-S3.star.bam
samtools sort -o sorted/G141-S3.star.sorted.bam star/G141-S3.star.bam
samtools index sorted/G141-S3.star.sorted.bam
qualimap bamqc -bam sorted/G141-S3.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S3/" -outfile "G141-S3.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S3.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S8_S8_L001_R1_001.fastq,fastq/G141-S8_S8_L002_R1_001.fastq,fastq/G141-S8_S8_L003_R1_001.fastq,fastq/G141-S8_S8_L004_R1_001.fastq fastq/G141-S8_S8_L001_R2_001.fastq,fastq/G141-S8_S8_L002_R2_001.fastq,fastq/G141-S8_S8_L003_R2_001.fastq,fastq/G141-S8_S8_L004_R2_001.fastq  --outFileNamePrefix star/G141-S8.star
samtools view -S -b star/G141-S8.starAligned.out.sam > star/G141-S8.star.bam
samtools sort -o sorted/G141-S8.star.sorted.bam star/G141-S8.star.bam
samtools index sorted/G141-S8.star.sorted.bam
qualimap bamqc -bam sorted/G141-S8.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S8/" -outfile "G141-S8.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S8.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S9_S9_L001_R1_001.fastq,fastq/G141-S9_S9_L002_R1_001.fastq,fastq/G141-S9_S9_L003_R1_001.fastq,fastq/G141-S9_S9_L004_R1_001.fastq fastq/G141-S9_S9_L001_R2_001.fastq,fastq/G141-S9_S9_L002_R2_001.fastq,fastq/G141-S9_S9_L003_R2_001.fastq,fastq/G141-S9_S9_L004_R2_001.fastq  --outFileNamePrefix star/G141-S9.star
samtools view -S -b star/G141-S9.starAligned.out.sam > star/G141-S9.star.bam
samtools sort -o sorted/G141-S9.star.sorted.bam star/G141-S9.star.bam
samtools index sorted/G141-S9.star.sorted.bam
qualimap bamqc -bam sorted/G141-S9.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S9/" -outfile "G141-S9.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S9.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S7_S7_L001_R1_001.fastq,fastq/G141-S7_S7_L002_R1_001.fastq,fastq/G141-S7_S7_L003_R1_001.fastq,fastq/G141-S7_S7_L004_R1_001.fastq fastq/G141-S7_S7_L001_R2_001.fastq,fastq/G141-S7_S7_L002_R2_001.fastq,fastq/G141-S7_S7_L003_R2_001.fastq,fastq/G141-S7_S7_L004_R2_001.fastq  --outFileNamePrefix star/G141-S7.star
samtools view -S -b star/G141-S7.starAligned.out.sam > star/G141-S7.star.bam
samtools sort -o sorted/G141-S7.star.sorted.bam star/G141-S7.star.bam
samtools index sorted/G141-S7.star.sorted.bam
qualimap bamqc -bam sorted/G141-S7.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S7/" -outfile "G141-S7.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S7.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S4_S4_L001_R1_001.fastq,fastq/G141-S4_S4_L002_R1_001.fastq,fastq/G141-S4_S4_L003_R1_001.fastq,fastq/G141-S4_S4_L004_R1_001.fastq fastq/G141-S4_S4_L001_R2_001.fastq,fastq/G141-S4_S4_L002_R2_001.fastq,fastq/G141-S4_S4_L003_R2_001.fastq,fastq/G141-S4_S4_L004_R2_001.fastq  --outFileNamePrefix star/G141-S4.star
samtools view -S -b star/G141-S4.starAligned.out.sam > star/G141-S4.star.bam
samtools sort -o sorted/G141-S4.star.sorted.bam star/G141-S4.star.bam
samtools index sorted/G141-S4.star.sorted.bam
qualimap bamqc -bam sorted/G141-S4.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S4/" -outfile "G141-S4.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S4.star.sorted-qualimap.log" 
STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn fastq/G141-S1_S1_L001_R1_001.fastq,fastq/G141-S1_S1_L002_R1_001.fastq,fastq/G141-S1_S1_L003_R1_001.fastq,fastq/G141-S1_S1_L004_R1_001.fastq fastq/G141-S1_S1_L001_R2_001.fastq,fastq/G141-S1_S1_L002_R2_001.fastq,fastq/G141-S1_S1_L003_R2_001.fastq,fastq/G141-S1_S1_L004_R2_001.fastq  --outFileNamePrefix star/G141-S1.star
samtools view -S -b star/G141-S1.starAligned.out.sam > star/G141-S1.star.bam
samtools sort -o sorted/G141-S1.star.sorted.bam star/G141-S1.star.bam
samtools index sorted/G141-S1.star.sorted.bam
qualimap bamqc -bam sorted/G141-S1.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/G141-S1/" -outfile "G141-S1.star.sorted" -outformat "HTML" --java-mem-size=4G >& "G141-S1.star.sorted-qualimap.log" 
