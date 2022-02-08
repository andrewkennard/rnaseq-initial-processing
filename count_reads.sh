module load RSEM/1.3.3
module unload python3/3.9.6
module load python3/3.8.2
module load bowtie/1.3.0

rsem-calculate-expression --paired-end -p 4 fastq/G141-S10_S10_L001_R1_001.fastq,fastq/G141-S10_S10_L002_R1_001.fastq,fastq/G141-S10_S10_L003_R1_001.fastq,fastq/G141-S10_S10_L004_R1_001.fastq fastq/G141-S10_S10_L001_R2_001.fastq,fastq/G141-S10_S10_L002_R2_001.fastq,fastq/G141-S10_S10_L003_R2_001.fastq,fastq/G141-S10_S10_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S10.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S11_S11_L001_R1_001.fastq,fastq/G141-S11_S11_L002_R1_001.fastq,fastq/G141-S11_S11_L003_R1_001.fastq,fastq/G141-S11_S11_L004_R1_001.fastq fastq/G141-S11_S11_L001_R2_001.fastq,fastq/G141-S11_S11_L002_R2_001.fastq,fastq/G141-S11_S11_L003_R2_001.fastq,fastq/G141-S11_S11_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S11.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S4_S4_L001_R1_001.fastq,fastq/G141-S4_S4_L002_R1_001.fastq,fastq/G141-S4_S4_L003_R1_001.fastq,fastq/G141-S4_S4_L004_R1_001.fastq fastq/G141-S4_S4_L001_R2_001.fastq,fastq/G141-S4_S4_L002_R2_001.fastq,fastq/G141-S4_S4_L003_R2_001.fastq,fastq/G141-S4_S4_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S4.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S5_S5_L001_R1_001.fastq,fastq/G141-S5_S5_L002_R1_001.fastq,fastq/G141-S5_S5_L003_R1_001.fastq,fastq/G141-S5_S5_L004_R1_001.fastq fastq/G141-S5_S5_L001_R2_001.fastq,fastq/G141-S5_S5_L002_R2_001.fastq,fastq/G141-S5_S5_L003_R2_001.fastq,fastq/G141-S5_S5_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S5.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S6_S6_L001_R1_001.fastq,fastq/G141-S6_S6_L002_R1_001.fastq,fastq/G141-S6_S6_L003_R1_001.fastq,fastq/G141-S6_S6_L004_R1_001.fastq fastq/G141-S6_S6_L001_R2_001.fastq,fastq/G141-S6_S6_L002_R2_001.fastq,fastq/G141-S6_S6_L003_R2_001.fastq,fastq/G141-S6_S6_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S6.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S7_S7_L001_R1_001.fastq,fastq/G141-S7_S7_L002_R1_001.fastq,fastq/G141-S7_S7_L003_R1_001.fastq,fastq/G141-S7_S7_L004_R1_001.fastq fastq/G141-S7_S7_L001_R2_001.fastq,fastq/G141-S7_S7_L002_R2_001.fastq,fastq/G141-S7_S7_L003_R2_001.fastq,fastq/G141-S7_S7_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S7.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S1_S1_L001_R1_001.fastq,fastq/G141-S1_S1_L002_R1_001.fastq,fastq/G141-S1_S1_L003_R1_001.fastq,fastq/G141-S1_S1_L004_R1_001.fastq fastq/G141-S1_S1_L001_R2_001.fastq,fastq/G141-S1_S1_L002_R2_001.fastq,fastq/G141-S1_S1_L003_R2_001.fastq,fastq/G141-S1_S1_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S1.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S2_S2_L001_R1_001.fastq,fastq/G141-S2_S2_L002_R1_001.fastq,fastq/G141-S2_S2_L003_R1_001.fastq,fastq/G141-S2_S2_L004_R1_001.fastq fastq/G141-S2_S2_L001_R2_001.fastq,fastq/G141-S2_S2_L002_R2_001.fastq,fastq/G141-S2_S2_L003_R2_001.fastq,fastq/G141-S2_S2_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S2.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S3_S3_L001_R1_001.fastq,fastq/G141-S3_S3_L002_R1_001.fastq,fastq/G141-S3_S3_L003_R1_001.fastq,fastq/G141-S3_S3_L004_R1_001.fastq fastq/G141-S3_S3_L001_R2_001.fastq,fastq/G141-S3_S3_L002_R2_001.fastq,fastq/G141-S3_S3_L003_R2_001.fastq,fastq/G141-S3_S3_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S3.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S8_S8_L001_R1_001.fastq,fastq/G141-S8_S8_L002_R1_001.fastq,fastq/G141-S8_S8_L003_R1_001.fastq,fastq/G141-S8_S8_L004_R1_001.fastq fastq/G141-S8_S8_L001_R2_001.fastq,fastq/G141-S8_S8_L002_R2_001.fastq,fastq/G141-S8_S8_L003_R2_001.fastq,fastq/G141-S8_S8_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S8.rsem
rsem-calculate-expression --paired-end -p 4 fastq/G141-S9_S9_L001_R1_001.fastq,fastq/G141-S9_S9_L002_R1_001.fastq,fastq/G141-S9_S9_L003_R1_001.fastq,fastq/G141-S9_S9_L004_R1_001.fastq fastq/G141-S9_S9_L001_R2_001.fastq,fastq/G141-S9_S9_L002_R2_001.fastq,fastq/G141-S9_S9_L003_R2_001.fastq,fastq/G141-S9_S9_L004_R2_001.fastq  ~/nl/genomes/rsem/Naegr1.rsem rsem/G141-S9.rsem