import glob
import sys
import os

NAME_TEMPLATE = "fastq/{0}_{1}_L00{2}_R{3}_001.fastq"
STAR_TEMPLATE = """STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1 --readFilesIn {0} --outFileNamePrefix star/{1}.star""" 
RSEM_TEMPLATE = """rsem-calculate-expression --paired-end -p 4 {0} ~/nl/genomes/rsem/Naegr1.rsem rsem/{1}.rsem"""


def make_mate_string(prefix):
    run_string = prefix.split("-")[1]
    mate_string = ""
    for R_n in range(1,3):
        files = []
        for L_n in range(1,5):
            files.append(NAME_TEMPLATE.format(prefix, run_string, L_n, R_n))
        mate_string += ",".join(files) + " "
    return mate_string


def get_unique_prefixes():
    filenames = glob.glob(sys.argv[1] + "/*.fastq")
    prefixes = set([os.path.basename(filename).split("_")[0] for filename in filenames])
    return prefixes


def star_output(prefixes):
    print("""module load star/2.7.6a
    module load samtools/1.9
    module load qualimap/2.2.1\n""")

    samtools_template = """samtools view -S -b star/{0}.starAligned.out.sam > star/{0}.star.bam
    samtools sort -o sorted/{0}.star.sorted.bam star/{0}.star.bam
    samtools index sorted/{0}.star.sorted.bam"""

    qualimap_template = """qualimap bamqc -bam sorted/{0}.star.sorted.bam -gff "../../genomes/Naegr1_best_models_gff3.gff" -outdir "./4_qualimap/{0}/" -outfile "{0}.star.sorted" -outformat "HTML" --java-mem-size=4G >& "{0}.star.sorted-qualimap.log" """
    for prefix in prefixes:
        mate_string = make_mate_string(prefix)
        print(STAR_TEMPLATE.format(mate_string, prefix))
        print(samtools_template.format(prefix))
        print(qualimap_template.format(prefix))
    

def rsem_output(prefixes):

    print("""module load RSEM/1.3.3
module unload python3/3.9.6
module load python3/3.8.2
module load bowtie/1.3.0\n""")
    for prefix in prefixes:
        mate_string = make_mate_string(prefix)
        print(RSEM_TEMPLATE.format(mate_string, prefix))

def main():

    prefixes = get_unique_prefixes()
    #star_output(prefix, mate_string)
    rsem_output(prefixes)

if __name__ == "__main__":
    main()

