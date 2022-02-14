import glob
import sys
import os
import stat

# Globals
INITIAL_FOLDER_LIST = ["run_commands", "star", "rsem", "qualimap", "sorted"]
STAR_TEMPLATE = """# run STAR
    STAR --runThreadN 4 --genomeDir ~/nl/genomes/Naegr1/star --readFilesCommand zcat --readFilesIn {mate_string} --outFileNamePrefix {out_dir}/star/{prefix}.star""" 
BSUB_TEMPLATE = """bsub -n 4 -R "span[hosts=1] rusage[mem=4096]" -W 240 -q short -o "queue_logfile.txt" -e "queue_errors.txt" """
RSEM_TEMPLATE = """# run RSEM
    rsem-calculate-expression --star --star-gzipped-read-file --no-bam-output --paired-end -p 4 {mate_string} ~/nl/genomes/Naegr1/rsem/Naegr1.rsem {out_dir}/rsem/{prefix}.rsem"""
SAMTOOLS_TEMPLATE = """# run samtools
    samtools view -S -b {out_dir}/star/{prefix}.starAligned.out.sam > {out_dir}/star/{prefix}.star.bam  
    samtools sort -o {out_dir}/sorted/{prefix}.star.sorted.bam {out_dir}/star/{prefix}.star.bam 
    samtools index {out_dir}/sorted/{prefix}.star.sorted.bam"""
QUALIMAP_TEMPLATE = """# run qualimap
    qualimap bamqc -bam {out_dir}/sorted/{prefix}.star.sorted.bam -gff "~/nl/genomes/Naegr1/Naegr1_best_models_gff3.gff" -outdir "{out_dir}/qualimap/{prefix}/" -outfile "{prefix}.star.sorted" -outformat "HTML" --java-mem-size=4G >& "{out_dir}/{prefix}.star.sorted-qualimap.log" """

PIPELINE_NAME = "process_{0}.sh"

OUT_DIR = os.path.normpath(sys.argv[1])
#Make sure fastq files are always stored in a subfolder named fastq
FASTQ_DIR = os.path.join(OUT_DIR, "fastq")


def make_mate_string(prefix):
    read1_list = []
    read2_list = []
    dirname_list = sorted(glob.glob(os.path.join(FASTQ_DIR, prefix + "_*/")))
    for dirname in dirname_list:
        read_names = sorted(glob.glob(os.path.join(dirname, "*.fastq.gz")))
        for k,read_list in enumerate([read1_list, read2_list]):
                                      read_list.append(read_names[k])
    mate_string = ",".join(read1_list) + " " + ",".join(read2_list)
    return mate_string


def get_unique_prefixes():
    filenames = glob.glob(os.path.join(FASTQ_DIR,  "*/"))
    prefixes = set([os.path.basename(os.path.dirname(filename)).split("_")[0] 
                        for filename in filenames])
    return prefixes


def make_folder(folder_name):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


def initialize_folders():
    for name in INITIAL_FOLDER_LIST:
        folder_name = os.path.join(OUT_DIR, name)
        make_folder(folder_name)


def generate_pipeline_script(prefix):
    mate_string = make_mate_string(prefix)
    star_string = STAR_TEMPLATE.format(mate_string=mate_string,
                                       prefix=prefix,
                                       out_dir=OUT_DIR)
    samtools_string = SAMTOOLS_TEMPLATE.format(prefix=prefix,
                                               out_dir=OUT_DIR)
    qualimap_string = QUALIMAP_TEMPLATE.format(prefix=prefix,
                                              out_dir=OUT_DIR)
    rsem_string = RSEM_TEMPLATE.format(mate_string=mate_string,
                                       prefix=prefix,
                                       out_dir=OUT_DIR)
    full_pipeline = "{0}\n\n{1}\n\n{2}\n\n{3}".format(star_string,
                                                  samtools_string,
                                                  qualimap_string)
                                                  rsem_string)
    return full_pipeline


def save_executable_pipeline(prefix, pipeline):
    pipeline_name = PIPELINE_NAME.format(prefix)
    pipeline_path = os.path.join(OUT_DIR, "run_commands", pipeline_name)
    with open(pipeline_path, "w") as f:
            f.write(pipeline)
    #Add executable permissions for user, group, others
    make_file_executable(pipeline_path)


def make_file_executable(file_name):
    #Add executable permissions for user, group, and others (like chmod +x)
    st = os.stat(file_name)
    os.chmod(file_name, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def save_executable_queue_request(prefixes):
    queue_request_path = os.path.join(OUT_DIR, "run_commands", "request_processing.sh")
    with open(queue_request_path, "w") as f:
        for prefix in prefixes:
            pipeline_path = os.path.join(OUT_DIR, "run_commands", PIPELINE_NAME.format(prefix))
            request = "{0} {1}\n".format(BSUB_TEMPLATE, pipeline_path)
            f.write(request)
    make_file_executable(queue_request_path)


def main():
    initialize_folders()
    prefixes = get_unique_prefixes()
    for prefix in prefixes:
        #Write a script for each prefix/dataset and make it executable
        full_pipeline = generate_pipeline_script(prefix)
        save_executable_pipeline(prefix, full_pipeline)
        print(full_pipeline)
    #Write a script to request a job for each dataset and make it executable
    save_executable_queue_request(prefixes)



if __name__ == "__main__":
    main()

