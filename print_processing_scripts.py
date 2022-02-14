import argparse
import glob
import sys
import os
import stat


parser = argparse.ArgumentParser(description="Write shell scripts to run initial processing on RNA-Seq data: read mapping, quality checks, and transcript quantification.")
parser.add_argument('-o',
                    '--out_dir',
                    type=str,
                    help='path to the project containing dataset to be processed. Datasets should be in a fastq subdirectory of this specified directory.')

# Globals
INITIAL_FOLDER_LIST = ["run_commands", "star", "rsem", "qualimap", "sorted"] #folders to generate

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


def get_fastq_dir(out_dir):
    """Return the path to the fastq subdirectory of the specified directory out_dir.
    
    Args:
        out_dir: directory

    Returns:
        path to fastq directory
    """
    return os.path.join(out_dir, "fastq")


def make_mate_string(prefix, out_dir):
    """Create the paired lists of fastq files to process for a given dataset

    Illumina datasets are downloaded from BaseSpace in a series of folders which are 
    specified by the sample name (which forms a prefix) and the lane number.
    Within each of these folders are the two "mate" fastq files with each end of
    the paired-end reads. For read mapping of a given dataset, we want to pass
    all the fastq files for a given sample (from different lanes) to STAR or RSEM. 
    
    The required format is two comma separated lists of paths to fastq files, with
    the read1 files in the first list and the read2 files in the second list. For
    example, for a dataset S1 spread over 3 lanes, we need the following 
    mate string to specify the files to read:
        
        S1_L001_R1.fastq,S1_L002_R1.fastq,S1_L003_R1.fastq S1_L001_R2.fastq,S1_L002_R2.fastq,S1_L003_R2.fastq
    
    where each filename above is a path to that file. 

    Args:
        prefix: A string that specifies the particular sample to generate a mate
          string for
        out_dir: The directory containing the fastq subdirectory of files to process

    Returns:
        The mate string, which consists of two comma-separated lists of paths to
        fastq files corresponding to the different "mate" fastq files, separated
        by a space
    """

    read1_list = []
    read2_list = []
    dirname_list = sorted(glob.glob(os.path.join(get_fastq_dir(out_dir), prefix + "_*/")))
    for dirname in dirname_list:
        read_names_1, read_names_2 = sorted(glob.glob(os.path.join(dirname, "*.fastq.gz")))
        read1_list.append(read_names_1)
        read2_list.append(read_names_2)
    mate_string = ",".join(read1_list) + " " + ",".join(read2_list)
    return mate_string


def make_file_executable(file_name):
    """Add executable permissions to a file (like chmod +x)

    Args:
        file_name: path of the file whose permissions to modify
    """

    st = os.stat(file_name)
    #Add executable permissions for user, group, and others (like chmod +x)
    os.chmod(file_name, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def initialize_folders(out_dir):
    """Create subfolders in the specified directory.

    Args:
        out_dir: directory to create the folders in
    """

    for name in INITIAL_FOLDER_LIST:
        folder_name = os.path.join(out_dir, name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)


def get_unique_prefixes(out_dir):
    """Get a list of unique prefixes, specifying datasets to process.

    Args:
        out_dir: folder containing the fastq subdirectory with files to process

    Returns:
        prefixes: a set of unique prefixes, which identify different samples
    """

    filenames = glob.glob(os.path.join(get_fastq_dir(out_dir),  "*/"))
    prefixes = set([os.path.basename(os.path.dirname(filename)).split("_")[0] 
                        for filename in filenames])
    return prefixes


def generate_pipeline_script(prefix, out_dir):
    """Write shell commands to process reads from a sample and make it executable.

    The processing steps that must be done on a sample can be done independently
    on each sample, so writing a shell script that will be executed in a different
    cluster job for each sample will make most efficient use of time. In conjunction
    with the global variables that define the content of each processing command, 
    this function specifies which processing steps will be carried out on each 
    sample, and in what order.

    Args:
        prefix: a filename prefix that specifies all files corresponding to this sample
        out_dir: folder containing fastq files for this sample (in a fastq subdirectory)

    Returns:
        full_pipeline: A string consisting of fully specified shell commands to
          process this sample
    """

    mate_string = make_mate_string(prefix, out_dir)
    star_string = STAR_TEMPLATE.format(mate_string=mate_string,
                                       prefix=prefix,
                                       out_dir=out_dir)
    samtools_string = SAMTOOLS_TEMPLATE.format(prefix=prefix,
                                               out_dir=out_dir)
    qualimap_string = QUALIMAP_TEMPLATE.format(prefix=prefix,
                                              out_dir=out_dir)
    rsem_string = RSEM_TEMPLATE.format(mate_string=mate_string,
                                       prefix=prefix,
                                       out_dir=out_dir)
    full_pipeline = "{0}\n\n{1}\n\n{2}\n\n{3}".format(star_string,
                                                  samtools_string,
                                                  qualimap_string,
                                                  rsem_string)
    return full_pipeline


def save_executable_pipeline(prefix, pipeline, out_dir):
    """Save the shell commands to process a sample to a .sh file and make it executable.

    This function takes in a string containing commands to run in the shell and 
    saves it to a .sh file, and makes that file executable. The .sh file is saved
    in the run_commands subdirectory of out_dir.

    Args:
        prefix: a filename prefix that specifies all files corresponding to this sample
        pipeline: the string of shell commands that will be run to process this sample
        out_dir: folder containing fastq files for this sample (in a fastq subdirectory)
          as well as a run_commands subdirectory to save the shell script
    """
    pipeline_path = os.path.join(out_dir, "run_commands",
            PIPELINE_NAME.format(prefix))
    with open(pipeline_path, "w") as f:
            f.write(pipeline)
    #Add executable permissions for user, group, others
    make_file_executable(pipeline_path)


def save_executable_queue_request(prefixes, out_dir):
    """Save a shell command that requets jobs to process each sample (and does chmod +x).

    Once all the shell scripts for each sample have been written, this function writes
    another shell script that will submit a job for each sample's shell script so that
    they can run in parallel. The BSUB_TEMPLATE global can be modified to adjust the
    settings for each job requests.

    This shell script will also be saved in the run_commands folder of out_dir, wiht the
    name "request_processing.sh". This function also makes request_processing.sh executable

    Args:
        prefixes: a set of unique filename prefixes that specify samples to process.
        out_dir: folder containing fastq files for all samples (in a fastq subdirectory)
          and a run_commands subdirectory to save the resulting shell script.
    """

    queue_request_path = os.path.join(out_dir, "run_commands", "request_processing.sh")
    request_lines = []
    for prefix in prefixes:
        pipeline_path = os.path.join(out_dir, "run_commands", PIPELINE_NAME.format(prefix))
        request_lines.append("{0} {1}".format(BSUB_TEMPLATE, pipeline_path))

    with open(queue_request_path, "w") as f:
        f.write("\n".join(request_lines))
    make_file_executable(queue_request_path)



def main():
    args = parser.parse_args()
    out_dir = os.path.normpath(args.out_dir)
    initialize_folders(out_dir)
    prefixes = get_unique_prefixes(out_dir)
    for prefix in prefixes:
        #Write a script for each prefix/dataset and make it executable
        full_pipeline = generate_pipeline_script(prefix, out_dir)
        save_executable_pipeline(prefix, full_pipeline, out_dir)
        print(full_pipeline)
    #Write a script to request a job for each dataset and make it executable
    save_executable_queue_request(prefixes, out_dir)



if __name__ == "__main__":
    main()

