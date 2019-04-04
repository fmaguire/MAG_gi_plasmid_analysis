import os
import shlex
import subprocess
import glob

if __name__ = "__main__":

    # make single fasta reference genomes for quast
SAMN02603599	NC_004851	plasmid	39	medium
    reference_genomes = '../mag_assessment/reference_genomes'
    if not os.path.exists(reference_genomes):
        os.mkdir(reference_genomes)

    for genome in glob.glob(seq_data_folder + "/*"):
        fastas = [shlex.quote(x) for x in glob.glob(genome + "/*/*/*.fasta")]

        subprocess.check_call(f"cat {' '.join(fastas)} > {reference_genomes + shlex.quote(os.basename(genome))}", shell=True)



