import os
import random
import copy
import shlex
import pandas as pd
import glob
import subprocess
from Bio import SeqIO

def calculate_effective_copy_number(row, df):
    if row['seq_type'] == 'chromosome':
        return row['copy_number']

    elif row['seq_type'] == 'plasmid':
        chr_row = df[(df['taxa'] == row['taxa']) & \
                (df['seq_type'] == 'chromosome')]
        return row['copy_number'] * chr_row.iloc[0]['copy_number']

    else:
        raise ValueError

def build_metagenome_fasta(seq_data_folder, metadata):
    # for each taxa in our seq folder amplify and write to a file
    for ix, seq in metadata.iterrows():
        print(f"Amplifying x{seq['actual_copy_number']}: {seq['accession']}")
        in_fp = os.path.join(seq_data_folder, seq['taxa'], seq['sample'],
            seq['seq_type'], seq['accession'] + ".fasta")

        sequences = []
        record_ix = 0
        for record in SeqIO.parse(in_fp, "fasta"):
            record_ix += 1
            for copy_num in range(seq['actual_copy_number']):
                copy_record = copy.deepcopy(record)
                copy_record.id = copy_record.id + f"_{record_ix}_{copy_num}"
                copy_record.description = ""
                sequences.append(copy_record)

        with open('metagenome_fasta.fna', 'a') as fh:
            SeqIO.write(sequences, fh, "fasta")
    return 'metagenome_fasta.fna'

def select_chromosome_relative_abundance():
    """
    lognormal taxa distribution (mean 1 sigma 2)
    """
    mu = 1
    sigma = 2

    abundance = int(random.lognormvariate(mu, sigma))

    if abundance < 1:
        abundance = 1

    return abundance, 'chromosome'


def select_plasmid_copy_number():
    """
    gamma distribution biased towards lower bound for each
    """
    alpha = 4
    beta = 1
    regimes = {'low': [1, 20, 1],
               'medium': [20, 100, 10],
               'high': [500, 1000, 150]}

    copy_number_regime = random.choice(['low', 'low', 'medium', 'medium', 'high'])
    minimum, maximum, scaling = regimes[copy_number_regime]

    copy_number = random.gammavariate(alpha, beta)
    copy_number = int(copy_number * scaling)

    if copy_number < minimum:
        copy_number = minimum
    elif copy_number > maximum:
        copy_number = maximum

    return copy_number, copy_number_regime


def add_copy_numbers(metadata_fp):

    simulation_metadata_fp = 'metadata_for_simulation.tsv'
    out_fh = open(simulation_metadata_fp, 'w')
    out_fh.write("taxa\tsample\taccession\tseq_type\tcopy_number\tcopy_regime\n")
    with open(metadata_fp) as fh:
        for line in fh:
            line = line.strip().split('\t')
            if line[3] == 'chromosome':
                relative_abundance, regime = select_chromosome_relative_abundance()
                line.append(str(relative_abundance))
                line.append(regime)
            elif line[3] == 'plasmid':
                copy_number, regime = select_plasmid_copy_number()
                line.append(str(copy_number))
                line.append(regime)

            line = "\t".join(line)
            out_fh.write(line + "\n")
    out_fh.close()
    return simulation_metadata_fp

if __name__ == '__main__':

    random.seed(42)

    # simulate copy numbers
    sim_metadata_fp = add_copy_numbers('../data/taxa_metadata.tsv')

    # calculate the effective copy number i.e. taxa abundance x plasmid copy number
    metadata = pd.read_csv(sim_metadata_fp, sep='\t')

    metadata['actual_copy_number'] = metadata.apply(\
            lambda x: calculate_effective_copy_number(x, metadata), axis=1)

    seq_data_folder = '../data/sequences'
    metagenome_fasta_fp = build_metagenome_fasta(seq_data_folder, metadata)

    # fragement length is a bit arbitrary
    subprocess.check_call(f'art_illumina --rndSeed 42 --seqSys MSv3 --in {metagenome_fasta_fp} --len 250 --fcov 5 --out metagenome --mflen 1000 --sdev 50',
                            shell=True)

