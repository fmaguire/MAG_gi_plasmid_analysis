import random
import pandas as pd

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

    copy_number_regime = random.choice(['low', 'medium', 'high'])
    minimum, maximum, scaling = regimes[copy_number_regime]

    copy_number = random.gammavariate(alpha, beta)
    copy_number = int(copy_number * scaling)

    if copy_number < minimum:
        copy_number = minimum
    elif copy_number > maximum:
        copy_number = maximum

    return copy_number, copy_number_regime


if __name__ == '__main__':

    out_fh = open('metadata_for_simulation.tsv', 'w')
    with open('../data/taxa_metadata.tsv') as fh:
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


