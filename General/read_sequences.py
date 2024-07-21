from General.constants import *

def read_sequences(file_path):

    mutreg_regions = []
    protein_names=[]

    # read  protein coding-sequences from the file path
    with open(file_path) as file:
        for line in file.readlines():
            p_name, mutreg_region = line.strip().split('\t')
            mutreg_regions.append(mutreg_region)
            protein_names.append(p_name)

    full_sequences = []

    # add constant upstream and downstream regions to each sequence
    for mutreg_nt in mutreg_regions:
        sequence = upstream_nt + mutreg_nt + downstream_nt
        full_sequences.append(sequence)


    return mutreg_regions,full_sequences,protein_names
