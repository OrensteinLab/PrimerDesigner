import json
import gurobipy as gp
import pandas as pd
import primer3 as p3
from Bio.Seq import Seq


# constant upstream and downstream regions for all proteins
# upstream_nt = 'GCTAGTGGTGCTAGCCCC'

upstream_nt = 'GCTAGTGGTGCTAGCCCCGCGAAATTAATACGACTCACTATAGGGTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT'
downstream_nt = 'GGAGGGTCTGGGGGAGGAGGCAGTGGCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCG'

# defines global cross-hybridization thershold
MAX_TM = 45

def calc_tm(seq1, seq2):
    return p3.bindings.calc_heterodimer(
        seq1, seq2,
        mv_conc=50.0, dv_conc=1.5,
        dntp_conc=0.6, dna_conc=50.0,
        temp_c=37.0, max_loop=30
    ).tm

def get_model(file_path='gurobi.json'):

    with open(file_path, 'r') as json_file:
        params = json.load(json_file)

    env = gp.Env(params=params)

    # Create the model within the Gurobi environment
    model = gp.Model('min-sum', env=env)

    return model


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

def revcomp(seq):
  return str(Seq(seq).reverse_complement())

def subsequences(sequence,primer_lmin,primer_lmax): #Generates all subsequences w/ all poss. start-stop pairs
  ls = []
  for j in range(primer_lmin, primer_lmax+1): #length
    for i in range(len(sequence)-j+1): #starting index
      start = i
      stop = i+j
      ls.append([sequence[start:stop], start, stop, stop-start])
  return ls