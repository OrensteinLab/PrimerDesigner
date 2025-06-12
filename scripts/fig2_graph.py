import random as rand

import pandas as pd
import primer3
import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq
import random

primer_lmin, primer_lmax = 18,30

random_data = {}
real_data={}

# SpAP protein
upstream_nt = 'GCTAGTGGTGCTAGCCCCGCGAAATTAATACGACTCACTATAGGGTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT'
mutreg_nt='ATGCAAAGCCCAGCACCTGCCGCAGCGCCTGCCCCTGCGGCACGTTCCATCGCAGCTACGCCTCCTAAACTGATCGTGGCAATTAGCGTGGACCAGTTTAGTGCAGACTTGTTCTCGGAGTATCGTCAATATTACACCGGAGGTTTAAAGCGTCTTACATCCGAAGGAGCTGTGTTCCCACGTGGTTATCAGAGTCATGCGGCAACAGAAACGTGTCCTGGTCACTCAACGATCCTGACAGGATCACGTCCGTCACGTACGGGTATTATCGCTAATAACTGGTTCGACTTGGACGCAAAGCGTGAGGATAAAAATCTGTACTGTGCTGAGGATGAATCCCAACCCGGTAGTTCGTCTGACAAGTACGAAGCTTCGCCACTGCACTTAAAGGTACCCACCCTGGGGGGACGCATGAAAGCCGCCAATCCTGCGACTCGTGTCGTCTCTGTTGCCGGCAAGGATCGCGCGGCCATTATGATGGGTGGCGCCACAGCGGATCAGGTCTGGTGGTTAGGGGGGCCTCAGGGGTATGTTTCGTATAAGGGTGTAGCGCCAACTCCCCTTGTAACACAGGTCAATCAGGCCTTTGCACAGCGCTTAGCTCAGCCGAACCCGGGATTTGAGTTGCCTGCTCAGTGCGTCAGCAAGGACTTTCCTGTTCAAGCGGGAAATCGCACAGTGGGTACCGGCCGCTTCGCCCGTGATGCTGGTGACTACAAAGGTTTTCGCATTTCCCCGGAGCAGGATGCTATGACGCTTGCATTCGCTGCCGCGGCCATTGAAAATATGCAATTAGGGAAGCAGGCCCAGACCGATATTATTAGCATTGGACTGAGCGCTACGGATTACGTGGGACACACCTTCGGCACGGAGGGTACGGAGAGTTGCATCCAAGTGGATCGTTTAGACACGGAGCTTGGTGCATTCTTTGATAAACTGGATAAGGATGGGATTGACTACGTAGTAGTGCTGACTGCAGATCATGGAGGACACGATCTGCCCGAACGTCATCGTATGAATGCCATGCCGATGGAACAGCGCGTAGACATGGCCCTGACACCTAAAGCTCTGAATGCTACCATCGCTGAGAAAGCTGGCCTTCCGGGCAAAAAGGTTATTTGGTCAGATGGACCTTCTGGCGATATTTACTATGATAAGGGCCTTACAGCCGCTCAACGTGCCCGTGTTGAAACCGAGGCGTTAAAATACTTGCGCGCGCATCCCCAAGTACAGACTGTATTCACTAAGGCGGAAATCGCGGCTACCCCTTCTCCGTCGGGACCACCTGAGAGCTGGAGTTTGATCCAGGAAGCTCGCGCGTCATTTTACCCGTCGCGCTCCGGGGACCTGTTACTTTTATTGAAACCTCGTGTGATGAGCATTCCTGAGCAAGCAGTCATGGGCTCGGTTGCAACCCATGGATCTCCATGGGATACGGATCGCCGTGTGCCTATCCTGTTTTGGCGCAAAGGTATGCAGCATTTCGAACAACCCTTAGGAGTAGAGACTGTTGATATTTTGCCCTCCTTGGCTGCACTTATTAAGCTTCCTGTTCCTAAGGATCAGATCGACGGCCGCTGTCTGGACTTGGTCGCCGGCAAGGATGATTCCTGTGCTGGACAGGGA'
downstream_nt = 'GGAGGGTCTGGGGGAGGAGGCAGTGGCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCG'

sequence_nt = upstream_nt + mutreg_nt + downstream_nt

def generate_random_dna_sequence(length):
    nucleotides = ['A', 'C', 'G', 'T']
    random_sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    return random_sequence

def revcomp(seq):
    return str(Seq(seq).reverse_complement())

def calc_tm(piece1, piece2):
    tm = primer3.bindings.calc_heterodimer(piece1,piece2, mv_conc=50.0, dv_conc=1.5,
                                            dntp_conc=0.6,
                                            dna_conc=50.0, temp_c=37.0, max_loop=30).tm
    return tm


def main():
    num_primers=0

    while num_primers<10000:

        length1 = rand.randint(primer_lmin, primer_lmax)
        length2 = rand.randint(primer_lmin, primer_lmax)

        piece1=generate_random_dna_sequence(length1)
        piece2=generate_random_dna_sequence(length2)

        num_primers+=1
        tm = calc_tm(piece1,piece2)
        random_data[(piece1, piece2)] = tm

    num_primers=0

    while num_primers<10000:
        length1 = rand.randint(primer_lmin, primer_lmax)
        length2 = rand.randint(primer_lmin, primer_lmax)

        start1 = rand.randint(0, len(sequence_nt) - length1)
        start2 = rand.randint(0, len(sequence_nt) - length2)

        piece1 = sequence_nt[start1: length1 + start1 + 1]
        piece2 = sequence_nt[start2: length2 + start2 + 1]

        # # ensure pieces dont overlap
        if start1 + length1 <= start2 or start2 + length2 <= start1:
            num_primers += 1
            tm = calc_tm(piece1, piece2)
            real_data[(piece1, piece2)] = tm

    # Create a figure and a subplot with 1 row and 2 columns
    fig, axs = plt.subplots(1, 2, figsize=(18, 7))
    plt.subplots_adjust(bottom=0.13)

    plt.subplots_adjust(left=0.05, right=0.95, wspace=0.2, bottom=0.1, top=0.95)

    # Plot the first histogram in the first subplot
    axs[0].hist(real_data.values(), bins=50,edgecolor='black')
    axs[0].set_xlabel("Melting temperature (tm)",fontsize=18)
    axs[0].set_ylabel("Frequency",fontsize=18)
    axs[0].axvline(x=45, color='red', linestyle='--', linewidth=2, alpha=0.5)
    axs[0].text(45, 400, '45°C   Cross-hybridization risk', color='black', fontsize=14, va='center', ha='left',rotation=270)
    axs[0].set_xlim(-70, 60)
    axs[0].set_ylim(0, 780)
    axs[0].set_title("Real-protein primer pairs",fontsize=22,fontweight='bold')
    axs[0].axvspan(45, axs[0].get_xlim()[1], facecolor='grey', alpha=0.1)


    # Plot the second histogram in the second subplot
    axs[1].hist(random_data.values(), bins=50, edgecolor='black', color='blue')
    axs[1].set_xlabel("Melting temperature (tm)", fontsize=18)
    axs[1].set_ylabel("Frequency",fontsize=18)
    axs[1].axvline(x=45, color='red', linestyle='--', linewidth=2,alpha=0.5)
    axs[1].text(45, 400, '45°C   Cross-hybridization risk', color='black', fontsize=14, va='center', ha='left',rotation=270)
    axs[1].set_xlim(-70, 60)
    axs[1].set_ylim(0, 780)
    axs[1].set_title("Random primer pairs",fontsize=22,fontweight='bold')
    axs[1].axvspan(45, axs[1].get_xlim()[1], facecolor='grey', alpha=0.1)

    fig.text(0.015, 0.92, "A", fontsize=25, va='center', ha='center',fontweight='bold')
    fig.text(0.51, 0.92, "B", fontsize=25, va='center', ha='center',fontweight='bold')

    # fig.text(0.14,0.88, "Random primers", fontsize=20, va='center', ha='center', fontweight='bold', color='black')
    # fig.text(0.645, 0.88, "Real-protein primers", fontsize=20, va='center', ha='center', fontweight='bold',color= 'black')

    for ax in axs:
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

    plt.subplots_adjust(top=0.88)

    # # Show the combined plot
    plt.savefig("../results/figure2.png",dpi=300)
    
    results_tm_pairs = pd.DataFrame()

    results_tm_pairs['random_tm'] =random_data.values()
    results_tm_pairs['real_tm'] = real_data.values()

    # results_tm_pairs.to_csv("../results/tm_pairs.csv")


if __name__ == "__main__":
    main()
