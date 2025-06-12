import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import matplotlib.lines as mlines

length_df = pd.read_csv('../results/PD_mul_var_different_lengths.csv')

# Calculate the number of primers
length_df['num_primers_ilp'] = length_df["ILP Solution"].str.count(r"\(") / 2
length_df['num_primers_greedy'] = length_df["Greedy Solution"].str.count(r"\(") / 2

num_primers_ilp=length_df['num_primers_ilp']
num_primers_greedy=length_df['num_primers_greedy']

# Extracting data from the DataFrame
lengths = length_df['seq length:']

ilp_cost = length_df['ILP Objective']

length_df['total ilp runtime'] = length_df['Time (ILP)'] + length_df['Time (Setup)'] + length_df["Time (Graph)"]

ilp_time = length_df['total ilp runtime']

greedy_time = length_df['Greedy Time'] + length_df["Time (Graph)"]

greedy_cost = length_df["Greedy Objective"]

# Creating a 1x2 subplot grid
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(23, 20))

# Plotting sequence length versus cost for ILP and Greedy on the first subplot
ax1.plot(lengths, ilp_cost, label='PrimerDesigner', marker='o', linestyle='-')
ax1.plot(lengths, greedy_cost, label='Greedy', marker='o', linestyle='-')
ax1.set_xlabel('Coding sequence length [nt]', fontsize=20)
ax1.set_ylabel('Cost', fontsize=20)
ax1.set_ylim(55,140)
# Create a custom legend entry
primer_pairs_entry = mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                                   markersize=5, label='(•) Number of primer pairs')
# Retrieve the current handles and labels from the plot
handles, labels = ax1.get_legend_handles_labels()

# Add your custom legend entry
handles.append(primer_pairs_entry)

# Create the legend with the updated handles
ax1.legend(handles=handles, fontsize=20)
ax1.tick_params(axis='y', labelsize=16)
ax1.tick_params(axis='x', labelsize=16)
ax1.grid(False)

ilp_color = ax1.get_lines()[0].get_color()
greedy_color = ax1.get_lines()[1].get_color()

counter = 0
for x, y, n in zip(lengths, ilp_cost, num_primers_ilp):
    if counter > 0:
        ax1.text(x + 6, y - 1.5, "(" + str(int(n)) + ")", ha='center', va='top', color=ilp_color, fontsize=15, zorder=10)
    else:
        ax1.text(x + 6, y - 1, "(" + str(int(n)) + ")", ha='center', va='top', color=ilp_color, fontsize=15,
                 zorder=10)
    counter += 1

counter = 0
for x, y, n in zip(lengths, greedy_cost, num_primers_greedy):
    counter += 1
    if counter < len(lengths):
        ax1.text(x - 35, y + 1, "(" + str(int(n)) + ")", ha='center', va='bottom', color=greedy_color, fontsize=15,
                 zorder=10)
    else:
        ax1.text(x - 8, y, "(" + str(int(n)) + ")", ha='center', va='bottom', color=greedy_color, fontsize=15,
                 zorder=10)

# Plotting runtime for PrimerDesigner vs Greedy on the second subplot
ax2.plot(lengths, ilp_time, label='PrimerDesigner', marker='s', linestyle='--')
ax2.plot(lengths, greedy_time, label='Greedy', marker='s', linestyle='--')
ax2.set_xlabel('Coding sequence length [nt]', fontsize=20)
ax2.set_ylabel('Runtime [seconds]', fontsize=20)
ax2.tick_params(axis='y', labelsize=16)
ax2.tick_params(axis='x', labelsize=16)
ax2.grid(False)

counter=0
for x, y, n in zip(lengths, ilp_time, num_primers_ilp):
    if counter==0:
        ax2.text(x + 42, y - 30, str(int(y)), ha='center', va='top', color=ilp_color, fontsize=15, zorder=10)
    else:
        ax2.text(x + 25, y - 35, str(int(y)), ha='center', va='top', color=ilp_color, fontsize=15, zorder=10)
    counter+=1

counter=0
for x, y, n in zip(lengths, greedy_time, num_primers_greedy):
    if counter==0:
        ax2.text(x - 15, y + 15, str(int(y)), ha='center', va='bottom', color=greedy_color, fontsize=15, zorder=10)
    else:
        ax2.text(x + 25, y + 15, str(int(y)), ha='center', va='bottom', color=greedy_color, fontsize=15, zorder=10)
    counter+=1


variant_df = pd.read_csv('../results/PD_mul_var_different_num_variants.csv')

# Extract the necessary data for plotting
num_proteins = variant_df['Num proteins:']
num_primers_ilp = variant_df["ILP Length"]
num_primers_greedy = variant_df["Greedy Length"]
ilp_objective = variant_df['ILP Objective']
greedy_objective = variant_df['Greedy Objective']
ilp_time = variant_df['Time (ILP)'] + variant_df['Time (Graph)'] + variant_df['Time (Setup)']
greedy_time = variant_df['Greedy Time'] + variant_df['Time (Graph)']

# Plotting ILP and Greedy Objective
ax3.plot(num_proteins, ilp_objective, label='PrimerDesigner', marker='o', linestyle='-')
ax3.plot(num_proteins[:2], greedy_objective[:2], label='Greedy', marker='o', linestyle='-')
ax3.set_xlabel('Number of variants of the same protein', fontsize=20)
ax3.set_ylabel('Cost', fontsize=20)
ax3.grid(False)
ax3.set_xticks(np.arange(2, 7, 1))
ax3.tick_params(axis='y', labelsize=16)
ax3.tick_params(axis='x', labelsize=16)

counter = 0
for x, y, n in zip(num_proteins, ilp_objective, num_primers_ilp):
    if counter == 0:
        ax3.text(x+0.1, y-5, "(" + str(int(n)) + ")", ha='center', va='top', color=ilp_color, fontsize=15, zorder=10)
    else:
        ax3.text(x+0.07, y - 10, "(" + str(int(n)) + ")", ha='center', va='top', color=ilp_color, fontsize=15,
                 zorder=10)
    counter += 1

counter = 0
for x, y, n in zip(num_proteins[:2], greedy_objective[:2], num_primers_greedy[:2]):
    counter += 1
    if counter < len(lengths):
        ax3.text(x+0.1, y + 18, "(" + str(int(n)) + ")", ha='center', va='bottom', color=greedy_color, fontsize=15,
                 zorder=10)
    else:
        ax3.text(x+0.1, y + 10, "(" + str(int(n)) + ")", ha='center', va='bottom', color=greedy_color, fontsize=15,
                 zorder=10)


# Plotting ILP and Greedy Time
ax4.plot(num_proteins, ilp_time, label='ILP Time', marker='s', linestyle='--')
ax4.plot(num_proteins[:2], greedy_time[:2], label='Greedy Time', marker='s', linestyle='--')
ax4.set_xticks(np.arange(2, 7, 1))
ax4.set_xlabel('Number of variants of the same protein', fontsize=20)
ax4.set_ylabel('Runtime [seconds]', fontsize=20)
ax4.grid(False)
ax4.tick_params(axis='y', labelsize=16)
ax4.tick_params(axis='x', labelsize=16)

for x, y, n in zip(num_proteins, ilp_time, num_primers_ilp):
    ax4.text(x, y - 50, str(int(y)), ha='center', va='top', color=ilp_color, fontsize=15, zorder=10)

for x, y, n in zip(num_proteins[:2], greedy_time[:2], num_primers_greedy[:2]):
    ax4.text(x, y + 15, str(int(y)), ha='center', va='bottom', color=greedy_color, fontsize=15)


protein_df = pd.read_csv("../results/PD_mul-var_10_proteins.csv")

proteins = protein_df['Protein']

protein_accession =protein_df['CCDS accession code']

protein_length = protein_df['Sequence length (nt)']

num_pairs_greedy =  protein_df['Greedy pairs']

num_pairs_ilp =  protein_df['PrimerDesigner Pairs']

cost_ilp =  protein_df['PrimerDesigner cost']

cost_greedy = protein_df['Greedy cost']


bar_width = 0.45

index = np.arange(len(proteins))


xtick_labels = [f"{protein}\n{length}" for protein, accession, length in
                zip(proteins, protein_accession, protein_length)]
accession_labels = [f"{accession}" for accession in protein_accession]

ax5.set_xticks(index + bar_width / 2)
ax5.set_xticklabels(xtick_labels, fontsize=17, ha='center', va='top')
ax5.tick_params(axis='y', labelsize=16)

bars_cost_ilp = ax5.bar(index, cost_ilp, bar_width, label='PrimerDesigner')
bars_cost_greedy = ax5.bar(index + bar_width, cost_greedy, bar_width, label='Greedy')

# Add num_pairs above each bar for Cost
for rect, label in zip(bars_cost_ilp, num_pairs_ilp):
    height = rect.get_height()
    ax5.text(
        rect.get_x() + rect.get_width() / 2, height, "(" + str(label) + ")", ha="center", va="bottom",
        color='black', fontsize=15)

for rect, label in zip(bars_cost_greedy, num_pairs_greedy):
    height = rect.get_height()
    ax5.text(
        rect.get_x() + rect.get_width() / 2, height, "(" + str(label) + ")", ha="center", va="bottom",
        color='black', fontsize=15)

for i, cost in enumerate(cost_ilp):
    ax5.text(i, cost / 2, f"{math.ceil(cost):.0f}", ha='center', va='top', color='white', rotation=90, fontsize=15)

for i, cost in enumerate(cost_greedy):
    ax5.text(i + bar_width, cost / 2, f"{math.ceil(cost):.0f}", ha='center', va='top', color='black', rotation=90,
             fontsize=15)

ax5.set_xlabel('Protein\nCoding sequence length [nt]', fontsize=20)
ax5.set_ylabel('Cost', fontsize=20)
ax5.set_xticks(index + bar_width / 2)
ax5.set_xticklabels(xtick_labels, fontsize=17, ha='center')
ax5.set_ylim(0.1, 700)

# Set y-axis to logarithmic scale
ax5.set_yscale('log')

ax5.set_ylim(0.1, 10000)

runtime_ilp = [1010.414613, 1395.048043, 1110.6774, 1097.257544, 1102.269087, 828.4370906, 681.3019264, 804.9078047,
               1263.522754, 2164.191019]
setup_ilp = [1970.94, 2849.570877, 2187.58386, 2305.969112, 2171.757641, 1676.056677, 1308.60551, 1533.514255,
             2581.172325, 3547.182489]

total_run_ilp = []
for i in range(len(runtime_ilp)):
    total_run_ilp.append(runtime_ilp[i] + setup_ilp[i])

runtime_greedy = [244, 349, 264, 265, 257, 187, 145, 176, 330, 427]

bars_run_ilp = ax6.bar(index, total_run_ilp, bar_width, label='PrimerDesigner')
bars_run_greedy = ax6.bar(index + bar_width, runtime_greedy, bar_width, label='Greedy')

for i, run in enumerate(total_run_ilp):
    ax6.text(i, run / 2, f"{run:.0f}", ha='center', va='top', color='white', rotation=90, fontsize=15)

for i, run in enumerate(runtime_greedy):
    ax6.text(i + bar_width, run / 2, f"{run:.0f}", ha='center', va='top', color='black', fontsize=15, rotation=90)

ax6.set_xlabel('Protein\nCoding sequence length [nt]', fontsize=19)
ax6.set_ylabel('Runtime [seconds]', fontsize=19)
ax6.set_xticks(index + bar_width / 2)
ax6.set_xticklabels(xtick_labels, fontsize=17, ha='center')
ax6.set_ylim(0.1, 3000)
ax6.tick_params(axis='y', labelsize=16)

fig.text(0.015, 0.98, "A", fontsize=25, va='center', ha='center', fontweight='bold')
fig.text(0.51, 0.98, "B", fontsize=25, va='center', ha='center', fontweight='bold')

# Annotate the figure with labels
fig.text(0.03, 0.68, "C", fontsize=25, va='center', ha='center', fontweight='bold')
fig.text(0.52, 0.68, "D", fontsize=25, va='center', ha='center', fontweight='bold')

fig.text(0.03, 0.36, "E", fontsize=25, va='center', ha='center', fontweight='bold')
fig.text(0.52, 0.36, "F", fontsize=25, va='center', ha='center', fontweight='bold')


# Set y-axis to logarithmic scale
ax6.set_yscale('log')
ax6.set_ylim(0.1, 20000)

plt.subplots_adjust(left=0.05, right=0.96, wspace=0.15, bottom=0.08, top=0.98)

plt.savefig("../results/figure4.png", dpi=300)
