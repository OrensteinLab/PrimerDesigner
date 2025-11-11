import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import matplotlib.lines as mlines

length_df = pd.read_csv('../New_results/PD-var-ILP-increasing-lengths.csv')


num_primers_ilp = length_df['ilp_path_length']/2
num_primers_greedy = length_df['greedy_path_length']/2

# Extracting data from the DataFrame
lengths = length_df['seq_length']

ilp_efficiency = length_df['ilp_objective'] / (num_primers_ilp*2)

ilp_time = length_df['ilp_optimize_time_sec'] + length_df['ilp_setup_time_sec']

greedy_time = length_df['greedy_time_sec']

greedy_efficiency = length_df["greedy_objective"] / (num_primers_greedy*2)

# Creating a 1x2 subplot grid
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(23, 20))

# Plotting sequence length versus cost for ILP and Greedy on the first subplot
ax1.plot(lengths, ilp_efficiency, label='PD-var-ILP', marker='o', linestyle='-')
ax1.plot(lengths, greedy_efficiency, label='Greedy', marker='o', linestyle='-')
ax1.set_xlabel('Coding sequence length [nt]', fontsize=20)
ax1.set_ylabel('Average primer efficiency', fontsize=20)

ax1.set_ylim(0.25,0.60)

# Create a custom legend entry
primer_pairs_entry = mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                                   markersize=5, label='(•) Number of primer pairs')
# Retrieve the current handles and labels from the plot
handles, labels = ax1.get_legend_handles_labels()

# Add your custom legend entry
handles.append(primer_pairs_entry)

# Create the legend with the updated handles
ax1.legend(handles=handles,loc="upper left", fontsize=20)
ax1.tick_params(axis='y', labelsize=16)
ax1.tick_params(axis='x', labelsize=16)
ax1.grid(False)

ilp_color = ax1.get_lines()[0].get_color()
greedy_color = ax1.get_lines()[1].get_color()


for x, y, n in zip(lengths, ilp_efficiency, num_primers_ilp):
    if n == 18:
        ax1.text(x-5, y + 0.03, "(" + str(int(n)) + ")", ha='center', va='top', color=ilp_color, fontsize=15,
                 zorder=10)
        continue

    ax1.text(x-5, y + 0.02, "(" + str(int(n)) + ")", ha='center', va='top', color=ilp_color, fontsize=15,
                 zorder=10)


for x, y, n in zip(lengths, greedy_efficiency, num_primers_greedy):

    if n == 9 or n==16:
        ax1.text(x + 10, y - 0.028, "(" + str(int(n)) + ")", ha='center', va='bottom', color=greedy_color, fontsize=15,
                 zorder=10)
        continue
    elif n == 12:
        ax1.text(x + 25, y - 0.028, "(" + str(int(n)) + ")", ha='center', va='bottom', color=greedy_color, fontsize=15,
                 zorder=10)
        continue


    ax1.text(x + 10, y - 0.022, "(" + str(int(n)) + ")", ha='center', va='bottom', color=greedy_color, fontsize=15,
                 zorder=10)

# Plotting runtime for PrimerDesigner vs Greedy on the second subplot
ax2.plot(lengths, ilp_time, label='PD-mul-var', marker='s', linestyle='--')
ax2.plot(lengths, greedy_time, label='Greedy', marker='s', linestyle='--')
ax2.set_xlabel('Coding sequence length [nt]', fontsize=20)
ax2.set_ylabel('Runtime [seconds]', fontsize=20)
ax2.tick_params(axis='y', labelsize=16)
ax2.tick_params(axis='x', labelsize=16)
ax2.grid(False)
ax2.set_ylim(-50,2100)

for x, y, n in zip(lengths, ilp_time, num_primers_ilp):

    ax2.text(x + 25, y - 14, str(int(y)), ha='center', va='top', color=ilp_color, fontsize=15, zorder=10)



for x, y, n in zip(lengths, greedy_time, num_primers_greedy):

    ax2.text(x + 25, y + 6, str(int(y)), ha='center', va='bottom', color=greedy_color, fontsize=15, zorder=10)


variant_df = pd.read_csv('../New_results/PD-var-ILP-increasing_variants.csv')

# Extract the necessary data for plotting
num_proteins = variant_df['num_variants']
num_primers_ilp = variant_df["ilp_path_length"]/2
num_primers_greedy = variant_df["greedy_path_length"]/2
ilp_objective = variant_df['ilp_objective'] / (num_primers_ilp*2)
greedy_objective = variant_df['greedy_objective'] / (num_primers_greedy*2)
ilp_time = variant_df['ilp_optimize_time_sec']  + variant_df['ilp_setup_time_sec']
greedy_time = variant_df['greedy_time_sec'] 

# Plotting ILP and Greedy Objective
ax3.plot(num_proteins, ilp_objective, label='PD-var-ILP', marker='o', linestyle='-')
last_feasible_idx = 2   # index for x=4 if your num_proteins = [2,3,4,5,6]
ax3.plot(
    num_proteins[:last_feasible_idx+1],
    greedy_objective[:last_feasible_idx+1],
    label='Greedy',
    marker='o',
    linestyle='-',
    color=greedy_color
)

# Draw dashed horizontal infeasible line
x_start = num_proteins[last_feasible_idx]          # 4
x_end   = num_proteins.iloc[-1]                    # 6
y_level = greedy_objective[last_feasible_idx]      # y-value at x=4

ax3.hlines(
    y=y_level,
    xmin=x_start,
    xmax=x_end,
    linestyle="dashed",
    color='black',
    linewidth=1.8
)

# Add "infeasible" text under dashed line
ax3.text(
    (x_start + x_end) / 2,
    y_level - 0.006,            # slightly below dashed line
    "Infeasible",
    ha='center',
    va='top',
    fontsize=15,
    color='black'
)
ax3.set_xlabel('Number of variants of the same protein', fontsize=20)
ax3.set_ylabel('Average primer efficiency', fontsize=20)
ax3.grid(False)
ax3.set_xticks(np.arange(2, 7, 1))
ax3.tick_params(axis='y', labelsize=16)
ax3.tick_params(axis='x', labelsize=16)
ax3.set_ylim(0.37, 0.52)
for x, y, n in zip(num_proteins, ilp_objective, num_primers_ilp):

    ax3.text(x, y+0.008, "(" + str(int(n)) + ")", ha='center', va='top', color=ilp_color, fontsize=15, zorder=10)

for x, y, n in zip(num_proteins[:3], greedy_objective[:3], num_primers_greedy[:3]):

    ax3.text(x, y - 0.012, "(" + str(int(n)) + ")", ha='center', va='bottom', color=greedy_color, fontsize=15,zorder=10)

# Plotting ILP and Greedy Time
ax4.plot(num_proteins, ilp_time, label='ILP Time', marker='s', linestyle='--')
ax4.plot(num_proteins[:3], greedy_time[:3], label='Greedy Time', marker='s', linestyle='--')
ax4.set_xticks(np.arange(2, 7, 1))
ax4.set_xlabel('Number of variants of the same protein', fontsize=20)
ax4.set_ylabel('Runtime [seconds]', fontsize=20)
ax4.grid(False)
ax4.tick_params(axis='y', labelsize=16)
ax4.tick_params(axis='x', labelsize=16)
ax4.set_ylim(-30, 2300)

for x, y, n in zip(num_proteins, ilp_time, num_primers_ilp):
    ax4.text(x, y + 100, str(int(y)), ha='center', va='top', color=ilp_color, fontsize=15, zorder=10)

for x, y, n in zip(num_proteins[:3], greedy_time[:3], num_primers_greedy[:3]):
    ax4.text(x, y + 25, str(int(y)), ha='center', va='bottom', color=greedy_color, fontsize=15)


protein_df = pd.read_csv("../New_results/PD-var-ILP-different-proteins.csv")

proteins = protein_df['protein_name']

protein_length = protein_df['seq_length']

num_pairs_greedy =  protein_df['greedy_path_length']//2

num_pairs_ilp =  protein_df['ilp_path_length']//2

ilp_efficiency =  protein_df['ilp_objective'] / (num_pairs_ilp*2)

greedy_efficiency = protein_df['greedy_objective'] / (num_pairs_greedy*2)

total_run_ilp =  protein_df['ilp_optimize_time_sec'] + protein_df['ilp_setup_time_sec']

runtime_greedy = protein_df['greedy_time_sec']

bar_width = 0.45

index = np.arange(len(proteins))


xtick_labels = [f"{protein}\n{length}" for protein, length in
                zip(proteins, protein_length)]

ax5.set_xticks(index + bar_width / 2)
ax5.set_xticklabels(xtick_labels, fontsize=17, ha='center', va='top')
ax5.tick_params(axis='y', labelsize=16)

bars_cost_ilp = ax5.bar(index, ilp_efficiency, bar_width, label='PD-var-ILP')
bars_cost_greedy = ax5.bar(index + bar_width, greedy_efficiency, bar_width, label='Greedy')

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

for i, cost in enumerate(ilp_efficiency):
    ax5.text(i, cost / 2, f"{round(cost,2):.2f}", ha='center', va='top', color='white', rotation=90, fontsize=15)

for i, cost in enumerate(greedy_efficiency):
    ax5.text(i + bar_width, cost / 2, f"{round(cost,2):.2f}", ha='center', va='top', color='black', rotation=90,
             fontsize=15)

ax5.set_xlabel('Protein\nCoding sequence length [nt]', fontsize=20)
ax5.set_ylabel('Average primer efficiency', fontsize=20)
ax5.set_xticks(index + bar_width / 2)
ax5.set_xticklabels(xtick_labels, fontsize=17, ha='center')

# Set y-axis to logarithmic scale

# ax5.set_ylim(0, 33)


bars_run_ilp = ax6.bar(index, total_run_ilp, bar_width, label='PD-var-ILP')
bars_run_greedy = ax6.bar(index + bar_width, runtime_greedy, bar_width, label='Greedy')

for i, run in enumerate(total_run_ilp):
    ax6.text(i, run / 2, f"{run:.0f}", ha='center', va='top', color='white', rotation=90, fontsize=15)

for i, run in enumerate(runtime_greedy):
    ax6.text(i + bar_width, run / 2, f"{run:.0f}", ha='center', va='top', color='black', fontsize=15, rotation=90)

ax6.set_xlabel('Protein\nCoding sequence length [nt]', fontsize=19)
ax6.set_ylabel('Runtime [seconds]', fontsize=19)
ax6.set_xticks(index + bar_width / 2)
ax6.set_xticklabels(xtick_labels, fontsize=17, ha='center')
# ax6.set_ylim(0.1, 3000)
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
ax6.set_ylim(0.1, 2500)

plt.subplots_adjust(left=0.05, right=0.96, wspace=0.15, bottom=0.08, top=0.98)

plt.savefig("../results/figure3.png", dpi=300)
