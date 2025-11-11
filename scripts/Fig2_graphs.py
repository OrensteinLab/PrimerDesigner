import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import ast
import matplotlib.lines as mlines

run_df = pd.read_csv("../New_results/PD-mul-ILP.csv")

num_proteins= run_df["num_proteins"]

greedy_efficiency = run_df['greedy_objective']

ilp_efficiency = run_df['ilp_objective']

ilp_time = list(run_df['ilp_optimize_time_sec'] + run_df['ilp_setup_time_sec'])

greedy_time = list(run_df['greedy_time_sec'])

forbidden_pairs = run_df['ilp_single_forbidden_cnt'] + run_df['ilp_inter_forbidden_cnt']

num_primers_ILP = run_df['ilp_path_length']

num_primers_greedy = run_df['greedy_path_length']

avg_ilp_efficiency = ilp_efficiency/ num_primers_ILP
avg_greedy_efficiency = greedy_efficiency/ num_primers_greedy

num_primer_pairs_ILP = num_primers_ILP/2

# Creating a 1x2 subplot grid
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 13))
offset = 0.1
offset_num_proteins = num_proteins + offset

# Plotting sequence length versus cost for ILP and Greedy on the first subplot
ax1.plot(offset_num_proteins, avg_ilp_efficiency, label='PD-mul-ILP', marker='o', linestyle='-')
ax1.plot(num_proteins, avg_greedy_efficiency, label='Greedy', marker='o', linestyle='-')
# ax1.set_xlabel('Number of proteins\n(Number of potential cross-hybridization risks)', fontsize=20)
ax1.set_ylabel('Average primer efficiency', fontsize=26)

# Create a custom legend entry
primer_pairs_entry = mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                                   markersize=5, label='(•) Number of primer pairs')
# Retrieve the current handles and labels from the plot
handles, labels = ax1.get_legend_handles_labels()

# Add your custom legend entry
handles.append(primer_pairs_entry)

# Create the legend with the updated handles
ax1.legend(handles=handles,loc="upper left", fontsize=23)
ax1.set_ylim(0.34,0.49)

# ax1.legend(fontsize=20)
ax1.tick_params(axis='y', labelsize=24)
ax1.tick_params(axis='x', labelsize=24)
ax1.grid(False)

ilp_color = ax1.get_lines()[0].get_color()
greedy_color = ax1.get_lines()[1].get_color()

# Plotting runtime for PrimerDesigner vs Greedy on the second subplot
ax2.plot(num_proteins, ilp_time, label='PD-mul-ILP', marker='s', linestyle='--')
ax2.plot(num_proteins, greedy_time, label='Greedy', marker='s', linestyle='--')
# ax2.set_xlabel('Number of proteins\n(Number of potential cross-hybridization risks)', fontsize=20)
ax2.set_ylabel('Runtime [seconds]', fontsize=26)
ax2.tick_params(axis='y', labelsize=24)
ax2.tick_params(axis='x', labelsize=24)
ax2.grid(False)
ax2.set_ylim(-200,4000)

fig.text(0.25, 0.03, 'Number of proteins\n(Number of potential cross-hybridization risks)', ha='center', fontsize=26)
fig.text(0.75, 0.03, 'Number of proteins\n(Number of potential cross-hybridization risks)', ha='center', fontsize=26)

vertical_offset = 0.0035  # Adjust this value as needed to position the text above the markers

# Adjusting text annotations for ILP cost
for x, y in zip(num_proteins, avg_ilp_efficiency):
    ax1.text(x-0.12, y + vertical_offset, str(round(y,2)), ha='center', va='bottom', color=greedy_color, fontsize=20, zorder=10)

# Adjusting text annotations for ILP cost
for x, y, n in zip(num_proteins, avg_ilp_efficiency, num_primer_pairs_ILP):
    ax1.text(x+0.25, y-vertical_offset, "("+str(int(n))+")", ha='center', va='top', color=greedy_color, fontsize=20, zorder=10)


vertical_offset = 75

# Adjusting text annotations for runtime
for x, y in zip(num_proteins, ilp_time):
    ax2.text(x, y + vertical_offset, str(int(y)), ha='center', va='bottom', color=ilp_color, fontsize=20, zorder=10)

for x, y in zip(num_proteins, greedy_time):
    ax2.text(x, y - vertical_offset, str(int(y)), ha='center', va='top', color=greedy_color, fontsize=20, zorder=10)

# Adding text annotations for the x-ticks with numbers in parentheses
for i, (np, tp) in enumerate(zip(num_proteins, forbidden_pairs)):
    ax1.text(np, -0.05, f"({int(tp)})", ha='center', va='top', fontsize=24, rotation=45, transform=ax1.get_xaxis_transform())
    ax2.text(np, -0.05, f"({int(tp)})", ha='center', va='top', fontsize=24, rotation=45, transform=ax2.get_xaxis_transform())


# Set the x-ticks to show the number of proteins
ax1.set_xticks(num_proteins)
ax2.set_xticks(num_proteins)

fig.text(0.03, 0.93, "A", fontsize=30, va='center', ha='center', fontweight='bold')
fig.text(0.52, 0.93, "B", fontsize=30, va='center', ha='center', fontweight='bold')

plt.subplots_adjust(bottom=0.23, wspace=0.16, top=0.9, left=0.05, right=0.95)

plt.savefig("../results/figure2.png",dpi=300)
