import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import ast
import matplotlib.lines as mlines

forbidden_pairs = pd.read_csv("../results/Non_relaxed_forbidden_pairs.csv")

forbidden_pairs["Total"] = forbidden_pairs['single forbidden pairs'] + forbidden_pairs['multi forbidden pairs']

run_df = pd.read_csv("../results/Non_relaxed_10_proteins.csv")

num_proteins= run_df["num proteins:"]

greedy_cost = run_df['Greedy Objective']

ilp_cost = run_df['ILP Objective']

run_df['total ilp runtime'] = run_df['Time (ILP)'] + run_df['Time (Setup)'] + run_df["Time (Graph)"]+ run_df["Pair calculation time"]

run_df['total greedy runtime'] = run_df["Time (Graph)"]+ run_df["Pair calculation time"] + run_df['Greedy Time']

ilp_time = list(run_df['total ilp runtime'])

greedy_time = list(run_df['total greedy runtime'])

primer_lists=run_df['Greedy Solution']

num_primers=[]

for primer_ls in primer_lists:
    primer_ls=ast.literal_eval(primer_ls)
    num_primers.append(sum(len(ls) for ls in primer_ls))

# Your 8 data points
x = list(range(2, 10))  # This represents the positions of your values, assuming they are evenly spaced
y = ilp_time[:8]

# Perform linear regression
slope, intercept, _, _ , _= stats.linregress(x, y)

# Extrapolate the 10th value
x_10 = 10
y_10 = (slope * x_10) + intercept

ilp_time[8]=y_10

y = greedy_time[:8]

# Perform linear regression
slope, intercept,  _, _ , _ = stats.linregress(x, y)

# Extrapolate the 10th value
x_10 = 10
y_10 = (slope * x_10) + intercept

greedy_time[8]=y_10

# Creating a 1x2 subplot grid
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 13))
offset = 0.1
offset_num_proteins = num_proteins + offset

# Plotting sequence length versus cost for ILP and Greedy on the first subplot
ax1.plot(offset_num_proteins, ilp_cost, label='PrimerDesigner', marker='o', linestyle='-')
ax1.plot(num_proteins, greedy_cost, label='Greedy', marker='o', linestyle='-')
# ax1.set_xlabel('Number of proteins\n(Number of potential cross-hybridization risks)', fontsize=20)
ax1.set_ylabel('Cost', fontsize=26)

# Create a custom legend entry
primer_pairs_entry = mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                                   markersize=5, label='(•) Number of primer pairs')
# Retrieve the current handles and labels from the plot
handles, labels = ax1.get_legend_handles_labels()

# Add your custom legend entry
handles.append(primer_pairs_entry)

# Create the legend with the updated handles
ax1.legend(handles=handles, fontsize=23)

# ax1.legend(fontsize=20)
ax1.tick_params(axis='y', labelsize=24)
ax1.tick_params(axis='x', labelsize=24)
ax1.grid(False)
ax1.set_ylim(240,450)

ilp_color = ax1.get_lines()[0].get_color()
greedy_color = ax1.get_lines()[1].get_color()

# Plotting runtime for PrimerDesigner vs Greedy on the second subplot
ax2.plot(num_proteins, ilp_time, label='PrimerDesigner', marker='s', linestyle='--')
ax2.plot(num_proteins, greedy_time, label='Greedy', marker='s', linestyle='--')
# ax2.set_xlabel('Number of proteins\n(Number of potential cross-hybridization risks)', fontsize=20)
ax2.set_ylabel('Runtime [seconds]', fontsize=26)
ax2.tick_params(axis='y', labelsize=24)
ax2.tick_params(axis='x', labelsize=24)
ax2.grid(False)
ax2.set_ylim(-500,27000)

fig.text(0.25, 0.04, 'Number of proteins\n(Number of potential cross-hybridization risks)', ha='center', fontsize=26)
fig.text(0.75, 0.04, 'Number of proteins\n(Number of potential cross-hybridization risks)', ha='center', fontsize=26)

vertical_offset = 7  # Adjust this value as needed to position the text above the markers

# Adjusting text annotations for ILP cost
for x, y in zip(num_proteins, ilp_cost):
    ax1.text(x-0.05, y + vertical_offset, str(int(y)), ha='center', va='bottom', color=greedy_color, fontsize=22, zorder=10)

# Adjusting text annotations for ILP cost
for x, y, n in zip(num_proteins, ilp_cost,num_primers):
    ax1.text(x+0.25, y -7, "("+str(int(n))+")", ha='center', va='top', color=greedy_color, fontsize=22, zorder=10)


vertical_offset = 600

# Adjusting text annotations for runtime
for x, y in zip(num_proteins, ilp_time):
    ax2.text(x, y + vertical_offset, str(int(y)), ha='center', va='bottom', color=ilp_color, fontsize=22, zorder=10)

for x, y in zip(num_proteins, greedy_time):
    ax2.text(x, y - vertical_offset, str(int(y)), ha='center', va='top', color=greedy_color, fontsize=22, zorder=10)

# Adding text annotations for the x-ticks with numbers in parentheses
for i, (np, tp) in enumerate(zip(num_proteins, forbidden_pairs['Total'])):
    ax1.text(np, -0.05, f"({int(tp)})", ha='center', va='top', fontsize=24, rotation=45, transform=ax1.get_xaxis_transform())
    ax2.text(np, -0.05, f"({int(tp)})", ha='center', va='top', fontsize=24, rotation=45, transform=ax2.get_xaxis_transform())


# Set the x-ticks to show the number of proteins
ax1.set_xticks(num_proteins)
ax2.set_xticks(num_proteins)

fig.text(0.03, 0.93, "A", fontsize=30, va='center', ha='center', fontweight='bold')
fig.text(0.52, 0.93, "B", fontsize=30, va='center', ha='center', fontweight='bold')

plt.subplots_adjust(bottom=0.23, wspace=0.16, top=0.9, left=0.05, right=0.95)

plt.savefig("../results/figure3.png",dpi=300)
