import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

# ===============================
# Figure 2 – Tm distribution test
# ===============================

tm_df = pd.read_csv('../results/tm_pairs.csv')

random_tm_values = tm_df['random_tm']
real_tm_values = tm_df['real_tm']

# Wilcoxon rank-sum test (Mann–Whitney U test) – Figure 2
tm_statistic, tm_p_value = mannwhitneyu(random_tm_values, real_tm_values, alternative='two-sided')

print(f"Figure 2 – Wilcoxon rank-sum test statistic: {tm_statistic}")
print(f"Figure 2 – P-value: {tm_p_value}")

# ===============================
# Figure 3B – Runtime comparison
# ===============================

runtime_df = pd.read_csv("../results/Non_relaxed_10_proteins.csv")

# Total runtimes
runtime_df['total_ilp_runtime'] = (
    runtime_df['Time (ILP)'] +
    runtime_df['Time (Setup)'] +
    runtime_df['Time (Graph)'] +
    runtime_df['Pair calculation time']
)

runtime_df['total_greedy_runtime'] = (
    runtime_df["Time (Graph)"] +
    runtime_df["Pair calculation time"] +
    runtime_df['Greedy Time']
)

ilp_runtimes = runtime_df['total_ilp_runtime']
greedy_runtimes = runtime_df['total_greedy_runtime']

# Wilcoxon signed-rank test – Figure 3B
runtime_statistic, runtime_p_value = wilcoxon(ilp_runtimes, greedy_runtimes)

print(f"Figure 3B – Wilcoxon Signed-Rank Test: Statistic={runtime_statistic}, P-value={runtime_p_value}")

# =============================================
# Figure 4 – Cost and Runtime by sequence length
# =============================================

length_df = pd.read_csv('../results/PD_mul_var_different_lengths.csv')

sequence_lengths = length_df['seq length:']
ilp_costs = length_df['ILP Objective']
greedy_costs = length_df['Greedy Objective']

length_df['total_ilp_runtime'] = (
    length_df['Time (ILP)'] +
    length_df['Time (Setup)'] +
    length_df['Time (Graph)']
)

ilp_runtimes = length_df['total_ilp_runtime']
greedy_runtimes = length_df['Greedy Time'] + length_df["Time (Graph)"]

# Figure 4A – Cost comparison
cost_statistic, cost_p_value = wilcoxon(ilp_costs, greedy_costs)
print(f"Figure 4A – Cost Comparison: Statistic={cost_statistic}, P-value={cost_p_value}")

# Figure 4B – Runtime comparison
runtime_statistic, runtime_p_value = wilcoxon(ilp_runtimes, greedy_runtimes)
print(f"Figure 4B – Runtime Comparison: Statistic={runtime_statistic}, P-value={runtime_p_value}")

# ===================================
# Figure 4E & 4F – Protein comparisons
# ===================================

protein_df = pd.read_csv('../results/PD_mul-var_10_proteins.csv')

# Figure 4E – Protein cost comparison
ilp_protein_costs = protein_df['PrimerDesigner cost']
greedy_protein_costs = protein_df['Greedy cost']

cost_statistic_protein, cost_p_value_protein = wilcoxon(
    ilp_protein_costs, greedy_protein_costs, zero_method='zsplit'
)
print(f"Figure 4E – Cost Comparison: Statistic={cost_statistic_protein}, P-value={cost_p_value_protein}")

# Figure 4F – Protein runtime comparison
ilp_protein_runtimes = protein_df['PrimerDesigner runtime [s]']
greedy_protein_runtimes = protein_df['Greedy runtime [s]']

runtime_statistic_protein, runtime_p_value_protein = wilcoxon(
    ilp_protein_runtimes, greedy_protein_runtimes, zero_method='zsplit'
)
print(f"Figure 4F – Runtime Comparison: Statistic={runtime_statistic_protein}, P-value={runtime_p_value_protein}")
