import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

# ======================================================
# Figure Supp. – Tm distribution comparison (Random vs Real)
# ======================================================
tm_df = pd.read_csv('../results/tm_pairs.csv')

tm_random = tm_df['random_tm']
tm_real = tm_df['real_tm']

# Mann–Whitney U test (independent samples)
u_statistic, p_value_tm = mannwhitneyu(
    tm_random, tm_real, alternative='two-sided'
)

print(f"Figure Supp. – Mann–Whitney U test statistic: {u_statistic}")
print(f"Figure Supp. – P-value: {p_value_tm:.4e}")


# ======================================================
# Figure 2B – Runtime comparison (ILP vs Greedy)
# ======================================================
runtime_df = pd.read_csv("../results/PD-mul-ILP.csv")

runtime_ilp = runtime_df['ilp_optimize_time_sec']+runtime_df['ilp_setup_time_sec']
runtime_greedy = runtime_df['greedy_time_sec']

# Wilcoxon signed-rank test (paired samples)
stat_runtime, p_value_runtime = wilcoxon(
    runtime_ilp, runtime_greedy, alternative='greater'
)

print(f"Figure 2B – Wilcoxon Signed-Rank Test (Runtime): Statistic={stat_runtime}, P-value={p_value_runtime:.4e}")


# ======================================================
# Figure 3 – Efficiency and Runtime by Sequence Length
# ======================================================
length_df = pd.read_csv('../results/PD-var-ILP-increasing-lengths.csv')

ilp_primers = length_df['ilp_path_length']
greedy_primers = length_df['greedy_path_length']

seq_lengths = length_df['seq_length']

eff_ilp = length_df['ilp_objective']/ilp_primers
eff_greedy = length_df['greedy_objective']/greedy_primers


runtime_ilp = length_df['ilp_optimize_time_sec'] + length_df['ilp_setup_time_sec']
runtime_greedy = length_df['greedy_time_sec']

# Figure 3A – Efficiency comparison
stat_efficiency, p_value_efficiency = wilcoxon(
    eff_ilp, eff_greedy, alternative='greater'
)
print(f"Figure 3A – Efficiency Comparison: Statistic={stat_efficiency}, P-value={p_value_efficiency:.4e}")

# Figure 3B – Runtime comparison
stat_runtime_length, p_value_runtime_length = wilcoxon(
    runtime_ilp, runtime_greedy, alternative='greater'
)
print(f"Figure 3B – Runtime Comparison: Statistic={stat_runtime_length}, P-value={p_value_runtime_length:.4e}")


# ======================================================
# Figure 3E & 3F – Protein-specific Comparisons
# ======================================================
protein_df = pd.read_csv('../results/PD-var-ILP-different-proteins.csv')

ilp_primers = protein_df['ilp_path_length']
greedy_primers = protein_df['greedy_path_length']

eff_ilp_protein = protein_df['ilp_objective']/ilp_primers
eff_greedy_protein = protein_df['greedy_objective']/greedy_primers

runtime_ilp_protein = protein_df['ilp_optimize_time_sec']+protein_df['ilp_setup_time_sec']
runtime_greedy_protein = protein_df['greedy_time_sec']

# Figure 3E – Efficiency comparison
stat_eff_protein, p_value_eff_protein = wilcoxon(
    eff_ilp_protein, eff_greedy_protein, alternative='greater'
)
print(f"Figure 3E – Protein Efficiency Comparison: Statistic={stat_eff_protein}, P-value={p_value_eff_protein:.4e}")

# Figure 3F – Runtime comparison
stat_runtime_protein, p_value_runtime_protein = wilcoxon(
    runtime_ilp_protein, runtime_greedy_protein, alternative='greater'
)
print(f"Figure 3F – Protein Runtime Comparison: Statistic={stat_runtime_protein}, P-value={p_value_runtime_protein:.4e}")