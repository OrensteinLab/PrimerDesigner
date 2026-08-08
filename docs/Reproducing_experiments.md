# Reproducing PrimerDesigner experiments

Run every command from the **repository root** with `primer_env` active.

```bash
git clone https://github.com/OrensteinLab/PrimerDesigner.git
cd PrimerDesigner
conda activate primer_env
```

If the environment does not exist yet: `conda env create -f environment.yml`.

ILP experiments require `gurobi.json` (see the main README). Scripts write new files to **`./Results/`**. Committed copies for comparison are in `Experiments/Results/`.

Method comparisons on SpAP and CVB3 are documented separately: [Running_comparisons.md](Running_comparisons.md). Optional PCR-efficiency scoring: [pcr_Efficiency.md](pcr_Efficiency.md).

---

## PD-mul-Greedy on 100 CCDS proteins

Scales the greedy multi-protein designer to 100 CCDS sequences and records total efficiency, runtime, and cross-hybridization retries.

- **Data:** 100 human CCDS coding sequences (`data/100_ccds_protein_sequences.txt`) with SpAP upstream/downstream flanks (`configs/SPAP_experiment.json`).
- **Output:** `Results/PD_mul_Greedy_summary.csv` (and per-protein metrics). Reference: `Experiments/Results/PD-mul-Greedy.csv`

```bash
python -m Experiments.mul_greedy_100_cds
```

---

## PD-mul-Greedy with increasing protein count

Runs greedy design on the first 2, 3, …, 10 non-homologous coding sequences to measure how cost grows with the number of proteins.

- **Data:** Ten non-homologous protein-coding sequences (`data/10_protein_coding_sequences.txt`); the script uses the first *k* sequences for *k* = 2…10, with SpAP flanks (`configs/SPAP_experiment.json`).

```bash
python -m Experiments.mul_greedy_ILP_comparison
```

---

## PD-mul-ILP scaling

Solves the multi-protein ILP on the same 1…10 sequence prefixes. Requires Gurobi. Memory is high (committed peak ~10 GB at 10 proteins).

- **Data:** The same ten non-homologous coding sequences (`data/10_protein_coding_sequences.txt`), prefixes of length 1…10, with flanks from the root `config.json`.

```bash
python -m Experiments.mul_ILP
```

---

## PD-var-ILP: increasing sequence length

Designs primers for 3 variants of SpAP while lengthening the coding sequence (226–1626 nt, step 100). Requires Gurobi.

- **Data:** The SpAP coding sequence (`data/SPAP_reference.fa`), truncated to each target length, with SpAP flanks (`configs/SPAP_experiment.json`).

```bash
python -m Experiments.var_ILP_length
```

---

## PD-var-ILP: increasing number of variants

Holds the full SpAP coding sequence fixed and increases the number of variants from 2 to 8. Requires Gurobi.

- **Data:** The full SpAP coding sequence (`data/SPAP_reference.fa`) with SpAP flanks (`configs/SPAP_experiment.json`); only the requested variant count changes.

```bash
python -m Experiments.var_ILP_num
```

---

## PD-var-ILP: different proteins

Runs the variant ILP once per sequence in the 10-protein set. Requires Gurobi.

- **Data:** Each of the ten non-homologous coding sequences in turn (`data/10_protein_coding_sequences.txt`), with flanks from the root `config.json`.
- **Output:** `Results/PD-var-ILP-different-proteins.csv`. Reference: `Experiments/Results/PD-var-ILP-different-proteins.csv`

```bash
python -m Experiments.var_ILP_proteins
```
