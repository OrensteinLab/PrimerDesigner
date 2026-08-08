# Reproducing PrimerDesigner experiments

Run every command from the **repository root** after activating `primer_env`.

```bash
git clone https://github.com/OrensteinLab/PrimerDesigner.git
cd PrimerDesigner
conda activate primer_env
```

If the environment does not exist yet, create it first (`conda env create -f environment.yml`).

Committed result tables (to compare against a rerun):

- `Experiments/Results/` — scaling / ILP vs greedy CSVs
- `Comparisons/Results/` — SpAP and CVB3 method comparisons

Experiment scripts write new files to **`./Results/`** (repository root), not into `Experiments/Results/`. Compare your new CSVs to the committed copies.

---

## 1. Which license do you need?

| Analysis | Script(s) | Gurobi `gurobi.json` |
|---|---|---|
| SpAP / CVB3 comparisons | `Comparisons/compare_SpAP.py`, `Comparisons/compare_CVB3.py` | No |
| 100-CDS greedy scaling | `Experiments/mul_greedy_100_cds.py` | No |
| Greedy vs increasing protein count | `Experiments/mul_greedy_ILP_comparison.py` | No |
| PD-mul-ILP scaling | `Experiments/mul_ILP.py` | **Yes** |
| PD-var-ILP (length / variants / proteins) | `Experiments/var_ILP_*.py` | **Yes** |
| PCR-efficiency scores (optional) | `docs/pcr_Efficiency.md` | No (separate Python 2 env) |

Copy `gurobi.json.example` to `gurobi.json` and fill in a [Gurobi WLS](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) academic license.

ILP jobs are large (multi-GB RAM; `PD-mul-ILP` on 10 proteins can take hours).

---

## 2. Method comparisons (SpAP and CVB3)

Competing primer sets are already in `Comparisons/Primers/`. You do **not** need to reinstall Olivar / PrimalScheme to reproduce the comparison tables.

```bash
python -m Comparisons.compare_SpAP
python -m Comparisons.compare_CVB3
```

**Writes** (under `Comparisons/Results/`):

| File | Contents |
|---|---|
| `SpAP_comparison.csv` | Runtime, memory, mean efficiency vs QuickChange / PrimalScheme / Olivar |
| `SPAP_primers.csv` | Primer sequences per method |
| `null_paths_primers_SPAP.csv` | 1000 random graph paths (null distribution) |
| `CVB3_comparison.csv` | Same metrics vs the published CVB3 primer set |
| `CVB3_primers.csv` | Primer sequences |
| `null_paths_primers_CVB3.csv` | Null paths |

Configs: `configs/SPAP_experiment.json`, `configs/CVB3_experiment.json`.

Optional PCR-efficiency rescoring with the external **pcrEfficiency** tool: [pcr_Efficiency.md](pcr_Efficiency.md).  
Helper scripts: `Comparisons/pcrEffficiency_scripts/`.

---

## 3. Scaling experiments

### 3.1 PD-mul-Greedy on 100 CCDS proteins

```bash
python -m Experiments.mul_greedy_100_cds
```

- Input: `data/100_ccds_protein_sequences.txt`
- Config: `configs/SPAP_experiment.json`
- Output: `Results/PD_mul_Greedy_summary.csv` (and per-protein metrics). Committed reference: `Experiments/Results/PD-mul-Greedy.csv`
- Runtime: ~1 hour in the committed run (~4400 s)

### 3.2 PD-mul-Greedy, increasing protein count

```bash
python -m Experiments.mul_greedy_ILP_comparison
```

- Input: first 2…10 sequences in `data/10_protein_coding_sequences.txt`
- Output: `Results/PD-mul-Greedy_ILP_comparison.csv`

### 3.3 PD-mul-ILP scaling

Requires Gurobi. High memory (committed peak ~10 GB at 10 proteins).

```bash
python -m Experiments.mul_ILP
```

- Output: `Results/PD-mul-ILP.csv`
- Committed reference: `Experiments/Results/PD-mul-ILP.csv` (that copy also has extra greedy columns from a combined run)

### 3.4 PD-var-ILP

Requires Gurobi. Uses the SpAP coding sequence unless noted.

```bash
# Increasing coding-sequence length (226–1626 nt, step 100)
python -m Experiments.var_ILP_length

# Increasing number of variants (2–8)
python -m Experiments.var_ILP_num

# One run per protein in the 10-protein set
python -m Experiments.var_ILP_proteins
```

| Script | Output CSV | Committed copy |
|---|---|---|
| `var_ILP_length.py` | `Results/PD-var-ILP-increasing-lengths.csv` | `Experiments/Results/PD-var-ILP-increasing-lengths.csv` |
| `var_ILP_num.py` | `Results/PD-var-ILP-increasing_variants.csv` | `Experiments/Results/PD-var-ILP-increasing_variants.csv` |
| `var_ILP_proteins.py` | `Results/PD-var-ILP-different-proteins.csv` | `Experiments/Results/PD-var-ILP-different-proteins.csv` |

`var_ILP_proteins.py` uses root `config.json`. Length and variant sweeps use `configs/SPAP_experiment.json`.

---

## 4. Data files used in the paper

| Path | Role |
|---|---|
| `data/SPAP_reference.fa` / `.txt` | SpAP coding sequence |
| `data/CVB3_reference.fa` / `.txt` | CVB3 coding sequence |
| `data/10_protein_coding_sequences.txt` | 10 non-homologous CDSs (`name<TAB>sequence`) |
| `data/100_ccds_protein_sequences.txt` | 100 CCDS proteins for greedy scaling |
| `configs/SPAP_experiment.json` | Upstream/downstream flanks for SpAP experiments |
| `configs/CVB3_experiment.json` | Flanks for CVB3 |
| `Comparisons/Primers/` | Published / competitor primer BED/CSV inputs |

Input sequence files are **one protein per line**, name and DNA sequence separated by a **tab**.
