
# PrimerDesigner

## Introduction

This repository contains the source code for PrimerDesigner, a tool designed to find the most efficient primer set with complete coverage and no cross hybridization risk for protein synthesis using PCR. It accompanies the paper titled "PrimerDesigner: Designing efficient primers for synthesizing large protein libraries without cross-hybridization".

![alt text](https://github.com/OrensteinLab/PrimerDesigner/blob/main/docs/primer_design_illustration.png)

**License:** MIT (see [LICENSE](LICENSE)).  
**Contact:** jonathan.mandl2@gmail.com

---

## Creating the Conda environment

Clone the repository, then recreate the environment from `environment.yml`:

```bash
git clone https://github.com/OrensteinLab/PrimerDesigner.git
cd PrimerDesigner
conda env create -f environment.yml
conda activate primer_env
```

Run all later commands from this repository root with `primer_env` active.

`gurobipy` is installed with the environment. A Gurobi **license file is only required** for `PD-var-ILP` and `PD-mul-ILP`. `PD-single-LPath` and `PD-mul-Greedy` run without `gurobi.json`.

---

## Quickstart (no Gurobi)

Run **PD-single-LPath** on the SpAP coding sequence included in the repository (`data/SPAP_reference.txt`). 

```bash
python ./tool.py \
  --version PD-single-LPath \
  --file_path data/SPAP_reference.txt \
  --config config.json \
  --output quickstart_output
```

---

## Gurobi license (ILP versions only)

1. Request a free academic Web License Service (WLS) key:  
   https://www.gurobi.com/academia/academic-program-and-licenses/
2. Copy the template and fill in your credentials:

```bash
cp gurobi.json.example gurobi.json
```

```json
{
  "WLSACCESSID": "XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX",
  "WLSSECRET": "XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX",
  "LICENSEID": 12345
}
```

`gurobi.json` is gitignored — do not commit it.

---

## Setting up your own sequences

Create a text file of protein names and DNA coding sequences (**one protein per line, tab-separated**):

```text
SHP2	ATGACATCGCGGAGATGGTTTCACCCAAATATCACTGGTGTGGAGGCAGAAAACCTACTGTTGACAAGAGGAGT....
CXAR	ATGGCGCTCCTGCTGTGCTTCGTGCTCCTGTGCGGAGTAGTGGATTTCGCCAGAAGTTTGAGTATCACTACTCC....
```

Variant-specific formats are described in the version documentation below.

Create a `config.json` (or reuse the repo file) with constant upstream and downstream flanks added around each mutagenized coding sequence:

```json
{
  "upstream_nt": "GCTAGTGGTGCTAGCCCCGCGAAATTAATACGACTCACTATAGGGTCTAGA....",
  "downstream_nt": "GGAGGGTCTGGGGGAGGAGGCAGTGGCATGGTGAGCAAGGGCGAGGAGC...."
}
```

The checked-in `config.json` also sets `max_tm` (heterodimer Tm cutoff used when detecting cross-hybridization; default 45 °C).

---

## Versions

| Application | Version | Solver |
|---|---|---|
| One protein-coding sequence | **PD-single-LPath** | None (graph longest path) |
| Several variants of the same protein | **PD-var-ILP** | Gurobi |
| Several non-homologous proteins (fast) | **PD-mul-Greedy** | None |
| Several non-homologous proteins (optimal) | **PD-mul-ILP** | Gurobi |

Input and output formats for each version:

- [docs/PD-single-LPath.md](docs/PD-single-LPath.md)
- [docs/PD-var-ILP.md](docs/PD-var-ILP.md)
- [docs/PD-mul-Greedy.md](docs/PD-mul-Greedy.md)
- [docs/PD-mul-ILP.md](docs/PD-mul-ILP.md)

---

## Running PrimerDesigner

```bash
python ./tool.py --version <version> --file_path <file-path> --config config.json --output <output-folder>
```

- **file_path**: protein coding-sequence file  
- **version**: `PD-single-LPath`, `PD-mul-ILP`, `PD-mul-Greedy`, or `PD-var-ILP` (default: `PD-single-LPath`)  
- **output**: folder for CSV outputs  
- **config**: JSON with upstream and downstream regions  

### Parameters

- **primer_lmin**, **primer_lmax**  
  Minimum and maximum primer lengths.  
  *Default:* 18, 30  

- **oligo_lmin**, **oligo_lmax**  
  Minimum and maximum oligonucleotide lengths.  
  *Default:* 195, 205  

- **overlap_lmin**, **overlap_lmax**  
  Minimum and maximum overlap length between oligonucleotides.  
  *Default:* 45, 50  

- **allowed_overlap**  
  Allowed overlap between primer pairs.  
  *Default:* 6  

- **num_proteins**  
  Number of variants of the same sequence.  
  *Used for:* `PD-var-ILP` only.  
  *Default:* 3  

- **apply_threshold**  
  Boolean flag for applying primer quality threshold.  
  *Default:* False  

- **min_gc**, **max_gc**  
  Minimum and maximum GC content threshold (in %).  
  *Default:* 40, 60  

- **min_tm**, **max_tm**  
  Minimum and maximum melting temperature (Tm) thresholds (in °C).  
  *Default:* 58, 65  

- **max_difference**  
  Maximum allowed Tm difference between forward and reverse primers in a pair.  
  *Default:* 3  

- **merge_bins**  
  Boolean flag for merging bins corresponding to identical non-overlapping sequences.  
  *Used for:* `PD-var-ILP`.  
  *Default:* False  

---

## Reproducing paper experiments

Full commands, expected output files, and runtimes:

**[docs/Reproducing_experiments.md](docs/Reproducing_experiments.md)**

Also:

- Comparison workflow (short): [docs/Running_comparisons.md](docs/Running_comparisons.md)
- PCR efficiency with the external pcrEfficiency model: [docs/pcr_Efficiency.md](docs/pcr_Efficiency.md)

Committed result tables (for checking a rerun) are in `Experiments/Results/` and `Comparisons/Results/`.
