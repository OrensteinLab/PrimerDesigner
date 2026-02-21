# Reproducing PCR Efficiency Scores using pcrEfficiency

Some experiments in this project report PCR efficiency scores computed using the external tool **pcrEfficiency**, developed and maintained by its original authors:

https://github.com/maxibor/pcrEfficiency

Because `pcrEfficiency` requires a legacy Python 2 environment, it is not included as a direct dependency of PrimerDesigner.

---

## Step 1 — Clone pcrEfficiency

```bash
git clone https://github.com/maxibor/pcrEfficiency.git
cd pcrEfficiency
```

Create and activate a Python 2.7 environment:

```bash
conda create -n pcrEfficiency_py2 python=2.7.15
conda activate pcrEfficiency_py2
```

Install dependencies according to the pcrEfficiency repository documentation.

---

## Step 2 — Export PrimerDesigner outputs

PrimerDesigner produces required input files under:

```
Comparisons/Results/
```

These include:

- Primer sets:
  - `CVB3_primers.csv`
  - `SpAP_primers.csv`

- Null distribution primers:
  - `null_paths_primers_CVB3.csv`
  - `null_paths_primers_SPAP.csv`

- Reference sequences:
  - `CVB3_reference.fa`
  - `SPAP_reference.fa`

---

## Step 3 — Copy helper scripts

Copy the helper scripts from:

```
Comparisons/pcrEfficiency_scripts/
```

These scripts:

- Convert PrimerDesigner outputs to pcrEfficiency format
- Run efficiency predictions
- Export efficiency scores

Place them inside the `pcrEfficiency` directory.

---

## Step 4 — Run efficiency prediction

To compute PCR efficiency scores for the selected primers and compare PrimerDesigner with competing methods, run the following command:

```bash
python comparison.py   --primers_csv CVB3_primers.csv   --template_fasta CVB3_reference.fa
```
This will compute PCR efficiency scores for all primer pairs listed in CVB3_primers.csv.

Null distribution:

To compute PCR efficiency scores for the random primer designs (null distribution), run:
```bash
python null_dist.py   --null_csv null_paths_primers_CVB3.csv   --template_fasta CVB3_reference.fa   --out_csv null_efficiency_scores.csv
```
This will:
	•	evaluate PCR efficiency for each sampled null primer path
	•	save the results to null_efficiency_scores.csv
  
---

## Notes

- PrimerDesigner runs entirely in Python 3.
- Only the optional efficiency evaluation requires Python 2.
- Results should match those reported in the manuscript when using the same inputs.
