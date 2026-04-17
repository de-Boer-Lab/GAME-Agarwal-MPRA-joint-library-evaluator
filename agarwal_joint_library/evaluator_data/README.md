# Agarwal 2025 Joint Library — Evaluator Data

This directory contains the input data files required to run the Agarwal 2025 Joint Library MPRA Evaluator and the code to create them. The files were derived from the supplementary data of Agarwal et al. 2025 (https://www.nature.com/articles/s41586-024-08430-9) and processed into the formats expected by the evaluator scripts.


## File Descriptions

### `56k_measured_final_sequence_file.xlsx`
The primary input to the evaluator. Contains ~56,000 regulatory element sequences (200 nt each) for which complete expression measurements exist across all three cell types (WTC11, K562, HepG2).

| Column | Description |
|---|---|
| `name` | Unique sequence identifier |
| `sequence_trimmed` | 200 nt regulatory element (15 nt adaptors removed from both ends) |
| *(other columns)* | Additional metadata carried over from the original design file |

### `56k_measured_file.xlsx`
Contains the measured expression values for the sequences in `56k_measured_final_sequence_file.xlsx`. Used by `evaluator_metrics_calculator.py` to compute Pearson *r* between predicted and measured expression.

| Column | Description |
|---|---|
| `name` | Unique sequence identifier (matches `56k_measured_final_sequence_file.xlsx`) |
| `WTC11 [log2(rna/dna)]` | Measured expression in WTC11 iPSCs |
| `K562 [log2(rna/dna)]` | Measured expression in K562 lymphoblasts |
| `HepG2 [log2(rna/dna)]` | Measured expression in HepG2 hepatocytes |

Only sequences with complete measurements across **all three cell types** are retained (rows with any NaN values are dropped).

### `upstream_downstream_backbone.txt`
A tab-separated file containing the plasmid backbone sequences that flank each regulatory element, and the promoter coordinates used as `prediction_ranges` in the request payload. Read by `data_loader.py` at runtime.

| Row | Description |
|---|---|
| `upstream_padded` | Upstream flanking sequence prepended to each element |
| `downstream` | Downstream flanking sequence appended to each element |
| `promoter_coordinates` | JSON-encoded `[start, end]` coordinates of the promoter within the full construct, shared across all sequences |

---

## How the Files Were Generated

The script `make_evaluator_data.py` is used to prepare the Evaluator's data.

1. **Download raw data** from the paper supplement:
   - `2023-03-03628C-Table S11 - jointdata.xlsx` — measured expression values across all cell types
   - `2023-03-03628C-Table S10 - joint lib design.xlsx` — full joint library sequence design file

Note: `2023-03-03628C-Table S10 - joint lib design.xlsx` table was filtered to only keep probes with associated measured values in `2023-03-03628C-Table S11 - jointdata.xlsx`. 

2. **Filter to complete measurements** — rows in the measured values file with any missing cell-type measurement are dropped, ensuring all three cell types are present for every retained sequence.

3. **Filter sequences** — only sequences present in the filtered measured values file are kept.

4. **Trim adaptors** — the 15 nt 5′ and 15 nt 3′ adaptors are clipped from each 230 nt sequence, yielding the 200 nt regulatory element. All trimmed sequences are verified to be exactly 200 nt (the adapter sequences are included in the backbone).

5. **Save outputs** — the processed sequences and measurements are saved as `56k_measured_final_sequence_file.xlsx` and `56k_measured_file.xlsx` respectively.

## Creating the plasmid backbone

The `upstream_downstream_backbone.txt` file was generated from the full lentiviral MPRA plasmid sequence using the script `make_backbone.py`. It provides the biological flanking context that surrounds each regulatory element in the actual MPRA assay, and defines the coordinates of the reporter gene used as the `prediction_ranges` in each request.

The input is `Lenti_mpra.txt` — the full sequence of the lentiviral MPRA vector. In the plasmid, the regulatory element cloning site is represented as a run of 200+ ambiguous bases (`N`s), which serves as the delimiter between the upstream and downstream flanking sequences.

### How `upstream_downstream_backbone.txt` was generated

1. **Split on the N-run** — the plasmid sequence is split at the first run of 200 or more `N` bases. Everything before the N-run becomes `upstream_padded`; everything after becomes `downstream`. The 200bp Ns are a placeholder for any probe. 

2. **Locate the TSS** — the transcription start site of the reporter gene (EGFP) is identified by searching for the sequence `ATGGTGAGCAAGGG` (the start codon of the fluorescent reporter). Its position within the full plasmid sequence is recorded as `start_tss`.

3. **Define promoter coordinates** — the reporter gene body is defined as a 720 bp window starting at `start_tss`: `[start_tss, start_tss + 720]`. These coordinates are stored as `promoter_coordinates` and are used as the shared `prediction_ranges` for all sequences in the evaluator request.

### Output format

The resulting file is a tab-separated table with one row per entry:

| Row | `sequence` column |
|---|---|
| `upstream_padded` | Full upstream flanking sequence (plasmid sequence before the N-run) |
| `downstream` | Full downstream flanking sequence (plasmid sequence after the N-run) |
| `promoter_coordinates` | JSON-encoded `[start, end]` coordinates of the reporter gene (e.g. `[start_tss, start_tss + 720]`) |
