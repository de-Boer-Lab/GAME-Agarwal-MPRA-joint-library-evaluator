# GAME Agarwal MPRA Joint Library Evaluator

Agarwal et al. (https://www.nature.com/articles/s41586-024-08430-9) MPRA evaluator first requests log-scale predictions for ~56k probes across WTC11, K562, and HepG2 cells. It then calculates Pearson correlations between predicted and measured $(\log_2(rna/dna))$. It also assesses cell-type specificity by correlating predicted vs. observed expression differences.

---

## How It Works

The evaluator:

1. Loads the measured MPRA sequences from `56k_measured_final_sequence_file.xlsx` and the plasmid backbone from `upstream_downstream_backbone.txt`
2. Creates JSON to send to Predictor
3. Negotiates serialization format (JSON or MsgPack) with the Predictor
4. POSTs all sequences to the predictor's `/predict` endpoint in a single request, across three prediction tasks (WTC11, K562, HepG2)
5. Computes Pearson *r* between predicted expression (log scale) and measured log₂(RNA/DNA) values for each cell type
6. Computes cell-type specificity metrics across all three cell-type pairs
7. Saves raw predictions and a final evaluation summary CSV

Additional details about the Evaluator's data and the sequence file used can be found in the evaluator_data folder.
---

### Run Evaluator using Apptainer container

Download the container and evaluator data from Zenodo: `<LINK>`

```bash
apptainer run --containall \
  -B /absolute/path/to/evaluator_data:/evaluator_data \
  -B /absolute/path/to/predictions:/predictions \
  agarwal_evaluator.sif <predictor_ip> <predictor_port> /predictions
```

---

## Arguments

| Argument | Description |
|---|---|
| `predictor_ip` | IP address or hostname of the predictor server |
| `predictor_port` | Port the predictor is listening on |
| `output_dir` | Directory where prediction JSONs and metric CSVs are written |

---

## Request Structure

A single POST request is sent containing all ~56,000 sequences. The backbone upstream/downstream flanking sequences and the promoter coordinates (used as `prediction_ranges`) are read from `upstream_downstream_backbone.txt`.

Each request payload follows this structure:

```json
{
  "readout": "point",
  "prediction_tasks": [
    {
      "name": "agarwal_joint_lib_wtc11",
      "type": "expression",
      "cell_type": "induced pluripotent stem cells (iPS cells; WTC11)",
      "scale": "log",
      "species": "homo_sapiens"
    },
    {
      "name": "agarwal_joint_lib_k562",
      "type": "expression",
      "cell_type": "lymphoblasts (K562)",
      "scale": "log",
      "species": "homo_sapiens"
    },
    {
      "name": "agarwal_joint_lib_hepg2",
      "type": "expression",
      "cell_type": "human hepatocytes (HepG2)",
      "scale": "log",
      "species": "homo_sapiens"
    }
  ],
  "upstream_seq": "<plasmid upstream flanking sequence>",
  "downstream_seq": "<plasmid downstream flanking sequence>",
  "sequences": {
    "<seq_id>": "<200 nt sequence>",
    ...
  },
  "prediction_ranges": {
    "<seq_id>": [start, end],
    ...
  }
}
```

**Sequence structure:** Each sequence is 200nt long and corresponds to the MPRA probes. The `prediction_ranges` coordinates correspond to the promoter region within the backbone, as defined in `upstream_downstream_backbone.txt`, and are shared across all sequences.

---

## Outputs

All outputs are written to `output_dir`:

| File | Description |
|---|---|
| `agarwal_2025_joint_lib_<timestamp>_predictions_from_<predictor_name>.json` | Raw predictions returned by the predictor |
| `evaluation_summary_agarwal_2025_joint_lib_<timestamp>.csv` | Pearson *r* per cell type and cell-type specificity metrics, appended across runs |

**Metrics computed:**

The measurements for each probe are stored in `56k_measured_file.xlsx` in the `evaluator_data` directory.

- **Pearson *r*** — predicted vs. measured log₂(RNA/DNA) for each of the three cell types (WTC11, K562, HepG2)
- **Cell-type specificity** — Pearson *r* of predicted vs. measured differential expression across all three cell-type pairs:
  - HepG2 vs. WTC11
  - HepG2 vs. K562
  - WTC11 vs. K562

---

## Directory Structure

```
Agarwal_2025/
├── evaluator_RestAPI.py                # Main entry point
├── config.py                           # Settings (name, paths, formats, retries)
├── data_loader.py                      # Loads Excel + backbone, builds request payload
├── evaluator_content_handler.py        # Format negotiation, HTTP POST, deserialization
├── evaluator_metrics_calculator.py     # Pearson r, cell-type specificity, and CSV output
├── agarwal_evaluator.def               # Apptainer container definition
└── evaluator_data/
    └── 2023-03-03628-s5/
        ├── 56k_measured_final_sequence_file.xlsx
        └── upstream_downstream_backbone.txt
```

---
