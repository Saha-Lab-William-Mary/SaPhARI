# SaPhARI parameter determination, validation, and downstream search workflow

This directory contains the workflow used to derive, validate, score, and apply SaPhARI family parameters from a curated set of satellite genomes.

The curated reference dataset used for rule generation consisted of known satellite genomes that were first annotated with SaPhARI against a DIAMOND database built from the RefSeq archaeal, bacterial, and viral protein snapshot dated **November 18, 2025**. These annotated reference genomes were then used to derive family-level rules, validate the resulting parameter set, and score recovery performance. The finalized rule set was subsequently applied to a larger collection of new genomes for downstream satellite discovery.

## Directory structure

The GitHub repository preserves the workflow scripts, archived reference data, accession list, and key output products needed to document the analysis:

```text
saphari_pub/
├── all_accessions.txt
├── annotated_satellites.zip
├── compile_flagged_results.py
├── extract_clusters.py
├── score.py
├── output/
│   ├── compiled_flags/
│   ├── extraction/
│   │   ├── Families_final.json
│   │   ├── family_rules_best.tsv
│   │   ├── family_rules_readable.txt
│   │   ├── cluster_backmap.tsv
│   │   └── pairwise_confusion.tsv
│   ├── test_json/
│   └── new_genomes_fasta/
│       ├── 001/
│       │   └── RESULTS/
│       ├── 002/
│       │   └── RESULTS/
│       ├── 003/
│       │   └── RESULTS/
│       ├── 004/
│       │   └── RESULTS/
│       ├── 005/
│       │   └── RESULTS/
│       ├── 006/
│       │   └── RESULTS/
│       ├── 007/
│       │   └── RESULTS/
│       ├── 008/
│       │   └── RESULTS/
│       └── 009/
│           └── RESULTS/
└── slurm/
```

The full local workflow also included large input directories that are not retained in full in GitHub because of file size:

```text
annotated/         curated annotated reference genomes used for rule generation and validation
assembled_genomes/ larger genome collection used for downstream satellite searches
```

In the public repository, the curated reference set is provided as `annotated_satellites.zip`, while the downstream new-genome search is represented by `all_accessions.txt` together with the retained `output/new_genomes_fasta/*/RESULTS/` directories.

## Workflow summary

### 1. Annotated reference genomes

`annotated_satellites.zip` contains the curated satellite genome dataset used for rule generation and validation. In the original local workflow, these genomes were organized by family and stored as DIAMOND-filtered TSV annotation files produced by SaPhARI. These annotated family-organized inputs served as the reference dataset for parameter extraction.

### 2. Parameter extraction

`extract_clusters.py` is used to infer candidate family rules from the annotated reference set. This step generates a SaPhARI-compatible `Families_final.json` file together with supporting diagnostic outputs written to `output/extraction/`. The extraction outputs include:

* `Families_final.json`
* `family_rules_best.tsv`
* `family_rules_readable.txt`
* `cluster_backmap.tsv`
* `pairwise_confusion.tsv`

### 3. Validation

The extracted `Families_final.json` file is then used to run SaPhARI back on the same curated annotated reference set as a structured validation step. These validation runs are stored in `output/test_json/` as per-family output directories containing SaPhARI `RESULTS/` outputs.

### 4. Scoring

`score.py` is used to score the SaPhARI validation output. This step generates per-file results, family-level summaries, and overall summary statistics, providing the quantitative evaluation of the derived family rules.

### 5. Downstream search on new genomes

After validation and scoring, the finalized family rule set is applied to a larger collection of new genomes for downstream satellite discovery. The accessions for these genomes are listed in `all_accessions.txt`. Downstream SaPhARI runs are stored in batched directories under `output/new_genomes_fasta/`, with each batch retaining its corresponding `RESULTS/` folder.

## Database snapshot

The DIAMOND database used for SaPhARI annotation in this workflow was `SaPhARI_11_18_25.dmnd`, corresponding to a RefSeq archaeal, bacterial, and viral protein snapshot formatted for DIAMOND on **November 18, 2025**.

## Reproducibility

The `slurm/` directory contains the batch scripts used to run parameter extraction, validation, scoring, and downstream genome searches on the cluster. Together with the archived reference dataset, the retained workflow outputs, and the specified DIAMOND database snapshot, these scripts document how the study was executed and provide the structure needed to recreate the analysis when the corresponding input datasets are available.
