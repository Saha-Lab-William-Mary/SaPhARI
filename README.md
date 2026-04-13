# SaPhARI — Satellite Phage Algorithmic Recognition and Interpretation

**SaPhARI** (Satellite Phage Algorithmic Recognition and Interpretation) is a flexible bioinformatics tool for identifying satellite phages in genomic and metagenomic data. Rather than relying solely on sequence conservation, SaPhARI detects co-localized clusters of hallmark genes based on annotated function and positional context.

The pipeline begins with ORF prediction and annotation, followed by a fast homology search using DIAMOND or BLASTn. Detected proteins are passed through a rule-based filter where users can define required, synonym, and forbidden proteins, along with intergenic distance, span, and minimum protein count constraints. These customizable rules enable both targeted searches and exploratory discovery of novel or lineage-specific satellite families, even when involving hypothetical proteins.

SaPhARI produces easily-readable TSV reports and is distributed as a fully containerized package (Docker/Apptainer/Singularity) for reproducibility across local, cloud, and HPC environments.

> Developed for the 2024 William & Mary iGEM Team, SaPhARI was nominated for **Best Software Tool** at the iGEM Grand Jamboree.  
> Project page: [https://2024.igem.wiki/william-and-mary/software](https://2024.igem.wiki/william-and-mary/software)

## Installation

SaPhARI can be installed in multiple ways depending on your environment:

1.  **Docker (recommended):** Build the container image using the included `Dockerfile`.
    
2.  **Apptainer/Singularity:** Convert the Docker image to an `.sif` file for HPC clusters.
    
3.  **Conda environment:** Use the provided `environment.linux.yml` file to create a conda environment with all dependencies.
    
4.  **Source install (developers):** Clone the repository and install with pip (`pip install -e .`).
    

### Example (Docker)

```bash
git clone https://github.com/Saha-Lab-William-Mary/SaPhARI.git
cd SaPhARI
docker build -t saphari:0.1.0 -f containers/Dockerfile .
```

### Example (Conda)

```bash
git clone https://github.com/Saha-Lab-William-Mary/SaPhARI.git
cd SaPhARI
conda env create -f environment.linux.yml
conda activate saphari
saphari --help
```

## Quick Start

```bash
# Run an example
mkdir -p out

docker run --rm \
  -u $(id -u):$(id -g) \
  -v "$PWD":/work -v "$PWD/out":/out -w /work \
  saphari:0.1.0 \
  saphari --workflow n_t_a \
    --data_dir tests/ex1/genomes \
    --db tests/ex1/test.dmnd \
    --families_file src/saphari/data/Families.json \
    --out_dir /out/run1 \
    --cpus 4 --memory "4 GB"

```

## Key Features

-   **Annotation + Family Search:** ORF prediction → DIAMOND/BLASTn → rule-based family classification.
    
-   **Flexible Rules:** Required, forbidden, and synonym proteins; max intergenic distances; cluster span; minimum match count across proteins.
    
-   **Multiple Modes:** Full pipeline (annotation+search) or reuse existing CDS/TSVs for reanalysis.
    
-   **Containerized & Reproducible:** Identical runs on laptops, servers, or HPC.
    
-   **HPC-Ready:** Scales via Nextflow.
    

## CLI Usage

```text
usage: saphari [-h] [--code_dir CODE_DIR]
               [--workflow {n_t_a,m_n_t_a,a_t_a,n_b_a,m_n_b_a}] [--CDS]
               [--reanalyze] --db DB [--ndb NDB]
               [--families_file FAMILIES_FILE] [--search SEARCH [SEARCH ...]]
               [--output_type {pfpf,pf}] --data_dir DATA_DIR --out_dir OUT_DIR
               [--dim_e DIM_E] [--dim_m DIM_M] [--dim_c DIM_C] [--dim_p DIM_P]
               [--nb_e NB_E] [--nb_m NB_M] [--nb_c NB_C] [--nb_p NB_P]
               [--cpus CPUS] [--memory MEMORY] [--profile PROFILE]
               [--container CONTAINER]

```


#### Annotation Pipeline
Options controlling the initial annotation step and database selection.
- `--workflow` (`-w`): Selects the Nextflow workflow to execute (`n_t_a`, `m_n_t_a`, `a_t_a`, `n_b_a`, `m_n_b_a`). Required unless `--CDS` or `--reanalyze` is specified.  
- `--CDS`: Skips annotation and treats input files as GenBank CDS (protein sequences).  
- `--reanalyze`: Re-runs the family search using existing annotated TSV files without repeating annotation.  
- `--db` (`-d`): Path to DIAMOND protein database (`.dmnd`) 
- `--ndb` (`-n`): Path to nucleotide BLAST database, used only for BLASTn-based workflows.  
- `--code_dir` (`-c`): Path to the SaPhARI source directory containing `annotation/main.nf`. Defaults to the installed package location. (developer focused)


#### Family Search
Options defining which families to detect and the output structure.
- `--families_file` (`-F`): Path to a custom `Families.json`. If not provided, the packaged default is used.  
- `--search` (`-f`): One or more family names to restrict the search. Defaults to all families in the JSON file.  
- `--output_type` (`-t`): Output format:  
  - `pfpf`: per-file-per-family (default)  
  - `pf`: per-family aggregated across all files.


#### Input / Output
Options specifying input data location and output directory.
- `--data_dir` (`-i`): Required. Directory containing input genomes (FASTA/GenBank) or annotated TSVs if `--reanalyze` is used.  
- `--out_dir` (`-o`): Required. Output directory. Subdirectories `annotated/`, `RESULTS/`, and ... will be created.


#### Annotation Tuning
Parameters for homology search sensitivity and thresholds.
- **DIAMOND:** `--dim_e`, `--dim_m`, `--dim_c`, `--dim_p`  
- **BLASTn:** `--nb_e`, `--nb_m`, `--nb_c`, `--nb_p`  

These correspond to e-value cutoffs, minimum matches, coverage, and percent identity. Defaults are provided; values can be adjusted to balance sensitivity and specificity.

#### Runtime / Resources
Options controlling computational resources and execution environment.
- `--cpus`: Number of CPU cores allocated (default: 4).  
- `--memory`: Memory limit (default: `8 GB`).  
- `--profile`: Nextflow profile to use (default: `local`).  
- `--container`: Docker or Singularity/Apptainer container image to execute within (optional).



## Modes and Workflows

### Modes

-   **Default (full pipeline):** Annotates FASTA input then performs family search.
    
-   `--CDS`: Use GenBank CDS (AA). Skips annotation.
    
-   `--reanalyze`: Reuse annotated TSVs. Skips annotation.
    

### Workflows


Workflow   | Input Type          | Typical Length      | Method
-----------|---------------------|----------------------|------------------------------------------------------
n_t_a      | Nucleotide FASTA    | > 20 kb              | Prodigal ORFs → DIAMOND search (general use)
m_n_t_a    | Nucleotide FASTA    | 5–20 kb              | Prodigal (meta mode) → DIAMOND (short/fragmented contigs)
a_t_a      | Protein FASTA (CDS) | n/a                  | Use provided Amino Acid sequences → DIAMOND (skips ORF prediction)
n_b_a      | Nucleotide FASTA    | > 20 kb              | Direct BLASTn search (useful for high identity at nucleotide level)
m_n_b_a    | Nucleotide FASTA    | 5–20 kb              | BLASTn tuned for short contigs (e.g. prophage or viral bins)


## Families JSON Schema

The **Families JSON** file defines the rules that SaPhARI uses to detect specific satellite phage families. Each object corresponds to a family and specifies its required hallmark proteins, optional synonyms, and any disqualifying conditions.

| Field             | Type                                | Meaning                                                                 |
|-------------------|--------------------------------------|-------------------------------------------------------------------------|
| `title`           | string                              | Family name                                                             |
| `proteins`        | list[string or list[string]]        | Required proteins; synonym groups allowed                              |
| `length`          | int                      | Max span                                      |
| `minnumber`       | int                                 | Minimum number of proteins from proteins required (used as a heuristic)                      |
| `forbidden`       | list[string] (optional)             | Disqualifying proteins                                                 |
| `specific_distances` | list[{pair, distance}] (optional) | Max distance between specific pairs; `pair` may contain synonym lists  |


**Example:**

```json
{
  "title": "PICI",
  "proteins": ["integrase", "primase", ["major head", "major capsid"]],
  "length": 15000,
  "minnumber": 3,
  "forbidden": ["PhoQ"],
  "specific_distances": [
    {"pair": ["integrase", "icd-like"], "distance": 5000}
  ]
}

```
## Database Considerations  

The satellite phage families (**PICI**, **CFPICI**, **EPIP**, **P4**, **MFMR**, and **PHANIE**) defined in `Families.json` require a database that actually contains these proteins.  

### Recommended Databases  
- **RefSeq proteins (Bacteria, Archaea, Viruses)**: Deep coverage of microbial and viral proteins.  
- **Custom-augmented RefSeq**: RefSeq plus selected satellite phage proteins (e.g., PLE proteins) for better sensitivity.  

Large protein databases can exceed **100–400 GB** once indexed with DIAMOND. Make sure your system has sufficient storage and memory before running searches.  

### Build a Custom Database  
```bash
diamond makedb --in protein_fasta.faa --db db_name.dmnd --threads num_of_threads
````

* `protein_fasta.faa` → multi-FASTA of amino acid sequences
* `db_name.dmnd` → name of the output database
* `num_of_threads` → CPUs to allocate during indexing

See the [DIAMOND documentation](https://github.com/bbuchfink/diamond) for more options.

## Output Structure
SaPhARI creates a structured output directory. This ensures that intermediate files (annotation, BLAST results, etc.) are separated from the final per-family reports.
```
<out_dir>/
├── amino_acid_orf/ 
├── nucleotide_orf/
├── protein_DIAMOND_BLAST/
├── BLASTn/
├── annotated/
├── RESULTS/
│   └── <FAMILY>/
│       ├── <sample>_<FAMILY>_output.tsv
├── work/
├── .nextflow/
├── .nextflow.log

```

## Example Per-Family Output

A typical per-family TSV (here: `test_genome_PICI_output.tsv`) is written to:

```
<out_dir>/RESULTS/<FAMILY>/<sample>_<FAMILY>_output.tsv
```

Each row corresponds to an ORF (open reading frame) within the detected satellite region, annotated with positional and BLAST/DIAMOND information.


### Core Region Example

Below is an excerpt from the **Core** section of `PICI REGION 1 - test_genome`:

| Region                     | Section | ORF No. | Start    | Stop     | Strand | Protein ID      | %ID   | %Cov  | Blast e-value | Protein Name                                                                 |
|-----------------------------|---------|---------|----------|----------|--------|-----------------|-------|-------|---------------|-------------------------------------------------------------------------------|
| PICI REGION 1 - test_genome | Core    | 1308    | 1377930  | 1379159  | +      | CORE            | 100.0 | 99.8  | 6.49e-308     | PROTEIN: WP_000085269.1 MULTISPECIES: site-specific recombinase [Enterobacteriaceae] |
| PICI REGION 1 - test_genome | Core    | 1309    | 1379766  | 1381052  | +      | CORE            | 55.9  | 80.2  | 3.09e-117     | PROTEIN: WP_327058574.1 host cell division inhibitor [Enterobacteriaceae]            |
| PICI REGION 1 - test_genome | Core    | 1310    | 1381049  | 1381279  | +      | WP_000076671.1  | 100.0 | 98.7  | 4.57e-51      | MULTISPECIES: hypothetical protein [Enterobacteriaceae]                             |
| PICI REGION 1 - test_genome | Core    | 1311    | 1381269  | 1381490  | +      | WP_001372127.1  | 100.0 | 98.6  | 1.76e-49      | MULTISPECIES: hypothetical protein [Enterobacteriaceae]                             |
| PICI REGION 1 - test_genome | Core    | 1312    | 1381435  | 1381704  | +      | AAN79960.1      | 98.7  | 87.8  | 1.26e-52      | Hypothetical protein c1491 [Escherichia coli CFT073]                               |
| PICI REGION 1 - test_genome | Core    | 1313    | 1381706  | 1381939  | +      | WP_001204962.1  | 100.0 | 98.7  | 1.46e-52      | MULTISPECIES: hypothetical protein [Enterobacteriaceae]                             |
| PICI REGION 1 - test_genome | Core    | 1314    | 1381945  | 1382244  | +      | WP_000770152.1  | 100.0 | 99.0  | 1.72e-70      | DUF4222 domain-containing protein [Escherichia coli]                               |
| PICI REGION 1 - test_genome | Core    | 1315    | 1382241  | 1383995  | +      | CORE            | 100.0 | 99.8  | 0.0           | PROTEIN: WP_000761835.1 phage/plasmid primase, replication protein                  |
| PICI REGION 1 - test_genome | Core    | 1316    | 1384284  | 1384541  | +      | AAN79965.1      | 100.0 | 98.8  | 1.96e-56      | Unknown protein encoded by prophage [Escherichia coli CFT073]                       |
| PICI REGION 1 - test_genome | Core    | 1317    | 1384538  | 1384948  | +      | WP_000126698.1  | 100.0 | 99.3  | 2.27e-97      | MULTISPECIES: single-stranded DNA-binding protein [Enterobacteriaceae]              |


### Understanding the Columns

- **Region**: Locus detected as a candidate satellite (e.g. *PICI REGION 1*).  

- **Section**: Subdivision of the locus for contextual interpretation:  
  - *Upstream* – genes flanking the family region  
  - *Core* – hallmark gene cluster that defines the family  
  - *Downstream* – genes flanking the family region  
  
- **ORF No.**: Identifier of the open reading frame from Prodigal.  
- **Start / Stop**: Genomic coordinates (bp).  
- **Strand**: `+` (forward) or `–` (reverse).  
- **Protein ID**: Best BLAST/DIAMOND hit accession. If reported as `CORE`, it indicates a **rule-defining protein** required by the family definition (`families.json`).  
- **%ID / %Cov**: Percent identity and coverage of the alignment.  
- **Blast e-value**: Statistical significance of the match. Lower is stronger.  
- **Protein Name**: Functional annotation of the top hit.  


## Example Runs

### Full Pipeline

```bash
saphari --workflow n_t_a --data_dir tests/ex1/genomes --db tests/ex1/test.dmnd \
  --families_file src/saphari/data/Families.json --out_dir out/run1 --cpus 4

```

### Reanalysis

```bash
saphari --reanalyze --data_dir out/run1/annotated --db tests/ex1/test.dmnd \
  --out_dir out/search_only --search PICI cfPICI --output_type pf

```

### CDS Mode

```bash
saphari --CDS --data_dir tests/ex1/cds_genbank --db tests/ex1/test.dmnd \
  --out_dir out/cds_run --cpus 8

```


## HPC Usage

### Singularity

```bash
export SINGULARITYENV_NXF_HOME=$PWD/out
singularity exec --home $PWD/out -B $PWD:/work -B $PWD/out:/out saphari.sif \
  saphari --workflow n_t_a --data_dir /work/tests/ex1/genomes \
    --db /work/tests/ex1/test.dmnd --out_dir /out/run1

```

### SLURM

```bash
#!/bin/bash -l
#SBATCH --job-name=saphari
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=job.%j.out

module load singularity
export SINGULARITYENV_NXF_HOME=/out

singularity exec --home $PWD/out -B $PWD:/work -B $PWD/out:/out saphari.sif \
  saphari --workflow n_t_a --data_dir /work/tests/ex1/genomes \
    --db /work/tests/ex1/test.dmnd --out_dir /out/run1 --cpus $SLURM_CPUS_PER_TASK

```


## Troubleshooting

-   **No input files found** → Check `--data_dir` and bind mounts.
    
-   **No TSVs produced** → Ensure annotation finished, or provide the right `--reanalyze` path.
    
-   **Permissions errors** → Add `-u $(id -u):$(id -g)` in Docker.
    
-   **Slow performance** → Increase `--cpus`/`--memory`; place DBs on SSD/local disk.
    

## Feedback & Issues

Please use the repository’s **Issues** page to report problems or suggest features.

## Contributing

We welcome contributions! Use GitHub Flow:

1.  **Fork** the repo, create a feature branch.
    
2.  Implement changes with clear commits.
    
3.  Open a **Pull Request** against `master` with details.
    
4.  We’ll review, provide feedback, and merge when ready.
    


## License

MIT License

## Authors and Acknowledgements

Author: Namit Nallapaneni

Thank you to Emma Holley, Eric Walter, Margaret Saha, Eric Bradley, and all of our project sponsors for their advice in developing this tool.
