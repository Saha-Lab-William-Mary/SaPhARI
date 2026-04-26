SaPhARI blastn validation workflow

Purpose:
Extract SaPhARI-called nucleotide CORE regions, group them by predicted family,
blastn them against same-family satellite nucleotide reference databases, and
summarize nucleotide-level support.

Scripts:
01_extract_core_region_fastas.py
  Extracts core-only nucleotide FASTAs from SaPhARI result TSVs.
  Default input:
    output/new_genomes_fasta_reanalyzed
  Default output:
    output/core_region_fastas

02_prep_core_family_query_fastas.sh
  Concatenates individual core FASTAs into one query FASTA per family.
  Output:
    output/blast_core_queries/<FAMILY>_core_queries.fasta

03_build_family_blastdb.sh
  Builds blastn nucleotide databases from:
    /sciclone/scr10/ncnallapaneni/satellite_fasta_db/<FAMILY>_refs.fasta
  Output DB prefix:
    /sciclone/scr10/ncnallapaneni/satellite_blastdb/<FAMILY>

04_run_core_family_blastn.slurm
  Slurm array blastn job.
  Queries:
    output/blast_core_queries/<FAMILY>_core_queries.fasta
  DBs:
    /sciclone/scr10/ncnallapaneni/satellite_blastdb/<FAMILY>
  Output:
    output/blast_core_results/<FAMILY>_core_vs_refs.tsv

05_summarize_core_family_wide_support.py
  Summarizes blastn support for each core query region.
  Output:
    output/blast_core_family_support


