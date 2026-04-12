#!/usr/bin/env python3

from pathlib import Path
import csv
import argparse
from collections import defaultdict

def accession_from_filename(filename: str, family: str) -> str:
    suffix = f"_{family}_output.tsv"
    if filename.endswith(suffix):
        return filename[:-len(suffix)]
    return filename.replace("_output.tsv", "")

def main():
    parser = argparse.ArgumentParser(
        description="Collect SaPhARI result files above a size threshold and summarize flags by accession/family."
    )
    parser.add_argument(
        "root",
        help="Root folder like /sciclone/home/ncnallapaneni/scr10/saphari_pub/output/new_genomes_fasta"
    )
    parser.add_argument(
        "--min-bytes",
        type=int,
        default=89,
        help="Minimum file size to count as a real hit. Default: 89"
    )
    parser.add_argument(
        "--outdir",
        default="compiled_flags",
        help="Output directory for summary files"
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    family_hits = defaultdict(list)
    accession_to_families = defaultdict(set)
    accession_to_files = defaultdict(list)

    # expected structure: root/<batch>/RESULTS/<FAMILY>/*.tsv
    for batch_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        results_dir = batch_dir / "RESULTS"
        if not results_dir.is_dir():
            continue

        for family_dir in sorted(p for p in results_dir.iterdir() if p.is_dir()):
            family = family_dir.name

            for tsv in sorted(family_dir.glob("*_output.tsv")):
                size = tsv.stat().st_size
                if size < args.min_bytes:
                    continue

                acc = accession_from_filename(tsv.name, family)

                family_hits[family].append({
                    "accession": acc,
                    "batch": batch_dir.name,
                    "family": family,
                    "file": str(tsv),
                    "size_bytes": size,
                })

                accession_to_families[acc].add(family)
                accession_to_files[acc].append({
                    "family": family,
                    "batch": batch_dir.name,
                    "file": str(tsv),
                    "size_bytes": size,
                })

    # 1) one row per flagged file
    all_flagged_file = outdir / "all_flagged_files.tsv"
    with all_flagged_file.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "family", "batch", "size_bytes", "file"])
        for family in sorted(family_hits):
            for row in sorted(family_hits[family], key=lambda x: (x["accession"], x["batch"], x["file"])):
                writer.writerow([
                    row["accession"],
                    row["family"],
                    row["batch"],
                    row["size_bytes"],
                    row["file"]
                ])

    # 2) one row per accession with all families
    accession_summary_file = outdir / "accession_family_summary.tsv"
    with accession_summary_file.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "n_families", "families"])
        for acc in sorted(accession_to_families):
            fams = sorted(accession_to_families[acc])
            writer.writerow([acc, len(fams), ",".join(fams)])

    # 3) only multi-family hits
    multi_file = outdir / "multi_family_hits.tsv"
    with multi_file.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "n_families", "families"])
        for acc in sorted(accession_to_families):
            fams = sorted(accession_to_families[acc])
            if len(fams) > 1:
                writer.writerow([acc, len(fams), ",".join(fams)])

    # 4) family-specific accession lists
    family_lists_dir = outdir / "family_lists"
    family_lists_dir.mkdir(exist_ok=True)

    for family in sorted(family_hits):
        fam_file = family_lists_dir / f"{family}_accessions.txt"
        accessions = sorted({row["accession"] for row in family_hits[family]})
        with fam_file.open("w") as f:
            for acc in accessions:
                f.write(acc + "\n")

    # 5) detailed per-accession breakdown
    detailed_file = outdir / "accession_detailed_breakdown.tsv"
    with detailed_file.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "family", "batch", "size_bytes", "file"])
        for acc in sorted(accession_to_files):
            for row in sorted(accession_to_files[acc], key=lambda x: (x["family"], x["batch"], x["file"])):
                writer.writerow([
                    acc,
                    row["family"],
                    row["batch"],
                    row["size_bytes"],
                    row["file"]
                ])

    # 6) quick stats
    stats_file = outdir / "summary_stats.txt"
    with stats_file.open("w") as f:
        total_flagged_files = sum(len(v) for v in family_hits.values())
        total_unique_accessions = len(accession_to_families)
        total_multi = sum(1 for v in accession_to_families.values() if len(v) > 1)

        f.write(f"Root: {root}\n")
        f.write(f"Min bytes: {args.min_bytes}\n")
        f.write(f"Total flagged files: {total_flagged_files}\n")
        f.write(f"Total unique flagged accessions: {total_unique_accessions}\n")
        f.write(f"Accessions flagged in multiple families: {total_multi}\n\n")

        f.write("Per-family unique accession counts:\n")
        for family in sorted(family_hits):
            count = len({row['accession'] for row in family_hits[family]})
            f.write(f"{family}\t{count}\n")

    print(f"Done. Output written to: {outdir}")
    print(f"Main files:")
    print(f"  {all_flagged_file}")
    print(f"  {accession_summary_file}")
    print(f"  {multi_file}")
    print(f"  {detailed_file}")
    print(f"  {stats_file}")

if __name__ == "__main__":
    main()
