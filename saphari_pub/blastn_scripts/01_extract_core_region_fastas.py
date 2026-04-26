#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

FAMILY_NAMES = {"CFPICI", "PICI", "EPIP", "MFMR", "P4", "PHANIE"}

FASTA_SUFFIXES = (
    ".fasta", ".fa", ".fna", ".fas",
    ".fasta.gz", ".fa.gz", ".fna.gz", ".fas.gz"
)


def open_maybe_gzip(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def sanitize_name(s: str) -> str:
    s = s.strip()
    s = s.replace(" ", "_")
    s = s.replace(" - ", "_")
    s = s.replace("/", "_")
    s = re.sub(r"[^A-Za-z0-9_.|=-]+", "_", s)
    return s


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def parse_fasta(path: Path):
    header = None
    chunks = []

    with open_maybe_gzip(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

        if header is not None:
            yield header, "".join(chunks)


def fasta_stem(path: Path) -> str:
    name = path.name
    lower = name.lower()
    for suf in FASTA_SUFFIXES:
        if lower.endswith(suf):
            return name[:-len(suf)]
    return path.stem


def build_genome_index(genomes_root: Path) -> dict[str, Path]:
    index = {}

    for path in genomes_root.rglob("*"):
        if not path.is_file():
            continue

        lower = path.name.lower()
        if not any(lower.endswith(suf) for suf in FASTA_SUFFIXES):
            continue

        index[fasta_stem(path)] = path

    return index


def find_genome_file_for_accession(accession: str, genome_index: dict[str, Path]) -> Path | None:
    if accession in genome_index:
        return genome_index[accession]

    candidates = []
    for stem, path in genome_index.items():
        if stem.startswith(accession):
            candidates.append(path)
        elif path.name.startswith(accession + "."):
            candidates.append(path)

    if len(candidates) == 1:
        return candidates[0]

    return None


def record_hint_from_region_rows(region_rows: list[dict]) -> str:
    for row in region_rows:
        raw = ""
        for key in row:
            if key and ("ORF" in key or key.lower().startswith("query") or key.lower().startswith("sequence")):
                raw = row.get(key, "")
                break

        if raw:
            m = re.match(r"^(contig_[0-9]+)_", raw)
            if m:
                return m.group(1)

    return ""


def choose_record(records: list[tuple[str, str]], accession: str, region_rows: list[dict]) -> tuple[str, str]:
    acc_lower = accession.lower()

    accession_matches = [(h, s) for h, s in records if acc_lower in h.lower()]
    if len(accession_matches) == 1:
        return accession_matches[0]
    if len(accession_matches) > 1:
        return max(accession_matches, key=lambda x: len(x[1]))

    hint = record_hint_from_region_rows(region_rows)
    if hint:
        hint_lower = hint.lower()
        hint_matches = [(h, s) for h, s in records if hint_lower in h.lower()]
        if len(hint_matches) == 1:
            return hint_matches[0]
        if len(hint_matches) > 1:
            return max(hint_matches, key=lambda x: len(x[1]))

    return max(records, key=lambda x: len(x[1]))


def get_region_number(region_name: str) -> str:
    m = re.search(r"REGION\s+(\d+)", region_name, re.IGNORECASE)
    if m:
        return m.group(1)
    return sanitize_name(region_name)


def accession_from_tsv_name(tsv_file: Path, family: str) -> str:
    suffix = f"_{family}_output.tsv"
    name = tsv_file.name
    if name.endswith(suffix):
        return name[:-len(suffix)]
    return name.replace("_output.tsv", "")


def main():
    parser = argparse.ArgumentParser(
        description="Extract core-only nucleotide FASTAs from SaPhARI result TSVs."
    )
    parser.add_argument(
        "--project-root",
        default="/sciclone/scr10/ncnallapaneni/saphari_pub",
        help="SaPhARI project root."
    )
    parser.add_argument(
        "--results-root",
        default=None,
        help="SaPhARI results root. Default: <project-root>/output/new_genomes_fasta_reanalyzed"
    )
    parser.add_argument(
        "--genomes-root",
        default=None,
        help="Genome FASTA root. Default: <project-root>/assembled_genomes"
    )
    parser.add_argument(
        "--out-root",
        default=None,
        help="Output root. Default: <project-root>/output/core_region_fastas"
    )
    parser.add_argument(
        "--min-bytes",
        type=int,
        default=100,
        help="Skip result TSVs smaller than this many bytes."
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    results_root = Path(args.results_root).resolve() if args.results_root else project_root / "output" / "new_genomes_fasta_reanalyzed"
    genomes_root = Path(args.genomes_root).resolve() if args.genomes_root else project_root / "assembled_genomes"
    out_root = Path(args.out_root).resolve() if args.out_root else project_root / "output" / "core_region_fastas"

    out_root.mkdir(parents=True, exist_ok=True)
    manifest_path = out_root / "core_region_manifest.tsv"

    genome_index = build_genome_index(genomes_root)
    manifest_rows = []

    for batch_dir in sorted(results_root.iterdir()):
        if not batch_dir.is_dir():
            continue

        batch = batch_dir.name
        results_dir = batch_dir / "RESULTS"
        if not results_dir.exists():
            continue

        for family_dir in sorted(results_dir.iterdir()):
            if not family_dir.is_dir():
                continue

            family = family_dir.name
            if family not in FAMILY_NAMES:
                continue

            family_out = out_root / family
            family_out.mkdir(parents=True, exist_ok=True)

            for tsv_file in sorted(family_dir.glob("*_output.tsv")):
                if tsv_file.stat().st_size <= args.min_bytes:
                    continue

                accession = accession_from_tsv_name(tsv_file, family)

                with open(tsv_file, "r", newline="") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    rows = list(reader)

                if not rows:
                    continue

                grouped = defaultdict(list)
                for row in rows:
                    region_name = row.get("Region", "").strip()
                    if region_name:
                        grouped[region_name].append(row)

                genome_path = find_genome_file_for_accession(accession, genome_index)

                for region_name, region_rows in grouped.items():
                    region_num = get_region_number(region_name)

                    base_manifest = {
                        "batch": batch,
                        "family": family,
                        "accession": accession,
                        "region_name": region_name,
                        "region_number": region_num,
                        "core_start": "",
                        "core_end": "",
                        "strand": "",
                        "n_orfs": len(region_rows),
                        "n_core_orfs": sum(
                            1 for r in region_rows
                            if r.get("Section", "").strip().lower() == "core"
                        ),
                        "genome_fasta": str(genome_path) if genome_path else "",
                        "output_fasta": "",
                        "status": "",
                        "source_tsv": str(tsv_file),
                    }

                    if genome_path is None:
                        base_manifest["status"] = "genome_not_found"
                        manifest_rows.append(base_manifest)
                        continue

                    records = list(parse_fasta(genome_path))
                    if not records:
                        base_manifest["status"] = "empty_genome_fasta"
                        manifest_rows.append(base_manifest)
                        continue

                    core_coords = []
                    strands = []

                    for row in region_rows:
                        if row.get("Section", "").strip().lower() != "core":
                            continue

                        try:
                            start = int(row["Start"])
                            stop = int(row["Stop"])
                        except Exception:
                            continue

                        lo = min(start, stop)
                        hi = max(start, stop)
                        core_coords.append((lo, hi))

                        strand = str(row.get("Strand", "")).strip()
                        if strand in {"1", "-1"}:
                            strands.append(strand)

                    if not core_coords:
                        base_manifest["status"] = "no_core_coordinates"
                        manifest_rows.append(base_manifest)
                        continue

                    chosen_header, genome_seq = choose_record(records, accession, region_rows)
                    genome_len = len(genome_seq)

                    core_start = min(x[0] for x in core_coords)
                    core_end = max(x[1] for x in core_coords)

                    strand_set = set(strands)
                    if len(strand_set) == 1:
                        region_strand = next(iter(strand_set))
                    elif len(strand_set) > 1:
                        region_strand = "mixed"
                    else:
                        region_strand = ""

                    clipped_start = max(1, core_start)
                    clipped_end = min(genome_len, core_end)

                    base_manifest["core_start"] = core_start
                    base_manifest["core_end"] = core_end
                    base_manifest["strand"] = region_strand

                    if clipped_start > clipped_end:
                        base_manifest["status"] = "invalid_core_span_after_clipping"
                        manifest_rows.append(base_manifest)
                        continue

                    seq = genome_seq[clipped_start - 1:clipped_end]

                    if region_strand == "-1":
                        seq = reverse_complement(seq)

                    out_name = f"{accession}_R{region_num}_core.fasta"
                    out_path = family_out / out_name

                    fasta_header = (
                        f">{accession}|{family}|R{region_num}|core"
                        f"|batch={batch}"
                        f"|record={sanitize_name(chosen_header)}"
                        f"|start={clipped_start}|end={clipped_end}"
                        f"|strand={region_strand}"
                    )

                    with open(out_path, "w") as out_fh:
                        out_fh.write(fasta_header + "\n")
                        out_fh.write(wrap_fasta(seq) + "\n")

                    base_manifest["output_fasta"] = str(out_path)
                    base_manifest["status"] = "ok"
                    manifest_rows.append(base_manifest)

    manifest_rows.sort(
        key=lambda r: (
            r["family"],
            r["accession"],
            str(r["region_number"]),
            str(r["core_start"])
        )
    )

    fieldnames = [
        "batch",
        "family",
        "accession",
        "region_name",
        "region_number",
        "core_start",
        "core_end",
        "strand",
        "n_orfs",
        "n_core_orfs",
        "genome_fasta",
        "output_fasta",
        "status",
        "source_tsv",
    ]

    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    ok_count = sum(1 for r in manifest_rows if r["status"] == "ok")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Extracted core regions: {ok_count}")
    print(f"Total manifest rows: {len(manifest_rows)}")


if __name__ == "__main__":
    main()
