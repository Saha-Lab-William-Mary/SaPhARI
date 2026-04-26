#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

BLAST_DIR = Path("/sciclone/scr10/ncnallapaneni/saphari_pub/output/blast_core_results")
OUT_DIR = Path("/sciclone/scr10/ncnallapaneni/saphari_pub/output/blast_core_family_support")
OUT_DIR.mkdir(parents=True, exist_ok=True)

FAMILY_FILES = {
    "CFPICI": BLAST_DIR / "CFPICI_core_vs_refs.tsv",
    "PICI": BLAST_DIR / "PICI_core_vs_refs.tsv",
    "EPIP": BLAST_DIR / "EPIP_core_vs_refs.tsv",
    "MFMR": BLAST_DIR / "MFMR_core_vs_refs.tsv",
    "P4": BLAST_DIR / "P4_core_vs_refs.tsv",
    "PHANIE": BLAST_DIR / "PHANIE_core_vs_refs.tsv",
}


def merge_intervals(intervals):
    if not intervals:
        return []

    norm = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [list(norm[0])]

    for start, end in norm[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start, end])

    return [(s, e) for s, e in merged]


def interval_len(intervals):
    return sum(e - s + 1 for s, e in intervals)


def weighted_identity_from_hsps(rows):
    total_len = sum(r["length"] for r in rows)
    if total_len == 0:
        return 0.0
    return sum(r["pident"] * r["length"] for r in rows) / total_len


def parse_query_metadata(qseqid):
    parts = qseqid.split("|")
    meta = {
        "accession": "",
        "family": "",
        "region": "",
        "core_label": "",
        "batch": "",
        "record": "",
        "start": None,
        "end": None,
        "strand": "",
    }

    if len(parts) >= 4:
        meta["accession"] = parts[0]
        meta["family"] = parts[1]
        meta["region"] = parts[2]
        meta["core_label"] = parts[3]

    for part in parts[4:]:
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        if k == "batch":
            meta["batch"] = v
        elif k == "record":
            meta["record"] = v
        elif k == "start":
            try:
                meta["start"] = int(v)
            except ValueError:
                pass
        elif k == "end":
            try:
                meta["end"] = int(v)
            except ValueError:
                pass
        elif k == "strand":
            meta["strand"] = v

    return meta


def load_rows(path):
    rows = []
    if not path.exists() or path.stat().st_size == 0:
        return rows

    with open(path, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            rows.append({
                "qseqid": row[0],
                "sseqid": row[1],
                "pident": float(row[2]),
                "length": int(row[3]),
                "qstart": int(row[6]),
                "qend": int(row[7]),
                "sstart": int(row[8]),
                "send": int(row[9]),
                "evalue": row[10],
                "bitscore": float(row[11]),
                "qlen": int(row[12]),
                "slen": int(row[13]),
            })
    return rows


def summarize_family(family, blast_path):
    rows = load_rows(blast_path)
    if not rows:
        return [], None

    by_query = defaultdict(list)
    for r in rows:
        by_query[r["qseqid"]].append(r)

    query_summaries = []

    for qseqid, qrows in by_query.items():
        meta = parse_query_metadata(qseqid)
        qlen = qrows[0]["qlen"]

        all_query_intervals = [(r["qstart"], r["qend"]) for r in qrows]
        merged_query = merge_intervals(all_query_intervals)

        supported_bp = interval_len(merged_query)
        core_family_coverage = supported_bp / qlen if qlen else 0.0
        weighted_pident = weighted_identity_from_hsps(qrows)

        by_ref = defaultdict(list)
        for r in qrows:
            by_ref[r["sseqid"]].append(r)

        ref_contribs = []
        for ref, rlist in by_ref.items():
            ref_q_intervals = [(r["qstart"], r["qend"]) for r in rlist]
            ref_merged = merge_intervals(ref_q_intervals)
            ref_supported_bp = interval_len(ref_merged)
            ref_weighted_pident = weighted_identity_from_hsps(rlist)
            bitscore_sum = sum(r["bitscore"] for r in rlist)

            subj_intervals = [(r["sstart"], r["send"]) for r in rlist]
            subj_merged = merge_intervals(subj_intervals)
            subj_supported_bp = interval_len(subj_merged)
            slen = rlist[0]["slen"]
            subj_cov = subj_supported_bp / slen if slen else 0.0

            ref_contribs.append({
                "reference": ref,
                "aligned_query_bp": ref_supported_bp,
                "query_fraction": ref_supported_bp / qlen if qlen else 0.0,
                "aligned_subject_bp": subj_supported_bp,
                "subject_coverage": subj_cov,
                "weighted_pident": ref_weighted_pident,
                "bitscore_sum": bitscore_sum,
                "n_hsps": len(rlist),
            })

        ref_contribs.sort(
            key=lambda x: (
                x["aligned_query_bp"],
                x["subject_coverage"],
                x["bitscore_sum"],
                x["weighted_pident"],
            ),
            reverse=True
        )

        top_ref = ref_contribs[0]["reference"] if ref_contribs else ""
        top_ref_bp = ref_contribs[0]["aligned_query_bp"] if ref_contribs else 0
        top_ref_fraction = ref_contribs[0]["query_fraction"] if ref_contribs else 0.0
        top_ref_pident = ref_contribs[0]["weighted_pident"] if ref_contribs else 0.0

        query_summaries.append({
            "query_id": qseqid,
            "accession": meta["accession"],
            "family": meta["family"] or family,
            "region": meta["region"],
            "batch": meta["batch"],
            "record": meta["record"],
            "query_length": qlen,
            "family_supported_bp_core_region": supported_bp,
            "core_region_family_coverage": round(core_family_coverage, 4),
            "family_supported_weighted_pident": round(weighted_pident, 3),
            "top_contributing_reference": top_ref,
            "top_reference_supported_bp": top_ref_bp,
            "top_reference_query_fraction": round(top_ref_fraction, 4),
            "top_reference_weighted_pident": round(top_ref_pident, 3),
            "n_references_hit": len(ref_contribs),
            "n_hsps_total": len(qrows),
        })

        contrib_path = OUT_DIR / f"{family}_core_reference_contributions.tsv"
        write_header = not contrib_path.exists()
        with open(contrib_path, "a", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "query_id",
                    "family",
                    "reference",
                    "aligned_query_bp",
                    "query_fraction",
                    "aligned_subject_bp",
                    "subject_coverage",
                    "weighted_pident",
                    "bitscore_sum",
                    "n_hsps",
                ],
                delimiter="\t",
            )
            if write_header:
                writer.writeheader()

            for rc in ref_contribs:
                writer.writerow({
                    "query_id": qseqid,
                    "family": meta["family"] or family,
                    "reference": rc["reference"],
                    "aligned_query_bp": rc["aligned_query_bp"],
                    "query_fraction": round(rc["query_fraction"], 4),
                    "aligned_subject_bp": rc["aligned_subject_bp"],
                    "subject_coverage": round(rc["subject_coverage"], 4),
                    "weighted_pident": round(rc["weighted_pident"], 3),
                    "bitscore_sum": round(rc["bitscore_sum"], 1),
                    "n_hsps": rc["n_hsps"],
                })

    query_summaries.sort(key=lambda r: (r["family"], r["accession"], r["region"]))

    fam_summary = {
        "family": family,
        "n_queries": len(query_summaries),
        "mean_core_region_family_coverage": round(mean(r["core_region_family_coverage"] for r in query_summaries), 4),
        "median_core_region_family_coverage": round(median(r["core_region_family_coverage"] for r in query_summaries), 4),
        "mean_family_supported_weighted_pident": round(mean(r["family_supported_weighted_pident"] for r in query_summaries), 3),
        "median_family_supported_weighted_pident": round(median(r["family_supported_weighted_pident"] for r in query_summaries), 3),
        "mean_n_references_hit": round(mean(r["n_references_hit"] for r in query_summaries), 2),
        "median_n_references_hit": round(median(r["n_references_hit"] for r in query_summaries), 2),
    }

    return query_summaries, fam_summary


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    for family in FAMILY_FILES:
        p = OUT_DIR / f"{family}_core_reference_contributions.tsv"
        if p.exists():
            p.unlink()

    all_queries = []
    family_summaries = []

    for family, blast_path in FAMILY_FILES.items():
        qrows, fsum = summarize_family(family, blast_path)

        if qrows:
            write_tsv(
                OUT_DIR / f"{family}_core_query_support.tsv",
                qrows,
                [
                    "query_id",
                    "accession",
                    "family",
                    "region",
                    "batch",
                    "record",
                    "query_length",
                    "family_supported_bp_core_region",
                    "core_region_family_coverage",
                    "family_supported_weighted_pident",
                    "top_contributing_reference",
                    "top_reference_supported_bp",
                    "top_reference_query_fraction",
                    "top_reference_weighted_pident",
                    "n_references_hit",
                    "n_hsps_total",
                ],
            )
            all_queries.extend(qrows)

        if fsum:
            family_summaries.append(fsum)

    if all_queries:
        write_tsv(
            OUT_DIR / "all_core_query_support.tsv",
            all_queries,
            [
                "query_id",
                "accession",
                "family",
                "region",
                "batch",
                "record",
                "query_length",
                "family_supported_bp_core_region",
                "core_region_family_coverage",
                "family_supported_weighted_pident",
                "top_contributing_reference",
                "top_reference_supported_bp",
                "top_reference_query_fraction",
                "top_reference_weighted_pident",
                "n_references_hit",
                "n_hsps_total",
            ],
        )

    if family_summaries:
        write_tsv(
            OUT_DIR / "core_family_summary.tsv",
            sorted(family_summaries, key=lambda r: r["family"]),
            [
                "family",
                "n_queries",
                "mean_core_region_family_coverage",
                "median_core_region_family_coverage",
                "mean_family_supported_weighted_pident",
                "median_family_supported_weighted_pident",
                "mean_n_references_hit",
                "median_n_references_hit",
            ],
        )

    print(f"Wrote outputs to: {OUT_DIR}")


if __name__ == "__main__":
    main()
