#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path
from statistics import mean, median

PROJECT = Path("/sciclone/scr10/ncnallapaneni/saphari_pub")

MANIFEST = PROJECT / "output/core_region_fastas/core_region_manifest.tsv"
SUPPORT = PROJECT / "output/blast_core_family_support/all_core_query_support.tsv"
OUTDIR = PROJECT / "output/blast_core_family_support"
OUTDIR.mkdir(parents=True, exist_ok=True)

FAMILIES = ["CFPICI", "PICI", "EPIP", "MFMR", "P4", "PHANIE"]


def pct(x):
    if x is None:
        return ""
    return f"{x * 100:.2f}%"


def safe_mean(vals):
    return mean(vals) if vals else None


def safe_median(vals):
    return median(vals) if vals else None


def fmt_num(x, nd=2):
    if x is None:
        return "—"
    return f"{x:.{nd}f}"


def fmt_int(x):
    if x is None:
        return "—"
    return str(int(round(x)))


def load_manifest():
    rows = []

    with MANIFEST.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            if r.get("status") != "ok":
                continue

            region = f"R{r['region_number']}"
            key = (r["accession"], r["family"], region)

            try:
                core_start = int(r["core_start"])
                core_end = int(r["core_end"])
                core_len = abs(core_end - core_start) + 1
            except Exception:
                core_start = None
                core_end = None
                core_len = None

            rows.append({
                "key": key,
                "family": r["family"],
                "accession": r["accession"],
                "region": region,
                "core_start": core_start,
                "core_end": core_end,
                "core_len_bp": core_len,
                "n_core_orfs": int(r["n_core_orfs"]) if r["n_core_orfs"] else None,
                "output_fasta": r["output_fasta"],
                "source_tsv": r["source_tsv"],
            })

    return rows


def load_support():
    support = {}

    if not SUPPORT.exists():
        return support

    with SUPPORT.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            key = (r["accession"], r["family"], r["region"])
            support[key] = {
                "query_length": int(r["query_length"]),
                "supported_bp": int(r["family_supported_bp_core_region"]),
                "coverage": float(r["core_region_family_coverage"]),
                "weighted_identity": float(r["family_supported_weighted_pident"]),
                "top_reference": r["top_contributing_reference"],
                "top_reference_fraction": float(r["top_reference_query_fraction"]),
                "top_reference_identity": float(r["top_reference_weighted_pident"]),
                "n_references_hit": int(r["n_references_hit"]),
                "n_hsps_total": int(r["n_hsps_total"]),
                "query_id": r["query_id"],
            }

    return support


def coverage_bin(cov):
    if cov == 0:
        return "0%"
    if cov < 0.10:
        return ">0-10%"
    if cov < 0.25:
        return "10-25%"
    if cov < 0.50:
        return "25-50%"
    if cov < 0.75:
        return "50-75%"
    return "75-100%"


def write_tsv(path, rows, fieldnames):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})


def summarize_family(fam, rows):
    extracted_n = len(rows)
    supported_rows = [r for r in rows if r["had_blastn_support"] == 1]
    supported_n = len(supported_rows)

    core_lens = [r["core_len_bp"] for r in rows if r["core_len_bp"] is not None]
    supported_bp_all = [r["supported_bp"] for r in rows]
    supported_bp_supported_only = [r["supported_bp"] for r in supported_rows]

    coverage_all_zeroes = [r["core_region_family_coverage"] for r in rows]
    coverage_supported = [r["core_region_family_coverage"] for r in supported_rows]

    identities_supported = [
        r["family_supported_weighted_identity"]
        for r in supported_rows
        if r["family_supported_weighted_identity"] != ""
    ]

    refs_supported = [r["n_references_hit"] for r in supported_rows]
    hsps_supported = [r["n_hsps_total"] for r in supported_rows]

    total_core_bp = sum(core_lens)
    total_supported_bp = sum(r["supported_bp"] for r in rows)

    family_wide_bp_coverage = (total_supported_bp / total_core_bp) if total_core_bp else None

    n_cov_ge_10 = sum(1 for x in coverage_all_zeroes if x >= 0.10)
    n_cov_ge_25 = sum(1 for x in coverage_all_zeroes if x >= 0.25)
    n_cov_ge_50 = sum(1 for x in coverage_all_zeroes if x >= 0.50)
    n_cov_ge_75 = sum(1 for x in coverage_all_zeroes if x >= 0.75)

    n_cov50_id90 = sum(
        1 for r in supported_rows
        if r["core_region_family_coverage"] >= 0.50
        and r["family_supported_weighted_identity"] != ""
        and r["family_supported_weighted_identity"] >= 90
    )

    return {
        "family": fam,
        "extracted_core_regions": extracted_n,
        "core_regions_with_blastn_support": supported_n,
        "blastn_support_rate": pct(supported_n / extracted_n) if extracted_n else "—",

        "total_core_bp_tested": total_core_bp,
        "total_blastn_supported_bp": total_supported_bp,
        "family_wide_bp_coverage": pct(family_wide_bp_coverage) if family_wide_bp_coverage is not None else "—",

        "mean_core_region_len_bp": fmt_num(safe_mean(core_lens), 1),
        "median_core_region_len_bp": fmt_num(safe_median(core_lens), 1),

        "mean_supported_bp_all_regions_zero_if_no_hit": fmt_num(safe_mean(supported_bp_all), 1),
        "median_supported_bp_all_regions_zero_if_no_hit": fmt_num(safe_median(supported_bp_all), 1),

        "mean_supported_bp_supported_regions_only": fmt_num(safe_mean(supported_bp_supported_only), 1),
        "median_supported_bp_supported_regions_only": fmt_num(safe_median(supported_bp_supported_only), 1),

        "mean_coverage_all_regions_zero_if_no_hit": pct(safe_mean(coverage_all_zeroes)) if coverage_all_zeroes else "—",
        "median_coverage_all_regions_zero_if_no_hit": pct(safe_median(coverage_all_zeroes)) if coverage_all_zeroes else "—",

        "mean_coverage_supported_regions_only": pct(safe_mean(coverage_supported)) if coverage_supported else "—",
        "median_coverage_supported_regions_only": pct(safe_median(coverage_supported)) if coverage_supported else "—",

        "mean_weighted_identity_supported_regions_only": f"{safe_mean(identities_supported):.2f}%" if identities_supported else "—",
        "median_weighted_identity_supported_regions_only": f"{safe_median(identities_supported):.2f}%" if identities_supported else "—",

        "mean_n_references_hit_supported_regions_only": fmt_num(safe_mean(refs_supported), 2),
        "median_n_references_hit_supported_regions_only": fmt_num(safe_median(refs_supported), 2),

        "mean_n_hsps_supported_regions_only": fmt_num(safe_mean(hsps_supported), 2),
        "median_n_hsps_supported_regions_only": fmt_num(safe_median(hsps_supported), 2),

        "n_regions_coverage_ge_10pct": n_cov_ge_10,
        "n_regions_coverage_ge_25pct": n_cov_ge_25,
        "n_regions_coverage_ge_50pct": n_cov_ge_50,
        "n_regions_coverage_ge_75pct": n_cov_ge_75,
        "n_regions_coverage_ge_50pct_and_identity_ge_90pct": n_cov50_id90,
    }


def main():
    manifest_rows = load_manifest()
    support = load_support()

    joined = []
    for r in manifest_rows:
        b = support.get(r["key"])

        if b is None:
            cov = 0.0
            ident = ""
            supported = 0
            supported_bp = 0
            query_length = r["core_len_bp"] if r["core_len_bp"] is not None else ""
            top_reference = ""
            top_reference_fraction = ""
            top_reference_identity = ""
            n_references_hit = 0
            n_hsps_total = 0
            query_id = ""
        else:
            cov = b["coverage"]
            ident = b["weighted_identity"]
            supported = 1
            supported_bp = b["supported_bp"]
            query_length = b["query_length"]
            top_reference = b["top_reference"]
            top_reference_fraction = b["top_reference_fraction"]
            top_reference_identity = b["top_reference_identity"]
            n_references_hit = b["n_references_hit"]
            n_hsps_total = b["n_hsps_total"]
            query_id = b["query_id"]

        joined.append({
            "family": r["family"],
            "accession": r["accession"],
            "region": r["region"],
            "core_start": r["core_start"],
            "core_end": r["core_end"],
            "core_len_bp": r["core_len_bp"],
            "n_core_orfs": r["n_core_orfs"],
            "had_blastn_support": supported,
            "query_length": query_length,
            "supported_bp": supported_bp,
            "core_region_family_coverage": cov,
            "core_region_family_coverage_percent": pct(cov),
            "coverage_bin": coverage_bin(cov),
            "family_supported_weighted_identity": ident,
            "family_supported_weighted_identity_percent": "" if ident == "" else f"{ident:.2f}%",
            "top_reference": top_reference,
            "top_reference_query_fraction": top_reference_fraction,
            "top_reference_query_fraction_percent": "" if top_reference_fraction == "" else pct(top_reference_fraction),
            "top_reference_weighted_identity": top_reference_identity,
            "n_references_hit": n_references_hit,
            "n_hsps_total": n_hsps_total,
            "query_id": query_id,
            "output_fasta": r["output_fasta"],
            "source_tsv": r["source_tsv"],
        })

    region_fields = [
        "family",
        "accession",
        "region",
        "core_start",
        "core_end",
        "core_len_bp",
        "n_core_orfs",
        "had_blastn_support",
        "query_length",
        "supported_bp",
        "core_region_family_coverage",
        "core_region_family_coverage_percent",
        "coverage_bin",
        "family_supported_weighted_identity",
        "family_supported_weighted_identity_percent",
        "top_reference",
        "top_reference_query_fraction",
        "top_reference_query_fraction_percent",
        "top_reference_weighted_identity",
        "n_references_hit",
        "n_hsps_total",
        "query_id",
        "output_fasta",
        "source_tsv",
    ]

    write_tsv(
        OUTDIR / "core_region_blastn_region_level_with_zeroes.tsv",
        joined,
        region_fields,
    )

    summary_rows = []
    bin_rows = []

    for fam in FAMILIES:
        rows = [r for r in joined if r["family"] == fam]
        summary_rows.append(summarize_family(fam, rows))

        extracted_n = len(rows)
        bins = {
            "0%": 0,
            ">0-10%": 0,
            "10-25%": 0,
            "25-50%": 0,
            "50-75%": 0,
            "75-100%": 0,
        }

        for r in rows:
            bins[coverage_bin(r["core_region_family_coverage"])] += 1

        for bname, count in bins.items():
            bin_rows.append({
                "family": fam,
                "coverage_bin": bname,
                "n_regions": count,
                "percent_of_extracted_regions": pct(count / extracted_n) if extracted_n else "—",
            })

    summary_rows.append(summarize_family("Total", joined))

    summary_fields = [
        "family",
        "extracted_core_regions",
        "core_regions_with_blastn_support",
        "blastn_support_rate",

        "total_core_bp_tested",
        "total_blastn_supported_bp",
        "family_wide_bp_coverage",

        "mean_core_region_len_bp",
        "median_core_region_len_bp",

        "mean_supported_bp_all_regions_zero_if_no_hit",
        "median_supported_bp_all_regions_zero_if_no_hit",

        "mean_supported_bp_supported_regions_only",
        "median_supported_bp_supported_regions_only",

        "mean_coverage_all_regions_zero_if_no_hit",
        "median_coverage_all_regions_zero_if_no_hit",

        "mean_coverage_supported_regions_only",
        "median_coverage_supported_regions_only",

        "mean_weighted_identity_supported_regions_only",
        "median_weighted_identity_supported_regions_only",

        "mean_n_references_hit_supported_regions_only",
        "median_n_references_hit_supported_regions_only",
        "mean_n_hsps_supported_regions_only",
        "median_n_hsps_supported_regions_only",

        "n_regions_coverage_ge_10pct",
        "n_regions_coverage_ge_25pct",
        "n_regions_coverage_ge_50pct",
        "n_regions_coverage_ge_75pct",
        "n_regions_coverage_ge_50pct_and_identity_ge_90pct",
    ]

    bin_fields = [
        "family",
        "coverage_bin",
        "n_regions",
        "percent_of_extracted_regions",
    ]

    write_tsv(OUTDIR / "core_blastn_summary_metrics_expanded.tsv", summary_rows, summary_fields)
    write_tsv(OUTDIR / "core_blastn_coverage_bins.tsv", bin_rows, bin_fields)

    print(f"Wrote: {OUTDIR / 'core_region_blastn_region_level_with_zeroes.tsv'}")
    print(f"Wrote: {OUTDIR / 'core_blastn_summary_metrics_expanded.tsv'}")
    print(f"Wrote: {OUTDIR / 'core_blastn_coverage_bins.tsv'}")


if __name__ == "__main__":
    main()
