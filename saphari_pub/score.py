#!/usr/bin/env python3
"""Score SaPhARI outputs and optionally write clean TSV summaries.


"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional


EXCLUDED_PREFIXES = ("PLE",)


def normalize_true_family(out_dir_name: str) -> Optional[str]:
    """Convert out_<family> folder name into the expected family label.

    Example
    --------
    out_CFPICI   -> CFPICI
    """
    if not out_dir_name.startswith("out_"):
        return None

    fam = out_dir_name[len("out_") :]
    if fam.startswith(EXCLUDED_PREFIXES):
        return None
    return fam


def count_lines(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in fh)


def resolve_default_summary_paths(
    write_details: Optional[str],
    write_family_summary: Optional[str],
    write_overall_summary: Optional[str],
) -> tuple[Optional[Path], Optional[Path], Optional[Path]]:
    """Resolve output paths for detail/family/overall TSVs.

    If --write-details is supplied but summary paths are omitted, summary TSVs are
    automatically written in the same directory as the details file.
    """
    details_path = Path(write_details).resolve() if write_details else None
    family_path = Path(write_family_summary).resolve() if write_family_summary else None
    overall_path = Path(write_overall_summary).resolve() if write_overall_summary else None

    if details_path is not None:
        base_dir = details_path.parent
        if family_path is None:
            family_path = base_dir / "family_summary.tsv"
        if overall_path is None:
            overall_path = base_dir / "overall_summary.tsv"

    return details_path, family_path, overall_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Score SaPhARI outputs for FP/FN using out_<true family> and "
            "RESULTS/<pred family>, with optional clean TSV summary exports."
        )
    )
    parser.add_argument(
        "root",
        help="Root directory containing out_* folders",
    )
    parser.add_argument(
        "--write-details",
        default=None,
        help="Optional path to write per-file details TSV",
    )
    parser.add_argument(
        "--write-family-summary",
        default=None,
        help=(
            "Optional path to write per-family summary TSV. If omitted but "
            "--write-details is provided, defaults to family_summary.tsv in the "
            "same directory."
        ),
    )
    parser.add_argument(
        "--write-overall-summary",
        default=None,
        help=(
            "Optional path to write one-row overall summary TSV. If omitted but "
            "--write-details is provided, defaults to overall_summary.tsv in the "
            "same directory."
        ),
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        raise SystemExit(f"ERROR: root directory does not exist: {root}")

    details_path, family_summary_path, overall_summary_path = resolve_default_summary_paths(
        args.write_details,
        args.write_family_summary,
        args.write_overall_summary,
    )

    total_files_seen = 0
    tp = 0
    fp = 0
    fn = 0
    tn = 0

    family_stats: dict[str, Counter] = defaultdict(Counter)
    detail_rows: list[dict[str, str | int]] = []

    for out_dir in sorted(p for p in root.iterdir() if p.is_dir() and p.name.startswith("out_")):
        expected_family = normalize_true_family(out_dir.name)
        if expected_family is None:
            continue

        results_dir = out_dir / "RESULTS"
        if not results_dir.is_dir():
            print(f"WARNING: missing RESULTS dir in {out_dir}")
            continue

        found_expected_file = False

        for pred_dir in sorted(p for p in results_dir.iterdir() if p.is_dir()):
            predicted_family = pred_dir.name

            for tsv in sorted(pred_dir.glob("*_output.tsv")):
                total_files_seen += 1
                nlines = count_lines(tsv)
                positive = nlines > 1
                is_expected = predicted_family == expected_family

                if is_expected:
                    found_expected_file = True
                    if positive:
                        tp += 1
                        family_stats[expected_family]["TP"] += 1
                        call_type = "TP"
                    else:
                        fn += 1
                        family_stats[expected_family]["FN"] += 1
                        call_type = "FN"
                else:
                    if positive:
                        fp += 1
                        family_stats[expected_family]["FP"] += 1
                        family_stats[predicted_family]["FP_as_predicted"] += 1
                        call_type = "FP"
                    else:
                        tn += 1
                        family_stats[expected_family]["TN"] += 1
                        call_type = "TN"

                detail_rows.append(
                    {
                        "true_input_folder": out_dir.name,
                        "expected_family": expected_family,
                        "predicted_family": predicted_family,
                        "file": str(tsv),
                        "line_count": nlines,
                        "positive_call": int(positive),
                        "classification": call_type,
                    }
                )

        if not found_expected_file:
            print(
                f"WARNING: no predicted folder matching expected family '{expected_family}' inside {results_dir}"
            )

    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0

    print("\n=== OVERALL SUMMARY ===")
    print(f"Root scored: {root}")
    print(f"Total output TSVs seen: {total_files_seen}")
    print(f"TP: {tp}")
    print(f"FP: {fp}")
    print(f"FN: {fn}")
    print(f"TN: {tn}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1:        {f1:.4f}")

    print("\n=== PER EXPECTED FAMILY ===")
    print("Family\tTP\tFP\tFN\tTN\tPrecision\tRecall")

    family_summary_rows: list[dict[str, str | int | float]] = []
    for fam in sorted(family_stats):
        fam_tp = family_stats[fam]["TP"]
        fam_fp = family_stats[fam]["FP"]
        fam_fn = family_stats[fam]["FN"]
        fam_tn = family_stats[fam]["TN"]
        fam_prec = fam_tp / (fam_tp + fam_fp) if (fam_tp + fam_fp) else 0.0
        fam_rec = fam_tp / (fam_tp + fam_fn) if (fam_tp + fam_fn) else 0.0
        print(f"{fam}\t{fam_tp}\t{fam_fp}\t{fam_fn}\t{fam_tn}\t{fam_prec:.4f}\t{fam_rec:.4f}")
        family_summary_rows.append(
            {
                "family": fam,
                "TP": fam_tp,
                "FP": fam_fp,
                "FN": fam_fn,
                "TN": fam_tn,
                "precision": round(fam_prec, 6),
                "recall": round(fam_rec, 6),
            }
        )

    overall_summary_row = {
        "root": str(root),
        "total_output_tsvs_seen": total_files_seen,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
    }

    if details_path:
        details_path.parent.mkdir(parents=True, exist_ok=True)
        with details_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "true_input_folder",
                    "expected_family",
                    "predicted_family",
                    "file",
                    "line_count",
                    "positive_call",
                    "classification",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(detail_rows)
        print(f"\nWrote details to: {details_path}")

    if family_summary_path:
        family_summary_path.parent.mkdir(parents=True, exist_ok=True)
        with family_summary_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["family", "TP", "FP", "FN", "TN", "precision", "recall"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(family_summary_rows)
        print(f"Wrote family summary to: {family_summary_path}")

    if overall_summary_path:
        overall_summary_path.parent.mkdir(parents=True, exist_ok=True)
        with overall_summary_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "root",
                    "total_output_tsvs_seen",
                    "TP",
                    "FP",
                    "FN",
                    "TN",
                    "precision",
                    "recall",
                    "f1",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerow(overall_summary_row)
        print(f"Wrote overall summary to: {overall_summary_path}")


if __name__ == "__main__":
    main()
