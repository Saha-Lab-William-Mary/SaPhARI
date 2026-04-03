#!/usr/bin/env python3
"""
SaPhARI CLI
"""
from __future__ import annotations

import os
import sys
import glob
import json
import argparse
import subprocess
from typing import Dict, List, Any, Tuple

# Package-relative imports (works when installed as a package or run from src)
try:
    from .satellite import Satellite  # type: ignore
    from .convert_output import process_file, export_regions_to_tsv  # type: ignore
except Exception:  # pragma: no cover - fallback for flat module layout
    from satellite import Satellite  # type: ignore
    from convert_output import process_file, export_regions_to_tsv  # type: ignore


# ----------------------------- helpers -------------------------------------

def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def resolve_code_root(arg_code_dir: str) -> str:
    """Return a directory that contains 'annotation/main.nf'.
    Prefers an installed package path; falls back to user-supplied or nearby.
    """
    candidates: List[str] = []

    # 1) If a code_dir was provided, try it first
    if arg_code_dir:
        candidates.append(os.path.abspath(arg_code_dir))

    # 2) Environment override
    env_cd = os.environ.get("SAPHARI_CODE_DIR", "")
    if env_cd:
        candidates.append(os.path.abspath(env_cd))

    # 3) Installed package directory via importlib.resources
    try:
        import importlib.resources as ir  # py3.9+
        pkg_root = os.path.abspath(str(ir.files("saphari")))  # type: ignore[arg-type]
        candidates.append(pkg_root)
    except Exception:
        pass

    # 4) Relative to this file
    here = os.path.abspath(os.path.dirname(__file__))
    candidates.extend([
        here,
        os.path.join(here, "saphari"),
        os.path.join(os.path.dirname(here), "saphari"),
    ])

    for cand in candidates:
        nf = os.path.join(cand, "annotation", "main.nf")
        if os.path.isfile(nf):
            return cand

    sys.stderr.write(
        "Error: Could not locate 'annotation/main.nf'. Try passing --code_dir pointing to the\n"
        "installed 'saphari' package directory (the one that contains 'annotation/').\n"
        f"Tried: {candidates}\n"
    )
    sys.exit(1)


def load_families_from_json(json_path: str) -> List[Dict[str, Any]]:
    """Load Families.json and validate that it is a LIST of family objects.
    We do NOT coerce shapes here; normalization happens right before add_family.
    """
    with open(json_path, "r") as fp:
        data = json.load(fp)
    if not isinstance(data, list):
        raise ValueError("Families.json must be a list of family objects.")
    for i, fam in enumerate(data):
        if not isinstance(fam, dict) or "title" not in fam:
            raise ValueError(f"Family index {i} must be an object with a 'title' field.")
    return data


def resolve_default_families_json() -> str:
    """Return the installed Families.json path using importlib.resources."""
    try:
        import importlib.resources as ir
        path = ir.files("saphari").joinpath("data").joinpath("Families.json")  # type: ignore[arg-type]
        return os.fspath(path)
    except Exception:
        # Fallback: relative to this file
        here = os.path.abspath(os.path.dirname(__file__))
        alt = os.path.join(here, "data", "Families.json")
        if os.path.isfile(alt):
            return alt
        # One more fallback: sibling saphari/data
        alt2 = os.path.join(os.path.dirname(here), "saphari", "data", "Families.json")
        if os.path.isfile(alt2):
            return alt2
        raise FileNotFoundError("Could not locate default Families.json in installed package.")


def build_title_map(families_list: List[Dict[str, Any]]) -> Dict[str, str]:
    """Map lowercase title -> canonical title."""
    return {fam["title"].lower(): fam["title"] for fam in families_list}


def normalize_search_list(user_search: List[str], title_map: Dict[str, str]) -> List[str]:
    out: List[str] = []
    for name in user_search:
        key = name.lower()
        if key in title_map:
            out.append(title_map[key])
        else:
            # allow partial case-insensitive match on substring
            candidates = [v for k, v in title_map.items() if key in k]
            if candidates:
                out.extend(sorted(set(candidates)))
            else:
                sys.stderr.write(f"[WARN] Requested family '{name}' not found in Families.json. Skipping.\n")
    # de-dupe preserving order
    seen = set()
    deduped: List[str] = []
    for x in out:
        if x not in seen:
            seen.add(x)
            deduped.append(x)
    return deduped


def _flatten_proteins(proteins: List[Any]) -> List[str]:
    """Flatten synonym groups (lists) into individual string tokens.
    De-duplicate while preserving order.
    """
    flat: List[str] = []
    for p in proteins:
        if isinstance(p, list):
            flat.extend([str(x) for x in p])
        else:
            flat.append(str(p))
    seen: set = set()
    out: List[str] = []
    for tok in flat:
        if tok not in seen:
            seen.add(tok)
            out.append(tok)
    return out


def _expand_specific_distances(sd_list: List[Dict[str, Any]]) -> Dict[Tuple[str, str], int]:
    """Expand pairs with synonym lists into all string x string pairs.
    If the same pair appears multiple times, keep the max distance.
    """
    out: Dict[Tuple[str, str], int] = {}
    for item in sd_list or []:
        pair = item.get("pair")
        dist = int(item.get("distance", 0))
        if not pair or len(pair) != 2:
            continue
        a, b = pair
        A = a if isinstance(a, list) else [a]
        B = b if isinstance(b, list) else [b]
        for aa in A:
            for bb in B:
                key = (str(aa), str(bb))
                prev = out.get(key, 0)
                if dist > prev:
                    out[key] = dist
    return out


def post_process_results(out_dir: str) -> None:
    """Convert RESULTS/*/*_output.txt -> *.tsv and REMOVE the .txt after success.
    Does NOT touch annotated TSVs.
    """
    results_root = os.path.join(out_dir, "RESULTS")
    if not os.path.isdir(results_root):
        return
    for family_name in os.listdir(results_root):
        fam_dir = os.path.join(results_root, family_name)
        if not os.path.isdir(fam_dir):
            continue
        for fname in os.listdir(fam_dir):
            if not fname.endswith("_output.txt"):
                continue
            txt_path = os.path.join(fam_dir, fname)
            try:
                regions = process_file(txt_path)
                tsv_path = os.path.splitext(txt_path)[0] + ".tsv"
                export_regions_to_tsv(regions, tsv_path)
                os.remove(txt_path)
                print(f"Converted {txt_path} -> {tsv_path}")
            except Exception as e:
                print(f"[WARNING] Failed converting {txt_path}: {e}")


# ------------------------------ CLI ----------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="saphari",
        description=(
            "Satellite Phage Algorithmic Recognition and Interpretation (SaPhARI)\n\n"
            "Process nucleotide sequences (FASTA/GenBank), run annotation (Prodigal/DIAMOND/BLASTn), "
            "and search for putative satellite-phage gene families."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Annotation pipeline params
    g = parser.add_argument_group("Annotation Pipeline")
    g.add_argument("--code_dir", "-c", default="",
                   help="Path that contains 'annotation/main.nf'. Defaults to the installed package dir.")
    g.add_argument("--workflow", "-w",
                   choices=["n_t_a", "m_n_t_a", "a_t_a", "n_b_a", "m_n_b_a"],
                   help="Nextflow workflow to run (omit if --CDS or --reanalyze).");
    g.add_argument("--CDS", action="store_true",
                   help="Skip annotation: treat input as GenBank CDS (AA).");
    g.add_argument("--reanalyze", action="store_true",
                   help="Skip annotation: re-run family search on existing annotated TSVs.");
    g.add_argument("--db", "-d", required=True,
                   help="Path to DIAMOND protein DB (.dmnd) or BLAST+ DB (for n_b_a workflows).");
    g.add_argument("--ndb", "-n", default="",
                   help="Path to nucleotide BLAST database (if using BLASTn workflows).");

    # Family search
    s = parser.add_argument_group("Family Search")
    s.add_argument("--families_file", "-F",
                   help="Optional JSON file defining families (overrides built-in).");
    s.add_argument("--search", "-f", nargs="+",
                   help="List of family names to search for. Default = all in JSON.");
    s.add_argument("--output_type", "-t",
                   type=lambda v: v.strip().lower(),
                   choices=["pfpf", "pf"], default="pfpf",
                   help="pfpf: per-file-per-family; pf: per-family aggregated.");

    # I/O
    io = parser.add_argument_group("Input / Output")
    io.add_argument("--data_dir", "-i", required=True,
                    help="Input directory: genomes for annotation or annotated TSVs if --reanalyze.");
    io.add_argument("--out_dir", "-o", required=True,
                    help="Output directory (annotated/, RESULTS/, work/).");

    # Tuning
    tune = parser.add_argument_group("Annotation Tuning")
    tune.add_argument("--dim_e", default="0.0001")
    tune.add_argument("--dim_m", default="10")
    tune.add_argument("--dim_c", default="20")
    tune.add_argument("--dim_p", default="20")
    tune.add_argument("--nb_e", default="0.0001")
    tune.add_argument("--nb_m", default="5")
    tune.add_argument("--nb_c", default="20")
    tune.add_argument("--nb_p", default="20")

    # Runtime
    r = parser.add_argument_group("Runtime / Resources")
    r.add_argument("--cpus", type=int, default=4)
    r.add_argument("--memory", default="8 GB")
    r.add_argument("--profile", default="local")
    r.add_argument("--container", default="");

    args = parser.parse_args()

    # Validate input/output
    data_dir = os.path.abspath(args.data_dir)
    if not os.path.isdir(data_dir):
        sys.stderr.write(f"Error: --data_dir '{data_dir}' not found.\n")
        sys.exit(1)
    out_dir = ensure_dir(os.path.abspath(args.out_dir))

    # Resolve code root for annotation
    CODE_ROOT = resolve_code_root(args.code_dir)
    ANNOTATION_DIR = os.path.join(CODE_ROOT, "annotation")

    # Families
    if args.families_file:
        if not os.path.isfile(args.families_file):
            sys.stderr.write(f"Error: families file '{args.families_file}' not found.\n")
            sys.exit(1)
        families_list = load_families_from_json(args.families_file)
    else:
        try:
            default_fam_json = resolve_default_families_json()
            families_list = load_families_from_json(default_fam_json)
        except Exception as e:
            sys.stderr.write(f"Error: Could not load default Families.json: {e}\n")
            sys.exit(1)

    title_map = build_title_map(families_list)

    # Initialize Satellite and add families with FLATTENED shapes (Option A)
    sat = Satellite()
    for fam in families_list:
        sat.add_family(
            title=fam["title"],
            proteins=_flatten_proteins(fam.get("proteins", [])),  # strings only
            length=fam.get("length"),
            minnumber=fam.get("minnumber"),
            forbidden=fam.get("forbidden", []),
            specific_distances=_expand_specific_distances(fam.get("specific_distances", [])),  # {(a,b): dist}
        )

    if args.search:
        families_to_search = normalize_search_list(args.search, title_map)
        if not families_to_search:
            sys.stderr.write("Error: none of the requested families matched the JSON titles.\n")
            sys.exit(1)
    else:
        families_to_search = [fam["title"] for fam in families_list]

    # -------------------------- Modes --------------------------------------
    # 1) CDS mode (AA GenBank)
    if args.CDS:
        seq_patterns = ["*.gbk", "*.gb", "*.gbff"]
        sequence_files = sorted(sum((glob.glob(os.path.join(data_dir, p)) for p in seq_patterns), []))
        if not sequence_files:
            sys.stderr.write("Error: No GenBank CDS files (*.gbk|*.gb|*.gbff) found in --data_dir.\n")
            sys.exit(1)
        for gb_file in sequence_files:
            try:
                sat.process_files(
                    file_list=[gb_file],
                    use_cds=True,
                    output_type=args.output_type,
                    inputpath=data_dir,
                    familiestosearchfor=families_to_search,
                    outDir=out_dir,
                )
            except Exception as e:
                sys.stderr.write(f"[WARN] Skipping CDS file '{gb_file}': {e}\n")
        post_process_results(out_dir)
        sys.exit(0)

    # 2) Reanalyze mode (use existing annotated TSVs)
    if args.reanalyze:
        annotated_files = sorted(glob.glob(os.path.join(data_dir, "**", "*.tsv"), recursive=True))
        if annotated_files:
            try:
                tsv_inputpath = os.path.commonpath(annotated_files)
            except ValueError:
                tsv_inputpath = os.path.dirname(annotated_files[0])
        else:
            # Also try out_dir/annotated
            annotated_folder = os.path.join(out_dir, "annotated")
            annotated_files = sorted(glob.glob(os.path.join(annotated_folder, "*.tsv")))
            tsv_inputpath = annotated_folder if annotated_files else None

        if not annotated_files or not tsv_inputpath:
            sys.stderr.write("Error: No annotated TSV files found in --data_dir or out_dir/annotated.\n")
            sys.exit(1)

        sat.process_files(
            file_list=annotated_files,
            use_cds=False,
            output_type=args.output_type,
            inputpath=tsv_inputpath,
            familiestosearchfor=families_to_search,
            outDir=out_dir,
        )
        post_process_results(out_dir)
        sys.exit(0)

    # 3) Full pipeline (Nextflow + Satellite)
    if not args.workflow:
        sys.stderr.write("Error: --workflow is required unless --CDS or --reanalyze.\n")
        sys.exit(1)

    # Sanity-check inputs exist
    if not glob.glob(os.path.join(data_dir, "*")):
        sys.stderr.write("Error: No input files found in --data_dir.\n")
        sys.exit(1)

    nf_script = os.path.join(ANNOTATION_DIR, "main.nf")
    if not os.path.isfile(nf_script):
        sys.stderr.write(f"Error: Nextflow script not found at '{nf_script}'.\n")
        sys.exit(1)

    # Nextflow execution in a writable cwd with explicit work dir & NXF_HOME
    work_dir = ensure_dir(os.path.join(out_dir, "work"))
    env = os.environ.copy()
    env.setdefault("NXF_HOME", out_dir)

    nf_cmd = [
        "nextflow", "run", nf_script,
        "--data_dir", data_dir,
        "--sequences", os.path.join(data_dir, "*"),
        "--code_dir", CODE_ROOT,
        "--out_dir", out_dir,
        "--pdbpath", os.path.abspath(args.db),
        "--ndbpath", os.path.abspath(args.ndb),
        "--dimevalue", str(args.dim_e),
        "--dimmatchamount", str(args.dim_m),
        "--dimcov", str(args.dim_c),
        "--dimpi", str(args.dim_p),
        "--nbevalue", str(args.nb_e),
        "--nbmatchamount", str(args.nb_m),
        "--nbcov", str(args.nb_c),
        "--nbpi", str(args.nb_p),
        "--cpus", str(args.cpus),
        "--memory", str(args.memory),
        "--container", str(args.container),
        "-entry", args.workflow,
        "-work-dir", work_dir,
    ]
    if args.profile:
        nf_cmd.extend(["-profile", args.profile])

    print("Running Nextflow annotation pipeline with command:")
    print(" ".join(nf_cmd), "\n")

    try:
        proc = subprocess.Popen(
            nf_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=out_dir,  # ensure .nextflow/* and history can be created
            env=env,
        )
    except FileNotFoundError as e:
        sys.stderr.write(f"Error: {e}\nMake sure 'nextflow' is installed and on PATH.\n")
        sys.exit(1)

    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(line.decode("utf-8"))
    proc.wait()
    if proc.returncode != 0:
        sys.stderr.write(f"Error: Nextflow exited with code {proc.returncode}.\n")
        sys.exit(proc.returncode)

    # Collect annotated outputs generated by Nextflow
    annotated_folder = os.path.join(out_dir, "annotated")
    annotated_files = sorted(glob.glob(os.path.join(annotated_folder, "*.tsv")))
    if not annotated_files:
        sys.stderr.write(
            f"Error: No annotated TSV files found in '{annotated_folder}'.\n"
            "Check Nextflow logs and ensure the input glob is correct.\n"
        )
        sys.exit(1)

    # Satellite search against the annotated folder (correct inputpath)
    sat.process_files(
        file_list=annotated_files,
        use_cds=False,
        output_type=args.output_type,
        inputpath=annotated_folder,
        familiestosearchfor=families_to_search,
        outDir=out_dir,
    )

    # Convert RESULTS text outputs -> TSV
    post_process_results(out_dir)
    sys.exit(0)


if __name__ == "__main__":
    main()
