#!/usr/bin/env python3
"""In-memory beam-search optimizer for SaPhARI family-rule determination.

This standalone script infers SaPhARI-compatible family rules directly from
annotated TSV collections. Each family is represented by a rule consisting of:

- a positive marker list,
- a minimum marker count (``minnumber``), and
- an optional forbidden-marker list.

The optimizer parses the annotated training set once, evaluates candidate rules
in memory, and searches the resulting rule space with a family-level beam
search. The implementation is designed for parameter determination workflows in
which the primary objective is to reduce one-vs-rest false positives and false
negatives while preserving interpretable, JSON-ready rules.

Algorithmic features
--------------------
- single-pass parsing of annotated TSV inputs;
- direct optimization in ``Families.json`` rule space;
- exact in-memory TP/FP/FN/TN scoring for candidate rules;
- rival-aware penalties, including configurable pairwise confusion weights;
- joint optimization of positive markers, ``minnumber``, and forbidden markers;
- optional freezing of already well-performing families to focus search effort
  on difficult families.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from statistics import median
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from statistics import median
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Set, Tuple


# ---------------------------------------------------------------------------
# Canonicalization backbone
# ---------------------------------------------------------------------------
NORMALIZE_MAP = {
    "hypothetical protein": "hypothetical",
    "helix-turn-helix domain-containing protein": "htH",
    "helix-turn-helix transcriptional regulator": "htH",
    "HTH domain-containing protein": "htH",

    "tyrosine-type recombinase/integrase": "integrase",
    "site-specific recombinase": "integrase",
    "site-specific integrase": "integrase",
    "serine recombinase": "integrase",
    "recombinase": "integrase",
    "integrase": "integrase",
    "phage integrase": "integrase",

    "excisionase": "excisionase",
    "excisionase and transcriptional regulator": "excisionase",

    "phage major capsid protein": "capsid",
    "phage capsid protein": "capsid",
    "major capsid protein": "capsid",
    "capsid protein": "capsid",
    "major head protein": "major_head",

    "phage portal protein": "portal",
    "portal protein": "portal",

    "phage tail protein": "tail",
    "tail protein": "tail",

    "phage tail tube protein": "tail_tube",
    "tail tube protein": "tail_tube",

    "phage tail fiber protein": "tail_fiber",
    "tail fiber protein": "tail_fiber",

    "phage tape measure protein": "tape_measure",
    "tape measure protein": "tape_measure",

    "terminase large subunit": "terminase_large",
    "terminase small subunit": "terminase_small",
    "phage terminase small subunit P27 family": "terminase_small_p27",
    "P27 family phage terminase small subunit": "terminase_small_p27",

    "DNA helicase": "helicase",
    "ATP-dependent helicase": "helicase",
    "DEAD/DEAH box helicase": "helicase",

    "DNA primase": "primase",
    "DNA primase family protein": "primase",
    "primase": "primase",

    "single-stranded DNA-binding protein": "ssDNA_binding",
    "single stranded DNA-binding protein": "ssDNA_binding",
    "ssDNA-binding protein": "ssDNA_binding",

    "replication protein": "replication",
    "replication initiation protein": "replication",
    "replication protein RepA": "replication",
    "replication initiator protein A": "replication",

    "transcriptional regulator": "transcription_regulator",
}

STOPWORDS = {
    "protein", "domain", "containing", "like", "family", "subunit", "putative",
    "phage", "bacteriophage", "type", "dependent", "large", "small",
    "dna", "rna", "atp", "binding", "partial", "related", "predicted",
    "multispecies",
}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SampleRecord:
    family: str
    sample_id: str
    path: str
    present: frozenset[str]


@dataclass
class FamilyRule:
    family: str
    markers: List[str] = field(default_factory=list)
    minnumber: int = 3
    forbidden_cids: List[str] = field(default_factory=list)
    forbidden_names: List[str] = field(default_factory=list)
    ordered_pool: List[str] = field(default_factory=list)
    singleton_family: bool = False
    objective: float = float("-inf")
    metrics: Dict[str, float] = field(default_factory=dict)
    pairwise_fp_sources: Dict[str, int] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Parsing helpers (kept compatible with your old builder)
# ---------------------------------------------------------------------------
def strip_taxon_brackets(subject: str) -> str:
    if not subject:
        return ""
    return re.sub(r"\s+\[[^\]]+\]\s*$", "", subject).strip()


def parse_accession(subject_no_taxon: str) -> str:
    """
    Extract the leading protein accession from an annotation subject.
    Examples: WP_..., YP_..., NP_..., XP_...
    """
    if not subject_no_taxon:
        return ""
    m = re.match(r"^((?:WP|YP|NP|XP)_[0-9]+\.[0-9]+)\b", subject_no_taxon.strip())
    return m.group(1) if m else ""

def parse_raw_name(subject_no_taxon: str) -> str:
    """
    Remove the leading accession and annotation metadata prefixes,
    leaving the broad protein name.
    """
    if not subject_no_taxon:
        return ""
    s = subject_no_taxon.strip()
    s = re.sub(r"^(?:WP|YP|NP|XP)_[0-9]+\.[0-9]+\s+", "", s)
    s = re.sub(r"^MULTISPECIES:\s*", "", s, flags=re.IGNORECASE)
    return s.strip()

def normalize_function_label(raw_name: str) -> str:
    if not raw_name:
        return ""

    key = raw_name.strip()
    key = re.sub(r"^MULTISPECIES:\s*", "", key, flags=re.IGNORECASE)

    if key in NORMALIZE_MAP:
        return NORMALIZE_MAP[key]
    for k, v in NORMALIZE_MAP.items():
        if k.lower() == key.lower():
            return v

    s = key.lower()
    if "hypothetical" in s:
        return ""

    toks = re.split(r"[^a-z0-9]+", s)
    toks = [t for t in toks if t and t not in STOPWORDS]
    if not toks:
        return ""
    return "_".join(toks)


def iter_subject_column(tsv_path: str) -> Iterable[str]:
    with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            yield cols[1]


def is_hypothetical(name: str) -> bool:
    """
    Return True for generic hypothetical protein annotations.
    These should not become broad functional markers.
    """
    if not name:
        return False
    s = name.lower()
    return "hypothetical protein" in s or s.strip() == "hypothetical"


def read_sample(tsv_path: str) -> List[Tuple[str, str, str, str]]:
    """
    Returns tuples: (cluster_id, raw_name, evidence_subject_no_taxon, accession)
      * hypothetical -> HYP:<accession>
      * else         -> FUNC:<normalized_function>
    """
    hits: List[Tuple[str, str, str, str]] = []
    for subj in iter_subject_column(tsv_path):
        ev = strip_taxon_brackets(subj)
        acc = parse_accession(ev)
        raw = parse_raw_name(ev)

        if not ev or not acc:
            continue

        if is_hypothetical(raw):
            cluster_id = f"HYP:{acc}"
            hits.append((cluster_id, raw.strip(), ev, acc))
        else:
            func = normalize_function_label(raw)
            if not func:
                continue
            cluster_id = f"FUNC:{func}"
            hits.append((cluster_id, raw.strip(), ev, acc))
    return hits


def discover_families(annotated_root: str) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for name in sorted(os.listdir(annotated_root)):
        fam_path = os.path.join(annotated_root, name)
        ann_dir = os.path.join(fam_path, "annotated")
        if os.path.isdir(fam_path) and os.path.isdir(ann_dir):
            out.append((name, ann_dir))
    return out


# ---------------------------------------------------------------------------
# Dataset loading
# ---------------------------------------------------------------------------
def prevalence(count: int, denom: int) -> float:
    return count / denom if denom > 0 else 0.0


def choose_representative_raw_name(raw_names: Set[str]) -> str:
    if not raw_names:
        return ""
    return sorted(raw_names, key=lambda s: (len(s), s.lower()))[0]


def clean_broad_marker_name(raw_name: str) -> str:
    """
    Clean a broad protein label for use as a SaPhARI marker entry.
    This keeps the searchable biological phrase, not the normalized FUNC key.
    """
    if not raw_name:
        return ""
    s = raw_name.strip()
    s = re.sub(r"^MULTISPECIES:\s*", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def extract_accessions_from_evidence(evidence_subjects: Set[str]) -> List[str]:
    """
    Pull accession-only marker entries from evidence strings.
    This avoids redundant entries like:
      WP_214396676.1 phage head completion protein
    and instead keeps:
      WP_214396676.1
    """
    accs: List[str] = []
    for ev in evidence_subjects:
        m = re.match(r"^((?:WP|YP|NP|XP)_[0-9]+\.[0-9]+)\b", ev.strip())
        if m:
            accs.append(m.group(1))
    return sorted(set(accs))


def dedupe_keep_order(items: List[str]) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    for x in items:
        x = x.strip()
        if not x or x in seen:
            continue
        seen.add(x)
        out.append(x)
    return out


def build_json_group(cid: str, raw_names: Set[str], evidence_subjects: Set[str]) -> Optional[List[str]]:
    """
    Build JSON marker entries optimized for current SaPhARI search behavior.

    Current SaPhARI flattens nested JSON marker groups before searching.
    Therefore, this output intentionally separates:
      1. broad functional marker names
      2. accession-only reference-specific markers

    For functionally named proteins:
      - keep one clean broad marker name
      - keep accession-only markers

    For hypothetical proteins:
      - keep accession-only markers only
      - do NOT emit broad "hypothetical protein"
    """
    accessions = extract_accessions_from_evidence(evidence_subjects)

    if cid.startswith("FUNC:"):
        broad_candidates = [
            clean_broad_marker_name(x)
            for x in sorted(raw_names, key=lambda s: (len(s), s.lower()))
            if x
        ]

        broad_candidates = [
            x for x in broad_candidates
            if x and not is_hypothetical(x)
        ]

        items: List[str] = []

        if broad_candidates:
            items.append(broad_candidates[0])

        items.extend(accessions)

        out = dedupe_keep_order(items)
        return out if out else None

    if cid.startswith("HYP:"):
        out = dedupe_keep_order(accessions)
        return out if out else None

    return None


def cluster_id_to_forbidden_name(
    cid: str,
    cluster_to_raw: Dict[str, Set[str]],
    cluster_to_evs: Dict[str, Set[str]],
    forbid_hyp_as_accession: bool = False,
) -> Optional[str]:
    if cid.startswith("FUNC:"):
        rep = choose_representative_raw_name(cluster_to_raw.get(cid, set()))
        return rep or None
    if cid.startswith("HYP:") and forbid_hyp_as_accession:
        return cid.split("HYP:", 1)[1].strip() or None
    return None


def load_dataset(annotated_root: str):
    fam_dirs = discover_families(annotated_root)
    if not fam_dirs:
        raise SystemExit(f"No family folders found under: {annotated_root}")

    family_to_samples: Dict[str, List[SampleRecord]] = defaultdict(list)
    family_counts: Dict[str, Counter] = defaultdict(Counter)
    family_union: Dict[str, Set[str]] = defaultdict(set)
    family_to_paths: Dict[str, List[str]] = defaultdict(list)
    global_counts: Counter = Counter()
    cluster_to_raw: DefaultDict[str, Set[str]] = defaultdict(set)
    cluster_to_evs: DefaultDict[str, Set[str]] = defaultdict(set)
    all_samples: List[SampleRecord] = []

    for fam, ann_dir in fam_dirs:
        tsvs = sorted(glob.glob(os.path.join(ann_dir, "*.tsv")))
        for tsv in tsvs:
            hits = read_sample(tsv)
            present = frozenset(cid for cid, *_ in hits)
            rec = SampleRecord(family=fam, sample_id=os.path.basename(tsv), path=tsv, present=present)
            family_to_samples[fam].append(rec)
            family_to_paths[fam].append(tsv)
            all_samples.append(rec)

            for cid, raw, ev, _acc in hits:
                if raw:
                    cluster_to_raw[cid].add(raw)
                if ev:
                    cluster_to_evs[cid].add(ev)

            for cid in present:
                family_counts[fam][cid] += 1
                family_union[fam].add(cid)
                global_counts[cid] += 1

    family_sizes = {fam: len(samples) for fam, samples in family_to_samples.items()}
    total_samples = len(all_samples)
    if total_samples == 0:
        raise SystemExit("No TSV samples found across any families.")

    return {
        "family_to_samples": family_to_samples,
        "family_counts": family_counts,
        "family_union": family_union,
        "family_to_paths": family_to_paths,
        "family_sizes": family_sizes,
        "global_counts": global_counts,
        "cluster_to_raw": cluster_to_raw,
        "cluster_to_evs": cluster_to_evs,
        "all_samples": all_samples,
        "total_samples": total_samples,
        "families": sorted(family_to_samples.keys()),
    }


# ---------------------------------------------------------------------------
# Evaluation utilities
# ---------------------------------------------------------------------------
def safe_div(a: float, b: float) -> float:
    return a / b if b else 0.0


def count_hits(sample: SampleRecord, markers: Set[str]) -> int:
    return len(sample.present & markers)


def sample_predicted_positive(sample: SampleRecord, markers: Set[str], m_required: int, forbidden: Set[str]) -> bool:
    if forbidden and (sample.present & forbidden):
        return False
    return count_hits(sample, markers) >= m_required


def evaluate_rule(
    family: str,
    markers: Sequence[str],
    m_required: int,
    forbidden_cids: Sequence[str],
    all_samples: Sequence[SampleRecord],
) -> Dict[str, object]:
    marker_set = set(markers)
    forbidden_set = set(forbidden_cids)

    tp = fp = fn = tn = 0
    pairwise_fp_sources: Counter = Counter()
    pairwise_fn_targets: Counter = Counter()
    true_support_counts: List[int] = []
    out_support_counts: List[int] = []
    fp_samples_by_source: DefaultDict[str, List[SampleRecord]] = defaultdict(list)

    for sample in all_samples:
        hits = count_hits(sample, marker_set)
        pred = sample_predicted_positive(sample, marker_set, m_required, forbidden_set)
        truth = sample.family == family

        if truth:
            true_support_counts.append(hits)
        else:
            out_support_counts.append(hits)

        if pred and truth:
            tp += 1
        elif pred and not truth:
            fp += 1
            pairwise_fp_sources[sample.family] += 1
            fp_samples_by_source[sample.family].append(sample)
        elif (not pred) and truth:
            fn += 1
            pairwise_fn_targets[sample.family] += 1
        else:
            tn += 1

    prec = safe_div(tp, tp + fp)
    rec = safe_div(tp, tp + fn)
    spec = safe_div(tn, tn + fp)
    fpr = safe_div(fp, fp + tn)
    f1 = safe_div(2 * prec * rec, prec + rec) if (prec + rec) else 0.0
    num = (tp * tn) - (fp * fn)
    den = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    mcc = num / math.sqrt(den) if den > 0 else 0.0

    return {
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "precision": prec,
        "recall": rec,
        "specificity": spec,
        "fpr": fpr,
        "f1": f1,
        "mcc": mcc,
        "true_support_counts": true_support_counts,
        "out_support_counts": out_support_counts,
        "pairwise_fp_sources": dict(pairwise_fp_sources),
        "fp_samples_by_source": dict(fp_samples_by_source),
    }
@dataclass(order=True)
class CandidateState:
    """Candidate family rule tracked within the beam search."""
    sort_key: Tuple = field(init=False, repr=False)
    objective: float
    family: str = field(compare=False)
    markers: Tuple[str, ...] = field(compare=False)
    minnumber: int = field(compare=False)
    forbids: Tuple[str, ...] = field(compare=False)
    metrics: Dict[str, float] = field(compare=False, default_factory=dict)
    pairwise_fp_sources: Dict[str, int] = field(compare=False, default_factory=dict)
    notes: List[str] = field(compare=False, default_factory=list)

    def __post_init__(self) -> None:
        self.sort_key = (
            self.objective,
            self.metrics.get("recall", 0.0),
            self.metrics.get("precision", 0.0),
            self.metrics.get("specificity", 0.0),
            self.minnumber,
            len(self.markers),
            -self.metrics.get("fp", 0),
            -self.metrics.get("fn", 0),
            self.metrics.get("mcc", 0.0),
        )


# ---------------------------------------------------------------------------
# General utilities
# ---------------------------------------------------------------------------
def safe_div(a: float, b: float) -> float:
    """Return ``a / b`` and guard against division by zero."""
    return a / b if b else 0.0


def quantile(values: Sequence[int], q: float) -> float:
    """Compute a linear-interpolated empirical quantile."""
    if not values:
        return 0.0
    xs = sorted(values)
    if len(xs) == 1:
        return float(xs[0])
    q = max(0.0, min(1.0, q))
    pos = (len(xs) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(xs[lo])
    frac = pos - lo
    return xs[lo] * (1.0 - frac) + xs[hi] * frac


def unique_keep_order(items: Iterable[str]) -> List[str]:
    """Deduplicate an iterable while preserving first-occurrence order."""
    out: List[str] = []
    seen: Set[str] = set()
    for x in items:
        if x in seen:
            continue
        seen.add(x)
        out.append(x)
    return out


def parse_pair_weights(raw: str) -> Dict[Tuple[str, str], float]:
    """Parse pairwise confusion weights of the form ``SRC->DST=value``."""
    out: Dict[Tuple[str, str], float] = {}
    if not raw:
        return out
    for part in raw.split(","):
        part = part.strip()
        if not part:
            continue
        try:
            lhs, weight_str = part.split("=")
            src, dst = lhs.split("->")
            out[(src.strip(), dst.strip())] = float(weight_str)
        except Exception as exc:
            raise SystemExit(f"Could not parse --pair-weights entry '{part}': {exc}")
    return out


# ---------------------------------------------------------------------------
# Rival discovery and pairwise confusion accounting
# ---------------------------------------------------------------------------
def compute_confusion_matrix(rules: Dict[str, CandidateState], dataset: Dict[str, object]) -> Dict[str, Dict[str, int]]:
    """Count cross-family false positives induced by the current rule set."""
    out: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    all_samples: List[SampleRecord] = dataset["all_samples"]
    for target_family, rule in rules.items():
        if not rule.markers or rule.minnumber < 3:
            continue
        marker_set = set(rule.markers)
        forbid_set = set(rule.forbids)
        for sample in all_samples:
            if sample.family == target_family:
                continue
            if sample_predicted_positive(sample, marker_set, rule.minnumber, forbid_set):
                out[target_family][sample.family] += 1
    return {fam: dict(srcs) for fam, srcs in out.items()}


def initial_overlap_rivals(family: str, dataset: Dict[str, object], top_n: int) -> List[str]:
    """Seed a rival list using feature overlap before confusion-driven refinement."""
    fam_union = dataset["family_union"].get(family, set())
    if not fam_union:
        return []
    n_f = dataset["family_sizes"].get(family, 0)
    if n_f <= 0:
        return []
    fam_counts = dataset["family_counts"]
    fam_sizes = dataset["family_sizes"]

    scored: List[Tuple[float, str]] = []
    for other in dataset["families"]:
        if other == family:
            continue
        other_union = dataset["family_union"].get(other, set())
        if not other_union:
            continue
        inter = fam_union & other_union
        if not inter:
            continue
        n_o = fam_sizes.get(other, 0)
        if n_o <= 0:
            continue
        overlap = 0.0
        for cid in inter:
            p_f = prevalence(fam_counts[family].get(cid, 0), n_f)
            p_o = prevalence(fam_counts[other].get(cid, 0), n_o)
            overlap += min(p_f, p_o)
        score = overlap / max(1.0, math.sqrt(len(fam_union) * len(other_union)))
        scored.append((score, other))
    scored.sort(key=lambda kv: (-kv[0], kv[1]))
    return [fam for score, fam in scored[:top_n] if score > 0.0]


def merged_rivals(
    family: str,
    current_rule: Optional[CandidateState],
    confusion: Dict[str, Dict[str, int]],
    initial_rivals: Dict[str, List[str]],
    top_n: int,
    forced_pairs: Dict[Tuple[str, str], float],
) -> List[str]:
    """Combine confusion-derived, overlap-derived, and forced rival families."""
    sources = confusion.get(family, {})
    ranked = [src for src, cnt in sorted(sources.items(), key=lambda kv: (-kv[1], kv[0])) if cnt > 0]
    seed = ranked[:top_n] + initial_rivals.get(family, [])
    # Current rule might already know a useful confusion source.
    if current_rule is not None:
        seed += sorted(current_rule.pairwise_fp_sources, key=lambda s: (-current_rule.pairwise_fp_sources[s], s))[:top_n]
    # Force special boundary rivals.
    for (src, dst), _w in forced_pairs.items():
        if dst == family:
            seed.append(src)
        if src == family:
            seed.append(dst)
    return unique_keep_order(seed)[:max(top_n, 1) * 2]


# ---------------------------------------------------------------------------
# Marker statistics and candidate-pool construction
# ---------------------------------------------------------------------------
def cluster_stats_for_family(family: str, dataset: Dict[str, object], rivals: Sequence[str]) -> Dict[str, Dict[str, float]]:
    """Summarize within-family, outside-family, and rival prevalence per marker."""
    fam_sizes: Dict[str, int] = dataset["family_sizes"]
    fam_counts: Dict[str, Counter] = dataset["family_counts"]
    global_counts: Counter = dataset["global_counts"]
    total_samples: int = dataset["total_samples"]
    fam_union: Dict[str, Set[str]] = dataset["family_union"]

    n_in = fam_sizes.get(family, 0)
    n_out = total_samples - n_in
    counts_in = fam_counts[family]

    stats: Dict[str, Dict[str, float]] = {}
    for cid in fam_union.get(family, set()):
        c_in = counts_in.get(cid, 0)
        p_in = prevalence(c_in, n_in)
        c_global = global_counts.get(cid, 0)
        c_out = max(0, c_global - c_in)
        p_out = prevalence(c_out, n_out) if n_out > 0 else 0.0

        rival_vals: List[float] = []
        for rival in rivals:
            n_r = fam_sizes.get(rival, 0)
            if n_r <= 0:
                continue
            rival_vals.append(prevalence(fam_counts[rival].get(cid, 0), n_r))
        max_rival = max(rival_vals) if rival_vals else 0.0
        mean_rival = sum(rival_vals) / len(rival_vals) if rival_vals else 0.0
        stats[cid] = {
            "p_in": p_in,
            "p_out": p_out,
            "max_rival_prev": max_rival,
            "mean_rival_prev": mean_rival,
        }
    return stats


def build_positive_pool(
    family: str,
    dataset: Dict[str, object],
    rivals: Sequence[str],
    args: argparse.Namespace,
) -> List[str]:
    """Assemble a ranked candidate pool of positive family markers."""
    family_samples: List[SampleRecord] = dataset["family_to_samples"][family]
    other_samples: List[SampleRecord] = [s for s in dataset["all_samples"] if s.family != family]
    rival_set = set(rivals)
    rival_samples = [s for s in other_samples if s.family in rival_set]
    nonrival_samples = [s for s in other_samples if s.family not in rival_set]
    stats = cluster_stats_for_family(family, dataset, rivals)
    if not stats:
        return []

    candidates = list(stats.keys())
    if len(family_samples) > 1:
        filtered = [
            cid for cid in candidates
            if stats[cid]["p_in"] >= args.min_marker_support or stats[cid]["p_out"] <= args.max_marker_out
        ]
        if filtered:
            candidates = filtered

    fam_contains = {cid: [i for i, s in enumerate(family_samples) if cid in s.present] for cid in candidates}
    rival_contains = {cid: [i for i, s in enumerate(rival_samples) if cid in s.present] for cid in candidates}
    other_contains = {cid: [i for i, s in enumerate(nonrival_samples) if cid in s.present] for cid in candidates}

    true_cover = [0] * len(family_samples)
    rival_cover = [0] * len(rival_samples)
    other_cover = [0] * len(nonrival_samples)

    ordered: List[str] = []
    remaining: Set[str] = set(candidates)
    max_pool = min(args.positive_pool_cap, max(len(candidates), args.m_min))

    while remaining and len(ordered) < max_pool:
        best_cid = None
        best_score = float("-inf")
        for cid in remaining:
            st = stats[cid]
            marginal_true = sum(1.0 / (1.0 + true_cover[i]) for i in fam_contains.get(cid, []))
            marginal_rival = sum(1.0 / (1.0 + rival_cover[i]) for i in rival_contains.get(cid, []))
            marginal_other = sum(1.0 / (1.0 + other_cover[i]) for i in other_contains.get(cid, []))
            if len(family_samples) == 1:
                score = (
                    args.singleton_disc_weight * (1.0 - st["p_out"])
                    - args.singleton_rival_weight * st["max_rival_prev"]
                    + args.singleton_true_weight * marginal_true
                    - args.singleton_other_weight * marginal_other
                )
            else:
                score = (
                    args.cover_true_weight * marginal_true
                    + args.support_in_weight * st["p_in"]
                    - args.support_out_weight * st["p_out"]
                    - args.rival_prev_weight * st["max_rival_prev"]
                    - args.cover_rival_weight * marginal_rival
                    - args.cover_other_weight * marginal_other
                )
            if score > best_score:
                best_score = score
                best_cid = cid
        if best_cid is None:
            break
        ordered.append(best_cid)
        remaining.remove(best_cid)
        for i in fam_contains.get(best_cid, []):
            true_cover[i] += 1
        for i in rival_contains.get(best_cid, []):
            rival_cover[i] += 1
        for i in other_contains.get(best_cid, []):
            other_cover[i] += 1
        if len(ordered) >= max(args.m_min, 3) and best_score < args.marker_stop_score:
            break

    if len(ordered) < max(args.m_min, 3):
        leftovers = sorted(
            (cid for cid in candidates if cid not in ordered),
            key=lambda cid: (
                stats[cid]["p_out"],
                stats[cid]["max_rival_prev"],
                -stats[cid]["p_in"],
                cid,
            )
        )
        for cid in leftovers:
            if len(ordered) >= max(args.m_min, 3):
                break
            ordered.append(cid)

    return ordered


def build_static_forbid_pool(
    family: str,
    rivals: Sequence[str],
    dataset: Dict[str, object],
    args: argparse.Namespace,
    exclude_markers: Optional[Set[str]] = None,
) -> List[str]:
    """Rank forbidden-marker candidates using prevalence and rival enrichment."""
    exclude_markers = exclude_markers or set()
    fam_sizes = dataset["family_sizes"]
    fam_counts = dataset["family_counts"]
    total_samples = dataset["total_samples"]
    global_counts = dataset["global_counts"]
    n_in = fam_sizes.get(family, 0)
    if n_in <= 0:
        return []
    n_out = total_samples - n_in
    search_space: Set[str] = set()
    if rivals:
        for rival in rivals:
            search_space.update(dataset["family_union"].get(rival, set()))
    else:
        for other in dataset["families"]:
            if other != family:
                search_space.update(dataset["family_union"].get(other, set()))

    ranked: List[Tuple[float, str]] = []
    for cid in search_space:
        if cid in exclude_markers:
            continue
        if (not args.forbid_include_hyp) and cid.startswith("HYP:"):
            continue
        p_target = prevalence(fam_counts[family].get(cid, 0), n_in)
        if p_target > args.forbid_target_prev_max:
            continue
        rival_vals = []
        for rival in rivals:
            n_r = fam_sizes.get(rival, 0)
            if n_r <= 0:
                continue
            rival_vals.append(prevalence(fam_counts[rival].get(cid, 0), n_r))
        max_rival = max(rival_vals) if rival_vals else 0.0
        mean_rival = sum(rival_vals) / len(rival_vals) if rival_vals else 0.0
        c_global = global_counts.get(cid, 0)
        c_out = max(0, c_global - fam_counts[family].get(cid, 0))
        p_out = prevalence(c_out, n_out) if n_out > 0 else 0.0
        if max_rival < args.forbid_rival_prev_min and p_out < args.forbid_out_prev_min:
            continue
        score = (
            args.forbid_rival_weight * max_rival
            + args.forbid_mean_rival_weight * mean_rival
            + args.forbid_out_weight * p_out
            - args.forbid_target_weight * p_target
        )
        ranked.append((score, cid))
    ranked.sort(key=lambda kv: (-kv[0], kv[1]))
    return [cid for _score, cid in ranked[:args.forbid_pool_cap]]


def build_fp_forbid_pool(
    family: str,
    current_rule: Optional[CandidateState],
    dataset: Dict[str, object],
    args: argparse.Namespace,
    exclude_markers: Optional[Set[str]] = None,
) -> List[str]:
    """Rank forbidden markers using the current family false positives."""
    exclude_markers = exclude_markers or set()
    if current_rule is None or not current_rule.markers:
        return []

    metrics = evaluate_rule(family, current_rule.markers, current_rule.minnumber, current_rule.forbids, dataset["all_samples"])
    fp_samples_by_source = metrics.get("fp_samples_by_source", {})
    fam_sizes = dataset["family_sizes"]
    fam_counts = dataset["family_counts"]
    n_in = fam_sizes.get(family, 0)
    if n_in <= 0:
        return []

    candidates: Set[str] = set()
    fp_samples: List[SampleRecord] = []
    for src, samples in fp_samples_by_source.items():
        fp_samples.extend(samples)
    for sample in fp_samples:
        candidates.update(sample.present)

    ranked: List[Tuple[float, str]] = []
    fp_n = len(fp_samples)
    for cid in candidates:
        if cid in exclude_markers:
            continue
        if (not args.forbid_include_hyp) and cid.startswith("HYP:"):
            continue
        p_target = prevalence(fam_counts[family].get(cid, 0), n_in)
        if p_target > args.forbid_target_prev_max:
            continue
        fp_prev = safe_div(sum(1 for s in fp_samples if cid in s.present), fp_n)
        if fp_prev < args.forbid_fp_prev_min:
            continue
        score = args.forbid_fp_prev_weight * fp_prev - args.forbid_target_weight * p_target
        ranked.append((score, cid))
    ranked.sort(key=lambda kv: (-kv[0], kv[1]))
    return [cid for _score, cid in ranked[:args.forbid_pool_cap]]


# ---------------------------------------------------------------------------
# Objective and rule evaluation
# ---------------------------------------------------------------------------
def family_objective(
    family: str,
    metrics: Dict[str, object],
    k: int,
    m: int,
    forbids: Sequence[str],
    dataset: Dict[str, object],
    args: argparse.Namespace,
    target_recall: float,
    forced_pair_weights: Dict[Tuple[str, str], float],
) -> float:
    """Compute the family-level optimization objective for a candidate rule."""
    n_in = dataset["family_sizes"].get(family, 0)
    n_out = dataset["total_samples"] - n_in
    rec = float(metrics["recall"])
    prec = float(metrics["precision"])
    spec = float(metrics["specificity"])
    f1 = float(metrics["f1"])
    fp = int(metrics["fp"])
    fn = int(metrics["fn"])
    true_support = list(metrics.get("true_support_counts", []))

    fp_rate = safe_div(fp, n_out)
    fn_rate = safe_div(fn, n_in)

    pairwise_pen = 0.0
    for src, cnt in metrics.get("pairwise_fp_sources", {}).items():
        src_n = dataset["family_sizes"].get(src, 0)
        base_rate = cnt / src_n if src_n > 0 else 0.0
        weight = forced_pair_weights.get((src, family), args.default_pairwise_weight)
        if family == src:
            continue
        pairwise_pen += weight * base_rate

    ratio = safe_div(m, max(1, k))
    loose_pen = max(0.0, args.loose_ratio_target - ratio)
    q_low = quantile(true_support, args.support_quantile_low) if true_support else float(m)
    strict_pen = max(0.0, (m - q_low) / max(1.0, k))
    recall_shortfall = max(0.0, target_recall - rec)

    score = (
        args.obj_recall_weight * rec
        + args.obj_precision_weight * prec
        + args.obj_specificity_weight * spec
        + args.obj_f1_weight * f1
        + args.obj_m_abs_weight * m
        + args.obj_m_ratio_weight * ratio
        + args.obj_k_weight * math.log1p(k)
        + args.obj_forbid_count_weight * math.log1p(len(forbids))
        - args.obj_fp_weight * fp_rate
        - args.obj_fn_weight * fn_rate
        - args.obj_pairwise_weight * pairwise_pen
        - args.obj_loose_shape_weight * loose_pen
        - args.obj_strict_shape_weight * strict_pen
        - args.obj_recall_shortfall_weight * (recall_shortfall ** 2)
    )
    return score


def optimize_m_for_rule(
    family: str,
    markers: Sequence[str],
    forbids: Sequence[str],
    dataset: Dict[str, object],
    args: argparse.Namespace,
    target_recall: float,
    forced_pair_weights: Dict[Tuple[str, str], float],
    note_prefix: Optional[str] = None,
) -> CandidateState:
    """Select the best ``minnumber`` for a fixed marker and forbidden set."""
    if len(markers) < args.m_min:
        return CandidateState(
            objective=float("-inf"),
            family=family,
            markers=tuple(markers),
            minnumber=args.m_min,
            forbids=tuple(forbids),
            metrics={"precision": 0.0, "recall": 0.0, "specificity": 1.0, "f1": 0.0, "mcc": 0.0, "tp": 0, "fp": 0, "fn": dataset["family_sizes"].get(family, 0), "tn": dataset["total_samples"] - dataset["family_sizes"].get(family, 0)},
            pairwise_fp_sources={},
            notes=[f"too_few_markers<{args.m_min}"] if note_prefix is None else [note_prefix, f"too_few_markers<{args.m_min}"]
        )

    best: Optional[CandidateState] = None
    max_m = min(args.m_max, len(markers))
    for m in range(args.m_min, max_m + 1):
        mets = evaluate_rule(family, markers, m, forbids, dataset["all_samples"])
        score = family_objective(family, mets, len(markers), m, forbids, dataset, args, target_recall, forced_pair_weights)
        cand = CandidateState(
            objective=score,
            family=family,
            markers=tuple(markers),
            minnumber=m,
            forbids=tuple(forbids),
            metrics={
                "precision": float(mets["precision"]),
                "recall": float(mets["recall"]),
                "specificity": float(mets["specificity"]),
                "fpr": float(mets["fpr"]),
                "f1": float(mets["f1"]),
                "mcc": float(mets["mcc"]),
                "tp": int(mets["tp"]), "fp": int(mets["fp"]), "fn": int(mets["fn"]), "tn": int(mets["tn"]),
                "support_min": min(mets["true_support_counts"]) if mets["true_support_counts"] else 0,
                "support_median": float(median(mets["true_support_counts"])) if mets["true_support_counts"] else 0.0,
                "support_max": max(mets["true_support_counts"]) if mets["true_support_counts"] else 0,
            },
            pairwise_fp_sources=dict(mets.get("pairwise_fp_sources", {})),
            notes=[] if note_prefix is None else [note_prefix],
        )
        if best is None or cand.sort_key > best.sort_key:
            best = cand
    assert best is not None
    return best


# ---------------------------------------------------------------------------
# Seed building and beam search
# ---------------------------------------------------------------------------
def build_seed_rules(
    family: str,
    positive_pool: Sequence[str],
    forbid_pool: Sequence[str],
    dataset: Dict[str, object],
    args: argparse.Namespace,
    target_recall: float,
    forced_pair_weights: Dict[Tuple[str, str], float],
    current_rule: Optional[CandidateState] = None,
) -> List[CandidateState]:
    """Generate the initial beam-search seed states for a family."""
    seeds: List[CandidateState] = []
    prefix_sizes = sorted({
        max(args.m_min, 3),
        min(len(positive_pool), max(args.m_min, 3) + 2),
        min(len(positive_pool), max(args.m_min, 3) + 4),
        min(len(positive_pool), 10),
        min(len(positive_pool), 14),
        min(len(positive_pool), args.k_max),
    })
    for k in prefix_sizes:
        if k < args.m_min:
            continue
        seeds.append(optimize_m_for_rule(
            family, list(positive_pool[:k]), [], dataset, args, target_recall, forced_pair_weights, note_prefix=f"seed_prefix_k={k}"
        ))
        # Try first few forbids jointly from the beginning.
        for fcount in range(1, min(args.seed_forbid_size, len(forbid_pool)) + 1):
            seeds.append(optimize_m_for_rule(
                family, list(positive_pool[:k]), list(forbid_pool[:fcount]), dataset, args, target_recall, forced_pair_weights,
                note_prefix=f"seed_prefix_k={k}_forb={fcount}"
            ))
    if current_rule is not None and current_rule.markers:
        seeds.append(current_rule)
    # Dedup by state.
    uniq: Dict[Tuple[Tuple[str, ...], Tuple[str, ...], int], CandidateState] = {}
    for cand in seeds:
        key = (tuple(cand.markers), tuple(cand.forbids), cand.minnumber)
        if key not in uniq or cand.sort_key > uniq[key].sort_key:
            uniq[key] = cand
    ranked = sorted(uniq.values(), reverse=True)
    return ranked[:args.beam_width]


def beam_search_family(
    family: str,
    current_rule: Optional[CandidateState],
    frozen_rules: Dict[str, CandidateState],
    dataset: Dict[str, object],
    rivals: Sequence[str],
    args: argparse.Namespace,
    target_recall: float,
    forced_pair_weights: Dict[Tuple[str, str], float],
    hard_family: bool,
) -> CandidateState:
    """Optimize one family rule with a local beam search."""
    positive_pool = build_positive_pool(family, dataset, rivals, args)
    if len(positive_pool) < args.m_min:
        # Hard fallback for tiny families.
        fallback_markers = list(positive_pool[:max(args.m_min, min(len(positive_pool), args.k_max))])
        if len(fallback_markers) < args.m_min:
            # last resort: take all observed family-union markers ranked by p_out asc.
            fam_union = list(dataset["family_union"].get(family, set()))
            fam_union.sort(key=lambda cid: (
                prevalence(dataset["global_counts"].get(cid, 0) - dataset["family_counts"][family].get(cid, 0), max(1, dataset["total_samples"] - dataset["family_sizes"][family])),
                cid,
            ))
            fallback_markers = fam_union[:args.m_min]
        cand = optimize_m_for_rule(family, fallback_markers, [], dataset, args, target_recall, forced_pair_weights, note_prefix="fallback_tiny_family")
        cand.notes.append("low_confidence_fallback")
        return cand

    static_forbids = build_static_forbid_pool(family, rivals, dataset, args, exclude_markers=set(positive_pool))
    fp_forbids = build_fp_forbid_pool(family, current_rule, dataset, args, exclude_markers=set(positive_pool))
    forbid_pool = unique_keep_order(list(fp_forbids) + list(static_forbids))[:args.forbid_pool_cap]

    beam = build_seed_rules(
        family, positive_pool, forbid_pool, dataset, args, target_recall, forced_pair_weights, current_rule=current_rule
    )
    if not beam:
        return optimize_m_for_rule(family, list(positive_pool[:max(args.m_min, 3)]), [], dataset, args, target_recall, forced_pair_weights)

    iters = args.hard_beam_iters if hard_family else args.easy_beam_iters
    pos_ops = args.hard_neighbor_add if hard_family else args.easy_neighbor_add
    forbid_ops = args.hard_forbid_add if hard_family else args.easy_forbid_add

    best_overall = max(beam)
    for _ in range(iters):
        candidates: Dict[Tuple[Tuple[str, ...], Tuple[str, ...]], CandidateState] = {}
        for state in beam:
            marker_set = set(state.markers)
            forbid_set = set(state.forbids)

            def push(markers: List[str], forbids: List[str], note: str) -> None:
                if len(markers) < args.m_min:
                    return
                markers = unique_keep_order(markers)[:args.k_max]
                forbids = unique_keep_order(f for f in forbids if f not in markers)[:args.max_forbids]
                key = (tuple(sorted(markers)), tuple(sorted(forbids)))
                cand = optimize_m_for_rule(family, markers, forbids, dataset, args, target_recall, forced_pair_weights, note_prefix=note)
                prev = candidates.get(key)
                if prev is None or cand.sort_key > prev.sort_key:
                    candidates[key] = cand

            # keep current
            push(list(state.markers), list(state.forbids), "keep")

            # add positive markers
            adds = [cid for cid in positive_pool if cid not in marker_set][:pos_ops]
            for cid in adds:
                push(list(state.markers) + [cid], list(state.forbids), f"add_marker:{cid}")

            # remove weakest markers (toward tail of ordered pool)
            if len(state.markers) > args.m_min:
                removable = list(state.markers[-min(args.remove_marker_cap, len(state.markers)):])
                for cid in removable:
                    push([x for x in state.markers if x != cid], list(state.forbids), f"remove_marker:{cid}")

            # swap in top pool candidates for tail markers
            if len(state.markers) > args.m_min:
                removable = list(state.markers[-min(args.swap_marker_cap, len(state.markers)):])
                for rem in removable:
                    for add in adds[:args.swap_add_cap]:
                        if add == rem:
                            continue
                        push([x for x in state.markers if x != rem] + [add], list(state.forbids), f"swap_marker:{rem}->{add}")

            # add forbids
            forbid_adds = [cid for cid in forbid_pool if cid not in forbid_set and cid not in marker_set][:forbid_ops]
            for cid in forbid_adds:
                push(list(state.markers), list(state.forbids) + [cid], f"add_forbid:{cid}")

            # remove forbids
            for cid in list(state.forbids)[-min(args.remove_forbid_cap, len(state.forbids)):]:
                push(list(state.markers), [x for x in state.forbids if x != cid], f"remove_forbid:{cid}")

            # swap forbids
            for rem in list(state.forbids)[-min(args.swap_forbid_cap, len(state.forbids)):]:
                for add in forbid_adds[:args.swap_add_cap]:
                    if add == rem:
                        continue
                    push(list(state.markers), [x for x in state.forbids if x != rem] + [add], f"swap_forbid:{rem}->{add}")

        if not candidates:
            break
        ranked = sorted(candidates.values(), reverse=True)
        beam = ranked[:args.beam_width]
        if beam[0].sort_key > best_overall.sort_key:
            best_overall = beam[0]

    # Final stronger forbid refinement around the best rule.
    refined_forbid_pool = unique_keep_order(
        build_fp_forbid_pool(family, best_overall, dataset, args, exclude_markers=set(best_overall.markers)) +
        build_static_forbid_pool(family, rivals, dataset, args, exclude_markers=set(best_overall.markers))
    )[:args.final_forbid_pool_cap]
    cur_best = best_overall
    improved = True
    while improved and len(cur_best.forbids) < args.max_forbids:
        improved = False
        for cid in refined_forbid_pool:
            if cid in cur_best.forbids or cid in cur_best.markers:
                continue
            cand = optimize_m_for_rule(
                family,
                list(cur_best.markers),
                list(cur_best.forbids) + [cid],
                dataset,
                args,
                target_recall,
                forced_pair_weights,
                note_prefix=f"greedy_forbid:{cid}",
            )
            if cand.sort_key > cur_best.sort_key:
                cur_best = cand
                improved = True
                break

    # Refit positives with forbids fixed using a narrower search over pool prefixes + swaps.
    final_candidates: List[CandidateState] = [cur_best]
    prefix_sizes = sorted(set([max(args.m_min, 3), len(cur_best.markers), min(args.k_max, len(cur_best.markers) + 2), min(args.k_max, 10), min(args.k_max, 14)]))
    for k in prefix_sizes:
        if k < args.m_min or k > len(positive_pool):
            continue
        final_candidates.append(optimize_m_for_rule(
            family, list(positive_pool[:k]), list(cur_best.forbids), dataset, args, target_recall, forced_pair_weights,
            note_prefix=f"refit_prefix_k={k}"
        ))
    # Also keep a version with the strongest current markers first.
    for k in range(max(args.m_min, len(cur_best.markers) - 2), min(args.k_max, len(cur_best.markers) + 2) + 1):
        if k < args.m_min or k > len(positive_pool):
            continue
        pool_slice = unique_keep_order(list(cur_best.markers[:min(len(cur_best.markers), k)]) + [cid for cid in positive_pool if cid not in set(cur_best.markers)])[:k]
        final_candidates.append(optimize_m_for_rule(
            family, pool_slice, list(cur_best.forbids), dataset, args, target_recall, forced_pair_weights,
            note_prefix=f"refit_mix_k={k}"
        ))

    best_final = max(final_candidates)
    return best_final


# ---------------------------------------------------------------------------
# Top-level optimization
# ---------------------------------------------------------------------------
def initialize_rules(dataset: Dict[str, object], args: argparse.Namespace, forced_pair_weights: Dict[Tuple[str, str], float]) -> Tuple[Dict[str, CandidateState], Dict[str, List[str]]]:
    """Construct initial family rules and overlap-based rival seeds."""
    initial_rivals = {fam: initial_overlap_rivals(fam, dataset, args.initial_rival_top_n) for fam in dataset["families"]}
    rules: Dict[str, CandidateState] = {}
    for family in dataset["families"]:
        rivals = initial_rivals.get(family, [])
        positive_pool = build_positive_pool(family, dataset, rivals, args)
        target_recall = args.singleton_target_recall if dataset["family_sizes"].get(family, 0) == 1 else args.target_recall
        static_forbids = build_static_forbid_pool(family, rivals, dataset, args, exclude_markers=set(positive_pool))
        seeds = build_seed_rules(family, positive_pool, static_forbids, dataset, args, target_recall, forced_pair_weights)
        if seeds:
            rules[family] = max(seeds)
        else:
            rules[family] = optimize_m_for_rule(family, list(positive_pool[:max(args.m_min, 3)]), [], dataset, args, target_recall, forced_pair_weights, note_prefix="init_fallback")
    return rules, initial_rivals


def optimize_rules(dataset: Dict[str, object], args: argparse.Namespace) -> Dict[str, CandidateState]:
    """Run global optimization passes and return the final family rules."""
    forced_pair_weights = parse_pair_weights(args.pair_weights)
    rules, initial_rivals = initialize_rules(dataset, args, forced_pair_weights)

    for pass_idx in range(args.global_passes):
        confusion = compute_confusion_matrix(rules, dataset)
        new_rules: Dict[str, CandidateState] = {}
        for family in dataset["families"]:
            current = rules.get(family)
            rivals = merged_rivals(family, current, confusion, initial_rivals, args.rival_top_n, forced_pair_weights)
            target_recall = args.singleton_target_recall if dataset["family_sizes"].get(family, 0) == 1 else args.target_recall
            hard_family = True
            if current is not None:
                if current.metrics.get("precision", 0.0) >= args.freeze_precision and current.metrics.get("recall", 0.0) >= args.freeze_recall:
                    hard_family = False
            if (family in args.always_optimize_families):
                hard_family = True
            if (family in args.never_freeze_families):
                hard_family = True

            if not hard_family and current is not None:
                frozen = current
                frozen.notes = list(frozen.notes) + [f"pass{pass_idx+1}:frozen"]
                new_rules[family] = frozen
                continue

            optimized = beam_search_family(
                family=family,
                current_rule=current,
                frozen_rules=new_rules,
                dataset=dataset,
                rivals=rivals,
                args=args,
                target_recall=target_recall,
                forced_pair_weights=forced_pair_weights,
                hard_family=hard_family,
            )
            optimized.notes = list(optimized.notes) + [f"pass{pass_idx+1}:rivals={','.join(rivals) if rivals else 'none'}"]
            new_rules[family] = optimized
        rules = new_rules
    return rules


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
def write_outputs(rules: Dict[str, CandidateState], dataset: Dict[str, object], args: argparse.Namespace) -> None:
    """Write the optimized rule set and supporting audit files."""
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    out_json = os.path.join(out_dir, args.emit_families_json)
    out_best = os.path.join(out_dir, "family_rules_best.tsv")
    out_readable = os.path.join(out_dir, "family_rules_readable.txt")
    out_conf = os.path.join(out_dir, "pairwise_confusion.tsv")
    out_backmap = os.path.join(out_dir, "cluster_backmap.tsv")

    cluster_to_raw = dataset["cluster_to_raw"]
    cluster_to_evs = dataset["cluster_to_evs"]

    families_json: List[Dict[str, object]] = []
    for family in dataset["families"]:
        rule = rules[family]
        groups: List[List[str]] = []
        for cid in rule.markers:
            grp = build_json_group(cid, cluster_to_raw.get(cid, set()), cluster_to_evs.get(cid, set()))
            if grp:
                groups.append(grp)
        if not groups:
            # Last resort JSON salvage: use any marker with available evidence.
            salvage = []
            for cid in dataset["family_union"].get(family, set()):
                grp = build_json_group(cid, cluster_to_raw.get(cid, set()), cluster_to_evs.get(cid, set()))
                if grp:
                    salvage.append(grp)
                if len(salvage) >= max(args.m_min, 3):
                    break
            groups = salvage
        if not groups:
            continue
        forbidden_names: List[str] = []
        for cid in rule.forbids:
            nm = cluster_id_to_forbidden_name(cid, cluster_to_raw, cluster_to_evs, forbid_hyp_as_accession=args.forbid_hyp_as_accession)
            if nm:
                forbidden_names.append(nm)
        obj: Dict[str, object] = {
            "title": family,
            "proteins": groups,
            "length": args.default_length,
            "minnumber": int(rule.minnumber),
        }
        if forbidden_names:
            obj["forbidden"] = unique_keep_order(forbidden_names)
        families_json.append(obj)

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(families_json, f, indent=2)

    with open(out_best, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "family", "n_samples", "objective", "k", "m",
            "precision", "recall", "specificity", "fpr", "f1", "mcc",
            "tp", "fp", "fn", "tn",
            "support_min", "support_median", "support_max",
            "pairwise_fp_sources", "forbidden_count", "forbidden_cids", "markers", "notes",
        ])
        for family in dataset["families"]:
            rule = rules[family]
            mets = rule.metrics
            w.writerow([
                family, dataset["family_sizes"].get(family, 0), f"{rule.objective:.6f}",
                len(rule.markers), rule.minnumber,
                f"{mets.get('precision', 0.0):.6f}",
                f"{mets.get('recall', 0.0):.6f}",
                f"{mets.get('specificity', 0.0):.6f}",
                f"{mets.get('fpr', 0.0):.6f}",
                f"{mets.get('f1', 0.0):.6f}",
                f"{mets.get('mcc', 0.0):.6f}",
                mets.get("tp", 0), mets.get("fp", 0), mets.get("fn", 0), mets.get("tn", 0),
                mets.get("support_min", 0), mets.get("support_median", 0.0), mets.get("support_max", 0),
                json.dumps(rule.pairwise_fp_sources, sort_keys=True),
                len(rule.forbids), ", ".join(rule.forbids), ", ".join(rule.markers), " | ".join(rule.notes),
            ])

    with open(out_readable, "w", encoding="utf-8") as f:
        for family in dataset["families"]:
            rule = rules[family]
            mets = rule.metrics
            f.write(f"=== {family} ===\n")
            f.write(f"n_samples={dataset['family_sizes'].get(family, 0)}  objective={rule.objective:.6f}  k={len(rule.markers)}  m={rule.minnumber}\n")
            f.write(
                f"precision={mets.get('precision', 0.0):.4f}  recall={mets.get('recall', 0.0):.4f}  specificity={mets.get('specificity', 0.0):.4f}  fpr={mets.get('fpr', 0.0):.4f}  f1={mets.get('f1', 0.0):.4f}  mcc={mets.get('mcc', 0.0):.4f}\n"
            )
            f.write(f"TP={mets.get('tp', 0)} FP={mets.get('fp', 0)} FN={mets.get('fn', 0)} TN={mets.get('tn', 0)}\n")
            if rule.pairwise_fp_sources:
                f.write(f"pairwise_fp_sources: {json.dumps(rule.pairwise_fp_sources, sort_keys=True)}\n")
            if rule.notes:
                f.write(f"notes: {' | '.join(rule.notes)}\n")
            f.write("markers:\n")
            n_in = dataset["family_sizes"].get(family, 0)
            n_out = dataset["total_samples"] - n_in
            for cid in rule.markers:
                c_in = dataset["family_counts"][family].get(cid, 0)
                c_global = dataset["global_counts"].get(cid, 0)
                c_out = max(0, c_global - c_in)
                p_in = prevalence(c_in, n_in)
                p_out = prevalence(c_out, n_out)
                f.write(f"  - {cid} (p_in={p_in:.3f}, p_out={p_out:.3f})\n")
            if rule.forbids:
                f.write("forbidden:\n")
                for cid in rule.forbids:
                    nm = cluster_id_to_forbidden_name(cid, cluster_to_raw, cluster_to_evs, forbid_hyp_as_accession=args.forbid_hyp_as_accession)
                    f.write(f"  - {cid}")
                    if nm:
                        f.write(f" => {nm}")
                    f.write("\n")
            f.write("\n")

    confusion = compute_confusion_matrix(rules, dataset)
    with open(out_conf, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["target_family", "source_family", "count"])
        for target in sorted(confusion):
            for source, count in sorted(confusion[target].items(), key=lambda kv: (-kv[1], kv[0])):
                w.writerow([target, source, count])

    with open(out_backmap, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["family", "cluster_id", "count_in_family", "n_in", "p_in", "count_outside_family", "n_out", "p_out", "raw_names", "evidence_subjects"])
        for family in dataset["families"]:
            n_in = dataset["family_sizes"].get(family, 0)
            n_out = dataset["total_samples"] - n_in
            for cid in sorted(dataset["family_union"].get(family, set())):
                c_in = dataset["family_counts"][family].get(cid, 0)
                c_global = dataset["global_counts"].get(cid, 0)
                c_out = max(0, c_global - c_in)
                p_in = prevalence(c_in, n_in)
                p_out = prevalence(c_out, n_out)
                raw_join = " | ".join(sorted(cluster_to_raw.get(cid, set())))
                ev_join = " | ".join(sorted(cluster_to_evs.get(cid, set())))
                w.writerow([family, cid, c_in, n_in, f"{p_in:.4f}", c_out, n_out, f"{p_out:.4f}", raw_join, ev_join])

    print("Wrote:")
    print(" -", out_json)
    print(" -", out_best)
    print(" -", out_readable)
    print(" -", out_conf)
    print(" -", out_backmap)


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------
def parse_family_list(raw: str) -> Set[str]:
    """Parse a comma-delimited family list argument."""
    if not raw:
        return set()
    return {x.strip() for x in raw.split(",") if x.strip()}


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface for rule optimization."""
    ap = argparse.ArgumentParser(description="Fast in-memory SaPhARI rule optimizer with beam search and forbidden-first refinement")
    ap.add_argument("--annotated-root", default="annotated")
    ap.add_argument("--out-dir", default="cluster_params/output_beam")
    ap.add_argument("--emit-families-json", default="Families_beam.json")
    ap.add_argument("--default-length", type=int, default=15000)

    # Pools / signatures
    ap.add_argument("--k-max", type=int, default=30)
    ap.add_argument("--m-min", type=int, default=3)
    ap.add_argument("--m-max", type=int, default=12)
    ap.add_argument("--positive-pool-cap", type=int, default=40)
    ap.add_argument("--forbid-pool-cap", type=int, default=40)
    ap.add_argument("--seed-forbid-size", type=int, default=2)
    ap.add_argument("--min-marker-support", type=float, default=0.10)
    ap.add_argument("--max-marker-out", type=float, default=0.85)
    ap.add_argument("--marker-stop-score", type=float, default=-0.20)

    # Positive pool weights
    ap.add_argument("--cover-true-weight", type=float, default=1.8)
    ap.add_argument("--support-in-weight", type=float, default=1.6)
    ap.add_argument("--support-out-weight", type=float, default=1.5)
    ap.add_argument("--rival-prev-weight", type=float, default=2.4)
    ap.add_argument("--cover-rival-weight", type=float, default=1.3)
    ap.add_argument("--cover-other-weight", type=float, default=0.7)

    # Singleton priors
    ap.add_argument("--singleton-disc-weight", type=float, default=3.0)
    ap.add_argument("--singleton-rival-weight", type=float, default=2.4)
    ap.add_argument("--singleton-true-weight", type=float, default=1.2)
    ap.add_argument("--singleton-other-weight", type=float, default=0.8)
    ap.add_argument("--singleton-target-recall", type=float, default=0.70)

    # Objective
    ap.add_argument("--target-recall", type=float, default=0.75)
    ap.add_argument("--obj-recall-weight", type=float, default=4.5)
    ap.add_argument("--obj-precision-weight", type=float, default=2.4)
    ap.add_argument("--obj-specificity-weight", type=float, default=2.2)
    ap.add_argument("--obj-f1-weight", type=float, default=1.0)
    ap.add_argument("--obj-m-abs-weight", type=float, default=0.08)
    ap.add_argument("--obj-m-ratio-weight", type=float, default=0.8)
    ap.add_argument("--obj-k-weight", type=float, default=0.15)
    ap.add_argument("--obj-forbid-count-weight", type=float, default=0.12)
    ap.add_argument("--obj-fp-weight", type=float, default=3.1)
    ap.add_argument("--obj-fn-weight", type=float, default=4.6)
    ap.add_argument("--obj-pairwise-weight", type=float, default=2.8)
    ap.add_argument("--obj-loose-shape-weight", type=float, default=1.5)
    ap.add_argument("--obj-strict-shape-weight", type=float, default=1.4)
    ap.add_argument("--obj-recall-shortfall-weight", type=float, default=10.0)
    ap.add_argument("--support-quantile-low", type=float, default=0.15)
    ap.add_argument("--loose-ratio-target", type=float, default=0.28)
    ap.add_argument("--default-pairwise-weight", type=float, default=1.0)
    ap.add_argument("--pair-weights", default="CFPICI->PICI=5.0,PICI->CFPICI=5.0")

    # Rivals / passes / freezing
    ap.add_argument("--initial-rival-top-n", type=int, default=3)
    ap.add_argument("--rival-top-n", type=int, default=3)
    ap.add_argument("--global-passes", type=int, default=4)
    ap.add_argument("--freeze-precision", type=float, default=0.95)
    ap.add_argument("--freeze-recall", type=float, default=0.95)
    ap.add_argument("--always-optimize-families", default="CFPICI,PICI,VEIME_2,VEIME_3,VEIME_4,VEIME_6,VEIME_8,VEIME_9,VEIME_11,VEIME_M1")
    ap.add_argument("--never-freeze-families", default="CFPICI,PICI,VEIME_11,VEIME_M1")

    # Beam search
    ap.add_argument("--beam-width", type=int, default=10)
    ap.add_argument("--hard-beam-iters", type=int, default=14)
    ap.add_argument("--easy-beam-iters", type=int, default=6)
    ap.add_argument("--hard-neighbor-add", type=int, default=8)
    ap.add_argument("--easy-neighbor-add", type=int, default=4)
    ap.add_argument("--hard-forbid-add", type=int, default=8)
    ap.add_argument("--easy-forbid-add", type=int, default=4)
    ap.add_argument("--remove-marker-cap", type=int, default=6)
    ap.add_argument("--remove-forbid-cap", type=int, default=4)
    ap.add_argument("--swap-marker-cap", type=int, default=4)
    ap.add_argument("--swap-forbid-cap", type=int, default=2)
    ap.add_argument("--swap-add-cap", type=int, default=3)

    # Forbids
    ap.add_argument("--max-forbids", type=int, default=10)
    ap.add_argument("--final-forbid-pool-cap", type=int, default=25)
    ap.add_argument("--forbid-fp-prev-min", type=float, default=0.25)
    ap.add_argument("--forbid-target-prev-max", type=float, default=0.20)
    ap.add_argument("--forbid-rival-prev-min", type=float, default=0.35)
    ap.add_argument("--forbid-out-prev-min", type=float, default=0.35)
    ap.add_argument("--forbid-fp-prev-weight", type=float, default=3.8)
    ap.add_argument("--forbid-rival-weight", type=float, default=2.6)
    ap.add_argument("--forbid-mean-rival-weight", type=float, default=1.1)
    ap.add_argument("--forbid-out-weight", type=float, default=0.9)
    ap.add_argument("--forbid-target-weight", type=float, default=3.2)
    ap.add_argument("--forbid-include-hyp", action="store_true")
    ap.add_argument("--forbid-hyp-as-accession", action="store_true")
    return ap


def main() -> None:
    """Entry point for command-line execution."""
    args = build_parser().parse_args()
    if args.m_min < 3:
        raise SystemExit("--m-min must be at least 3")
    if args.m_max < args.m_min:
        raise SystemExit("--m-max must be >= --m-min")
    if args.k_max < args.m_min:
        raise SystemExit("--k-max must be >= --m-min")
    args.always_optimize_families = parse_family_list(args.always_optimize_families)
    args.never_freeze_families = parse_family_list(args.never_freeze_families)

    dataset = load_dataset(os.path.abspath(args.annotated_root))
    rules = optimize_rules(dataset, args)
    write_outputs(rules, dataset, args)


if __name__ == "__main__":
    main()
