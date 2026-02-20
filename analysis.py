"""
Détection de motifs dans une séquence ADN.
• Nettoyage robuste (suppression de tout caractère non-ACGT[N])
• Détection regex (motifs exacts)
• Détection heuristique (contenu GC/AT, homopolymères, etc.)
"""

import re
from typing import List, Dict
from collections import defaultdict

from motifs_data import (
    MOTIFS,
    SYNERGY_RULES,
    CATEGORY_WEIGHTS,
    MOTIF_WEIGHTS,
    reload_config_if_changed,
)
import math

# ----- Constantes -----
WINDOW = 100
IUPAC = {
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "N": "[ACGT]",
}

# Hardcoded motif weights (learned from 2026-01-20 dataset).
_HARD_CODED_MOTIF_WEIGHTS_RAW = {
    "AAV ITR (3ƒ?ý ou 5ƒ?ý)": 0.891953,
    "DnaA box (E. coli)": 0.142225,
    "G-quadruplex motif": 0.068802,
    "HomopolymÇùre long": -0.489931,
    "Iteron 9-mer": 0.703687,
    "ORI thermosensible (mutation GƒÅ'A)": 0.327131,
    "Palindrome long": -0.062260,
    "Pause GC-riche": 0.401954,
    "Retroviral LTR (type A, 3ƒ?ý & 5ƒ?ý)": 0.624350,
    "RÇ¸gion haute AT": -2.270334,
    "RÇ¸gion haute GC": -0.721605,
    "Skew GC fort": 1.458753,
    "Tandem repeat long": 1.081868,
    "Terminateur intrinsÇùque": 0.049241,
    "i-motif motif": -0.255025,
}


def _norm_key(name: str) -> str:
    # Normalize names to avoid encoding mismatch.
    return re.sub(r"[^A-Za-z0-9]+", "", name).lower()


HARD_CODED_MOTIF_WEIGHTS = {
    _norm_key(k): v for k, v in _HARD_CODED_MOTIF_WEIGHTS_RAW.items()
}

# ------------------------------------------------------------------
# Nettoyage robuste
def _clean_sequence(raw: str, allow_n: bool = True) -> str:
    """Supprime entêtes FASTA, espaces/retours & caractères non-ACGT(N)."""
    seq_lines = [l for l in raw.splitlines() if not l.startswith(">")]
    seq = "".join(seq_lines).upper()
    alphabet = "ACGTN" if allow_n else "ACGT"
    return re.sub(f"[^{alphabet}]", "", seq)

# ------------------------------------------------------------------
# Fonctions utilitaires
def _iupac_to_regex(pat: str) -> str:
    regex = pat.upper()
    for code, repl in IUPAC.items():
        regex = regex.replace(code, repl)
    return regex

def _windows(seq: str, size: int):
    for i in range(0, len(seq) - size + 1):
        yield i, seq[i : i + size]

def _gc_content(fragment: str) -> float:
    gc = fragment.count("G") + fragment.count("C")
    return gc / len(fragment)

def _at_content(fragment: str) -> float:
    at = fragment.count("A") + fragment.count("T")
    return at / len(fragment)

def _reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]

# ------------------------------------------------------------------
# Détections heuristiques
def _detect_high_gc(seq: str) -> List[int]:
    return [i + 1 for i, frag in _windows(seq, WINDOW) if _gc_content(frag) >= 0.80]

def _detect_high_at(seq: str) -> List[int]:
    return [i + 1 for i, frag in _windows(seq, WINDOW) if _at_content(frag) >= 0.75]

def _detect_homopolymer(seq: str) -> List[int]:
    return [m.start() + 1 for m in re.finditer(r"(A{8,}|C{8,}|G{8,}|T{8,})", seq)]

def _detect_tandem(seq: str) -> List[int]:
    return [m.start() + 1 for m in re.finditer(r"([ACGT]{2,6})(\1){5,}", seq)]

def _detect_palindrome(seq: str) -> List[int]:
    pos = []
    for length in range(10, 31):         # bras = 10–30 pb
        total = length * 2
        for i in range(len(seq) - total + 1):
            left = seq[i : i + length]
            right = seq[i + length : i + total]
            if left == _reverse_complement(right):
                pos.append(i + 1)
    return pos

def _detect_gc_skew(seq: str) -> List[int]:
    positions = []
    for i, frag in _windows(seq, WINDOW):
        g = frag.count("G")
        c = frag.count("C")
        if g + c == 0:
            continue
        ratio = g / (g + c)
        if ratio >= 0.70 or ratio <= 0.30:
            positions.append(i + 1)
    return positions

HEURISTIC_FUNCS = {
    "Région haute GC": _detect_high_gc,
    "Région haute AT": _detect_high_at,
    "Homopolymère long": _detect_homopolymer,
    "Tandem repeat long": _detect_tandem,
    "Palindrome long": _detect_palindrome,
    "Skew GC fort": _detect_gc_skew,
}

# ------------------------------------------------------------------
# Utilities for merging window-based hits into contiguous basepair segments.
def _merge_window_hits_to_segments(window_starts: List[int], window_size: int):
    """Merge consecutive window starts (1-based) into segments.

    Rule: consecutive window starts (i, i+1, i+2...) form a run.
    Segment start_bp = first_window_start
    Segment end_bp = last_window_start + window_size - 1
    Returns list of dicts: {start, end, length_bp}
    """
    if not window_starts:
        return []
    sorted_starts = sorted(window_starts)
    segments = []
    run_start = sorted_starts[0]
    prev = run_start
    for s in sorted_starts[1:]:
        if s == prev + 1:
            prev = s
            continue
        # close run
        seg = {"start": run_start, "end": prev + window_size - 1}
        seg["length_bp"] = seg["end"] - seg["start"] + 1
        segments.append(seg)
        run_start = s
        prev = s
    # final run
    seg = {"start": run_start, "end": prev + window_size - 1}
    seg["length_bp"] = seg["end"] - seg["start"] + 1
    segments.append(seg)
    return segments


def _segment_coverage(segments: List[dict], seq_len: int):
    """Compute total_bp_covered and coverage_fraction for segments.

    segments: list of {start,end,length_bp,...}
    """
    if not segments:
        return {"total_bp_covered": 0, "coverage_fraction": 0.0}
    # Merge overlapping segments conservatively (though segments from windows shouldn't overlap)
    ranges = sorted([(s["start"], s["end"]) for s in segments])
    merged = []
    cur_s, cur_e = ranges[0]
    for s, e in ranges[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    total = sum(e - s + 1 for s, e in merged)
    return {"total_bp_covered": total, "coverage_fraction": total / max(1, seq_len)}


# Default threshold used when motif-specific value is missing
DEFAULT_HEURISTIC_THRESHOLD = 0.05

def _scan_window_metric(seq: str, metric: str):
    """Scan windows producing list of (start_1based, value).

    metric: 'gc', 'at', 'skew_dev'
    """
    out = []
    for i, frag in _windows(seq, WINDOW):
        if metric == "gc":
            val = _gc_content(frag)
            if val >= 0.80:
                out.append((i + 1, val))
        elif metric == "at":
            val = _at_content(frag)
            if val >= 0.75:
                out.append((i + 1, val))
        elif metric == "skew_dev":
            g = frag.count("G")
            c = frag.count("C")
            if g + c == 0:
                continue
            ratio = g / (g + c)
            deviation = abs(ratio - 0.5)
            if ratio >= 0.70 or ratio <= 0.30:
                out.append((i + 1, deviation))
    return out


def analyze_sequence_v2(raw_sequence: str) -> Dict:
    """Analyse renvoyant le nouveau schéma (v2) décrit dans la tâche.
    - Non-heuristic motifs keep `positions` list (unchanged)
    - Window-based heuristics return merged `segments` and coverage summary
    - Risk scoring by category and total, with simple synergy rules

    This function preserves the exact regex matching behaviour for non-heuristic motifs.
    """
    reload_config_if_changed()
    seq = _clean_sequence(raw_sequence)
    seq_len = len(seq)

    motifs_out = []

    # 1) regex motifs (non-heuristic)
    for motif in [m for m in MOTIFS if not m["heuristic"]]:
        regex = _iupac_to_regex(motif["pattern"]) if motif.get("pattern") else None
        positions = []
        if regex:
            try:
                compiled = re.compile(regex)
            except re.error:
                compiled = None
            if compiled:
                positions = [m.start() + 1 for m in compiled.finditer(seq)]
        if positions:
            motifs_out.append(
                dict(
                    name=motif["name"],
                    category=motif["category"],
                    score=motif["score"],
                    impact=motif.get("description", ""),
                    heuristic=False,
                    positions=positions,
                )
            )

    # 2) heuristics: handle window-based heuristics specially to produce segments
    # window-based metrics mapping
    window_based = {
        "Région haute GC": "gc",
        "Région haute AT": "at",
        "Skew GC fort": "skew_dev",
    }

    # collect heuristic results
    for motif in [m for m in MOTIFS if m["heuristic"]]:
        name = motif["name"]
        if name in window_based:
            metric = window_based[name]
            # allow per-motif window size if present
            window_size = motif.get("window_size", WINDOW)
            # if not default window, temporarily override WINDOW for scanning
            if window_size != WINDOW:
                # create local scanning using the provided window size
                hits = []
                for i in range(0, len(seq) - window_size + 1):
                    frag = seq[i : i + window_size]
                    if metric == "gc":
                        val = _gc_content(frag)
                        if val >= 0.80:
                            hits.append((i + 1, val))
                    elif metric == "at":
                        val = _at_content(frag)
                        if val >= 0.75:
                            hits.append((i + 1, val))
                    elif metric == "skew_dev":
                        g = frag.count("G")
                        c = frag.count("C")
                        if g + c == 0:
                            continue
                        ratio = g / (g + c)
                        deviation = abs(ratio - 0.5)
                        if ratio >= 0.70 or ratio <= 0.30:
                            hits.append((i + 1, deviation))
            else:
                hits = _scan_window_metric(seq, metric)
            # hits: list of (start, value)
            starts = [s for s, v in hits]
            segments = _merge_window_hits_to_segments(starts, WINDOW)
            # attach max_value per segment if available
            if hits and segments:
                # build a map of start->value for quick lookup
                val_map = {s: v for s, v in hits}
                for seg in segments:
                    seg_starts = list(range(seg["start"], seg["end"] - WINDOW + 2))
                    # compute max value among contributing window starts
                    max_val = 0.0
                    for s in seg_starts:
                        v = val_map.get(s)
                        if v is not None and v > max_val:
                            max_val = v
                    if max_val:
                        seg["max_value"] = max_val
            coverage = _segment_coverage(segments, seq_len)
            # use motif-specific threshold if present
            threshold = motif.get("threshold_fraction", motif.get("threshold", motif.get("threshold_fraction", None)))
            if threshold is None:
                threshold = motif.get("threshold_fraction", motif.get("threshold", DEFAULT_HEURISTIC_THRESHOLD))
            # store the threshold used for traceability
            motifs_out.append(
                dict(
                    name=name,
                    category=motif["category"],
                    score=motif["score"],
                    impact=motif.get("description", ""),
                    heuristic=True,
                    segments=segments,
                    count_segments=len(segments),
                    total_bp_covered=coverage["total_bp_covered"],
                    coverage_fraction=coverage["coverage_fraction"],
                    threshold_fraction=threshold,
                )
            )
        else:
            # non-window heuristic: reuse existing funcs for exact positions
            func = HEURISTIC_FUNCS.get(name)
            positions = func(seq) if func else []
            if positions:
                motifs_out.append(
                    dict(
                        name=name,
                        category=motif["category"],
                        score=motif["score"],
                        impact=motif.get("description", ""),
                        heuristic=True,
                        positions=positions,
                    )
                )

    # Risk scoring
    risk_by_category = defaultdict(float)

    # contributions per motif to allow sorting heuristics by contribution
    motif_contributions = {}
    base_contributions = {}

    for m in motifs_out:
        base = m.get("score", 0.0)
        cat = m.get("category", "heuristic")
        if m.get("heuristic"):
            # heuristics: weight by coverage fraction for window-based, else by presence
            cov = m.get("coverage_fraction")
            if cov is not None:
                thresh = m.get("threshold_fraction", DEFAULT_HEURISTIC_THRESHOLD)
                contrib = base * min(1.0, cov / thresh)
            else:
                # non-window heuristic -> presence counts once
                contrib = base
        else:
            # non-heuristic: presence counts once (prevents "more hits = worse" bias)
            contrib = base
        base_contrib = contrib
        base_contributions[m["name"]] = base_contrib
        motif_weight = MOTIF_WEIGHTS.get(m.get("name"), 1.0)
        contrib *= motif_weight
        risk_by_category[cat] += contrib
        motif_contributions[m["name"]] = contrib

    # Synergy rules (data-driven via SYNERGY_RULES from motifs_data)
    synergy_flags = []
    notes = []
    names_present = {m["name"] for m in motifs_out}
    for rule in SYNERGY_RULES:
        required = set(rule.get("requires", []))
        if required.issubset(names_present):
            cat = rule.get("category")
            if cat and cat in risk_by_category:
                if "multiplier" in rule:
                    before = risk_by_category.get(cat, 0.0)
                    risk_by_category[cat] = before * rule["multiplier"]
                elif "add" in rule:
                    risk_by_category[cat] = risk_by_category.get(cat, 0.0) + rule["add"]
            synergy_flags.append(rule.get("flag"))
            if rule.get("note"):
                notes.append(rule.get("note"))

    weighted_by_category = {
        cat: score * CATEGORY_WEIGHTS.get(cat, 1.0)
        for cat, score in risk_by_category.items()
    }
    total_risk = sum(weighted_by_category.values())

    # Hardcoded risk (independent from config/venv)
    risk_new_by_motif = {}
    for name, base_contrib in base_contributions.items():
        weight = HARD_CODED_MOTIF_WEIGHTS.get(_norm_key(name), 0.0)
        if base_contrib == 0.0 and weight == 0.0:
            continue
        risk_new_by_motif[name] = base_contrib * weight
    risk_new_total = sum(risk_new_by_motif.values())

    # Sort motifs: non-heuristic by score desc, heuristics by contribution desc
    def sort_key(m):
        if not m.get("heuristic"):
            return (0, -m.get("score", 0.0), m.get("positions", [0])[0])
        return (1, -motif_contributions.get(m.get("name"), 0.0))

    motifs_out.sort(key=sort_key)

    return {
        "sequence_length": seq_len,
        "motifs": motifs_out,
        "risk": {
            "total": total_risk,
            "by_category": dict(risk_by_category),
            "by_category_weighted": weighted_by_category,
            "weights": dict(CATEGORY_WEIGHTS),
            "motif_weights": dict(MOTIF_WEIGHTS),
            "notes": notes,
            "synergy_flags": synergy_flags,
        },
        "risk_new": {
            "version": "v2-hardcoded",
            "total": risk_new_total,
            "by_motif": risk_new_by_motif,
        },
    }

# ------------------------------------------------------------------
def analyze_sequence(raw_sequence: str) -> List[Dict]:
    """Analyse et retourne une liste de motifs trouvés."""
    seq = _clean_sequence(raw_sequence)

    results = []

    # 1) motifs regex
    for motif in [m for m in MOTIFS if not m["heuristic"]]:
        regex = _iupac_to_regex(motif["pattern"])
        try:
            compiled = re.compile(regex)
        except re.error:
            continue
        positions = [m.start() + 1 for m in compiled.finditer(seq)]
        if positions:
            results.append(
                dict(
                    name=motif["name"],
                    positions=positions,
                    score=motif["score"],
                    impact=motif["description"],
                    category=motif["category"],
                    heuristic=False,
                )
            )

    # 2) heuristiques
    for motif in [m for m in MOTIFS if m["heuristic"]]:
        func = HEURISTIC_FUNCS[motif["name"]]
        positions = func(seq)
        if positions:
            results.append(
                dict(
                    name=motif["name"],
                    positions=positions,
                    score=motif["score"],
                    impact=motif["description"],
                    category=motif["category"],
                    heuristic=True,
                )
            )

    # Tri par score DESC puis première position ASC
    results.sort(key=lambda r: (-r["score"], r["positions"][0]))
    return results
