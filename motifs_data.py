"""
Catalogue codé en dur de motifs problématiques.
Chaque entrée possède :
  name, pattern, description, score, heuristic (bool), category
"""

# On tente d’utiliser le moteur « regex » (optionnel). Aucune tolérance
# aux mismatches n’est utilisée ; si 'regex' n’est pas présent, on bascule
# sur le module standard 're'.
try:
    import regex as re          # pip install regex
except ImportError:
    import re

MOTIFS_BY_CATEGORY = {
    # ------------------------------------------------------------------
    "replication": [
        {
            "name": "DnaA box (E. coli)",
            "pattern": "TTATCCACA",
            "description": "Site de liaison DnaA – initiation de la réplication.",
            "score": 5.0,
            "heuristic": False,
        },
        {
            "name": "Iteron 9-mer",
            "pattern": r"ATGCAAAT[AT]",
            "description": "Séquence répétée régulant la réplication plasmidique.",
            "score": 4.0,
            "heuristic": False,
        },
        {
            "name": "Tus/Ter site",
            "pattern": r"AGTATGTTGTAA[ACGT]{7}C",
            "description": "Site Ter se liant à Tus – blocage de la fourche.",
            "score": 3.5,
            "heuristic": False,
        },
        {
            "name": "Telomeric repeat (TTAGGG)n",
            "pattern": r"(TTAGGG){4,}",
            "description": "Répétitions télomériques formant G4 aux extrémités.",
            "score": 3.0,
            "heuristic": False,
        },
        # ------------------------------------------------------------------
        # ITR des plasmides AAV : couvre les orientations 3′ ET 5′
        {
            "name": "AAV ITR (3′ ou 5′)",
            "pattern": (
                r"(?:"
                # 3′ ITR : cœur 130 nt + queue variable 0-15 nt
                r"AGGAACCCCTAGTGATGGAGTTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCC"
                r"GGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAGCGAGCGAGCGCGCAG"
                r"[ACGT]{0,15}"
                r"|"
                # 5′ ITR : tête variable 0-26 nt + cœur reverse-complement
                r"[ACGT]{0,26}"
                r"CTGCGCGCTCGCTCGCTCACTGAGGCCGCCCGGGCAAAGCCCGGGCGTCGGGCGACCTTTGGTC"
                r"GCCCGGCCTCAGTGAGCGAGCGAGCGCGCAGAGAGGGAGTGGCCAACTCCATCACTAGGGGTTCCT"
                r")"
            ),
            "description": (
                "Inverted terminal repeat d’AAV ; couvre les orientations 3′ et 5′ "
                "sur un seul brin (queues variables : 0-15 nt ou 0-26 nt)."
            ),
            "score": 4.8,
            "heuristic": False,
        },
        {
            "name": "Microsatellite (CAG)n",
            "pattern": r"(CAG){6,}",
            "description": "Triplet répétitif propice au glissement de polymérase.",
            "score": 3.0,
            "heuristic": False,
        },
        {
            "name": "ORI thermosensible (mutation G→A)",
            "pattern": "CGGCTACACTAGAAGAACAGTATTTGGT",
            "description": (
                "Mutation ponctuelle dans l'origine de réplication : G→A (thermosensible). "
                "Augmente la réplication car le complexe initateur se lie moins fortement."
            ),
            "score": 4.3,
            "heuristic": False,
        },
    ],
    # ------------------------------------------------------------------
    "structure": [
        {
            "name": "G-quadruplex motif",
            "pattern": r"G{3}[ACGT]{1,7}G{3}[ACGT]{1,7}G{3}[ACGT]{1,7}G{3}",
            "description": "Séquence G-riche formant un quadruplex de guanine.",
            "score": 4.5,
            "heuristic": False,
        },
        {
            "name": "i-motif motif",
            "pattern": r"C{3}[ACGT]{1,7}C{3}[ACGT]{1,7}C{3}[ACGT]{1,7}C{3}",
            "description": "Séquence C-riche formant un i-motif.",
            "score": 3.5,
            "heuristic": False,
        },
        {
            "name": "Z-DNA (CG)n",
            "pattern": r"(CG){6,}",
            "description": "Répétition purine/pyrimidine adoptant la conformation Z.",
            "score": 3.0,
            "heuristic": False,
        },
        {
            "name": "Triplex (H-DNA) purine",
            "pattern": r"(GA){10,}",
            "description": (
                "Répétition miroir homopurine formant un triplex intramoléculaire."
            ),
            "score": 3.0,
            "heuristic": False,
        },
    ],
    # ------------------------------------------------------------------
    "transcription": [
        # ------------------------------------------------------------------
        # Retroviral LTR – type A (HIV-1 dérivé)
        {
            "name": "Retroviral LTR (type A, 3′ & 5′)",
            "pattern": (
                r"[ACGT]{0,60}"
                r"GGKYYYYYYKGKTWRRMCMRRWYYKRRSCYKGGRRSYYYYYKGSYWAMYWRGGRAMCCMM"
                r"YKSYTWARSCYYMAWWAARSYTKSCYTKRRKKSYTYMARKWRKKKKKKSCCSKYYKKTKKKKKRMY"
                r"YYKGKWAMYWRRRRWYCCYYMRRMCCYTTTWRKYMRKKKKGRAAAWYYYYWRSMR"
                r"[ACGT]{0,5}"
            ),
            "description": (
                "LTR rétroviral commun aux plasmides testés (famille HIV-1). "
                "Couvre les 3′-LTR complètes et les 5′-LTR tronquées sur un seul brin "
                "(zones variables : tête 0-60 nt, queue 0-5 nt)."
            ),
            "score": 4.8,
            "heuristic": False,
        },
        # ------------------------------------------------------------------
        {
            "name": "Terminateur intrinsèque",
            "pattern": r"GC{2,3}[ACGT]{1,5}GC{2,3}T{4,}",
            "description": "Hairpin GC-riche + queue poly-T – décroche l’ARN Pol.",
            "score": 4.0,
            "heuristic": False,
        },
        {
            "name": "Pause GC-riche",
            "pattern": r"[GC]{20,}",
            "description": "Segment très GC induisant une pause de RNAP.",
            "score": 2.5,
            "heuristic": False,
        },
    ],
    # ------------------------------------------------------------------
    "heuristic": [
        {
            "name": "Région haute GC",
            "pattern": None,
            "description": "Fenêtre de 100 pb avec > 80 % GC.",
            "score": 2.5,
            # window-based heuristic parameters (optional):
            "window_size": 100,
            # coverage fraction at which the heuristic gives full contribution
            "threshold_fraction": 0.05,
            "heuristic": True,
        },
        {
            "name": "Région haute AT",
            "pattern": None,
            "description": "Fenêtre de 100 pb avec > 75 % AT.",
            "score": 2.0,
            "window_size": 100,
            "threshold_fraction": 0.05,
            "heuristic": True,
        },
        {
            "name": "Homopolymère long",
            "pattern": None,
            "description": "≥ 8 nucléotides identiques consécutifs.",
            "score": 2.0,
            "heuristic": True,
        },
        {
            "name": "Tandem repeat long",
            "pattern": None,
            "description": "Unité 2-6 pb répétée > 6 fois.",
            "score": 2.5,
            "heuristic": True,
        },
        {
            "name": "Palindrome long",
            "pattern": None,
            "description": "Répétition inversée totale > 20 pb.",
            "score": 3.0,
            "heuristic": True,
        },
        {
            "name": "Skew GC fort",
            "pattern": None,
            "description": "Skew G/C ≥ 70 % sur 100 pb.",
            "score": 2.5,
            "window_size": 100,
            "threshold_fraction": 0.05,
            "heuristic": True,
        },
    ],
}

# ----------------------------------------------------------------------
# Liste aplatie exportée
MOTIFS: list[dict] = []
for cat, lst in MOTIFS_BY_CATEGORY.items():
    for m in lst:
        m["category"] = cat
        MOTIFS.append(m)

# ----------------------------------------------------------------------
# Synergy rules: data-driven list of rules applied after motif collection.
# Each rule defines `requires` (list of motif names), the `category` affected,
# and either a `multiplier` or an `add` value. `flag` is the identifier returned
# in analysis risk.synergy_flags for transparency.
SYNERGY_RULES = [
    {
        "flag": "AAV_ITR+ORI_ts",
        "requires": ["AAV ITR (3′ ou 5′)", "ORI thermosensible (mutation G→A)"],
        "category": "replication",
        "multiplier": 1.25,
        "note": "AAV ITR and thermosensitive ORI co-occurred; replication risk amplified.",
    }
]

# Default category weights (can be overridden via config).
# Default calibrated weights (can be overridden via config).
CATEGORY_WEIGHTS = {
    "replication": 0.092389,
    "structure": 0.044266,
    "transcription": 0.000866,
    "heuristic": 0.052381,
}

# Optional per-motif weights learned from data (empty by default).
MOTIF_WEIGHTS = {}


# Optional external configuration loader. If a file `config/motifs_config.json`
# exists it can override per-motif parameters (e.g. `threshold_fraction`,
# `window_size`) and replace or extend `SYNERGY_RULES` without editing this file.
from pathlib import Path
import json

_CFG_DIR = Path(__file__).parent / "config"
_cfg_mtime = None


def _resolve_cfg_path() -> Path | None:
    json_path = _CFG_DIR / "motifs_config.json"
    if json_path.exists():
        return json_path
    return None


def _load_config(path: Path) -> None:
    global SYNERGY_RULES
    try:
        cfg = json.loads(path.read_text()) or {}
        # Apply motif overrides
        motifs_cfg = cfg.get("motifs", {})
        if motifs_cfg:
            # motifs_cfg expected: mapping motif_name -> {key: value}
            name_to_motif = {m["name"]: m for m in MOTIFS}
            for nm, overrides in motifs_cfg.items():
                if nm in name_to_motif:
                    name_to_motif[nm].update(overrides)
        # Replace or extend synergy rules
        if "synergy_rules" in cfg:
            if isinstance(cfg["synergy_rules"], list):
                SYNERGY_RULES = cfg["synergy_rules"]
        # Optional category weights for risk aggregation
        cat_weights = cfg.get("category_weights")
        if isinstance(cat_weights, dict):
            for k, v in cat_weights.items():
                try:
                    CATEGORY_WEIGHTS[k] = float(v)
                except Exception:
                    continue
        # Optional per-motif weights
        motif_weights = cfg.get("motif_weights")
        if isinstance(motif_weights, dict):
            for k, v in motif_weights.items():
                try:
                    MOTIF_WEIGHTS[k] = float(v)
                except Exception:
                    continue
    except Exception:
        # If parsing fails, ignore and keep built-in defaults
        pass


def reload_config_if_changed() -> None:
    """Reload config if file exists and changed on disk."""
    global _cfg_mtime
    path = _resolve_cfg_path()
    if not path:
        return
    try:
        mtime = path.stat().st_mtime
    except Exception:
        return
    if _cfg_mtime is None or mtime != _cfg_mtime:
        _cfg_mtime = mtime
        _load_config(path)


# Expose current config path for debug.
def get_active_config_path() -> str | None:
    path = _resolve_cfg_path()
    return str(path) if path else None


# initial load
_path = _resolve_cfg_path()
if _path:
    _load_config(_path)
