from dataclasses import dataclass
from bisect import bisect_left
from io import StringIO
import json
import os
from pathlib import Path
import random
from threading import Timer
from typing import Dict, List
import webbrowser

from Bio import SeqIO, motifs
from Bio.Seq import Seq
from flask import Flask, render_template, request

from analysis import analyze_sequence_v2

APP_NAME = "Promoter + DNA Motifs Dashboard"
BASE_DIR = Path(__file__).resolve().parent
PRESET_OPTIONS = ("custom", "strict", "balanced", "sensitive")
PRESET_VALUES = {
    "strict": {
        "min_spacing": "15",
        "max_spacing": "19",
        "min_norm_35": "0.60",
        "min_norm_10": "0.60",
        "min_combined": "0.70",
        "max_hits": "40",
        "circular_sequence": "on",
    },
    "balanced": {
        "min_spacing": "14",
        "max_spacing": "20",
        "min_norm_35": "0.42",
        "min_norm_10": "0.42",
        "min_combined": "0.52",
        "max_hits": "120",
        "circular_sequence": "on",
    },
    "sensitive": {
        "min_spacing": "12",
        "max_spacing": "23",
        "min_norm_35": "0.20",
        "min_norm_10": "0.20",
        "min_combined": "0.32",
        "max_hits": "350",
        "circular_sequence": "on",
    },
}
PRESET_SIGNIFICANCE_ALPHA = {"strict": 0.02, "balanced": 0.04, "sensitive": 0.12}
PRESET_BACKGROUND_SHUFFLES = {"strict": 2, "balanced": 2, "sensitive": 1}
PRESET_OVERLAP_DISTANCE = {"strict": 4, "balanced": 3, "sensitive": 2}

app = Flask(__name__)
app.config["SECRET_KEY"] = "dev"


@dataclass
class PromoterHit:
    strand: str
    box35_start: int
    box10_start: int
    spacing: int
    raw35: float
    raw10: float
    norm35: float
    norm10: float
    spacing_bonus: float
    combined: float
    promoter_seq: str
    p_value: float = 1.0


def load_jaspar(jaspar_path: Path):
    with open(jaspar_path, encoding="utf-8") as handle:
        record = motifs.parse(handle, "jaspar")
        return record[0]


def normalize_score(value: float, score_min: float, score_max: float) -> float:
    if score_max <= score_min:
        return 0.0
    normalized = (value - score_min) / (score_max - score_min)
    return max(0.0, min(1.0, float(normalized)))


def parse_float(name: str, value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid number for '{name}'.") from exc


def parse_int(name: str, value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid integer for '{name}'.") from exc


def parse_sequence_from_input(fasta_file, pasted_sequence: str) -> str:
    if fasta_file and fasta_file.filename:
        fasta_content = fasta_file.read().decode("utf-8", errors="ignore")
        records = list(SeqIO.parse(StringIO(fasta_content), "fasta"))
        if not records:
            raise ValueError("The uploaded FASTA file does not contain any sequence record.")
        sequence = "".join(str(record.seq) for record in records)
    elif pasted_sequence.strip():
        raw = pasted_sequence.strip()
        if raw.startswith(">"):
            records = list(SeqIO.parse(StringIO(raw), "fasta"))
            if not records:
                raise ValueError("The pasted FASTA content does not contain any sequence record.")
            sequence = "".join(str(record.seq) for record in records)
        else:
            sequence = "".join(raw.split())
    else:
        raise ValueError("Please upload a FASTA file or paste a DNA sequence.")

    sequence = sequence.upper()
    allowed = {"A", "C", "G", "T", "N"}
    invalid = sorted({char for char in sequence if char not in allowed})
    if invalid:
        invalid_str = ", ".join(invalid[:10])
        raise ValueError(f"Invalid DNA characters found: {invalid_str}")
    if not sequence:
        raise ValueError("Sequence is empty after parsing.")
    return sequence


def build_spacing_bonus(spacing: int, optimal_spacing: int, sigma: float) -> float:
    distance = spacing - optimal_spacing
    return float(2.718281828 ** (-(distance * distance) / (2.0 * sigma * sigma)))


def _scan_candidates_on_sequence(
    sequence: str,
    motif_35,
    motif_10,
    min_spacing: int,
    max_spacing: int,
    min_norm_35: float,
    min_norm_10: float,
    min_combined: float,
    circular_sequence: bool,
) -> List[PromoterHit]:
    seq_len = len(sequence)
    box35_len = len(motif_35)
    box10_len = len(motif_10)
    required_span = box35_len + max_spacing + box10_len
    pssm_35 = motif_35.pssm
    pssm_10 = motif_10.pssm
    sigma = max(1.0, (max_spacing - min_spacing) / 2.0)
    optimal_spacing = 17
    context_len = required_span + 2
    scanned = []
    dedupe = set()

    strands = {
        "+": sequence,
        "-": str(Seq(sequence).reverse_complement()),
    }

    for strand, strand_seq in strands.items():
        if circular_sequence:
            effective_seq = strand_seq + strand_seq[:context_len]
            start_limit = seq_len
        else:
            effective_seq = strand_seq
            start_limit = max(0, seq_len - box35_len + 1)

        for pos35 in range(start_limit):
            try:
                raw35 = float(pssm_35.calculate(effective_seq[pos35 : pos35 + box35_len]))
            except Exception:
                continue

            norm35 = normalize_score(raw35, pssm_35.min, pssm_35.max)
            if norm35 < min_norm_35:
                continue

            for spacing in range(min_spacing, max_spacing + 1):
                pos10 = pos35 + box35_len + spacing
                end10 = pos10 + box10_len

                if not circular_sequence and end10 > seq_len:
                    continue

                try:
                    raw10 = float(pssm_10.calculate(effective_seq[pos10:end10]))
                except Exception:
                    continue

                norm10 = normalize_score(raw10, pssm_10.min, pssm_10.max)
                if norm10 < min_norm_10:
                    continue

                spacing_bonus = build_spacing_bonus(spacing, optimal_spacing, sigma)
                combined = 0.43 * norm35 + 0.47 * norm10 + 0.10 * spacing_bonus
                if combined < min_combined:
                    continue

                if strand == "+":
                    box35_start = pos35 % seq_len
                    box10_start = pos10 % seq_len
                else:
                    box35_start = (seq_len - ((pos35 % seq_len) + box35_len)) % seq_len
                    box10_start = (seq_len - ((pos10 % seq_len) + box10_len)) % seq_len

                dedupe_key = (strand, box35_start, box10_start, spacing)
                if dedupe_key in dedupe:
                    continue
                dedupe.add(dedupe_key)

                scanned.append(
                    PromoterHit(
                        strand=strand,
                        box35_start=box35_start,
                        box10_start=box10_start,
                        spacing=spacing,
                        raw35=raw35,
                        raw10=raw10,
                        norm35=norm35,
                        norm10=norm10,
                        spacing_bonus=spacing_bonus,
                        combined=combined,
                        promoter_seq=effective_seq[pos35:end10],
                    )
                )

    return scanned


def _shuffle_sequence_preserving_composition(sequence: str, rng: random.Random) -> str:
    chars = list(sequence)
    rng.shuffle(chars)
    return "".join(chars)


def _estimate_background_scores(
    sequence: str,
    motif_35,
    motif_10,
    min_spacing: int,
    max_spacing: int,
    min_norm_35: float,
    min_norm_10: float,
    min_combined: float,
    circular_sequence: bool,
    shuffles: int,
) -> List[float]:
    rng = random.Random(17)
    scores: List[float] = []
    for _ in range(max(1, shuffles)):
        shuffled = _shuffle_sequence_preserving_composition(sequence, rng)
        null_hits = _scan_candidates_on_sequence(
            sequence=shuffled,
            motif_35=motif_35,
            motif_10=motif_10,
            min_spacing=min_spacing,
            max_spacing=max_spacing,
            min_norm_35=min_norm_35,
            min_norm_10=min_norm_10,
            min_combined=min_combined,
            circular_sequence=circular_sequence,
        )
        scores.extend(hit.combined for hit in null_hits)
    return scores


def _distance(pos_a: int, pos_b: int, seq_len: int, circular_sequence: bool) -> int:
    direct = abs(pos_a - pos_b)
    if not circular_sequence:
        return direct
    return min(direct, seq_len - direct)


def _suppress_overlapping_hits(
    hits: List[PromoterHit],
    seq_len: int,
    circular_sequence: bool,
    overlap_distance: int,
) -> List[PromoterHit]:
    if overlap_distance <= 0:
        return hits

    selected: List[PromoterHit] = []
    for candidate in sorted(hits, key=lambda hit: hit.combined, reverse=True):
        conflict = False
        for kept in selected:
            if kept.strand != candidate.strand:
                continue
            if _distance(kept.box35_start, candidate.box35_start, seq_len, circular_sequence) > overlap_distance:
                continue
            if _distance(kept.box10_start, candidate.box10_start, seq_len, circular_sequence) > overlap_distance:
                continue
            conflict = True
            break
        if not conflict:
            selected.append(candidate)
    return selected


def scan_promoters(
    sequence: str,
    motif_35,
    motif_10,
    min_spacing: int,
    max_spacing: int,
    min_norm_35: float,
    min_norm_10: float,
    min_combined: float,
    max_hits: int,
    circular_sequence: bool,
    sensitivity_profile: str = "custom",
) -> List[PromoterHit]:
    if min_spacing > max_spacing:
        raise ValueError("Minimum spacing cannot be greater than maximum spacing.")
    if not 0 <= min_norm_35 <= 1 or not 0 <= min_norm_10 <= 1 or not 0 <= min_combined <= 1:
        raise ValueError("Normalized thresholds must be between 0.0 and 1.0.")
    if max_hits < 1:
        raise ValueError("Maximum hits must be at least 1.")

    seq_len = len(sequence)
    box35_len = len(motif_35)
    box10_len = len(motif_10)
    required_span = box35_len + max_spacing + box10_len

    if not circular_sequence and seq_len < required_span:
        raise ValueError(f"Sequence too short for current settings. Need at least {required_span} bases.")

    scanned = _scan_candidates_on_sequence(
        sequence=sequence,
        motif_35=motif_35,
        motif_10=motif_10,
        min_spacing=min_spacing,
        max_spacing=max_spacing,
        min_norm_35=min_norm_35,
        min_norm_10=min_norm_10,
        min_combined=min_combined,
        circular_sequence=circular_sequence,
    )
    if not scanned:
        return []

    alpha = PRESET_SIGNIFICANCE_ALPHA.get(sensitivity_profile)
    if alpha is not None:
        background_scores = _estimate_background_scores(
            sequence=sequence,
            motif_35=motif_35,
            motif_10=motif_10,
            min_spacing=min_spacing,
            max_spacing=max_spacing,
            min_norm_35=min_norm_35,
            min_norm_10=min_norm_10,
            min_combined=min_combined,
            circular_sequence=circular_sequence,
            shuffles=PRESET_BACKGROUND_SHUFFLES.get(sensitivity_profile, 2),
        )
        if background_scores:
            sorted_background = sorted(background_scores)
            total = len(sorted_background)
            significant_hits = []
            for hit in scanned:
                left_idx = bisect_left(sorted_background, hit.combined)
                greater_or_equal = total - left_idx
                hit.p_value = (greater_or_equal + 1) / (total + 1)
                if hit.p_value <= alpha:
                    significant_hits.append(hit)
            scanned = significant_hits

    if not scanned:
        return []

    scanned.sort(key=lambda hit: hit.combined, reverse=True)
    overlap_distance = PRESET_OVERLAP_DISTANCE.get(sensitivity_profile, 0)
    non_overlapping = _suppress_overlapping_hits(
        hits=scanned,
        seq_len=seq_len,
        circular_sequence=circular_sequence,
        overlap_distance=overlap_distance,
    )
    return non_overlapping[:max_hits]


def generate_visualization(seq_len: int, hits: List[PromoterHit], box35_len: int, box10_len: int) -> List[str]:
    lines = []
    for idx, hit in enumerate(hits, start=1):
        line = ["."] * seq_len
        marker35 = "#" if hit.strand == "+" else "="
        marker10 = "*" if hit.strand == "+" else "+"

        for i in range(box35_len):
            line[(hit.box35_start + i) % seq_len] = marker35
        for i in range(box10_len):
            line[(hit.box10_start + i) % seq_len] = marker10

        summary = (
            f"[{idx}] strand:{hit.strand} -35@{hit.box35_start} -10@{hit.box10_start} "
            f"spacing:{hit.spacing} combined:{hit.combined:.3f}"
        )
        compact = "".join(line[:200]) + " ..." if seq_len > 200 else "".join(line)
        lines.append(f"{summary}\n{compact}")
    return lines


def apply_preset_to_form_values(form_values: Dict[str, str]) -> None:
    preset = form_values.get("preset", "balanced")
    if preset == "custom":
        return
    if preset not in PRESET_VALUES:
        form_values["preset"] = "balanced"
        preset = "balanced"
    form_values.update(PRESET_VALUES[preset])


def _build_promoter_report(hits: List[PromoterHit], max_hits: int) -> str:
    lines = []
    for i, hit in enumerate(hits, start=1):
        header = (
            f">promoter_{i} strand:{hit.strand} -35:{hit.box35_start} -10:{hit.box10_start} "
            f"spacing:{hit.spacing} score:{hit.combined:.3f} "
            f"n35:{hit.norm35:.3f} n10:{hit.norm10:.3f} p:{hit.p_value:.4g}"
        )
        lines.append(f"{header}\n{hit.promoter_seq}")

    summary = f"# detected_promoters:{len(hits)}"
    if len(hits) >= max_hits:
        summary += f" (capped_at_max_hits:{max_hits})"
    return summary + "\n" + "\n".join(lines)


@app.route("/", methods=["GET", "POST"])
def index():
    form_values: Dict[str, str] = {
        "preset": "balanced",
        "min_spacing": "14",
        "max_spacing": "20",
        "min_norm_35": "0.35",
        "min_norm_10": "0.35",
        "min_combined": "0.45",
        "max_hits": "100",
        "sequence_input": "",
        "circular_sequence": "on",
    }

    sequence_value = ""
    sequence_length = 0
    input_error = ""
    motif_error = ""
    promoter_error = ""
    motif_analysis = None
    promoter_result = ""
    promoter_matches = []

    if request.method == "POST":
        form_values.update(
            {
                "min_spacing": request.form.get("min_spacing", form_values["min_spacing"]),
                "max_spacing": request.form.get("max_spacing", form_values["max_spacing"]),
                "min_norm_35": request.form.get("min_norm_35", form_values["min_norm_35"]),
                "min_norm_10": request.form.get("min_norm_10", form_values["min_norm_10"]),
                "min_combined": request.form.get("min_combined", form_values["min_combined"]),
                "max_hits": request.form.get("max_hits", form_values["max_hits"]),
                "preset": request.form.get("preset", form_values["preset"]),
                "sequence_input": request.form.get("sequence_input", ""),
                "circular_sequence": "on" if request.form.get("circular_sequence") else "",
            }
        )
        if form_values["preset"] not in PRESET_OPTIONS:
            form_values["preset"] = "balanced"
        apply_preset_to_form_values(form_values)

        try:
            fasta_file = request.files.get("fasta_file")
            sequence_value = parse_sequence_from_input(fasta_file, form_values["sequence_input"])
            sequence_length = len(sequence_value)
        except Exception as exc:
            input_error = str(exc)

        if not input_error:
            try:
                motif_analysis = analyze_sequence_v2(sequence_value)
            except Exception as exc:
                motif_error = f"DNA motifs analysis failed: {exc}"

            try:
                min_spacing = parse_int("min_spacing", form_values["min_spacing"])
                max_spacing = parse_int("max_spacing", form_values["max_spacing"])
                min_norm_35 = parse_float("min_norm_35", form_values["min_norm_35"])
                min_norm_10 = parse_float("min_norm_10", form_values["min_norm_10"])
                min_combined = parse_float("min_combined", form_values["min_combined"])
                max_hits = parse_int("max_hits", form_values["max_hits"])
                circular_sequence = form_values["circular_sequence"] == "on"

                motif_35 = load_jaspar(BASE_DIR / "sigma70_box35_fixed.jaspar")
                motif_10 = load_jaspar(BASE_DIR / "sigma70_box10_fixed.jaspar")
                promoter_hits = scan_promoters(
                    sequence=sequence_value,
                    motif_35=motif_35,
                    motif_10=motif_10,
                    min_spacing=min_spacing,
                    max_spacing=max_spacing,
                    min_norm_35=min_norm_35,
                    min_norm_10=min_norm_10,
                    min_combined=min_combined,
                    max_hits=max_hits,
                    circular_sequence=circular_sequence,
                    sensitivity_profile=form_values["preset"],
                )

                if promoter_hits:
                    promoter_result = _build_promoter_report(promoter_hits, max_hits=max_hits)
                    promoter_matches = generate_visualization(
                        seq_len=sequence_length,
                        hits=promoter_hits[:50],
                        box35_len=len(motif_35),
                        box10_len=len(motif_10),
                    )
                else:
                    promoter_result = (
                        "No promoter found with current thresholds. "
                        "Try lower normalized thresholds (e.g. 0.25) or widen spacing range."
                    )
            except Exception as exc:
                promoter_error = f"Promoter scan failed: {exc}"

    return render_template(
        "index.html",
        app_name=APP_NAME,
        form_values=form_values,
        preset_options=PRESET_OPTIONS,
        sequence_value=sequence_value,
        sequence_length=sequence_length,
        input_error=input_error,
        motif_error=motif_error,
        promoter_error=promoter_error,
        motif_analysis=motif_analysis,
        promoter_result=promoter_result,
        promoter_matches=promoter_matches,
    )


@app.route("/admin/config", methods=["GET", "POST"])
def admin_config():
    cfg_path = BASE_DIR / "config" / "motifs_config.json"
    message = ""
    if request.method == "POST":
        content = request.form.get("config_json", "")
        try:
            parsed = json.loads(content) if content.strip() else {}
            cfg_path.parent.mkdir(parents=True, exist_ok=True)
            cfg_path.write_text(json.dumps(parsed, indent=2), encoding="utf-8")
            message = "Config saved."
        except Exception as exc:
            message = f"Invalid JSON: {exc}"
        current = cfg_path.read_text(encoding="utf-8") if cfg_path.exists() else "{}"
        return render_template("admin.html", config_text=current, message=message, app_name=APP_NAME)

    current = cfg_path.read_text(encoding="utf-8") if cfg_path.exists() else "{}"
    return render_template("admin.html", config_text=current, message=message, app_name=APP_NAME)


def _open_browser(host: str, port: int) -> None:
    webbrowser.open_new_tab(f"http://{host}:{port}/")


if __name__ == "__main__":
    host = os.environ.get("HOST", "127.0.0.1")
    port = int(os.environ.get("PORT", "5000"))
    debug_mode = os.environ.get("FLASK_DEBUG", "1") != "0"

    if os.environ.get("WERKZEUG_RUN_MAIN") == "true" or not debug_mode:
        Timer(1.0, _open_browser, args=(host, port)).start()

    app.run(host=host, port=port, debug=debug_mode, use_reloader=debug_mode)
