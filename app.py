from dataclasses import dataclass
from bisect import bisect_left
from datetime import datetime, timezone
from io import BytesIO, StringIO
import json
import os
from pathlib import Path
import random
from threading import Timer
from typing import Dict, List
import webbrowser

from Bio import SeqIO, motifs
from Bio.Seq import Seq
from flask import Flask, render_template, request, send_file

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


def _run_dashboard_analysis(sequence_value: str, form_values: Dict[str, str]) -> Dict[str, object]:
    motif_error = ""
    promoter_error = ""
    motif_analysis = None
    promoter_result = ""
    promoter_matches: List[str] = []
    promoter_hits: List[PromoterHit] = []
    sequence_length = len(sequence_value)

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

    return {
        "sequence_length": sequence_length,
        "motif_analysis": motif_analysis,
        "motif_error": motif_error,
        "promoter_error": promoter_error,
        "promoter_result": promoter_result,
        "promoter_matches": promoter_matches,
        "promoter_hits": promoter_hits,
    }


def _pdf_safe_text(value) -> str:
    text = str(value if value is not None else "")
    return text.encode("latin-1", "replace").decode("latin-1")


def _clip_text(value, max_len: int = 180) -> str:
    text = _pdf_safe_text(value).replace("\n", " ")
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def _build_copy_paste_table(headers: List[str], rows: List[List[str]]) -> str:
    safe_headers = [_pdf_safe_text(h) for h in headers]
    safe_rows = [[_pdf_safe_text(cell) for cell in row] for row in rows]
    widths = [len(h) for h in safe_headers]
    for row in safe_rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(cell))

    def build_row(cells: List[str]) -> str:
        return " | ".join(cells[i].ljust(widths[i]) for i in range(len(widths)))

    separator = "-+-".join("-" * w for w in widths)
    lines = [build_row(safe_headers), separator]
    for row in safe_rows:
        lines.append(build_row(row))
    return "\n".join(lines)


def _build_pdf_report(
    sequence_value: str,
    form_values: Dict[str, str],
    analysis_result: Dict[str, object],
    input_error: str = "",
) -> bytes:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import PageBreak, Paragraph, Preformatted, SimpleDocTemplate, Spacer, Table, TableStyle

    generated_on = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    motif_analysis = analysis_result.get("motif_analysis")
    promoter_hits: List[PromoterHit] = analysis_result.get("promoter_hits", [])  # type: ignore[assignment]
    motif_error = str(analysis_result.get("motif_error", ""))
    promoter_error = str(analysis_result.get("promoter_error", ""))

    buffer = BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        leftMargin=16 * mm,
        rightMargin=16 * mm,
        topMargin=16 * mm,
        bottomMargin=14 * mm,
        title="Promoter Motif Dashboard Report",
        author=APP_NAME,
    )

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        name="ReportTitle",
        parent=styles["Heading1"],
        fontName="Helvetica-Bold",
        fontSize=20,
        leading=24,
        textColor=colors.HexColor("#1b2a3a"),
        spaceAfter=6,
    )
    section_style = ParagraphStyle(
        name="SectionHeading",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=13,
        leading=16,
        textColor=colors.HexColor("#243f5a"),
        spaceBefore=8,
        spaceAfter=6,
    )
    body_style = ParagraphStyle(
        name="Body",
        parent=styles["BodyText"],
        fontName="Helvetica",
        fontSize=9,
        leading=12,
    )
    mono_style = ParagraphStyle(
        name="Mono",
        parent=styles["Code"],
        fontName="Courier",
        fontSize=7.6,
        leading=9.2,
    )

    story = []
    story.append(Paragraph("Promoter + DNA Motifs Dashboard Report", title_style))
    story.append(Paragraph(_pdf_safe_text(f"Generated on: {generated_on}"), body_style))
    story.append(Paragraph(_pdf_safe_text(f"Sequence length: {len(sequence_value)} bp"), body_style))
    story.append(
        Paragraph(
            _pdf_safe_text(
                f"Preset: {form_values['preset']} | Circular sequence: {'yes' if form_values['circular_sequence'] == 'on' else 'no'}"
            ),
            body_style,
        )
    )
    story.append(Spacer(1, 6))

    if input_error:
        story.append(Paragraph(_pdf_safe_text(f"Input parsing error: {input_error}"), body_style))

    settings_rows = [
        ["Min spacing", form_values["min_spacing"]],
        ["Max spacing", form_values["max_spacing"]],
        ["Min -35 normalized score", form_values["min_norm_35"]],
        ["Min -10 normalized score", form_values["min_norm_10"]],
        ["Min combined score", form_values["min_combined"]],
        ["Max promoter hits", form_values["max_hits"]],
    ]
    settings_table = Table([["Parameter", "Value"]] + settings_rows, colWidths=[70 * mm, 40 * mm])
    settings_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#d9e5f2")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#0f2236")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#95a8bc")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7fafc")]),
            ]
        )
    )
    story.append(Paragraph("Analysis Parameters", section_style))
    story.append(settings_table)
    story.append(Spacer(1, 8))

    summary_rows = [
        ["Detected promoters", str(len(promoter_hits))],
        [
            "Detected motifs",
            str(len(motif_analysis["motifs"])) if motif_analysis and motif_analysis.get("motifs") else "0",
        ],
        [
            "Risk total (config-driven)",
            f"{motif_analysis['risk']['total']:.6f}" if motif_analysis and motif_analysis.get("risk") else "n/a",
        ],
        [
            "Risk total (hardcoded v2)",
            f"{motif_analysis['risk_new']['total']:.6f}" if motif_analysis and motif_analysis.get("risk_new") else "n/a",
        ],
    ]
    summary_table = Table([["Metric", "Value"]] + summary_rows, colWidths=[80 * mm, 30 * mm])
    summary_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dae7da")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#103818")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#96b296")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7fbf7")]),
            ]
        )
    )
    story.append(Paragraph("Executive Summary", section_style))
    story.append(summary_table)
    story.append(Spacer(1, 8))

    story.append(
        Paragraph(
            _pdf_safe_text(
                "Interpretation guide: promoter combined score blends normalized -35 and -10 box quality with spacing quality. "
                "For motifs, non-heuristic motifs are pattern matches while heuristic motifs are signal-driven segments."
            ),
            body_style,
        )
    )

    story.append(Paragraph("DNA Motif Details", section_style))
    if motif_error:
        story.append(Paragraph(_pdf_safe_text(motif_error), body_style))
    elif motif_analysis and motif_analysis.get("motifs"):
        motif_table_data = [["Motif", "Category", "Type", "Score", "Details"]]
        for motif_item in motif_analysis["motifs"]:
            if motif_item.get("heuristic") and motif_item.get("segments") is not None:
                details = (
                    f"segments={motif_item.get('count_segments', 0)}, "
                    f"coverage={motif_item.get('coverage_fraction', 0.0) * 100:.2f}%"
                )
            elif motif_item.get("positions") is not None:
                positions = motif_item.get("positions", [])
                details = ", ".join(str(pos) for pos in positions[:12])
                if len(positions) > 12:
                    details += f" (+{len(positions) - 12} more)"
            else:
                details = "-"
            motif_table_data.append(
                [
                    _clip_text(motif_item.get("name", ""), 32),
                    _clip_text(motif_item.get("category", ""), 14),
                    "heuristic" if motif_item.get("heuristic") else "regex",
                    f"{motif_item.get('score', 0):.4g}",
                    _clip_text(details, 66),
                ]
            )
        motif_table = Table(
            motif_table_data,
            colWidths=[38 * mm, 24 * mm, 18 * mm, 18 * mm, 72 * mm],
            repeatRows=1,
        )
        motif_table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f3e6ce")),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                    ("FONTSIZE", (0, 0), (-1, -1), 8),
                    ("ALIGN", (3, 1), (3, -1), "RIGHT"),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#c4b79f")),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#fffaf2")]),
                ]
            )
        )
        story.append(motif_table)
    else:
        story.append(Paragraph("No motifs found.", body_style))

    story.append(Spacer(1, 8))
    story.append(Paragraph("Risk Model Details", section_style))
    if motif_analysis and motif_analysis.get("risk"):
        risk_data = motif_analysis["risk"]
        raw_by_category = risk_data.get("by_category", {})
        weighted_by_category = risk_data.get("by_category_weighted", {})
        risk_rows = []
        for category in sorted(set(raw_by_category.keys()) | set(weighted_by_category.keys())):
            raw_val = float(raw_by_category.get(category, 0.0))
            weighted_val = float(weighted_by_category.get(category, 0.0))
            risk_rows.append([_clip_text(category, 20), f"{raw_val:.6f}", f"{weighted_val:.6f}"])
        if risk_rows:
            risk_table = Table([["Category", "Raw score", "Weighted score"]] + risk_rows, colWidths=[52 * mm, 32 * mm, 34 * mm])
            risk_table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8e1f2")),
                        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
                        ("ALIGN", (1, 1), (-1, -1), "RIGHT"),
                        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#b8accc")),
                        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#faf8fd")]),
                    ]
                )
            )
            story.append(risk_table)

        risk_notes = risk_data.get("notes", [])
        synergy_flags = risk_data.get("synergy_flags", [])
        if synergy_flags:
            story.append(Spacer(1, 4))
            story.append(Paragraph(_pdf_safe_text("Triggered synergy flags: " + ", ".join(synergy_flags)), body_style))
        if risk_notes:
            story.append(Paragraph(_pdf_safe_text("Risk notes: " + " | ".join(risk_notes)), body_style))

        hardcoded = motif_analysis.get("risk_new", {}).get("by_motif", {})
        if hardcoded:
            story.append(Spacer(1, 6))
            story.append(Paragraph("Hardcoded v2 motif contributions", body_style))
            hardcoded_rows = []
            for motif_name, value in sorted(hardcoded.items(), key=lambda item: item[1], reverse=True):
                hardcoded_rows.append([_clip_text(motif_name, 58), f"{float(value):.6f}"])
            hardcoded_table = Table([["Motif", "Contribution"]] + hardcoded_rows, colWidths=[98 * mm, 20 * mm], repeatRows=1)
            hardcoded_table.setStyle(
                TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#efe8d8")),
                        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                        ("FONTSIZE", (0, 0), (-1, -1), 8),
                        ("ALIGN", (1, 1), (1, -1), "RIGHT"),
                        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#ccbfa8")),
                    ]
                )
            )
            story.append(hardcoded_table)
    else:
        story.append(Paragraph("Risk data unavailable due to analysis error.", body_style))

    story.append(Spacer(1, 8))
    story.append(Paragraph("Promoter Candidate Details", section_style))
    if promoter_error:
        story.append(Paragraph(_pdf_safe_text(promoter_error), body_style))
    elif promoter_hits:
        promoter_table_data = [
            [
                "Rank",
                "Strand",
                "-35 start",
                "-10 start",
                "Spacing",
                "Norm -35",
                "Norm -10",
                "Combined",
                "p-value",
                "Promoter sequence",
            ]
        ]
        for idx, hit in enumerate(promoter_hits, start=1):
            promoter_table_data.append(
                [
                    str(idx),
                    hit.strand,
                    str(hit.box35_start),
                    str(hit.box10_start),
                    str(hit.spacing),
                    f"{hit.norm35:.3f}",
                    f"{hit.norm10:.3f}",
                    f"{hit.combined:.3f}",
                    f"{hit.p_value:.4g}",
                    _clip_text(hit.promoter_seq, 36),
                ]
            )
        promoter_table = Table(
            promoter_table_data,
            colWidths=[10 * mm, 12 * mm, 16 * mm, 16 * mm, 13 * mm, 16 * mm, 16 * mm, 15 * mm, 16 * mm, 52 * mm],
            repeatRows=1,
        )
        promoter_table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dceaf1")),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                    ("FONTSIZE", (0, 0), (-1, -1), 7.6),
                    ("ALIGN", (0, 1), (8, -1), "CENTER"),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#9cb2c1")),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f8fbfd")]),
                ]
            )
        )
        story.append(promoter_table)
    else:
        story.append(Paragraph("No promoters found with the current thresholds.", body_style))

    story.append(PageBreak())
    story.append(Paragraph("Copy-Paste Tables", section_style))
    story.append(
        Paragraph(
            _pdf_safe_text(
                "These plain-text tables are included for direct copy and paste into spreadsheets or notebooks."
            ),
            body_style,
        )
    )
    story.append(Spacer(1, 4))

    if motif_analysis and motif_analysis.get("motifs"):
        motif_rows = []
        for motif_item in motif_analysis["motifs"]:
            if motif_item.get("segments") is not None:
                location = "; ".join(
                    f"[{seg.get('start')}-{seg.get('end')}] {seg.get('length_bp')}bp"
                    for seg in motif_item.get("segments", [])
                )
                if not location:
                    location = "segments: none"
            else:
                location = ",".join(str(pos) for pos in motif_item.get("positions", []))
            motif_rows.append(
                [
                    _clip_text(motif_item.get("name", ""), 42),
                    _clip_text(motif_item.get("category", ""), 16),
                    "heuristic" if motif_item.get("heuristic") else "regex",
                    _clip_text(location, 130),
                    str(motif_item.get("score", 0)),
                    _clip_text(motif_item.get("impact", ""), 82),
                ]
            )
        motif_copy_table = _build_copy_paste_table(
            ["motif", "category", "type", "locations", "score", "impact"],
            motif_rows,
        )
        story.append(Paragraph("Motif table", body_style))
        story.append(Preformatted(motif_copy_table, mono_style))
        story.append(Spacer(1, 8))

    if promoter_hits:
        promoter_rows = []
        for idx, hit in enumerate(promoter_hits, start=1):
            promoter_rows.append(
                [
                    str(idx),
                    hit.strand,
                    str(hit.box35_start),
                    str(hit.box10_start),
                    str(hit.spacing),
                    f"{hit.norm35:.4f}",
                    f"{hit.norm10:.4f}",
                    f"{hit.combined:.4f}",
                    f"{hit.p_value:.5g}",
                    _clip_text(hit.promoter_seq, 120),
                ]
            )
        promoter_copy_table = _build_copy_paste_table(
            ["rank", "strand", "box35_start", "box10_start", "spacing", "norm35", "norm10", "combined", "p_value", "promoter_seq"],
            promoter_rows,
        )
        story.append(Paragraph("Promoter table", body_style))
        story.append(Preformatted(promoter_copy_table, mono_style))

    story.append(PageBreak())
    story.append(Paragraph("Glossary", section_style))
    glossary_items = [
        ("-35 box", "Upstream promoter element recognized by sigma factor near ~35 bases before transcription start."),
        ("-10 box", "Promoter element close to ~10 bases before transcription start, often rich in A/T."),
        ("Normalized score", "Motif score scaled between motif-specific minimum and maximum to compare candidates."),
        ("Combined score", "Weighted score using normalized -35, normalized -10, and spacing bonus."),
        ("Spacing", "Distance in base pairs between the end of -35 box and start of -10 box."),
        ("p-value", "Empirical significance estimate against composition-preserving shuffled background."),
        ("Heuristic motif", "Pattern detected by sequence properties (for example GC-rich windows), not only regex."),
        ("Coverage fraction", "Fraction of the full sequence covered by merged heuristic segments."),
        ("Synergy rule", "Rule that modifies risk when specific motif combinations are present together."),
    ]
    for term, definition in glossary_items:
        story.append(Paragraph(_pdf_safe_text(f"<b>{term}</b>: {definition}"), body_style))
        story.append(Spacer(1, 2))

    def draw_footer(canvas, doc_obj):
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#5a6775"))
        canvas.drawRightString(A4[0] - 16 * mm, 8 * mm, f"Page {doc_obj.page}")
        canvas.drawString(16 * mm, 8 * mm, _pdf_safe_text(APP_NAME))

    doc.build(story, onFirstPage=draw_footer, onLaterPages=draw_footer)
    pdf_bytes = buffer.getvalue()
    buffer.close()
    return pdf_bytes


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
    report_ready = False

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
            analysis_result = _run_dashboard_analysis(sequence_value, form_values)
            motif_analysis = analysis_result["motif_analysis"]
            motif_error = str(analysis_result["motif_error"])
            promoter_error = str(analysis_result["promoter_error"])
            promoter_result = str(analysis_result["promoter_result"])
            promoter_matches = analysis_result["promoter_matches"]
            report_ready = True

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
        report_ready=report_ready,
    )


@app.route("/export/report.pdf", methods=["POST"])
def export_report_pdf():
    form_values: Dict[str, str] = {
        "preset": request.form.get("preset", "balanced"),
        "min_spacing": request.form.get("min_spacing", "14"),
        "max_spacing": request.form.get("max_spacing", "20"),
        "min_norm_35": request.form.get("min_norm_35", "0.35"),
        "min_norm_10": request.form.get("min_norm_10", "0.35"),
        "min_combined": request.form.get("min_combined", "0.45"),
        "max_hits": request.form.get("max_hits", "100"),
        "sequence_input": request.form.get("sequence_input", ""),
        "circular_sequence": "on" if request.form.get("circular_sequence") else "",
    }
    if form_values["preset"] not in PRESET_OPTIONS:
        form_values["preset"] = "balanced"
    apply_preset_to_form_values(form_values)

    input_error = ""
    sequence_value = ""
    sequence_for_report = request.form.get("sequence_for_report", "").strip()
    try:
        if sequence_for_report:
            sequence_value = parse_sequence_from_input(None, sequence_for_report)
        else:
            fasta_file = request.files.get("fasta_file")
            sequence_value = parse_sequence_from_input(fasta_file, form_values["sequence_input"])
    except Exception as exc:
        input_error = str(exc)

    if input_error:
        return f"Unable to export PDF report: {input_error}", 400

    analysis_result = _run_dashboard_analysis(sequence_value, form_values)
    try:
        pdf_bytes = _build_pdf_report(sequence_value, form_values, analysis_result, input_error=input_error)
    except ImportError as exc:
        return f"PDF export unavailable: {exc}. Install dependencies from requirements.txt.", 500
    except Exception as exc:
        return f"Failed to generate PDF report: {exc}", 500

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"promoter_motif_report_{stamp}.pdf"
    return send_file(
        BytesIO(pdf_bytes),
        mimetype="application/pdf",
        as_attachment=True,
        download_name=filename,
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
