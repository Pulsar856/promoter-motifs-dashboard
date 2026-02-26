# Promoter + DNA Motifs Dashboard

Single Flask app that merges:
- `promoter-finder-flask` (Sigma70 promoter scan)
- `dna_motifs_app` (motif detection + risk scoring)

## Features
- One dashboard form for DNA input:
  - Paste sequence/FASTA text
  - Upload FASTA file
- Promoter scan results (FASTA-like output + compact visualization)
- DNA motif and risk panel (v2 analysis schema)
- Professional PDF report export (detailed sections, wrapped print-safe tables, copy-paste tables, glossary)
- Browser auto-opens on app start
- Motif config editor at `/admin/config`

## Project Layout
- `app.py`: merged Flask backend + routes
- `analysis.py`, `motifs_data.py`: motif/risk engine
- `sigma70_box35_fixed.jaspar`, `sigma70_box10_fixed.jaspar`: promoter PWM files
- `config/motifs_config.json`: optional motif/risk overrides
- `templates/index.html`: unified dashboard UI
- `templates/admin.html`: config editor

## Run
```bash
.\setup_venv.ps1
.\run.ps1
```

Alternative (manual):
```bash
.venv\Scripts\activate
python app.py
```

`app.py` auto-switches to `.\.venv\Scripts\python.exe` when available, so launching `python app.py`
from a global interpreter no longer fails on missing local dependencies.

Default URL: `http://127.0.0.1:5000`
