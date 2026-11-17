# Three-drug FAERS ocular spectrum analysis

This directory contains a reproducible analysis of tadalafil, sildenafil, and vardenafil ocular reporting patterns in FAERS, 2013Q1-2025Q4.

## Run order

Use Python 3.12 with `pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, and `Pillow`.

```powershell
py -3 extract_all_PT.py
py -3 build_ocular_dictionary.py
py -3 analyze_ocular_PT.py
py -3 source_and_version_audit.py
py -3 build_ocular_background_and_sensitivity.py
py -3 independent_count_audit.py
py -3 finalize_study.py
```

The first extraction scripts depend on the local FAERS paths already recorded in the scripts. `finalize_study.py` reads only frozen intermediate results and can be rerun quickly.

## Principal outputs

- `results/ocular_PT_signal_stability.csv`: PT-level estimates and stability audits.
- `results/active_comparator_evaluable.csv`: direct comparisons with at least three events per drug.
- `results/active_comparator_BH_significant.csv`: multiplicity-adjusted direct differences.
- `results/final_QA.json`: machine-readable integrity checks.
- `figures/*.pdf`, `*.tiff`, and `*.png`: publication outputs.
- `manuscript_draft.md`: publication-ready English draft and figure legends.

## Interpretation boundary

Rows are FAERS reports, not verified independent patients. Source-family collapse is an audit of potential dependence, not patient deduplication. The ocular PT dictionary is project-specific and is not the official MedDRA Eye disorders SOC. Results describe reporting disproportionality and cannot estimate event frequency, establish causality, or provide a clinical safety ranking.
