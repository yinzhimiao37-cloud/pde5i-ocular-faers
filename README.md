# Three-drug PDE5 inhibitor ocular FAERS analysis

This repository contains the analysis code, frozen protocols, non-identifying summary tables, and publication figures for the comparison of tadalafil, sildenafil, and vardenafil in the FDA Adverse Event Reporting System (FAERS), covering 2013Q1–2025Q4.

## Scope and reproducibility

The repository is an auditable code and results archive. It does not redistribute FDA raw FAERS files, report-level identifiers, the local DuckDB database, or intermediate tables containing `PRIMARYID` or `CASEID`. The raw quarterly files must be downloaded directly from the FDA and kept locally.

The analysis scripts were run in the original project layout and expect an audited local input tree containing the FAERS working database and the mapped three-drug exposure tables. Because those inputs are not included here, this archive is not advertised as a zero-setup one-click rerun. The scripts, protocols, input manifest hashes, frozen ocular PT dictionary, summary outputs, and QA record are provided for review and independent audit.

## Analysis definition

- FAERS window: 2013Q1–2025Q4 (52 quarters).
- Primary role: `ROLE_COD=PS`; all-role analyses are sensitivity analyses.
- Drugs: tadalafil, sildenafil, and vardenafil using exact generic or declared ingredient-token mapping; fixed brand-name mappings are sensitivity analyses.
- Outcomes: a project-specific, semantically reviewed ocular preferred-term dictionary; it is not claimed to be the licensed MedDRA Eye Disorders SOC.
- Statistics: reporting odds ratios, proportional reporting ratios, BCPNN information components, exact intervals for sparse comparisons, multiplicity adjustment, direct active-comparator analyses, country influence analyses, and source-dependence audits.
- Unit of analysis: FAERS reports after CASEID version selection. Report counts are not confirmed independent patient counts.

## Layout

`review_20261013/mural_microenvironment_extension/FAERS_class_compare/` contains the three-drug comparison scripts, frozen protocol, summary tables, and figures. `analysis/faers_bayes.py` contains the audited BCPNN implementation. `review_20261014_protocol.md` records the ocular-spectrum extension protocol.

## Local run prerequisites

Python 3.12 is recommended. Install:

```powershell
python -m pip install -r requirements.txt
```

The full rerun additionally requires the original local FAERS extraction tree and the intermediate audited tables described in the scripts. Run from the repository root:

```powershell
python review_20261013/mural_microenvironment_extension/FAERS_class_compare/run_class_compare.py
python review_20261013/mural_microenvironment_extension/FAERS_class_compare/source_followup.py
```

The final plotting and QA script consumes the frozen intermediate results:

```powershell
python review_20261013/mural_microenvironment_extension/FAERS_class_compare/finalize_study.py
```

The protocol files are documentation, not executable code. Before rerunning with a different FAERS release, create a new output directory and record the release checksums; do not mix cached intermediate tables from different input windows.

## Main outputs

- `results/ocular_PT_primary_PS.csv`: complete primary-suspect PT-level table.
- `results/ocular_PT_primary_PS_signals.csv`: rows meeting the prespecified signal-priority rule.
- `results/active_comparator_evaluable.csv`: direct drug comparisons meeting the minimum event-count rule.
- `results/active_comparator_BH_significant.csv`: multiplicity-adjusted direct comparisons.
- `results/final_QA.json`: machine-readable integrity checks.
- `figures/`: publication figure exports.

## Data availability wording for the manuscript

The FAERS quarterly releases are publicly available from the U.S. Food and Drug Administration. The analysis code, frozen ocular preferred-term dictionary, non-identifying derived summary tables, protocols, and publication figures are archived in this repository. A permanent public DOI should be added after the repository is published or deposited in Zenodo.

## Interpretation boundary

These analyses identify reporting disproportionality and source-dependent reporting patterns. They do not estimate event frequency, establish causality, or rank clinical safety between drugs.
