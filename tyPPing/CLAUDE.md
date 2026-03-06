# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Python CLI reimplementation of the tyPPing R scripts for phage-plasmid (P-P) detection and typing. Uses HMM profile-based detection to classify genomes into 10 P-P types: AB_1, P1_1, P1_2, N15, SSU5_pHCM2, pMT1, pCAV, pSLy3, pKpn, cp32.

## Key Commands

```bash
# Install (editable for development)
pip install -e .

# Run all tests
pytest tests/ -v

# Run unit tests only
pytest tests/test_io.py -v

# Run integration tests only (compares output against R-generated expected files)
pytest tests/test_integration.py -v

# Run a single test
pytest tests/test_io.py::test_parse_hmmer_domtblout -v

# Manual verification against expected R output
typping run \
    -m test_example/protein_to_genome.tsv \
    -s test_example/genome_size.tsv \
    -i test_example/tyPPing_example_hmmsearch_output.tbl.out \
    -o /tmp/typping_test/ \
    --compositions tyPPing_input_data/Compositions_information_table.tsv \
    --profiles tyPPing_input_data/Profile_information_table.tsv

diff /tmp/typping_test/Final_prediction_table.tsv test_example/Final_prediction_table.tsv
diff /tmp/typping_test/All_hmm_hits_table.tsv test_example/All_hmm_hits_table.tsv
```

## Architecture

### Pipeline flow (`typping run`)

```
HMMER domtblout  ──→  io.py (parse + preprocess)  ──→  prediction.py  ──→  merge.py  ──→  Output TSVs
                                                         ├─ MinProteins branch
                                                         └─ Composition branch
```

**`io.py`** — Parses HMMER domtblout (22 fixed fields + description), reads TSV inputs. `preprocess_hmm_output()` computes profile coverage, joins with protein-to-genome mapping, keeps best hit per genome+profile, extracts `PP_category` from HMM profile names via regex `^(.+?)_pers_`.

**`prediction.py`** — Two parallel branches, each with complete/draft variants:
- **MinProteins**: filters by sequence score threshold, builds presence/absence matrix, computes weighted scores using `profile_scores_1` and `profile_scores_2`. Draft mode filters contigs with `score_0 >= 0.1` before aggregating to genome level.
- **Composition**: filters by `cov_profile >= 0.5`, joins with composition sets, requires exact match OR (hits >= size - tolerance AND hits >= 75% of size). Selects best composition per genome.

**`merge.py`** — Full-joins MinProteins and Composition results, assigns confidence (High/Medium/Low) based on: which branches predicted, genome size within type-specific ranges from `ALL_CUTOFFS`, and special deduplication logic for `CATEGORIES_SPECIAL` = {cp32, SSU5_pHCM2, pKpn, pSLy3}.

**`constants.py`** — `ALL_CUTOFFS` DataFrame with per-type thresholds. PP_category mapping dicts between short codes (AB, P11, P12, SSU5) and full names.

**`prepare.py`** — Generates protein_to_genome.tsv and genome_size.tsv from FASTA files (Prodigal naming convention).

**`hmmsearch.py`** — Subprocess wrapper for `hmmsearch`.

### Mode detection

Auto-detects complete vs. draft by checking if `contig_id` column exists in protein_to_genome file. Draft mode works at contig level first, then aggregates to genome.

### CLI subcommands

- `typping run` — Main prediction pipeline
- `typping prepare` — Generate input tables from FASTA files
- `typping hmmsearch` — Convenience wrapper for hmmsearch

## Critical Design Constraints

### R output compatibility

**This is the single most important constraint.** Python output must match the original R scripts exactly — integration tests compare TSV output cell-by-cell against R-generated expected files in `test_example/`.

Specific compatibility measures:
- `_r_round()` helper in `prediction.py` uses round-half-up to match R's `round()`, since Python's `round()` uses banker's rounding.
- `_format_output()` in `cli.py` converts NaN to `"NA"` strings and removes trailing `.0` from scores to match R output format.
- Integer columns are formatted as strings without decimals; score columns strip `.0` suffixes.

### After any change, always verify

```bash
pytest tests/test_integration.py -v
```

This runs `typping run` on the 20-genome test dataset and does cell-by-cell comparison of both output TSVs against R-generated expected files.

## Input Data

- `tyPPing_input_data/Profile_information_table.tsv` — 763 HMM profiles with sequence score thresholds and P-P/type scores
- `tyPPing_input_data/Compositions_information_table.tsv` — Composition sets (groups of profiles that define a P-P type)
- `test_example/` — 20-genome test dataset with expected output files
- `tyPPing_signature_profiles.zip` — HMM profiles (unzip before first use)

## Build System

- hatchling build backend, configured in `pyproject.toml`
- Python >= 3.10
- Dependencies: pandas, numpy, typer, biopython
- No linter configured
