# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bioinformatics tools for detection and typing of **phage-plasmids (P-Ps)** — mobile genetic elements with both phage and plasmid properties. Supports 10 P-P types: AB_1, P1_1, P1_2, N15, SSU5_pHCM2, pMT1, pCAV, pSLy3, pKpn, cp32.

Three detection methods:
1. **tyPPing** (primary) — Python CLI tool using HMM profile-based detection with MinProteins and Composition branches
2. **MM-GRC** — Random forest classification + weighted gene repertoire relatedness (wGRR), R scripts
3. **geNomad + vConTACT2** — Sequence classification and viral clustering pipeline

## Key Commands

### tyPPing setup and installation

```bash
cd tyPPing
pip install .          # Installs the `typping` CLI (Python >= 3.10)
unzip tyPPing_signature_profiles.zip  # First time only — extract HMM profiles
```

### Run tests

```bash
cd tyPPing
pip install pytest
pytest tests/ -v                    # All tests
pytest tests/test_io.py -v          # Unit tests only
pytest tests/test_integration.py -v # Integration tests (compare against R expected output)
```

### Test with example data (manual verification)

```bash
cd tyPPing
typping run \
    -m test_example/protein_to_genome.tsv \
    -s test_example/genome_size.tsv \
    -i test_example/tyPPing_example_hmmsearch_output.tbl.out \
    -o /tmp/typping_test/ \
    --compositions tyPPing_input_data/Compositions_information_table.tsv \
    --profiles tyPPing_input_data/Profile_information_table.tsv

# Verify output matches expected
diff /tmp/typping_test/Final_prediction_table.tsv test_example/Final_prediction_table.tsv
diff /tmp/typping_test/All_hmm_hits_table.tsv test_example/All_hmm_hits_table.tsv
```

### Docker

```bash
cd tyPPing
docker build -t typping .
docker run --rm -v /path/to/data:/data typping run \
    -m /data/protein_to_genome.tsv -s /data/genome_size.tsv \
    -i /data/hmm_output.tbl.out -o /data/output/ \
    --compositions /opt/tyPPing/tyPPing_input_data/Compositions_information_table.tsv \
    --profiles /opt/tyPPing/tyPPing_input_data/Profile_information_table.tsv
```

### MM-GRC (R-based, separate method)

```bash
R < MM-GRC/wGRR_MGEs.R --no-save
```

## Architecture

### tyPPing Python package (`tyPPing/src/typping/`)

The `typping` CLI (typer-based) has three subcommands: `run`, `prepare`, `hmmsearch`.

**Pipeline flow (`typping run`):**

1. **`io.py`** — Parses HMMER domtblout (22 fixed fields + description), reads TSV inputs (protein_to_genome, genome_size, compositions, profiles). `preprocess_hmm_output()` computes profile coverage, joins with protein-to-genome mapping, keeps best hit per genome+profile, extracts `PP_category` from HMM profile names via regex `^(.+?)_pers_`.

2. **`prediction.py`** — Two parallel branches, each with complete/draft variants:
   - **MinProteins**: filters hits by sequence score threshold, builds presence/absence matrix per genome (or per contig in draft mode), computes weighted scores using `profile_scores_1` and `profile_scores_2`. Draft mode filters contigs with `score_0 >= 0.1` before aggregating to genome level.
   - **Composition**: filters by `cov_profile >= 0.5`, joins with composition sets, requires exact match OR (hits >= size - tolerance AND hits >= 75% of size). Selects best composition per genome.

3. **`merge.py`** — Full-joins MinProteins and Composition results, assigns confidence (High/Medium/Low) based on: which branches predicted, genome size within type-specific ranges from `ALL_CUTOFFS`, and special deduplication logic for `CATEGORIES_SPECIAL` = {cp32, SSU5_pHCM2, pKpn, pSLy3}.

4. **`constants.py`** — `ALL_CUTOFFS` DataFrame with per-type thresholds (MinProteins_min_N, composition_tol_thr, size_min/max, pct10_mean_size). PP_category mapping dicts between short codes (AB, P11, P12, SSU5) and full names.

5. **`prepare.py`** — Generates protein_to_genome.tsv and genome_size.tsv from FASTA files (Prodigal naming convention).

6. **`hmmsearch.py`** — Subprocess wrapper for `hmmsearch`.

**Mode detection:** Auto-detects complete vs. draft by checking if `contig_id` column exists in protein_to_genome file. Draft mode works at contig level first, then aggregates to genome.

### Key design decisions

- Python output must match the original R scripts (`tyPPing.R` / `tyPPing_for_draft_genomes.R`) exactly — integration tests compare TSV output cell-by-cell against R-generated expected files.
- Uses `_r_round()` helper (round-half-up) to match R's rounding behavior, since Python's `round()` uses banker's rounding.
- Output formatting converts NaN to "NA" strings and removes trailing ".0" from scores to match R output format.
- The `_format_output()` function in `cli.py` handles this R-compatibility formatting before writing TSVs.

### Input data

- `tyPPing_input_data/Profile_information_table.tsv` — 763 HMM profiles with sequence score thresholds and P-P/type scores
- `tyPPing_input_data/Compositions_information_table.tsv` — Composition sets (groups of profiles that define a P-P type)
- `test_example/` — 20-genome test dataset with expected output files

### External data on Zenodo (doi:10.5281/zenodo.16616313)

HMM profiles (`tyPPing_signature_profiles.hmm`, `phage.hmm`), pre-trained random forest models, and figure reproduction data.

## Code Conventions

- Python package uses hatchling build system, typer for CLI, pandas/numpy for data processing
- R scripts in `MM-GRC/` and `Publication_related_data/` use `data.table::fread()` + tidyverse
- No linter configured; no formal R test framework — R validation via example datasets
- `Publication_related_data/` contains figure reproduction scripts (not part of the tool itself)
