# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Bioinformatics tools for detection and typing of **phage-plasmids (P-Ps)** — mobile genetic elements with both phage and plasmid properties. Accompanies the publication on efficient P-P detection. Supports 10 P-P types: AB_1, P1_1, P1_2, N15, SSU5_pHCM2, pMT1, pCAV, pSLy3, pKpn, cp32.

Three detection methods are provided:
1. **tyPPing** (primary) — HMM profile-based detection with MinProteins and Composition branches
2. **MM-GRC** — Random forest classification + weighted gene repertoire relatedness (wGRR)
3. **geNomad + vConTACT2** — Sequence classification and viral clustering pipeline

## Key Commands

### tyPPing (primary tool)

```bash
# Step 1: HMMER search (prerequisite)
hmmsearch -o tmp.out --domtblout hmm_output.tbl.out tyPPing_signature_profiles.hmm proteins.faa

# Step 2: Run tyPPing for complete genomes
Rscript tyPPing/tyPPing.R \
  -m protein_to_genome.tsv -s genome_size.tsv \
  -i hmm_output.tbl.out -o output_dir/ \
  --compositions tyPPing/tyPPing_input_data/Compositions_information_table.tsv \
  --profiles tyPPing/tyPPing_input_data/Profile_information_table.tsv

# For draft genomes / MAGs
Rscript tyPPing/tyPPing_for_draft_genomes.R \
  -m protein_to_genome.tsv -s contig_sizes.tsv \
  -i hmm_output.tbl.out -o output_dir/ \
  --compositions tyPPing/tyPPing_input_data/Compositions_information_table.tsv \
  --profiles tyPPing/tyPPing_input_data/Profile_information_table.tsv
```

### Test with example data

```bash
Rscript tyPPing/tyPPing.R \
  -m tyPPing/test_example/protein_to_genome.tsv \
  -s tyPPing/test_example/genome_size.tsv \
  -i tyPPing/test_example/tyPPing_example_hmmsearch_output.tbl.out \
  -o tyPPing/test_example/ \
  --compositions tyPPing/tyPPing_input_data/Compositions_information_table.tsv \
  --profiles tyPPing/tyPPing_input_data/Profile_information_table.tsv
```

Expected outputs for comparison: `tyPPing/test_example/Final_prediction_table.tsv` and `All_hmm_hits_table.tsv`.

### MM-GRC

```bash
R < MM-GRC/wGRR_MGEs.R --no-save
```

## Architecture

### tyPPing algorithm (tyPPing.R / tyPPing_for_draft_genomes.R)

- Parses HMMER domtblout → filters by score thresholds from `Profile_information_table.tsv`
- Runs two parallel branches:
  - **MinProteins**: score-based prediction requiring minimum number of signature proteins
  - **Composition**: checks for complete functional composition sets from `Compositions_information_table.tsv`
- Merges branch results → assigns confidence levels (High/Medium/Low) → outputs `Final_prediction_table.tsv`
- Scripts support both CLI mode (`argparse`) and interactive mode (edit `USER INPUT SECTION` at top)
- Branch functions use `_f` suffix naming (e.g., `MinProteins_score_based_prediction_f`)

### Input data requirements

- `protein_to_genome.tsv`: maps `protein_id` → `genome_id` (use `prepare_input_tables.R` for Prodigal-style headers)
- `genome_size.tsv`: maps `genome_id` → `size`
- HMM profiles: 763 signature profiles (from repo zip or Zenodo)

### External data on Zenodo (doi:10.5281/zenodo.16616313)

HMM profiles (`tyPPing_signature_profiles.hmm`, `phage.hmm`), pre-trained random forest models, and figure reproduction data.

## Code Conventions

- R is the primary language; uses `data.table::fread()` for I/O combined with `tidyverse` for transformations
- Section dividers: `################################################################################`
- Auto-install logic via `tryCatch` for required R packages
- No formal test framework — validation via example datasets and supplementary data comparison
- `Publication_related_data/` contains figure reproduction scripts (not part of the tool itself)
