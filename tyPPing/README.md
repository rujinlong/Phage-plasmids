# Introduction

**tyPPing** is a Python CLI tool designed for fast and accurate phage-plasmid (P-P) detection and typing. It uses a curated set of protein profiles trained on well-known P-P types to find patterns of conserved signature proteins. If a given sequence passes the thresholds of these patterns (MinProteins and Composition) and fits the typical size range of the corresponding P-P type, it is reported as a P-P with a confidence score. This allows users to quickly identify P-Ps, making tyPPing a practical tool for large-scale P-P screening and systematic classification.

Currently, tyPPing can detect the following P-P types:

**AB_1, P1_1, P1_2, N15, SSU5_pHCM2, pMT1, pCAV, pSLy3, pKpn, and cp32.**

It works with both **complete** and **draft genomes** using a single unified command with automatic mode detection.

**tyPPing's workflow**

<img width="2462" height="893" alt="image" src="https://github.com/user-attachments/assets/8b6af34b-997f-4303-bac7-fb638db0b14b" />


**Step 1. Protein search**

Use **HMMER** (`hmmsearch`) to compare your multi-protein FASTA file against the provided P-P signature profiles. tyPPing includes a convenience wrapper for this step (`typping hmmsearch`).

**Step 2. Run tyPPing**

tyPPing processes the `--domtblout` output from the HMM search using two approaches:

-   **MinProteins** branch counts how many highly conserved P-P proteins are detected by the profiles using the sequence score as a cutoff. Sequences are kept if they have at least a minimum type-specific number (=threshold).
-   **Composition** branch keeps elements only if they match unique sets of P-P proteins that are defined for a given P-P type. Proteins are considered as detected, if they cover at least 50% of the profile. We use type-specific variations of the composition sizes (incomplete sets) by allowing gaps (number of not detected proteins).

P-P type and confidence (**High**, **Medium**, or **Low**) to the prediction are assigned by considering MinProteins, Composition and genome size criteria. The results are summarized in `Final_prediction_table.tsv`.

The **`tyPPing/`** folder contains the CLI tool, input tables, HMMs, and an example data set.

```
tyPPing/
├── pyproject.toml                        # Python package configuration
├── src/typping/                          # Python source code
│   ├── __init__.py
│   ├── cli.py                            # CLI entry point (run, prepare, hmmsearch)
│   ├── constants.py                      # Cutoff tables and category mappings
│   ├── io.py                             # HMMER parser and TSV readers
│   ├── prediction.py                     # MinProteins + Composition branches
│   ├── merge.py                          # Result merging and confidence assignment
│   ├── prepare.py                        # FASTA → input tables
│   └── hmmsearch.py                      # hmmsearch subprocess wrapper
├── tests/                                # Test suite
│   ├── conftest.py
│   ├── test_io.py
│   └── test_integration.py
├── tyPPing_signature_profiles.zip        # Compressed HMM profiles (unzip before use)
├── tyPPing_input_data/
│   ├── Compositions_information_table.tsv
│   └── Profile_information_table.tsv
├── test_example/                         # Example input/output for 20 genomes
└── Dockerfile
```

# Usage

## Setup & installation

### Option A: Install from source (recommended)

```bash
# Clone or download the repository, then:
cd tyPPing
pip install .
```

This installs the `typping` command-line tool. Python 3.10 or higher is required.

### Option B: Docker

```bash
cd tyPPing
docker build -t typping .

# Run with Docker (mount your data directory)
docker run --rm -v /path/to/your/data:/data typping run \
    -m /data/protein_to_genome.tsv \
    -s /data/genome_size.tsv \
    -i /data/hmm_output.tbl.out \
    -o /data/output/ \
    --compositions /opt/tyPPing/tyPPing_input_data/Compositions_information_table.tsv \
    --profiles /opt/tyPPing/tyPPing_input_data/Profile_information_table.tsv
```

### Requirements

-   **Python** >= 3.10
-   [**HMMER**](https://github.com/EddyRivasLab/hmmer) must be installed and accessible via the command line (used for hmmsearch). Only needed for Step 1; not required if you already have the HMMER output.

### Dependencies (installed automatically)

-   pandas >= 2.0
-   numpy >= 1.24
-   typer >= 0.9
-   biopython >= 1.80

## CLI overview

After installation, tyPPing provides three subcommands:

```
typping run         # Main prediction pipeline
typping prepare     # Generate input tables from FASTA files
typping hmmsearch   # Convenience wrapper for hmmsearch
```

Use `typping --help` or `typping <subcommand> --help` for detailed usage information.

## Input data requirements

**1. For HMM search (Step 1):**

-   A multi-FASTA protein file (e.g., `proteins.faa`).

-   P-P HMMs `tyPPing_signature_profiles.hmm` (also available on [Zenodo repository](https://doi.org/10.5281/zenodo.16616313)).

**2. For tyPPing prediction (Step 2):**

-   HMMER output table (`--domtblout`) from Step 1.

-   tyPPing reference tables provided in the `tyPPing_input_data/` folder.

-   A protein-to-genome table: see format details below.

-   A genome/contig size table: see format details below.

## Running tyPPing

### Step 1. Search for P-P specific proteins in target genomes

Run `hmmsearch` to scan protein sequences against the P-P type-specific profile HMMs.

**Option A: Direct hmmsearch command:**

```bash
hmmsearch -o tmp.out --domtblout hmm_output.tbl.out \
    tyPPing_signature_profiles.hmm proteins.faa
```

**Option B: Using the tyPPing wrapper:**

```bash
typping hmmsearch \
    -f proteins.faa \
    -o output_dir/ \
    --hmm tyPPing_signature_profiles.hmm \
    -t 4
```

#### `typping hmmsearch` arguments

| Argument | Short | Required | Description |
|----------|-------|----------|-------------|
| `--fasta` | `-f` | Yes | Path to the protein FASTA file |
| `--outdir` | `-o` | Yes | Output directory (created if it does not exist) |
| `--hmm` | | Yes | Path to HMM profile database (`.hmm`) |
| `--threads` | `-t` | No | Number of CPU threads (default: 4) |

**Output:** `output_dir/hmmsearch_output.tbl.out`

### Step 2. Prepare the input files

tyPPing requires two input tables in TSV format.

#### For complete genomes

-   **Protein-to-genome** mapping file (`protein_to_genome.tsv`):

    | protein_id | genome_id |
    |------------|-----------|
    | NP_052607  | NC_002128 |
    | NP_052608  | NC_002128 |
    | NP_052609  | NC_002128 |
    | ...        | ...       |

-   **Genome-sizes** file (`genome_size.tsv`):

    | genome_id | size  |
    |-----------|-------|
    | NC_002128 | 92721 |
    | NC_005856 | 94800 |
    | ...       | ...   |

#### For draft genomes / MAGs

-   **Protein-to-genome** mapping file (`protein_to_genome.tsv`) — requires an additional `contig_id` column:

    | protein_id       | contig_id      | genome_id |
    |------------------|----------------|-----------|
    | 170D8_contig_1_1 | 170D8_contig_1 | 170D8     |
    | 170D8_contig_1_2 | 170D8_contig_1 | 170D8     |
    | 170D8_contig_1_3 | 170D8_contig_1 | 170D8     |

-   **Contig-sizes** file (`contig_sizes.tsv`):

    | contig_id       | size    | genome_id |
    |-----------------|---------|-----------|
    | 170D8_contig_1  | 5091073 | 170D8     |
    | 170D8_contig_13 | 198075  | 170D8     |
    | 170D8_contig_14 | 34813   | 170D8     |

#### Automatic table generation with `typping prepare`

If your protein sequences follow the default Prodigal naming convention (`genome_1`, `genome_2`, ...), you can generate these tables automatically:

```bash
# Generate protein_to_genome.tsv from a protein FASTA
typping prepare \
    -f proteins.faa \
    -o output_dir/ \
    --mode complete          # or --mode draft

# Also generate genome_size.tsv from genome FASTA files
typping prepare \
    -f proteins.faa \
    -o output_dir/ \
    --genome-fasta-dir /path/to/genome_fastas/ \
    --mode complete
```

#### `typping prepare` arguments

| Argument | Short | Required | Description |
|----------|-------|----------|-------------|
| `--fasta` | `-f` | Yes | Protein FASTA file |
| `--outdir` | `-o` | Yes | Output directory |
| `--genome-fasta-dir` | | No | Directory containing genome FASTA files (for genome_size.tsv) |
| `--mode` | | No | `complete` (default) or `draft` |
| `--sep` | | No | Separator between genome/contig ID and protein number (default: `_`) |
| `--fasta-pattern` | | No | Glob pattern for genome FASTA files (default: `*.fasta`) |

**Output:**
-   `output_dir/protein_to_genome.tsv`
-   `output_dir/genome_size.tsv` (only if `--genome-fasta-dir` is provided)

### Step 3. Run tyPPing prediction

```bash
typping run \
    -m protein_to_genome.tsv \
    -s genome_size.tsv \
    -i hmm_output.tbl.out \
    -o output_dir/ \
    --compositions tyPPing_input_data/Compositions_information_table.tsv \
    --profiles tyPPing_input_data/Profile_information_table.tsv
```

#### `typping run` arguments

| Argument | Short | Required | Description |
|----------|-------|----------|-------------|
| `--map` | `-m` | Yes | Protein-to-genome mapping file (TSV) |
| `--sizes` | `-s` | Yes | Genome/contig size file (TSV) |
| `--hmm-domtbl` | `-i` | Yes | HMMER domtblout output file (`*.tbl.out`) |
| `--outdir` | `-o` | Yes | Output directory (created automatically if needed) |
| `--compositions` | | Yes | Compositions information table |
| `--profiles` | | Yes | Profile information table |
| `--mode` | | No | `auto` (default), `complete`, or `draft` |

#### Mode detection

-   **`--mode auto`** (default): automatically detects the mode from the protein-to-genome file. If a `contig_id` column is present, draft mode is used; otherwise, complete mode.
-   **`--mode complete`**: forces complete genome mode (2-column protein-to-genome file).
-   **`--mode draft`**: forces draft genome mode (3-column protein-to-genome file with `contig_id`).

### Quick test with example data

```bash
cd tyPPing

typping run \
    -m test_example/protein_to_genome.tsv \
    -s test_example/genome_size.tsv \
    -i test_example/tyPPing_example_hmmsearch_output.tbl.out \
    -o /tmp/typping_test/ \
    --compositions tyPPing_input_data/Compositions_information_table.tsv \
    --profiles tyPPing_input_data/Profile_information_table.tsv
```

Compare your results against the expected output:

```bash
diff /tmp/typping_test/Final_prediction_table.tsv test_example/Final_prediction_table.tsv
diff /tmp/typping_test/All_hmm_hits_table.tsv test_example/All_hmm_hits_table.tsv
```

Both files should be identical.

## Standard output summary tables

tyPPing generates **two output tables**:

**1. `Final_prediction_table.tsv`**

This table contains the **final predictions** for elements detected by at least one branch (**MinProteins** and/or **Composition**) and includes confidence levels for each assignment.

<img width="1146" height="173" alt="prediction" src="https://github.com/user-attachments/assets/d38527cf-3786-4dac-adfc-30d7e657388e" />

**Column descriptions (complete mode):**

| Column | Description |
|--------|-------------|
| **Genome ID** | Unique identifier of the analyzed genome |
| **Genome size (bp)** | Total size of the genome, in base pairs |
| **P-P type** | Predicted phage-plasmid type assigned to the genome |
| **Confidence level** | Prediction confidence: **High**, **Medium**, or **Low** |
| **Predicted by** | Branch(es) that led to the prediction |
| **MinProteins cutoff** | Minimum required protein hits for a match to this P-P type |
| **MinProteins hits** | Number of signature proteins detected by MinProteins |
| **Composition** | Composition pattern ID matched for this P-P type |
| **Composition size** | Total number of signature profiles in the composition set |
| **Composition hits** | Number of composition-specific profiles found |

**2. `All_hmm_hits_table.tsv`**

This table lists **all HMM hits** (not just confident predictions). It is useful for inspection of intermediate data or investigating borderline cases.

<img width="1160" height="125" alt="all_hits" src="https://github.com/user-attachments/assets/c26bc977-cc86-4ada-b94e-b58e73815397" />

| Column | Description |
|--------|-------------|
| **Genome ID** | Unique identifier of the analyzed genome |
| **P-P type** | P-P type for which the hits are analyzed (can be several entries per genome) |
| **MinProteins cutoff** | Minimum required protein hits for a match |
| **MinProteins hits** | Number of conserved proteins detected by MinProteins |
| **Signature to all genes ratio** | Ratio of detected proteins to total proteins |
| **P-P score sum** | Sum of P-P scores for all matched signature profiles |
| **P-P type score sum** | Sum of P-P type scores for all matched signature profiles |
| **Composition** | ID of the matched set of proteins |
| **Composition size** | Number of P-P proteins in the defined composition set |
| **Tolerance to gaps** | Maximum number of missing proteins in the composition set |
| **Composition hits** | Number of matched profiles |
| **Genome size (bp)** | Total genome size, in base pairs |
| **Number of proteins** | Total number of encoded proteins in the genome |

## tyPPing for incomplete genomes (including MAGs)

tyPPing's scoring system works on incomplete genomes by aggregating counts across contigs within a genome. Conserved P-P proteins are counted across all contigs in a genome and tested to see if they collectively meet the MinProteins, composition, and size requirements. With this change, tyPPing detects P-Ps that are split across multiple contigs with a high density of conserved P-P proteins.

To use draft mode, either provide a protein-to-genome file with 3 columns (`protein_id`, `contig_id`, `genome_id`) and let auto-detection handle it, or explicitly pass `--mode draft`:

```bash
typping run \
    -m protein_to_genome.tsv \
    -s contig_sizes.tsv \
    -i hmm_output.tbl.out \
    -o output_dir/ \
    --compositions tyPPing_input_data/Compositions_information_table.tsv \
    --profiles tyPPing_input_data/Profile_information_table.tsv \
    --mode draft
```

**Draft mode output tables** include additional columns:

| Column | Description |
|--------|-------------|
| **MinProteins hits list** | Counts of conserved proteins per contig, separated by ";" |
| **MinProteins contigs list** | Contig IDs encoding matched MinProteins proteins, separated by ";" |
| **Composition contigs list** | Contig IDs that match Composition, separated by ";" |
| **Genome size MinProteins (bp)** | Combined size of contigs contributing to MinProteins hits |
| **Genome size Composition (bp)** | Combined size of contigs contributing to Composition hits |

**Comments on draft mode:**

-   This approach cannot detect multiple P-Ps of the same type within one single genome. We expect these events to be rare, since P-Ps of the same type are suggested to be incompatible (like plasmids).

-   It can detect P-Ps of different types in the same genome.

-   For MAGs, accurate binning is crucial. Mis-binned contigs may result in incorrect or chimeric P-P predictions.

## Full workflow example

Here is a complete end-to-end example starting from protein and genome FASTA files:

```bash
# 1. Unzip HMM profiles (first time only)
cd tyPPing
unzip tyPPing_signature_profiles.zip

# 2. Install tyPPing
pip install .

# 3. (Optional) Generate input tables from FASTA files
typping prepare \
    -f /path/to/proteins.faa \
    -o /path/to/output/ \
    --genome-fasta-dir /path/to/genome_fastas/ \
    --mode complete

# 4. Run hmmsearch
typping hmmsearch \
    -f /path/to/proteins.faa \
    -o /path/to/output/ \
    --hmm tyPPing_signature_profiles.hmm \
    -t 8

# 5. Run tyPPing prediction
typping run \
    -m /path/to/output/protein_to_genome.tsv \
    -s /path/to/output/genome_size.tsv \
    -i /path/to/output/hmmsearch_output.tbl.out \
    -o /path/to/output/ \
    --compositions tyPPing_input_data/Compositions_information_table.tsv \
    --profiles tyPPing_input_data/Profile_information_table.tsv

# 6. View results
cat /path/to/output/Final_prediction_table.tsv
```

## Performance notes

-   **Detection accuracy:**
    -   \>99% sensitivity (especially for cp32 P-Ps) and \>99% precision compared to MM-GRC.
-   **Running time:**
    -   The tyPPing prediction step typically completes in seconds. The HMM protein-to-profile comparison using HMMER is the bottleneck (~1h44m for >38,000 plasmids).

## Interpretation and recommendations

-   Prior to your analysis, please test the workflow using the small example dataset in `test_example/`. It includes 20 genomes, the required tables and the output files.

-   **High confidence** predictions are detected by **both branches** (*MinProteins* and *Composition*) and fit the expected **size** range. These are the most reliable predictions but may miss ~20% of P-Ps if only this category is considered.

-   **Medium confidence** predictions are detected by one or both branches and fit the size range with 10% tolerance. These often represent more divergent P-Ps (such as many cp32-like). Rare false positives (three elements in 05/23 dataset) were detected with Medium confidence.

-   **Low confidence** is assigned to elements with genome sizes significantly shorter or longer than expected. Short genomes may indicate degrading P-Ps, while unusually long ones may result from recombination. Although tyPPing reports these elements, they should not be considered reliable predictions without further manual inspection.

## Running tests

```bash
cd tyPPing
pip install pytest
pytest tests/ -v
```

# Citing tyPPing

-   If you use tyPPing, please cite the corresponding paper: **["Efficient detection and typing of phage-plasmids"](https://www.biorxiv.org/content/10.1101/2025.08.29.673033v1)**
