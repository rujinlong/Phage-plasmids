"""tyPPing CLI: phage-plasmid detection and typing."""

import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional

import typer

app = typer.Typer(
    name="typping",
    help="tyPPing: Bioinformatic tool for phage-plasmid (P-P) detection and typing.",
    add_completion=False,
)


def _predict_category(category, hmm_filtered, is_draft, minproteins_score,
                       minproteins_cutoff, compositions_df, size_df,
                       contig_to_genome=None):
    """Process a single PP category (worker function for parallel execution).

    Must be a top-level function to be picklable by ProcessPoolExecutor.
    """
    from .prediction import (
        composition_branch_complete,
        composition_branch_draft,
        minproteins_branch_complete,
        minproteins_branch_draft,
    )

    # SSU5 → SSU5_pHCM2
    if category == "SSU5":
        hmm_filtered = hmm_filtered.copy()
        hmm_filtered["PP_category"] = "SSU5_pHCM2"
        category = "SSU5_pHCM2"

    if is_draft:
        mp_result = minproteins_branch_draft(
            category, hmm_filtered, minproteins_score, minproteins_cutoff,
            size_df, contig_to_genome,
        )
        comp_result = composition_branch_draft(
            category, hmm_filtered, compositions_df, size_df,
        )
    else:
        mp_result = minproteins_branch_complete(
            category, hmm_filtered, minproteins_score, minproteins_cutoff, size_df,
        )
        comp_result = composition_branch_complete(
            category, hmm_filtered, compositions_df, size_df,
        )

    return category, mp_result, comp_result


@app.command()
def run(
    map_file: str = typer.Option(..., "-m", "--map", help="Protein-to-genome mapping file (TSV)"),
    sizes_file: str = typer.Option(..., "-s", "--sizes", help="Genome/contig size file (TSV)"),
    hmm_domtbl: str = typer.Option(..., "-i", "--hmm-domtbl", help="HMMER domtblout output file"),
    outdir: str = typer.Option(..., "-o", "--outdir", help="Output directory"),
    compositions: str = typer.Option(..., "--compositions", help="Compositions information table"),
    profiles: str = typer.Option(..., "--profiles", help="Profile information table"),
    mode: str = typer.Option("auto", "--mode", help="Mode: auto, complete, or draft"),
    threads: int = typer.Option(1, "-t", "--threads", help="Number of parallel workers for P-P type prediction (default: 1)"),
):
    """Run tyPPing prediction pipeline."""
    import numpy as np
    import pandas as pd

    from .constants import map_category_to_pptype
    from .io import (
        parse_hmmer_domtblout,
        preprocess_hmm_output,
        read_compositions,
        read_genome_size,
        read_profile_info,
        read_protein_to_genome,
    )
    from .merge import merge_predictions_complete, merge_predictions_draft

    start_time = time.time()

    typer.echo("\n-------------------  Welcome to tyPPing  ----------------------")
    typer.echo("-----------------------------START-----------------------------\n")

    # Validate input files
    for label, path in [("map", map_file), ("sizes", sizes_file), ("hmm_domtbl", hmm_domtbl),
                         ("compositions", compositions), ("profiles", profiles)]:
        if not Path(path).is_file():
            typer.echo(f"Error: {label} file not found: {path}", err=True)
            raise typer.Exit(1)

    # Create output directory
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # Detect mode from input file
    if mode == "auto":
        protein_to_genome = read_protein_to_genome(map_file)
        is_draft = "contig_id" in protein_to_genome.columns
    else:
        is_draft = mode == "draft"
        protein_to_genome = read_protein_to_genome(map_file)

    typer.echo(f"Mode: {'draft' if is_draft else 'complete'}")
    typer.echo(f"Protein to genome path: {map_file}")
    typer.echo(f"Genome size path:       {sizes_file}")
    typer.echo(f"HMM search output path: {hmm_domtbl}")
    typer.echo(f"Output directory path:  {outdir}")

    typer.echo("Reading and processing the input data...")

    # Read inputs
    genome_size = read_genome_size(sizes_file)
    compositions_df = read_compositions(compositions)
    minproteins_score, minproteins_cutoff = read_profile_info(profiles)

    if is_draft:
        # Draft mode: protein_to_genome has contig_id column
        n_protein_per_contig = protein_to_genome.groupby("contig_id").size().reset_index(name="n_protein")
        contig_size = genome_size.rename(columns={"genome_id": "contig_id"} if "genome_id" in genome_size.columns else {})
        if "contig_id" not in contig_size.columns:
            contig_size = contig_size.rename(columns={contig_size.columns[0]: "contig_id"})

        # Validate: check for duplicate column names after rename
        dup_cols = contig_size.columns[contig_size.columns.duplicated()].unique().tolist()
        if dup_cols:
            typer.echo(
                f"Error: Duplicate column names found in genome size file: {dup_cols}\n"
                f"Columns in file: {contig_size.columns.tolist()}\n"
                f"Please check your genome size file: {sizes_file}",
                err=True,
            )
            raise typer.Exit(1)

        # Validate: check for duplicate contig_id values
        dup_contigs = contig_size[contig_size["contig_id"].duplicated(keep=False)]
        if not dup_contigs.empty:
            n_dup = dup_contigs["contig_id"].nunique()
            examples = dup_contigs["contig_id"].unique()[:10].tolist()
            typer.echo(
                f"Error: Found {n_dup} duplicated contig_id(s) in genome size file: {sizes_file}\n"
                f"Examples: {examples}\n"
                f"Each contig_id must appear only once. Please deduplicate your input.",
                err=True,
            )
            raise typer.Exit(1)

        contig_size = contig_size.merge(n_protein_per_contig, on="contig_id", how="left").fillna(0)
        contig_size["n_protein"] = contig_size["n_protein"].astype(int)

        # Filter contigs: size <= 300000 and n_protein >= 3
        contig_size_filtered = contig_size[(contig_size["size"] <= 300000) & (contig_size["n_protein"] >= 3)]
        protein_to_genome = protein_to_genome.merge(contig_size_filtered[["contig_id"]], on="contig_id", how="inner")
        contig_to_genome = protein_to_genome[["contig_id", "genome_id"]].drop_duplicates()
    else:
        # Complete mode: compute n_protein per genome
        n_protein_per_genome = protein_to_genome.groupby("genome_id").size().reset_index(name="n_protein")
        genome_size = genome_size.merge(n_protein_per_genome, on="genome_id", how="left").fillna(0)
        genome_size["n_protein"] = genome_size["n_protein"].astype(int)

    typer.echo("Reading and processing the protein-to-profile HMM output file...")

    hmm_raw = parse_hmmer_domtblout(hmm_domtbl)
    hmm_output = preprocess_hmm_output(hmm_raw, protein_to_genome, is_draft=is_draft)

    typer.echo("Protein-to-profile HMM output file processed successfully")

    # Run prediction for each PP category
    typer.echo("P-P prediction for each P-P type has started...")

    pp_categories = hmm_output["PP_category"].dropna().unique()
    n_workers = min(threads, len(pp_categories))

    # Prepare shared arguments
    size_df = contig_size if is_draft else genome_size
    ctg_to_genome = contig_to_genome if is_draft else None

    mp_results = []
    comp_results = []

    if n_workers > 1:
        typer.echo(f"Using {n_workers} parallel workers for {len(pp_categories)} P-P types...")
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = []
            for category in pp_categories:
                hmm_filtered = hmm_output[hmm_output["PP_category"] == category].copy()
                future = executor.submit(
                    _predict_category, category, hmm_filtered, is_draft,
                    minproteins_score, minproteins_cutoff, compositions_df,
                    size_df, ctg_to_genome,
                )
                futures.append((category, future))

            for orig_category, future in futures:
                cat, mp_result, comp_result = future.result()
                typer.echo(f"Searching for  {map_category_to_pptype(orig_category)}")
                mp_results.append(mp_result)
                comp_results.append(comp_result)
    else:
        for category in pp_categories:
            mapped = map_category_to_pptype(category)
            typer.echo(f"Searching for  {mapped}")
            hmm_filtered = hmm_output[hmm_output["PP_category"] == category].copy()
            _, mp_result, comp_result = _predict_category(
                category, hmm_filtered, is_draft,
                minproteins_score, minproteins_cutoff, compositions_df,
                size_df, ctg_to_genome,
            )
            mp_results.append(mp_result)
            comp_results.append(comp_result)

    typer.echo("P-P prediction by P-P type finished successfully")

    # Merge results
    typer.echo("Summarizing the tyPPing predictions...")

    mp_all = pd.concat(mp_results, ignore_index=True) if mp_results else pd.DataFrame()
    comp_all = pd.concat(comp_results, ignore_index=True) if comp_results else pd.DataFrame()

    if is_draft:
        all_hits, final_pred = merge_predictions_draft(mp_all, comp_all)
    else:
        all_hits, final_pred = merge_predictions_complete(mp_all, comp_all, genome_size)

    # Write output
    typer.echo("Writing the results...")

    summary_path = Path(outdir) / "All_hmm_hits_table.tsv"
    prediction_path = Path(outdir) / "Final_prediction_table.tsv"

    # Format: convert numeric columns to match R output format, fill NaN with NA
    def _format_output(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        # Convert size/count columns to nullable int
        int_cols = [c for c in df.columns if any(k in c.lower() for k in
                    ["size", "proteins", "cutoff", "hits", "tolerance"])]
        for col in int_cols:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: str(int(x)) if pd.notna(x) and x == x else "NA")
        # Format score columns: remove trailing .0 (R prints 1 not 1.0)
        score_cols = [c for c in df.columns if "score" in c.lower() or "ratio" in c.lower()]
        for col in score_cols:
            if col in df.columns:
                def _fmt_score(x):
                    if pd.isna(x):
                        return "NA"
                    s = str(x)
                    if s.endswith(".0"):
                        s = s[:-2]
                    return s
                df[col] = df[col].apply(_fmt_score)
        # Fill remaining NaN with NA
        df = df.fillna("NA")
        return df

    all_hits = _format_output(all_hits)
    final_pred = _format_output(final_pred)

    all_hits.to_csv(summary_path, sep="\t", index=False)
    final_pred.to_csv(prediction_path, sep="\t", index=False)

    typer.echo(f"Results successfully written to:")
    typer.echo(f"- {summary_path}")
    typer.echo(f"- {prediction_path}")

    elapsed = time.time() - start_time
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    typer.echo(f"Running time:  {minutes} minutes {seconds} seconds")

    typer.echo("\n------------------------------END------------------------------")
    typer.echo("---------------------------------------------------------------\n")


@app.command()
def prepare(
    protein_fasta: str = typer.Option(..., "-f", "--fasta", help="Protein FASTA file"),
    outdir: str = typer.Option(..., "-o", "--outdir", help="Output directory"),
    genome_fasta_dir: Optional[str] = typer.Option(None, "--genome-fasta-dir", help="Directory with genome FASTA files"),
    mode: str = typer.Option("complete", "--mode", help="Mode: complete or draft"),
    sep: str = typer.Option("_", "--sep", help="Separator between genome/contig ID and protein number"),
    fasta_pattern: str = typer.Option("*.fasta", "--fasta-pattern", help="Glob pattern for genome FASTA files"),
):
    """Generate protein_to_genome.tsv and genome_size.tsv from FASTA files."""
    from .prepare import prepare_genome_size, prepare_protein_to_genome

    Path(outdir).mkdir(parents=True, exist_ok=True)

    typer.echo("Generating protein_to_genome.tsv...")
    p2g_path = prepare_protein_to_genome(protein_fasta, outdir, sep=sep, mode=mode)
    typer.echo(f"  Written: {p2g_path}")

    if genome_fasta_dir:
        typer.echo("Generating genome_size.tsv...")
        gs_path = prepare_genome_size(genome_fasta_dir, outdir, pattern=fasta_pattern)
        typer.echo(f"  Written: {gs_path}")

    typer.echo("Done.")


@app.command()
def hmmsearch(
    protein_fasta: str = typer.Option(..., "-f", "--fasta", help="Protein FASTA file"),
    outdir: str = typer.Option(..., "-o", "--outdir", help="Output directory"),
    hmm: str = typer.Option(..., "--hmm", help="HMM profile database (.hmm)"),
    threads: int = typer.Option(4, "-t", "--threads", help="Number of CPU threads"),
):
    """Run hmmsearch with tyPPing HMM profiles."""
    from .hmmsearch import run_hmmsearch

    typer.echo(f"Running hmmsearch with {threads} threads...")
    out_path = run_hmmsearch(protein_fasta, hmm, outdir, threads=threads)
    typer.echo(f"Output written to: {out_path}")


if __name__ == "__main__":
    app()
