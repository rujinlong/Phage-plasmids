"""Generate protein_to_genome.tsv and genome_size.tsv from FASTA files."""

import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def prepare_protein_to_genome(
    protein_fasta: str,
    output_dir: str,
    sep: str = "_",
    mode: str = "complete",
) -> Path:
    """Generate protein_to_genome.tsv from a protein FASTA file.

    Parameters
    ----------
    protein_fasta : str
        Path to multi-FASTA protein file (e.g. from Prodigal).
    output_dir : str
        Output directory.
    sep : str
        Separator between genome/contig ID and protein number.
    mode : str
        "complete" or "draft". Draft mode adds contig_id column.

    Returns
    -------
    Path
        Path to the output file.
    """
    records = []
    pattern = re.compile(rf"^(.+?){re.escape(sep)}(\d+)$")

    for record in SeqIO.parse(protein_fasta, "fasta"):
        protein_id = record.id.split("\t")[0]
        match = pattern.match(protein_id)
        if match:
            genome_or_contig_id = match.group(1)
        else:
            genome_or_contig_id = protein_id.rsplit(sep, 1)[0] if sep in protein_id else protein_id

        if mode == "draft":
            records.append({"protein_id": protein_id, "contig_id": genome_or_contig_id, "genome_id": genome_or_contig_id})
        else:
            records.append({"protein_id": protein_id, "genome_id": genome_or_contig_id})

    df = pd.DataFrame(records)
    out_path = Path(output_dir) / "protein_to_genome.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


def prepare_genome_size(
    genome_fasta_dir: str,
    output_dir: str,
    pattern: str = "*.fasta",
) -> Path:
    """Generate genome_size.tsv from genome FASTA files.

    Parameters
    ----------
    genome_fasta_dir : str
        Directory containing genome FASTA files.
    output_dir : str
        Output directory.
    pattern : str
        Glob pattern for FASTA files.

    Returns
    -------
    Path
        Path to the output file.
    """
    fasta_dir = Path(genome_fasta_dir)
    records = []

    for fasta_file in sorted(fasta_dir.glob(pattern)):
        total_length = sum(len(rec.seq) for rec in SeqIO.parse(str(fasta_file), "fasta"))
        genome_id = fasta_file.stem
        records.append({"genome_id": genome_id, "size": total_length})

    df = pd.DataFrame(records)
    out_path = Path(output_dir) / "genome_size.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    return out_path
