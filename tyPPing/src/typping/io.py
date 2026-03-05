"""I/O functions: parse HMMER domtblout and read TSV input files."""

import re

import numpy as np
import pandas as pd

from .constants import HMM_DOMTBL_COLUMNS, map_pptype_to_category


def parse_hmmer_domtblout(path: str) -> pd.DataFrame:
    """Parse HMMER domtblout file.

    The domtblout format has 22 fixed whitespace-delimited fields followed by
    a free-text description that may contain spaces. We split only the first
    22 fields and join the remainder as description.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(None, 22)
            if len(parts) < 22:
                continue
            # If no description field, add empty string
            if len(parts) == 22:
                parts.append("")
            rows.append(parts)

    df = pd.DataFrame(rows, columns=HMM_DOMTBL_COLUMNS)

    # Convert numeric columns
    numeric_cols = [
        "domain_len", "qlen",
        "sequence_evalue", "sequence_score", "sequence_bias",
        "domain_N", "domain_of",
        "domain_cevalue", "domain_ievalue", "domain_score", "domain_bias",
        "hmm_from", "hmm_to", "ali_from", "ali_to",
        "env_from", "env_to",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col])
    df["acc"] = pd.to_numeric(df["acc"])

    return df


def preprocess_hmm_output(hmm_raw: pd.DataFrame, protein_to_genome: pd.DataFrame,
                          is_draft: bool = False) -> pd.DataFrame:
    """Preprocess HMM search output: compute coverage, keep best hit per genome+profile.

    Parameters
    ----------
    hmm_raw : pd.DataFrame
        Raw parsed HMMER domtblout.
    protein_to_genome : pd.DataFrame
        Mapping of protein_id to genome_id (and contig_id for draft mode).
    is_draft : bool
        Whether this is draft genome mode.

    Returns
    -------
    pd.DataFrame
        Filtered HMM hits with PP_category assigned.
    """
    df = hmm_raw.copy()

    # Compute alignment coverage of the profile
    df["ali_length"] = (df["ali_to"] - df["ali_from"] + 1).abs()
    df["cov_profile"] = np.round(df["ali_length"] / df["qlen"], 3)

    # Select and rename columns
    df = df.rename(columns={"domain_name": "protein_id", "query_name": "hmm_profile"})
    df = df[["protein_id", "hmm_profile", "domain_ievalue", "cov_profile", "sequence_score"]]

    # Join with protein-to-genome mapping
    if is_draft:
        merge_cols = ["protein_id", "contig_id", "genome_id"]
    else:
        merge_cols = ["protein_id", "genome_id"]

    df = df.merge(protein_to_genome[merge_cols], on="protein_id", how="inner")

    # Keep only the best hit (lowest domain_ievalue) per genome_id + hmm_profile
    df = df.sort_values("domain_ievalue")
    df = df.groupby(["genome_id", "hmm_profile"], sort=False).first().reset_index()

    # Extract PP_category from hmm_profile name (everything before "_pers_")
    df["PP_category"] = df["hmm_profile"].str.extract(r"^(.+?)_pers_", expand=False)

    return df


def read_protein_to_genome(path: str) -> pd.DataFrame:
    """Read protein_to_genome TSV file.

    Returns DataFrame with columns: protein_id, genome_id (and contig_id for draft).
    """
    df = pd.read_csv(path, sep="\t")
    return df


def read_genome_size(path: str) -> pd.DataFrame:
    """Read genome_size TSV file. Returns DataFrame with columns: genome_id, size."""
    df = pd.read_csv(path, sep="\t")
    return df


def read_compositions(path: str) -> pd.DataFrame:
    """Read Compositions_information_table.tsv.

    Returns DataFrame with columns: PP_category, Composition, hmm_profile, composition_size.
    """
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "Composition ID": "Composition",
        "HMM profile in composition": "hmm_profile",
        "Composition size": "composition_size",
    })
    df["PP_category"] = df["P-P type"].map(map_pptype_to_category)
    return df


def read_profile_info(path: str):
    """Read Profile_information_table.tsv.

    Returns
    -------
    minproteins_score : pd.DataFrame
        Columns: hmm_profile, profile_scores_1, profile_scores_2, PP_category
    minproteins_cutoff : pd.DataFrame
        Columns: hmm_profile, sequence_score_thr
    """
    df = pd.read_csv(path, sep="\t")

    minproteins_score = df[["HMM profile ID", "P-P score", "P-P type score", "P-P type"]].copy()
    minproteins_score = minproteins_score.rename(columns={
        "HMM profile ID": "hmm_profile",
        "P-P score": "profile_scores_1",
        "P-P type score": "profile_scores_2",
    })
    minproteins_score["PP_category"] = minproteins_score["P-P type"].map(map_pptype_to_category)

    minproteins_cutoff = df[["HMM profile ID", "Sequence score threshold"]].copy()
    minproteins_cutoff = minproteins_cutoff.rename(columns={
        "HMM profile ID": "hmm_profile",
        "Sequence score threshold": "sequence_score_thr",
    })

    return minproteins_score, minproteins_cutoff
