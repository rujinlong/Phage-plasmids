"""Tests for io module."""

import pandas as pd

from typping.io import (
    parse_hmmer_domtblout,
    preprocess_hmm_output,
    read_compositions,
    read_genome_size,
    read_profile_info,
    read_protein_to_genome,
)


def test_parse_hmmer_domtblout(test_example_dir):
    df = parse_hmmer_domtblout(test_example_dir / "tyPPing_example_hmmsearch_output.tbl.out")
    assert not df.empty
    assert "domain_name" in df.columns
    assert "query_name" in df.columns
    assert df["sequence_score"].dtype == float


def test_read_protein_to_genome(test_example_dir):
    df = read_protein_to_genome(test_example_dir / "protein_to_genome.tsv")
    assert "protein_id" in df.columns
    assert "genome_id" in df.columns
    assert len(df) > 0


def test_read_genome_size(test_example_dir):
    df = read_genome_size(test_example_dir / "genome_size.tsv")
    assert "genome_id" in df.columns
    assert "size" in df.columns


def test_read_compositions(input_data_dir):
    df = read_compositions(input_data_dir / "Compositions_information_table.tsv")
    assert "PP_category" in df.columns
    assert "Composition" in df.columns
    assert "hmm_profile" in df.columns


def test_read_profile_info(input_data_dir):
    score, cutoff = read_profile_info(input_data_dir / "Profile_information_table.tsv")
    assert "hmm_profile" in score.columns
    assert "profile_scores_1" in score.columns
    assert "sequence_score_thr" in cutoff.columns


def test_preprocess_hmm_output(test_example_dir):
    hmm_raw = parse_hmmer_domtblout(test_example_dir / "tyPPing_example_hmmsearch_output.tbl.out")
    p2g = read_protein_to_genome(test_example_dir / "protein_to_genome.tsv")
    df = preprocess_hmm_output(hmm_raw, p2g, is_draft=False)
    assert "PP_category" in df.columns
    assert "cov_profile" in df.columns
    # Should have one hit per genome_id+hmm_profile
    assert df.groupby(["genome_id", "hmm_profile"]).size().max() == 1
