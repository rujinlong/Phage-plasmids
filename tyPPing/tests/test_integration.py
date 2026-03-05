"""Integration tests: validate Python output matches R expected output."""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from tests.conftest import INPUT_DATA_DIR, TEST_EXAMPLE_DIR


def test_complete_mode_matches_r_output():
    """Run tyPPing in complete mode on test data and compare with expected R output."""
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            "typping", "run",
            "-m", str(TEST_EXAMPLE_DIR / "protein_to_genome.tsv"),
            "-s", str(TEST_EXAMPLE_DIR / "genome_size.tsv"),
            "-i", str(TEST_EXAMPLE_DIR / "tyPPing_example_hmmsearch_output.tbl.out"),
            "-o", tmpdir,
            "--compositions", str(INPUT_DATA_DIR / "Compositions_information_table.tsv"),
            "--profiles", str(INPUT_DATA_DIR / "Profile_information_table.tsv"),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"typping run failed: {result.stderr}"

        # Compare Final_prediction_table.tsv
        expected_pred = pd.read_csv(TEST_EXAMPLE_DIR / "Final_prediction_table.tsv", sep="\t", dtype=str)
        actual_pred = pd.read_csv(Path(tmpdir) / "Final_prediction_table.tsv", sep="\t", dtype=str)

        # Sort both for comparison
        sort_cols = ["Genome ID", "P-P type"]
        expected_pred = expected_pred.sort_values(sort_cols).reset_index(drop=True)
        actual_pred = actual_pred.sort_values(sort_cols).reset_index(drop=True)

        pd.testing.assert_frame_equal(actual_pred, expected_pred)

        # Compare All_hmm_hits_table.tsv
        expected_hits = pd.read_csv(TEST_EXAMPLE_DIR / "All_hmm_hits_table.tsv", sep="\t", dtype=str)
        actual_hits = pd.read_csv(Path(tmpdir) / "All_hmm_hits_table.tsv", sep="\t", dtype=str)

        sort_cols_hits = ["Genome ID", "P-P type"]
        expected_hits = expected_hits.sort_values(sort_cols_hits).reset_index(drop=True)
        actual_hits = actual_hits.sort_values(sort_cols_hits).reset_index(drop=True)

        pd.testing.assert_frame_equal(actual_hits, expected_hits)


def test_final_prediction_genomes():
    """Verify specific genome predictions match expected results."""
    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            "typping", "run",
            "-m", str(TEST_EXAMPLE_DIR / "protein_to_genome.tsv"),
            "-s", str(TEST_EXAMPLE_DIR / "genome_size.tsv"),
            "-i", str(TEST_EXAMPLE_DIR / "tyPPing_example_hmmsearch_output.tbl.out"),
            "-o", tmpdir,
            "--compositions", str(INPUT_DATA_DIR / "Compositions_information_table.tsv"),
            "--profiles", str(INPUT_DATA_DIR / "Profile_information_table.tsv"),
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)

        pred = pd.read_csv(Path(tmpdir) / "Final_prediction_table.tsv", sep="\t")

        # Check known predictions
        ab1 = pred[pred["Genome ID"] == "NZ_CP107600"]
        assert len(ab1) == 1
        assert ab1.iloc[0]["P-P type"] == "AB_1"
        assert ab1.iloc[0]["Confidence level"] == "High"

        n15 = pred[pred["Genome ID"] == "NZ_CP064249"]
        assert len(n15) == 1
        assert n15.iloc[0]["P-P type"] == "N15"

        # cp32 with Composition-only prediction
        cp32_comp = pred[pred["Genome ID"] == "NZ_CP115597"]
        assert len(cp32_comp) == 1
        assert cp32_comp.iloc[0]["Predicted by"] == "Composition"
        assert cp32_comp.iloc[0]["Confidence level"] == "Medium"
