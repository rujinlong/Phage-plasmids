"""Convenience wrapper for running hmmsearch."""

import shutil
import subprocess
from pathlib import Path


def run_hmmsearch(
    protein_fasta: str,
    hmm_profiles: str,
    output_dir: str,
    threads: int = 4,
) -> Path:
    """Run hmmsearch and return path to domtblout output.

    Parameters
    ----------
    protein_fasta : str
        Path to protein FASTA file.
    hmm_profiles : str
        Path to HMM profile database (.hmm).
    output_dir : str
        Output directory.
    threads : int
        Number of CPU threads.

    Returns
    -------
    Path
        Path to the domtblout output file.
    """
    if shutil.which("hmmsearch") is None:
        raise RuntimeError("hmmsearch not found on PATH. Install HMMER first.")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    domtblout = out_dir / "hmmsearch_output.tbl.out"
    stdout_file = out_dir / "hmmsearch_output.out"

    cmd = [
        "hmmsearch",
        "--cpu", str(threads),
        "-o", str(stdout_file),
        "--domtblout", str(domtblout),
        str(hmm_profiles),
        str(protein_fasta),
    ]

    subprocess.run(cmd, check=True)
    return domtblout
