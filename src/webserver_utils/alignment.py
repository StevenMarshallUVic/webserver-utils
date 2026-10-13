"""Clustal sequence alignments."""

import logging
import shutil
import subprocess
from collections.abc import Collection
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(Path(__file__).name)


@dataclass(frozen=True)
class ClustalOmegaAligner:
    log_version: bool = True

    _PROGRAM_NAME = "clustalo"

    @property
    def is_installed(self):
        return shutil.which(self._PROGRAM_NAME) is not None

    @property
    def version(self):
        result = subprocess.run(
            [self._PROGRAM_NAME, "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()

    def __post_init__(self):
        if not self.is_installed:
            raise FileNotFoundError(
                f"Could not find {self._PROGRAM_NAME} in system path. "
                f"Are you sure it is installed?"
            )

    def write_clustal_omega_alignment(
            self,
            input_fasta: Path,
            output_aln: Path,
            force: bool = False,
    ) -> None:
        """

        Parameters
        ----------
        input_fasta
        output_aln
        force
            (Optional) When True, overwrites file at `output_aln` if it exists.

        Returns
        -------

        """

        if not input_fasta.is_file():
            raise FileNotFoundError(
                f"Could not find input FASTA at '{input_fasta}'."
            )
        if self.log_version:
            logger.info(
                "Performing alignment using Clustal Omega %s",
                self.version,
            )

        cmd = [
            self._PROGRAM_NAME,
            "-i", input_fasta,
            "-o", output_aln,
            "--outfmt=clu",
        ]
        if force:
            cmd.append("--force")
            if output_aln.is_file():
                logger.warning(
                    "Alignment already exists at %s, overwriting...",
                    output_aln
                )
        try:
            subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as e:
            e.add_note(f"Error running Clustal Omega: {e.stderr.decode()}")
            raise e

    def clustal_omega_align_fasta(
            self,
            input_fasta: Path,
    ) -> MultipleSeqAlignment:
        with TemporaryDirectory() as tmp_dir:
            alignment = Path(tmp_dir) / "alignment.afa"
            self.write_clustal_omega_alignment(
                input_fasta,
                alignment,
            )
            return AlignIO.read(alignment, "clustal")

    def clustal_omega_align_records(
            self,
            records: Collection[SeqRecord],
    ) -> MultipleSeqAlignment:
        if len(records) <= 1:
            raise ValueError(
                f"Requires at least two records to perform alignment, "
                f"got {len(records)}."
            )

        with TemporaryDirectory() as tmp_dir:
            fasta = Path(tmp_dir) / "input.fasta"
            SeqIO.write(records, fasta, "fasta-2line")
            return self.clustal_omega_align_fasta(fasta)
