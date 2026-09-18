import io
import tempfile
import zipfile
from dataclasses import dataclass
from enum import StrEnum
from functools import cached_property
from pathlib import Path

import pandas as pd
import gffpandas.gffpandas as gffpd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class SignalPeptideType(StrEnum):
    SIGNAL_PEPTIDE = "Signal Peptide (Sec/SPI)"
    LIPOPROTEIN_SIGNAL_PEPTIDE = "Lipoprotein signal peptide (Sec/SPII)"


@dataclass(frozen=True)
class SignalPeptide:
    """Predicted signal peptide from SignalP 6.0.

    Attributes
    ----------
    protein_id
        ID of protein.
    peptide_type
        Type of signal peptide.
    start_index
        1-based residue index of start of signal peptide.
    end_index
        1-based residue index of end of signal peptide.
    probability
        Confidence of signal peptide prediction (in range 0-1).
    cleaved_sequence
        Amino acid sequence with signal peptide cleaved.
    """

    protein_id: str
    peptide_type: SignalPeptideType
    start_index: int
    end_index: int
    probability: float
    cleaved_sequence: Seq


class SignalPResults:
    """Handles processing SignalP 6.0 results."""

    _results_zip: Path

    _output_gff3_name = "output.gff3"
    _processed_entries_fasta_name = "processed_entries.fasta"

    @cached_property
    def output_df(self) -> pd.DataFrame:
        """Output results table."""

        with tempfile.TemporaryDirectory() as temp_dir:
            with zipfile.ZipFile(self._results_zip, "r") as zf:
                zf.extract(self._output_gff3_name, temp_dir)
                extracted_gff3 = Path(temp_dir) / self._output_gff3_name
                return gffpd.read_gff3(extracted_gff3).df

    @cached_property
    def processed_entries(self) -> list[SeqRecord]:
        """Sequences with signal peptides removed."""

        with zipfile.ZipFile(self._results_zip, "r") as zf:
            with zf.open(
                    self._processed_entries_fasta_name,
                    "r"
            ) as fasta:
                return list(
                    SeqIO.parse(
                        io.TextIOWrapper(
                            fasta,
                            encoding="utf-8"
                        ),
                        "fasta"
                    )
                )

    @cached_property
    def signal_peptides(self) -> list[SignalPeptide]:
        """Signal peptides found in results."""

        df = self.output_df
        signal_peptides = []
        for record in self.processed_entries:
            protein_id = record.description
            protein_df = df.loc[df["seq_id"] == protein_id]
            if len(protein_df) == 0:
                raise ValueError(
                    f"Could not find protein in output files: '{protein_id}'."
                )
            if len(protein_df) > 1:
                raise ValueError(
                    f"Multiple proteins in output files with ID: "
                    f"'{protein_id}'."
                )

            protein = protein_df.iloc[0]
            cleaved_sequence = record.seq
            if not isinstance(cleaved_sequence, Seq):
                raise TypeError(
                    f"Expected cleaved sequence to have type of `Seq`. "
                    f"Got: {type(cleaved_sequence)}: {cleaved_sequence}."
                )

            signal_peptides.append(SignalPeptide(
                protein_id=protein_id,
                peptide_type=SignalPeptideType(protein["type"]),
                start_index=protein["start"],
                end_index=protein["end"],
                probability=protein["score"],
                cleaved_sequence=cleaved_sequence,
            ))

        return signal_peptides

    @cached_property
    def protein_id_to_signal_peptide(self) -> dict[str, SignalPeptide]:
        """Mapping between protein IDs and their signal peptide."""

        protein_id_to_signal_peptides = {}
        for peptide in self.signal_peptides:
            if peptide.protein_id in protein_id_to_signal_peptides:
                raise KeyError(
                    f"Multiple signal peptides with protein ID "
                    f"'{peptide.protein_id}'."
                )

            protein_id_to_signal_peptides[peptide.protein_id] = peptide
        return protein_id_to_signal_peptides

    def __init__(self, results_zip: Path):
        """Create instance.

        Parameters
        ----------
        results_zip
            Path to zip archive containing SignalP results.
        """

        if not results_zip.is_file():
            raise FileNotFoundError(
                f"Could not find SignalP results zip file at '{results_zip}'."
            )

        if results_zip.suffix != ".zip":
            raise TypeError(f"Expected zip file, got '{results_zip}'.")

        self._results_zip = results_zip
