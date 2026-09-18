import tempfile
import zipfile
from dataclasses import dataclass
from enum import StrEnum
from functools import cached_property
from pathlib import Path

import pandas as pd
import gffpandas.gffpandas as gffpd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class SignalPeptideType(StrEnum):
    # noinspection PyTypeChecker
    SIGNAL_PEPTIDE = ("signal_peptide", "Signal Peptide (Sec/SPI)")
    # noinspection PyTypeChecker
    LIPOPROTEIN_SIGNAL_PEPTIDE = (
        "lipoprotein_signal_peptide",
        "Lipoprotein signal peptide (Sec/SPII)"
    )

    def __new__(cls, value: str, display_name: str):
        obj = str.__new__(cls, value)
        obj._value_ = value
        obj.display_name = display_name
        return obj


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

    @classmethod
    def from_signalp_output_df_row(cls, row: pd.Series, **kwargs):
        return cls(
            protein_id=row["seq_id"],
            peptide_type=SignalPeptideType(row["type"]),
            start_index=row["start"],
            end_index=row["end"],
            probability=row["score"],
            **kwargs
        )


@dataclass(frozen=True)
class SignalPResults:
    """Handles processing SignalP 6.0 results.

    Attributes
    ----------
    output_df
        Output results table.
    processed_entries
        Sequences with signal peptides removed.
    """

    output_df: pd.DataFrame
    processed_entries: tuple[SeqRecord, ...]

    @cached_property
    def signal_peptides(self) -> list[SignalPeptide]:
        """Signal peptides found in results."""

        signal_peptides = []
        for record in self.processed_entries:
            protein_id = record.description
            protein_df = self.output_df.loc[
                self.output_df["seq_id"] == protein_id
            ]
            if len(protein_df) == 0:
                raise ValueError(
                    f"Could not find protein in output files: '{protein_id}'."
                )
            if len(protein_df) > 1:
                raise ValueError(
                    f"Multiple rows in output files with ID: "
                    f"'{protein_id}'."
                )

            cleaved_sequence = record.seq
            if not isinstance(cleaved_sequence, Seq):
                raise TypeError(
                    f"Expected cleaved sequence to have type of `Seq`. "
                    f"Got: {type(cleaved_sequence)}: {cleaved_sequence}."
                )

            signal_peptides.append(
                SignalPeptide.from_signalp_output_df_row(
                    row=protein_df.iloc[0],
                    cleaved_sequence=cleaved_sequence,
                )
            )

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

    @classmethod
    def from_signalp_results_files(
            cls,
            output_gff3: Path,
            processed_entries_fasta: Path,
    ):
        if not output_gff3.is_file():
            raise FileNotFoundError(
                f"Could not find output GFF3 file at '{output_gff3}'."
            )
        if not processed_entries_fasta.is_file():
            raise FileNotFoundError(
                f"Could not find processed entries FASTA at "
                f"'{processed_entries_fasta}'."
            )

        return cls(
            output_df=gffpd.read_gff3(output_gff3).df,
            processed_entries=tuple(list(
                SeqIO.parse(processed_entries_fasta, "fasta")
            ))
        )

    @classmethod
    def from_signalp_results_zip(
            cls,
            results_zip: Path,
            output_gff3_name="output.gff3",
            processed_entries_fasta_name="processed_entries.fasta",
    ):
        """Create instance.

        Parameters
        ----------
        results_zip
            Path to zip archive containing SignalP results.
        output_gff3_name
            Name of output gff3 file found in SignalP results zip.
        processed_entries_fasta_name
            Name of preprocessed entries FASTA found in SignalP results zip.
        """

        if not results_zip.is_file():
            raise FileNotFoundError(
                f"Could not find SignalP results zip file at '{results_zip}'."
            )

        if results_zip.suffix != ".zip":
            raise TypeError(f"Expected zip file, got '{results_zip}'.")


        with tempfile.TemporaryDirectory() as temp_dir:
            with zipfile.ZipFile(results_zip, "r") as zf:
                zf.extract(output_gff3_name, temp_dir)
                zf.extract(processed_entries_fasta_name, temp_dir)

                extracted_gff3 = Path(temp_dir) / output_gff3_name
                extracted_fasta = Path(temp_dir) / processed_entries_fasta_name
                return cls.from_signalp_results_files(
                    extracted_gff3,
                    extracted_fasta
                )
