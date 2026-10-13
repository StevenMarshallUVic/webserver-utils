"""Tests for alignments."""
from subprocess import CalledProcessError

import pytest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from webserver_utils.alignment import ClustalOmegaAligner

test_records: tuple[SeqRecord, ...] = tuple([
    SeqRecord(
        id="TP0021_var0001",
        name="",
        description="",
        seq=Seq("MDDARYAEWSASLVQLPDTHFFDLMRLYLGVLKTPFHKQRLVQQLSAFLQRKSIQNAVVQMLDELDLLFISVVMCVPRATLELLTIFF"),
    ),
    SeqRecord(
        id="TP0021_var0002",
        name="",
        description="",
        seq=Seq("MDDARYAEWSASLVQLPDTHFFDLMRLYLGVLKTPFHKQRLVQQLSAFLQRKSIQNAVVQMLDELDLLFISVVMCVPHATLELLTIFF"),
    ),
])


def test_clustal_omega_no_records():
    with pytest.raises(ValueError):
        ClustalOmegaAligner().align_records([])


def test_clustal_omega_one_record():
    with pytest.raises(ValueError):
        ClustalOmegaAligner().align_records([test_records[0]])


def test_clustal_omega_non_existent_fasta():
    with pytest.raises(FileNotFoundError):
        ClustalOmegaAligner().write_alignment(Path(), Path())


def test_clustal_omega_alignment_exists_no_force(tmp_path):
    fasta = tmp_path / "test.fasta"
    SeqIO.write(test_records, fasta, "fasta")

    alignment = tmp_path / "alignment.afa"
    ClustalOmegaAligner().write_alignment(fasta, alignment)
    with pytest.raises(
            CalledProcessError,
            match="FATAL: Cowardly refusing to overwrite already existing file"
    ):
        ClustalOmegaAligner().write_alignment(fasta, alignment)


def test_clustal_omega_write_alignment(tmp_path):
    fasta = tmp_path / "test.fasta"
    SeqIO.write(test_records, fasta, "fasta")

    alignment = tmp_path / "alignment.afa"
    ClustalOmegaAligner().write_alignment(fasta, alignment)

    assert alignment.is_file()
    assert len(list(SeqIO.parse(alignment, "clustal"))) == len(test_records)


def test_clustal_omega_align_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"
    SeqIO.write(test_records, fasta, "fasta")

    alignment = ClustalOmegaAligner().align_fasta(fasta)

    assert len(alignment) == len(test_records)


def test_clustal_omega_align_records():
    alignment = ClustalOmegaAligner().align_records(test_records)

    assert len(alignment) == len(test_records)
