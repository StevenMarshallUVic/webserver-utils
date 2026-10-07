"""Tests for BLAST.

TODO: Add tests for blastn.
TODO: Add tests for blastx.
TODO: Add tests for tblastn.
TODO: Add tests for tblastx.
"""

import pytest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from webserver_utils.blast import BlastSearcher, OutputFormat

database_records: tuple[SeqRecord, ...] = tuple([
    SeqRecord(
        id="sp|P00762|TRYP_PIG",
        name="",
        description="",
        seq=Seq("IVGGYTCGANTVPYQVSLNSGYHFCGGSLINSQWVVSAAHCYKSRIQVRLGEHNIKVLEGNEQFINAAKIIRHPKYNRDTLNNDIMLIKLSSPATLNSRVATVSLPRCCA"),
    ),
    SeqRecord(
        id="sp|P00698|LYSC_CHICK",
        name="",
        description="",
        seq=Seq("MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"),
    ),
])

query_records: tuple[SeqRecord, ...] = tuple([
    SeqRecord(
        id="Test_TRYP",
        name="",
        description="",
        seq=Seq("HFCGGSLINSQWVVSAAHCYKSRAQVRLGEHNIKVLEGNEQWINAAKIIRHPKYNRDTLNNDGMLIKLSSPATLNSRAATVSLPRCCA"),
    ),
    SeqRecord(
        id="Test_LYSC",
        name="",
        description="",
        seq=Seq("MRSLLILVLCFLPLAALGKVFGRCELGGMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCQIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAFRNRCKGTDVQAWIRGCRL"),
    ),
])


def test_no_records():
    with pytest.raises(ValueError):
        BlastSearcher(tuple())


def test_database_creation(tmp_path):
    searcher = BlastSearcher(
        database_records=database_records,
        intermediate_files_dir=tmp_path,
    )
    searcher._generate_database()
    assert searcher._database_dir.is_dir() and len(list(searcher._database_dir.iterdir())) > 0


def test_blastp_single_query(tmp_path):
    searcher = BlastSearcher(
        database_records=database_records,
        intermediate_files_dir=tmp_path,
    )
    results = searcher.blastp(
        query=query_records[0],
        output_format=OutputFormat.BLAST_XML
    )

    assert results is not None and len(list(results)) == 1


def test_blastp_single_query_not_searchio_parseable(tmp_path):
    searcher = BlastSearcher(
        database_records=database_records,
        intermediate_files_dir=tmp_path,
    )
    results = searcher.blastp(
        query=query_records[0],
        output_format=OutputFormat.CSV
    )

    assert results is None


def test_blastp_multi_query(tmp_path):
    searcher = BlastSearcher(
        database_records=database_records,
        intermediate_files_dir=tmp_path,
    )
    results = searcher.blastp(query=query_records)
    assert results is not None and len(list(results)) == 2
