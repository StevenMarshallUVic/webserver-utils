"""NCBI BLAST."""

import logging
import os
import subprocess
from dataclasses import dataclass
from enum import StrEnum, Enum
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Generator, Iterable

from Bio import SearchIO, SeqIO
# noinspection protected-member
from Bio.SearchIO._model.query import QueryResult
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(Path(__file__).name)


class InputType(StrEnum):
    FASTA = "fasta"
    BLASTDB = "blastdb"
    ASN1_TXT = "asn1_txt"
    ASN1_BIN = "asn1_bin"


class DatabaseType(StrEnum):
    NUCLEOTIDE = "nucl"
    PROTEIN = "prot"


class BlastProgram(StrEnum):
    BLASTN = "blastn"
    BLASTP = "blastp"
    BLASTX = "blastx"
    TBLASTN = "tblastn"
    TBLASTX = "tblastx"


class OutputFormat(Enum):
    PAIRWISE = (0, "blast-text")
    QUERY_ANCHORED_SHOW_IDENTITIES = (1, None)
    QUERY_ANCHORED_NO_IDENTITIES = (2, None)
    FLAT_QUERY_ANCHORED_SHOW_IDENTITIES = (3, None)
    FLAT_QUERY_ANCHORED_NO_IDENTITIES = (4, None)
    BLAST_XML = (5, "blast-xml")
    TABULAR = (6, "blast-tab")
    TABULAR_COMMENT_LINES = (7, "blast-tab")
    SEQALIGN_TEXT_ASN1 = (8, None)
    SEQALIGN_BINARY_ASN1 = (9, None)
    CSV = (10, None)
    BLAST_ARCHIVE_ASN1 = (11, None)
    SEQALIGN_JSON = (12, None)
    MULTI_FILE_BLAST_JSON = (13, None)
    MULTI_FILE_BLAST_XML2 = (14, "blast-xml")
    SINGLE_FILE_BLAST_JSON = (15, None)
    SINGLE_FILE_BLAST_XML2 = (16, "blast-xml")
    ORGANISM_REPORT = (18, None)
    CSV_WITH_HEADER_LINES = (20, None)

    def __init__(self, outfmt: int, searchio_format: str | None):
        self.outfmt = outfmt
        self.searchio_format = searchio_format


@dataclass(frozen=True)
class BlastSearcher:
    """Handles performing BLAST searches.

    Attributes
    ----------
    database_records
        Records to use in database.
    intermediate_files_dir
        (Optional) Directory to write intermediate files to.
        When not specified, intermediate files are written to a temporary
        directory that gets deleted following program execution.
    database_type
        (Optional) What type of database your `database_records` are in.
    database_title
        (Optional) Custom name for the BLAST database.
    database_prefix
        (Optional) Prefix to use for BLAST database file names.
    num_threads
        (Optional) Number of threads to use for performing BLAST searches.
        When not specified, uses all available CPUs.
    parse_seqids
        (Optional) Whether BLAST should strip input identifiers and
        restructure them into standard NCBI/UniProt identifiers,
        such as adding `sp|` prefixes before UniProt entry IDs.
    taxid
        (Optional) Taxonomy ID to assign to all sequences.
    """

    database_records: tuple[SeqRecord, ...]
    intermediate_files_dir: Path | None = None
    database_type: DatabaseType = DatabaseType.PROTEIN
    database_title: str = "BLAST Database"
    database_prefix: str = "database"
    num_threads: int = os.cpu_count() or 1
    parse_seqids: bool = False
    taxid: int | None = None
    verbose: bool = False

    _input_type: InputType = InputType.FASTA
    _temp_dir = TemporaryDirectory()

    @property
    def _temp_dir_path(self) -> Path:
        """Path to temp directory for this instance."""
        return Path(self._temp_dir.name)

    @property
    def _database_dir(self) -> Path:
        """Path to root directory of BLAST database."""
        return self.intermediate_files_dir / "database"

    @property
    def _database_path(self) -> Path:
        """Path to BLAST database, including prefix of database files."""
        return self._database_dir / self.database_prefix

    def __post_init__(self):
        if len(self.database_records) <= 0:
            raise ValueError("Must provide at least one database record.")

        if self.intermediate_files_dir is None:
            object.__setattr__(
                self,
                "intermediate_files_dir",
                self._temp_dir_path
            )

    def __del__(self):
        self._cleanup()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._cleanup()

        # Do not suppress exceptions
        return False

    def blastn(self, **kwargs):
        """Compare nucleotide queries against a nucleotide database.

        See `_blast()` method signature for arguments.
        """

        return self._blast(blast_program=BlastProgram.BLASTN, **kwargs)

    def blastp(self, **kwargs):
        """Compare protein queries against a protein database.

        See `_blast()` method signature for arguments.
        """

        return self._blast(blast_program=BlastProgram.BLASTP, **kwargs)

    def blastx(self, **kwargs):
        """Compare nucleotide queries (translated in 6 frames)
        against a protein database.

        See `_blast()` method signature for arguments.
        """
        return self._blast(blast_program=BlastProgram.BLASTX, **kwargs)

    def tblastn(self, **kwargs):
        """Compare protein queries against a nucleotide database
        (translated in 6 frames).

        See `_blast()` method signature for arguments.
        """
        return self._blast(blast_program=BlastProgram.TBLASTN, **kwargs)

    def tblastx(self, **kwargs):
        """Compare translated nucleotide queries against a translated
        nucleotide database.

        See `_blast()` method signature for arguments.
        """
        return self._blast(blast_program=BlastProgram.TBLASTX, **kwargs)

    def _generate_database(self) -> None:
        """Create BLAST database for `database_records`."""

        logger.debug(
            "Generating BLAST database with %s sequences named %s...",
            len(self.database_records),
            self.database_title,
        )
        if self._database_dir.is_dir():
            raise IsADirectoryError(
                f"Database already exists at '{self._database_dir}'."
            )

        database_fasta = self.intermediate_files_dir / "database.fasta"
        if not database_fasta.is_file():
            SeqIO.write(
                self.database_records,
                database_fasta,
                "fasta-2line"
            )

        self._database_dir.mkdir(parents=True)
        args = [
            "makeblastdb",
            "-in", database_fasta,
            "-input_type", self._input_type.value,
            "-dbtype", self.database_type.value,
            "-title", self.database_title,
            "-out", self._database_path,
        ]
        if self.parse_seqids:
            args.append("-parse_seqids")
        if self.taxid is not None:
            args.extend(["-taxid", str(self.taxid)])
        self._run_subprocess_for_args(args)

        if not self._database_dir.is_dir():
            raise NotADirectoryError(
                f"Could not find database at '{self._database_dir}'."
            )

        if len(list(self._database_dir.iterdir())) <= 0:
            raise FileNotFoundError("Database generation failed.")

    def _blast(
            self,
            blast_program: BlastProgram,
            query: SeqRecord | Iterable[SeqRecord] | Path,
            results_path: Path | None = None,
            output_format: OutputFormat = OutputFormat.BLAST_XML,
    ) -> Generator[QueryResult, None, None] | None:
        """Perform a query through a specified BLAST program.

        Parameters
        ----------
        blast_program
            What kind of BLAST search to perform.
        query
            The inputs to BLAST against the database.
        results_path
            (Optional) Where to write results to.
            If not specified, results are written to a temporary file and
            deleted following program execution.
        output_format
            What format the output should be formatted in.

        Returns
        -------
        Generator[QueryResult, None, None] | None
            Generator of results if `output_format` is supported by SearchIO,
            otherwise None.
        """

        if not self._database_dir.is_dir():
            self._generate_database()

        query_fasta: Path
        if isinstance(query, Path):
            query_fasta = query
        else:
            query_fasta = self.intermediate_files_dir / NamedTemporaryFile().name
            if query_fasta.is_file():
                raise FileExistsError(
                    f"Query FASTA already found at '{query_fasta}'."
                )

            SeqIO.write(query, query_fasta, "fasta-2line")

        if results_path is None:
            results_path = self.intermediate_files_dir / NamedTemporaryFile().name
            if results_path.is_file():
                raise FileExistsError(
                    f"BLAST results already found at '{results_path}'."
                )

        logger.debug(
            "Running %s for %s sequences against database '%s'...",
            blast_program.value,
            len(list(SeqIO.parse(query_fasta, "fasta"))),
            self.database_title,
        )

        args = [
            blast_program.value,
            "-query", query_fasta,
            "-db", self._database_path,
            "-out", results_path,
            "-outfmt", str(output_format.outfmt),
            "-num_threads", str(self.num_threads),
        ]
        self._run_subprocess_for_args(args)

        if output_format.searchio_format is None:
            logger.info(
                "%s output format cannot be parsed by SearchIO. "
                "Must specify `output_file` and manually handle "
                "downstream processing.",
                output_format.value,
            )
            return None

        return SearchIO.parse(results_path, output_format.searchio_format)

    def _run_subprocess_for_args(self, args) -> None:
        try:
            result = subprocess.run(
                args,
                capture_output=True,
                text=True,
                check=True
            )
            if self.verbose:
                if result.stdout:
                    print(f"STDOUT:\n{result.stdout}")
                if result.stderr:
                    print(f"STDERR:\n{result.stderr}")

        except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e}")
            print(f"STDOUT:\n{e.stdout}")
            print(f"STDERR:\n{e.stderr}")
            raise e

    def _cleanup(self):
        """Cleanup temp directory."""

        if hasattr(self, "_temp_dir"):
            self._temp_dir.cleanup()


if __name__ == "__main__":
    main()
