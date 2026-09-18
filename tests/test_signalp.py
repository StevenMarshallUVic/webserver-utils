from pathlib import Path

import pytest

from webserver_utils.signalp import SignalPeptideType, SignalPeptide, SignalPResults


class TestSignalPeptideType:
    def test_valid_signal_peptide_type(self):
        assert (
            SignalPeptideType("lipoprotein_signal_peptide")
            == SignalPeptideType.LIPOPROTEIN_SIGNAL_PEPTIDE
        )


    def test_invalid_signal_peptide_type(self):
        with pytest.raises(ValueError):
            SignalPeptideType("nonexistent_peptide_type")


    def test_signal_peptide_type_display_name(self):
        assert(
            SignalPeptideType.SIGNAL_PEPTIDE.display_name
            == "Signal Peptide (Sec/SPI)"
        )


class TestSignalPResults:
    @pytest.fixture
    def signalp_results_zip(self, test_data_dir):
        return test_data_dir / "signalp" / "output_all_results.zip"

    @pytest.fixture
    def signalp_results_tarball(self, test_data_dir):
        return test_data_dir / "signalp" / "output_all_results.tar.gz"

    def test_missing_results_zip(self, signalp_results_zip):
        with pytest.raises(FileNotFoundError):
            SignalPResults.from_signalp_results_zip(results_zip=Path())

    def test_incorrect_results_zip_suffix(self, signalp_results_tarball):
        with pytest.raises(TypeError):
            SignalPResults.from_signalp_results_zip(
                results_zip=signalp_results_tarball
            )

    def test_output_df_length(self, signalp_results_zip):
        assert len(
            SignalPResults.from_signalp_results_zip(
                signalp_results_zip
            ).output_df
        ) == 2

    def test_processed_entries_length(self, signalp_results_zip):
        assert len(
            SignalPResults.from_signalp_results_zip(
                signalp_results_zip
            ).processed_entries
        ) == 2

    def test_signal_peptides_length(self, signalp_results_zip):
        assert len(
            SignalPResults.from_signalp_results_zip(
                signalp_results_zip
            ).signal_peptides
        ) == 2

    def test_protein_id_to_signal_peptide_length(self, signalp_results_zip):
        assert len(
            SignalPResults.from_signalp_results_zip(signalp_results_zip).protein_id_to_signal_peptide
        ) == 2
