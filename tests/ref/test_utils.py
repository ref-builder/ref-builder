import pytest

from ref_builder.ncbi.utils import parse_refseq_comment
from ref_builder.utils import get_accession_key, pad_zeroes


class TestPadZeroes:
    def test_ok(self):
        """Test that integers are padded with zeroes to 8 digit strings."""
        assert pad_zeroes(112) == "00000112"
        assert pad_zeroes(1) == "00000001"
        assert pad_zeroes(87654321) == "87654321"

    @pytest.mark.parametrize("number", [100000000, 12345678901234567890])
    def test_too_big(self, number: int):
        """Test that integers larger than 8 digits are returned as strings."""
        with pytest.raises(ValueError) as e:
            pad_zeroes(number)

        assert str(e.value) == "Number is too large to pad"


class TestRefSeqCommentParser:
    """Test the parsing of standardized RefSeq comments."""

    @pytest.mark.parametrize(
        ("comment_string", "status", "predecessor_accession"),
        [
            (
                "PROVISIONAL REFSEQ: "
                "This record has not yet been subject to final NCBI review. "
                "The reference sequence is identical to EF546808,"
                "COMPLETENESS: full length.",
                "PROVISIONAL REFSEQ",
                "EF546808",
            )
        ],
    )
    def test_ok(self, comment_string: str, status: str, predecessor_accession: str):
        extracted_status, extracted_accession = parse_refseq_comment(comment_string)

        assert extracted_accession == predecessor_accession

        assert extracted_status == status

    def test_missing_comment(self):
        """Test that an exception is raised if no REFSEQ information can be parsed."""
        with pytest.raises(ValueError, match="Invalid RefSeq comment"):
            parse_refseq_comment("AGSDHJFKLDSFEU")


class TestGetAccessionKey:
    @pytest.mark.parametrize(
        ("string", "expected_accession"),
        [
            ("EF546808", "EF546808"),
            ("EF546808.1", "EF546808"),
            ("NC_123456", "NC_123456"),
            ("NC_123456.1", "NC_123456"),
        ],
    )
    def test_ok(self, string: str, expected_accession: str):
        """Test that keys are extracted from various valid accession strings."""
        assert get_accession_key(string) == expected_accession

    @pytest.mark.parametrize(
        "string", ["CLOCKTOWER", "AC269481.3.5", "123453.1", "123453"]
    )
    def test_bad_input(self, string: str):
        """Test that an exception is raised for various bad inputs."""
        with pytest.raises(ValueError):
            get_accession_key(string)
