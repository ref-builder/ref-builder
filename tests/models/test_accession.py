"""Tests for Accession model."""

import pytest

from ref_builder.models.accession import Accession


class TestAccession:
    """Test the Accession dataclass."""

    def test_from_string(self):
        """Test creating an Accession from a versioned string."""
        accession = Accession.from_string("MN908947.3")

        assert accession.key == "MN908947"
        assert accession.version == 3

    def test_from_string_refseq(self):
        """Test creating a RefSeq Accession from a versioned string."""
        accession = Accession.from_string("NC_045512.2")

        assert accession.key == "NC_045512"
        assert accession.version == 2

    def test_from_string_no_version(self):
        """Test that from_string fails when no version is provided."""
        with pytest.raises(ValueError, match="does not contain two parts"):
            Accession.from_string("MN908947")

    def test_from_string_invalid_version(self):
        """Test that from_string fails when version is not an integer."""
        with pytest.raises(ValueError, match="is not an integer"):
            Accession.from_string("MN908947.abc")

    def test_str(self):
        """Test string representation of Accession."""
        accession = Accession("MN908947", 3)

        assert str(accession) == "MN908947.3"

    def test_repr(self):
        """Test repr representation of Accession."""
        accession = Accession("NC_123456", 1)

        assert repr(accession) == "Accession(NC_123456, 1)"

    def test_equality(self):
        """Test equality comparison."""
        accession_1 = Accession("MN908947", 3)
        accession_2 = Accession("MN908947", 3)
        accession_3 = Accession("MN908947", 2)

        assert accession_1 == accession_2
        assert accession_1 != accession_3

    def test_hashable(self):
        """Test that Accession is hashable and works in sets/dicts."""
        accession_1 = Accession("NC_123456", 1)
        accession_2 = Accession("NC_123456", 1)
        accession_3 = Accession("NC_123456", 2)

        assert hash(accession_1) == hash(accession_2)
        assert hash(accession_1) != hash(accession_3)

        accession_set = {accession_1, accession_2, accession_3}
        assert len(accession_set) == 2

        accession_dict = {accession_1: "first", accession_3: "second"}
        assert accession_dict[accession_2] == "first"

    def test_sorting(self):
        """Test that accessions sort correctly by (key, version) tuple."""
        accessions = [
            Accession("NC_123456", 2),
            Accession("AB_999999", 1),
            Accession("NC_123456", 1),
            Accession("AB_111111", 3),
            Accession("AB_111111", 1),
        ]

        assert sorted(accessions) == [
            Accession("AB_111111", 1),
            Accession("AB_111111", 3),
            Accession("AB_999999", 1),
            Accession("NC_123456", 1),
            Accession("NC_123456", 2),
        ]


class TestFromString:
    """Test the from_string() class method."""

    def test_empty_string(self):
        """Test that from_string rejects empty string."""
        with pytest.raises(ValueError, match="cannot be empty or whitespace"):
            Accession.from_string("")

    def test_whitespace_only(self):
        """Test that from_string rejects whitespace-only string."""
        with pytest.raises(ValueError, match="cannot be empty or whitespace"):
            Accession.from_string("   ")

    def test_no_period(self):
        """Test that from_string rejects string without period."""
        with pytest.raises(
            ValueError, match="does not contain two parts delimited by a period"
        ):
            Accession.from_string("NC_123456")

    def test_multiple_periods(self):
        """Test that from_string rejects string with multiple periods."""
        with pytest.raises(ValueError, match="contains multiple periods"):
            Accession.from_string("NC_123456.1.2")

    def test_non_integer_version(self):
        """Test that from_string rejects non-integer version."""
        with pytest.raises(ValueError, match="is not an integer"):
            Accession.from_string("NC_123456.abc")


class TestIsRefSeq:
    """Test the is_refseq property."""

    @pytest.mark.parametrize(
        "accession_key",
        [
            "NC_010314",
            "NC_055390",
            "NC_003355",
            "NC_000001",
        ],
    )
    def test_refseq_6_digit(self, accession_key: str):
        """Test that 6-digit NC_ accessions are recognized as RefSeq."""
        accession = Accession(accession_key, 1)

        assert accession.is_refseq

    @pytest.mark.parametrize(
        "accession_key",
        [
            "NC_1234567",
            "NC_12345678",
        ],
    )
    def test_refseq_7_plus_digit(self, accession_key: str):
        """Test that 7+ digit NC_ accessions are recognized as RefSeq."""
        accession = Accession(accession_key, 1)

        assert accession.is_refseq

    @pytest.mark.parametrize(
        "accession_key",
        [
            "NC_ABC123",
            "NC_XYZ789",
            "NC_A1B2C3",
        ],
    )
    def test_refseq_alphanumeric(self, accession_key: str):
        """Test that alphanumeric NC_ accessions are recognized as RefSeq."""
        accession = Accession(accession_key, 1)

        assert accession.is_refseq

    @pytest.mark.parametrize(
        "accession_key",
        [
            "AB100000",
            "MN908947",
            "DQ178610",
            "MT240513",
        ],
    )
    def test_not_refseq(self, accession_key: str):
        """Test that GenBank accessions are not recognized as RefSeq."""
        accession = Accession(accession_key, 1)

        assert not accession.is_refseq

    @pytest.mark.parametrize(
        "accession_key",
        [
            "nc_010314",
            "Nc_010314",
            "NC010314",
            "NC_",
        ],
    )
    def test_invalid_refseq_format(self, accession_key: str):
        """Test that invalid NC_ formats are not recognized as RefSeq."""
        accession = Accession(accession_key, 1)

        assert not accession.is_refseq
