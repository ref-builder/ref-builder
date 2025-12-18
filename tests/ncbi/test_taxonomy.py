from datetime import UTC, datetime, timedelta
from pathlib import Path
from unittest.mock import MagicMock

from pytest_mock import MockerFixture

from ref_builder.ncbi.taxonomy import (
    _get_current_taxdump,
    _parse_plant_hosts,
    _parse_virus_species,
    fetch_plant_virus_taxids,
)


class TestFetchPlantVirusTaxids:
    """Test fetching plant virus taxids."""

    def test_ok(self, mocker: MockerFixture, tmp_path: Path):
        """Test successful fetch of plant virus taxids."""
        mocker.patch(
            "ref_builder.ncbi.taxonomy._taxdump_cache_path",
            tmp_path,
        )

        # Create a mock taxdump directory with timestamp
        timestamp = datetime.now(UTC).strftime("%Y-%m-%dT%H-%M-%S")
        taxdump_dir = tmp_path / f"taxdump_{timestamp}"
        taxdump_dir.mkdir()

        # nodes.dmp: taxid | parent_taxid | rank | ...
        # Structure: 10239 (Viruses) -> 100 (family) -> 200 (species, virus)
        #                            -> 300 (species, virus)
        #            1 (root) -> 400 (species, not a virus)
        nodes_content = (
            "1\t|\t1\t|\tno rank\t|\n"
            "10239\t|\t1\t|\tsuperkingdom\t|\n"
            "100\t|\t10239\t|\tfamily\t|\n"
            "200\t|\t100\t|\tspecies\t|\n"
            "300\t|\t100\t|\tspecies\t|\n"
            "400\t|\t1\t|\tspecies\t|\n"
        )
        (taxdump_dir / "nodes.dmp").write_text(nodes_content)

        # host.dmp: taxid | host_name | ...
        # 200 has plant host, 300 has mammal host
        host_content = (
            "200\t|\tplants\t|\n"
            "300\t|\tmammals\t|\n"
            "400\t|\tplants\t|\n"
        )
        (taxdump_dir / "host.dmp").write_text(host_content)

        result = fetch_plant_virus_taxids()

        # Should only include 200 (virus species with plant host)
        # Not 300 (virus but mammal host)
        # Not 400 (plant host but not a virus)
        assert result == {200}


class TestGetCurrentTaxdump:
    """Test cache staleness checking."""

    def test_ok(self, mocker: MockerFixture, tmp_path: Path):
        """Test finding a valid cached taxdump."""
        mocker.patch(
            "ref_builder.ncbi.taxonomy._taxdump_cache_path",
            tmp_path,
        )

        timestamp = datetime.now(UTC).strftime("%Y-%m-%dT%H-%M-%S")
        taxdump_dir = tmp_path / f"taxdump_{timestamp}"
        taxdump_dir.mkdir()

        result = _get_current_taxdump()

        assert result == taxdump_dir

    def test_stale(self, mocker: MockerFixture, tmp_path: Path):
        """Test that stale cache returns None."""
        mocker.patch(
            "ref_builder.ncbi.taxonomy._taxdump_cache_path",
            tmp_path,
        )

        # Create a taxdump directory with old timestamp
        old_time = datetime.now(UTC) - timedelta(days=8)
        timestamp = old_time.strftime("%Y-%m-%dT%H-%M-%S")
        taxdump_dir = tmp_path / f"taxdump_{timestamp}"
        taxdump_dir.mkdir()

        result = _get_current_taxdump()

        assert result is None

    def test_no_cache(self, mocker: MockerFixture, tmp_path: Path):
        """Test that missing cache returns None."""
        mocker.patch(
            "ref_builder.ncbi.taxonomy._taxdump_cache_path",
            tmp_path / "nonexistent",
        )

        result = _get_current_taxdump()

        assert result is None


class TestParseVirusSpecies:
    """Test parsing virus species from nodes.dmp."""

    def test_ok(self, tmp_path: Path):
        """Test parsing virus species taxids."""
        nodes_path = tmp_path / "nodes.dmp"
        nodes_content = (
            "1\t|\t1\t|\tno rank\t|\n"
            "10239\t|\t1\t|\tsuperkingdom\t|\n"
            "100\t|\t10239\t|\tfamily\t|\n"
            "200\t|\t100\t|\tspecies\t|\n"
            "201\t|\t200\t|\tstrain\t|\n"
            "300\t|\t100\t|\tgenus\t|\n"
            "400\t|\t1\t|\tspecies\t|\n"
        )
        nodes_path.write_text(nodes_content)

        result = _parse_virus_species(nodes_path)

        # Only species under Viruses (10239)
        assert result == {200}


class TestParsePlantHosts:
    """Test parsing plant hosts from host.dmp."""

    def test_ok(self, tmp_path: Path):
        """Test parsing plant host taxids."""
        host_path = tmp_path / "host.dmp"
        host_content = (
            "100\t|\tplants\t|\n"
            "200\t|\tmammals\t|\n"
            "300\t|\tPlants and fungi\t|\n"
        )
        host_path.write_text(host_content)

        result = _parse_plant_hosts(host_path)

        assert result == {100, 300}
