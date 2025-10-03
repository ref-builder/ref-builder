"""Repository services."""

from ref_builder.ncbi.client import NCBIClient
from ref_builder.repo import Repo


class Service:
    """Base class for repository services.

    Provides access to the repository and NCBI client.
    """

    def __init__(self, repo: Repo, ncbi_client: NCBIClient) -> None:
        """Initialize the service with a repository and NCBI client.

        :param repo: the repository instance
        :param ncbi_client: the NCBI client instance
        """
        self._repo = repo
        self._ncbi = ncbi_client

    @property
    def ncbi(self) -> NCBIClient:
        """The NCBI client instance."""
        return self._ncbi


__all__ = [
    "Service",
]
