"""Repository services."""

from ref_builder.ncbi.client import NCBIClientProtocol
from ref_builder.repo import Repo


class Service:
    """Base class for repository services.

    Provides access to the repository and NCBI client.
    """

    def __init__(self, repo: Repo, ncbi_client: NCBIClientProtocol) -> None:
        """Initialize the service with a repository and NCBI client.

        :param repo: the repository instance
        :param ncbi_client: the NCBI client instance
        """
        self._repo = repo
        self._ncbi = ncbi_client

    @property
    def ncbi(self) -> NCBIClientProtocol:
        """The NCBI client instance."""
        return self._ncbi


__all__ = [
    "Service",
]
