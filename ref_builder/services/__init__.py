"""Repository services."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ref_builder.ncbi.client import NCBIClientProtocol
from ref_builder.repo import Repo

if TYPE_CHECKING:
    from ref_builder.services.cls import Services


class Service:
    """Base class for repository services.

    Provides access to the repository and NCBI client.
    """

    def __init__(
        self,
        repo: Repo,
        ncbi_client: NCBIClientProtocol,
        services: Services,
    ) -> None:
        """Initialize the service with a repository and NCBI client.

        :param repo: the repository instance
        :param ncbi_client: the NCBI client instance
        :param services: services container for accessing other services
        """
        self._repo = repo
        self._ncbi = ncbi_client
        self._services = services

    @property
    def ncbi(self) -> NCBIClientProtocol:
        """The NCBI client instance."""
        return self._ncbi


__all__ = ["Service"]
