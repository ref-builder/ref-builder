from ref_builder.ncbi.client import NCBIClientProtocol
from ref_builder.repo import Repo
from ref_builder.services.isolate import IsolateService
from ref_builder.services.otu import OTUService


class Services:
    def __init__(self, repo: Repo, ncbi_client: NCBIClientProtocol) -> None:
        self.otu = OTUService(repo, ncbi_client)
        self.isolate = IsolateService(repo, ncbi_client)
