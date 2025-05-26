from ref_builder.resources import RepoOTU
from ref_builder.utils import IsolateName, IsolateNameType


def get_isolate_id_by_name(otu: RepoOTU, isolate_type: str, name_value: str):
    try:
        isolate_name = IsolateName(IsolateNameType(isolate_type), name_value)
    except ValueError:
        return None

    return otu.get_isolate_id_by_name(isolate_name)
