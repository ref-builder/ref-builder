class IsolateInconsistencyWarning(UserWarning):
    """Warn when an isolate contains both RefSeq and non-RefSeq accessions.

    All sequences in an isolate should be sourced from the same database.
    """
