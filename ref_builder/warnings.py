class IsolateInconsistencyWarning(UserWarning):
    """Warn when an isolate contains both RefSeq and non-RefSeq accessions.

    All sequences in an isolate should be sourced from the same database.
    """


class OTUDeletedWarning(UserWarning):
    """A warning raised when a previously deleted OTU is requested."""


class PlanWarning(UserWarning):
    """Warns when the plan does not follow best practices."""
