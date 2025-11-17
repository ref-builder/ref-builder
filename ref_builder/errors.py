from uuid import UUID


class HydrationError(Exception):
    """Raised when a hydration fails."""


class HydrationIsolateError(HydrationError):
    """Raised when there is a problem hydrating a isolate."""


class HydrationSequenceError(HydrationError):
    """Raised when there is a problem hydrating a sequence."""


class LockConflictError(Exception):
    """Raised when a lock is attempted on a locked repository."""

    def __init__(self):
        super().__init__("Repository is already locked by another process.")


class LockRequiredError(Exception):
    """Raised when an operation is attempted without a lock."""

    def __init__(self) -> None:
        super().__init__("Repository must be locked for this operation.")


class InvalidInputError(Exception):
    """Raised when a command input does not match expected paramters."""

    def __init__(self, message: str) -> None:
        self.message = message

        super().__init__(message)


class OTUDeletedError(Exception):
    """Raised when an OTU is deleted."""

    def __init__(self, otu_id: UUID):
        super().__init__(f"OTU {otu_id} has been marked for deletion.")


class OTUExistsError(Exception):
    """Raised when attempting to create an OTU with a taxonomy ID that already exists."""

    def __init__(self, taxid: int, otu_id: UUID):
        self.taxid = taxid
        self.otu_id = otu_id
        super().__init__(f"OTU already exists for taxonomy ID {taxid}")


class PlanCreationError(ValueError):
    """Raised when a plan cannot be created from provided parameters."""


class PlanValidationError(ValueError):
    """Raised when an isolate does not pass plan validation."""
