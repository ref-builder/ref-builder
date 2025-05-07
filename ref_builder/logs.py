import logging
import sys

import structlog


def configure_logger(verbosity: int, no_color: bool = False) -> None:
    """Configure structlog-based logging.

    :param verbosity: The verbosity level of the logger.
    :param no_color: Disable colored input, even if global settings allow it.
    """
    # Disable faker logging.
    logging.getLogger("faker").setLevel(logging.ERROR)

    processors = [
        structlog.processors.add_log_level,
        structlog.processors.TimeStamper(fmt="iso"),
    ]

    if verbosity == 0:
        level = logging.WARNING
    elif verbosity == 1:
        level = logging.INFO
    else:
        level = logging.DEBUG

        processors.append(
            structlog.processors.CallsiteParameterAdder(
                [
                    structlog.processors.CallsiteParameter.MODULE,
                    structlog.processors.CallsiteParameter.FUNC_NAME,
                ]
            )
        )

    if no_color:
        processors.append(structlog.processors.JSONRenderer())

    else:
        processors.append(structlog.dev.ConsoleRenderer())

    structlog.configure(
        logger_factory=structlog.PrintLoggerFactory(sys.stderr),
        processors=processors,
        wrapper_class=structlog.make_filtering_bound_logger(level),
    )
