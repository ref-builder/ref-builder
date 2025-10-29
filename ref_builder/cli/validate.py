"""Validation utilities for Click commands."""

from typing import TypeVar

import click
from click import Context, Parameter

T = TypeVar("T")


def validate_no_duplicate_accessions(
    _ctx: Context, _param: Parameter, value: list[str]
) -> list[str]:
    """Validate that a sequence of accessions does not contain duplicate values.

    Intended to be used as a callback for Click validation.

    :param _ctx: the Click context
    :param _param: the Click parameter
    :param value: the sequence to validate
    :return: the validated sequence
    """
    if len(value) != len(set(value)):
        raise click.BadParameter("Duplicate accessions are not allowed.")

    return value
