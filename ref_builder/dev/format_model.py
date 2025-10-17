"""Format Pydantic models as multi-line Python code."""

from enum import Enum
from typing import Any

from pydantic import BaseModel


def format_value(value: Any, indent: int = 0) -> str:
    """Format a value as Python code with proper indentation.

    Args:
        value: The value to format
        indent: Current indentation level

    Returns:
        Formatted Python code string

    """
    indent_str = "    " * indent
    next_indent_str = "    " * (indent + 1)

    if isinstance(value, BaseModel):
        # Format Pydantic model
        class_name = value.__class__.__name__
        fields = []

        for field_name, field_value in value.model_dump().items():
            formatted_value = format_value(getattr(value, field_name), indent + 1)

            # Check if the formatted value contains newlines
            if "\n" in formatted_value:
                fields.append(f"{next_indent_str}{field_name}={formatted_value}")
            else:
                fields.append(f"{next_indent_str}{field_name}={formatted_value}")

        if not fields:
            return f"{class_name}()"

        return f"{class_name}(\n" + ",\n".join(fields) + f"\n{indent_str})"

    if isinstance(value, list):
        if not value:
            return "[]"

        # Format list items
        items = []
        for item in value:
            formatted_item = format_value(item, indent + 1)
            items.append(f"{next_indent_str}{formatted_item}")

        return "[\n" + ",\n".join(items) + f"\n{indent_str}]"

    if isinstance(value, Enum):
        # Format enum as EnumClass.VALUE
        return f"{value.__class__.__name__}.{value.name}"

    if isinstance(value, str):
        # For very long strings (like sequences), keep them on one line but repr them
        return repr(value)

    if isinstance(value, (int, float, bool, type(None))):
        return repr(value)

    # Fallback to repr for other types
    return repr(value)


def format_model(model: BaseModel, indent: int = 0) -> str:
    """Format a Pydantic model as properly indented Python code.

    Args:
        model: The Pydantic model to format
        indent: Starting indentation level

    Returns:
        Multi-line Python code string

    """
    return format_value(model, indent)


def format_value_list(values: list[Any], indent: int = 0) -> str:
    """Format a list of values as properly indented Python code.

    Args:
        values: List of values to format
        indent: Starting indentation level

    Returns:
        Multi-line Python code string for the list

    """
    return format_value(values, indent)
