"""Format argparse options for the Sphinx command reference."""

import argparse
import re

from snipar.numrange import parseNumRange


def get_parser_doc(parser):
    """Render visible options and positional arguments as an RST definition list."""
    entries = []
    types = {str: "string", int: "integer", float: "float", parseNumRange: "integers or ranges"}
    for action in parser._actions:
        if action.help == argparse.SUPPRESS:
            continue

        names = action.option_strings or [action.dest]
        label = ", ".join(f"``{name}``" for name in names)
        if action.nargs == 0:
            details = ["flag"]
        else:
            value_type = types.get(action.type, "string")
            if action.nargs == "+":
                value_type = f"one or more {value_type} values"
            details = [value_type]
        if action.required:
            details.append("required" if action.option_strings else "required positional argument")
        if action.default not in (None, argparse.SUPPRESS) and action.nargs != 0:
            details.append(f"default={action.default}")

        description = " ".join((action.help or "").split())
        description = re.sub(r"(?<!\w)--[\w-]+", r"``\g<0>``", description)
        entries.append(f"    {label} : {', '.join(details)}\n        {description}")
    return "\n\n".join(entries)
