import argparse
from pathlib import Path
import re
import shutil
import subprocess

import pytest

from snipar.docgen import get_parser_doc
from snipar.numrange import parseNumRange


def test_command_reference_preserves_cli_names_and_defaults():
    parser = argparse.ArgumentParser()
    parser.add_argument("out", help="Output prefix.")
    parser.add_argument("--range", "--chr-range", type=parseNumRange, nargs="*", help="Chromosomes.")
    parser.add_argument("--count", type=int, default=0, help="Repeat\nwith --chr-range.")
    parser.add_argument("--hidden", help=argparse.SUPPRESS)

    reference = get_parser_doc(parser)
    assert "``out`` : string, required positional argument" in reference
    assert "``--range``, ``--chr-range`` : integers or ranges" in reference
    assert "``--count`` : integer, default=0\n        Repeat with ``--chr-range``." in reference
    assert "``-h``, ``--help`` : flag" in reference
    assert argparse.SUPPRESS not in reference
    assert "--hidden" not in reference


def test_correlation_reference_covers_runtime_flags():
    from snipar.scripts.correlate import build_parser

    reference = (Path(__file__).resolve().parents[1] / "docs/scripts.rst").read_text()
    options = reference.split(".. _correlate.py:", 1)[1].split("Args:", 1)[1].split("Results:", 1)[0]
    flags = {flag for action in build_parser()._actions for flag in action.option_strings}
    assert set(re.findall(r"``(-[\w-]+)", options)) == flags


@pytest.mark.skipif(not shutil.which("Rscript"), reason="Rscript is unavailable")
def test_binary_command_reference_covers_runtime_flags():
    root = Path(__file__).resolve().parents[1]
    reference = (root / "docs/scripts.rst").read_text()
    options = reference.split(".. _pgs.R:", 1)[1].split("Args:", 1)[1].split("Results:", 1)[0]
    help_text = subprocess.check_output(["Rscript", str(root / "snipar/scripts/pgs.R"), "--help"], text=True)
    runtime_flags = set(re.findall(r"--[\w-]+", help_text))
    assert set(re.findall(r"``(--[\w-]+)", options)) == runtime_flags
