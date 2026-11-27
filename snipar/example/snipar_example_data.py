#!/usr/bin/env python
"""Copies the bundled tutorial data into a new directory.

**Usage**

    ``snipar_example_data.py --dest example_data``

Args:
@parser@

Results:
    Example directory
        Genotypes, phenotypes, family information, and score weights for the tutorial.
"""
import os
import shutil
import snipar.example
import argparse

from snipar.docgen import get_parser_doc

parser = argparse.ArgumentParser()
parser.add_argument('--dest', type=str, default='example_data',
                    help='Directory to create for the example files')
__doc__ = __doc__.replace("@parser@", get_parser_doc(parser))


def load_example_data(dest="example_data"):
    data_dir_name = "example_data"
    dir_name = os.path.dirname(snipar.example.__file__)
    data_dir = os.path.join(dir_name, data_dir_name)
    shutil.copytree(data_dir, dest)

if __name__ == '__main__':
    args = parser.parse_args()
    load_example_data(args.dest)
