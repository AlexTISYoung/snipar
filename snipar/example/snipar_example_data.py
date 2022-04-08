#!/usr/bin/env python
import os
import shutil
import snipar.example
import argparse

def load_example_data(dest="example_data"):
    data_dir_name = "example_data"
    dir_name = os.path.dirname(snipar.example.__file__)
    data_dir = os.path.join(dir_name, data_dir_name)
    shutil.copytree(data_dir, dest)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dest',type=str,help='Directory to save example data to',default='example_data')
    args = parser.parse_args()
    load_example_data(args.dest)
