#!/usr/bin/env python

import argparse
import re
import numpy
import os
import glob
from subprocess import call

def get_args():
    parser = argparse.ArgumentParser(description='Analyze streamer_1d simulations', prog='analyze_sim_1d.py')
    parser.add_argument('name', type=str, help='simulation name')
    return parser.parse_args()

def get_file_number(sim_name, file_name):
    my_match = re.match(sim_name + r"_(\d+)\.txt$", file_name, re.I)
    if (my_match):
        return int(my_match.group(1))

if __name__ == '__main__':
    args = get_args()

    os.chdir("output")
    all_files = glob.glob(args.name + "*.txt")

    # Select numbered files, e.g., sim_part_12.txt etc
    sim_data_files = [ (get_file_number(args.name, f), f) for
                       f in all_files if get_file_number(args.name, f)]
    sim_data_files.sort(key=lambda x: x[0])
    n_files = len(sim_data_files)

    for ix in range(n_files):
        with open(sim_data_files[ix][1]) as f:
            print this_file
