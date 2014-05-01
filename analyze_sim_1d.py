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

if __name__ == '__main__':
    args = get_args()

    os.chdir("/output")
    all_files = glob.glob(args.name + "*.txt")

    # Select numbered files, e.g., sim_part_12.txt etc
    sim_data = [file for file in all_files
