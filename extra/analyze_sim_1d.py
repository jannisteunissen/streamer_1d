#!/usr/bin/env python

import argparse
import re
import numpy
import glob


def get_args():
    parser = argparse.ArgumentParser(prog='analyze_sim_1d.py')
    parser.add_argument('name', type=str, help='simulation name')
    parser.add_argument('-Ech', type=float,
                        default=1.0e4, help='E-field in channel')
    return parser.parse_args()


def get_file_number(sim_name, file_name):
    my_match = re.match(sim_name + r"_(\d+)\.txt$", file_name, re.I)
    if (my_match):
        return int(my_match.group(1))

if __name__ == '__main__':
    args = get_args()

    all_files = glob.glob(args.name + "*.txt")

    # Select numbered files, e.g., sim_part_12.txt etc
    sim_data_files = [ (get_file_number(args.name, f), f) for
                       f in all_files if get_file_number(args.name, f)]
    sim_data_files.sort(key=lambda x: x[0])
    n_files = len(sim_data_files)

    for ix in range(n_files):
        data = numpy.loadtxt(sim_data_files[ix][1])
        efield = data[:, 1]
        edens_ch = numpy.ma.masked_array(data[:, 2],
                                         mask=(efield > float(args.Ech)))

        print("file #" + str(sim_data_files[ix][0]) +
              " - n_ch: " + str(edens_ch.mean()))
