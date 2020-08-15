#! /usr/bin/env python3

"""
Parallel Pairwise Correlator
Greg Conan: conan@ohsu.edu
Created 2020-08-05
Updated 2020-08-05
"""

##################################
#
# Compare 2 groups of pconns subject-by-subject using .conc files by running
# multiple instance of pairwise_correlator.py in parallel
#
##################################

# Imports
import argparse
from conan_tools import (valid_whole_number,  valid_readable_file,
                         valid_conc_file, valid_output_dir,  validate)
import os
import subprocess
import sys
import time


# Functions


def main():

    # Collect parameters from user
    parser = argparse.ArgumentParser("Parallel pairwise pconn correlator")
    parser.add_argument(
        "-out", "--output", type=valid_output_dir
    )
    parser.add_argument(
        "-concs", "--concs-dir", type=valid_readable_file
    )
    parser.add_argument(
        "-range", "--conc-num-range", type=valid_whole_number, nargs=2
    )
    cli_args = parser.parse_args()
    
    # Run subprocess objects in parallel
    processes = []
    outlog = open("out_log.txt", "a+")
    conc_base = os.path.join(cli_args.concs_dir, "{}{}.conc")
    for conc_num in range(*cli_args.conc_num_range):
        cmd = ("python3", "pairwise_correlator.py",
               "--output", cli_args.output,
               "-conc1", conc_base.format("DCAN_to_DAIC_", conc_num),
               "-conc2", conc_base.format("DAIC_to_DCAN_", conc_num))
        print(cmd)
        processes.append(subprocess.Popen(cmd, stdout=outlog, bufsize=1,
                                          universal_newlines=True))
    outlog.close()
                

if __name__ == "__main__":
    main()

