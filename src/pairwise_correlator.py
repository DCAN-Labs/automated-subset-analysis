#! /usr/bin/env python3

"""
Pairwise Correlator
Greg Conan: conan@ohsu.edu
Created 2020-08-03
Updated 2020-08-03
"""

##################################
#
# Compare 2 groups of pconns subject-by-subject using .conc files
#
##################################

# Imports
import argparse
from conan_tools import (
    count_digits_of, get_cli_args, get_pwd, GP_MTR_FILE, load_matrix_from,
    other_group_n, now, read_file_into_list, spearman_rho, track_progress,
    update_progress, valid_conc_file, valid_output_dir, valid_readable_file,
    validate
)
import numpy as np
import os
import subprocess
import time

# Constants
PWD = get_pwd()


# Functions


def main():
    cli_args = get_cli_args("Pairwise connectivity matrix correlator", (
        "axis_font_size", GP_MTR_FILE.format(1), GP_MTR_FILE.format(2),
        "graph_title", "hide_legend", "marker_size", "output",
        "spearman_rho", "title_font_size", "y_range"
    ), PWD, validate_cli_args)

    # Import each subject's matrices from both .conc files and compare them
    all_IDs, matrix_paths = import_concs(cli_args)
    all_correls = correlate_pairwise(all_IDs, matrix_paths, cli_args)
    print(all_correls)


def validate_cli_args(cli_args, parser):

    # Matrix path .conc files are required
    if not (getattr(cli_args, "matrices_conc_1", None)
            and getattr(cli_args, "matrices_conc_2", None)):
        parser.error("Please include the --matrices-conc-1 and "
                     "--matrices-conc-2 arguments. Even though they are shown"
                     "as optional, this script requires them.")

    return cli_args


def import_concs(cli_args):

    # Local variables: String for error to show if .conc files' subjects do not
    # match exactly, and also collections of subject IDs/paths to return
    MATCH_ERR = ("Line {0} of {1} has subject {2}, but line {0} of {3} has "
                 "subject {4}. Both lines should have the same subject.")
    CONC_1 = os.path.basename(getattr(cli_args, GP_MTR_FILE.format(1)))
    CONC_2 = os.path.basename(getattr(cli_args, GP_MTR_FILE.format(2)))
    all_IDs = []
    matrix_paths = dict()

    # Import file paths from .conc files
    for i in (1, 2):
        matrix_paths[i] = read_file_into_list(getattr(cli_args,
                                                      GP_MTR_FILE.format(i)))

    # Ensure that both .conc files have the same number of matrix file paths
    if len(matrix_paths[1]) != len(matrix_paths[2]):
        raise ValueError("Both matrix .conc files must have exactly the same "
                         "number of file paths, 1 per line.")  
    
    # Ensure that both matrix file path lists have the same subjects
    for i in range(len(matrix_paths[1])):
        id1 = get_subject_id_from(matrix_paths[1][i])
        id2 = get_subject_id_from(matrix_paths[2][i])
        if id1 == id2:
            all_IDs.append(id1)
        else:
            raise ValueError(MATCH_ERR.format(i, CONC_1, id1, CONC_2, id2))
                                                
    return all_IDs, matrix_paths



def get_subject_id_from(path):
    """
    Accepts the path to a subject's matrix file, and gets that subject's ID
    :param path: String which is a valid path to a CIFTI2 matrix file
    :return: String which is the subject ID which was in path
    """
    sub_id_pos = path.find("INV")  # Where in path does ID start
    id_end_pos = sub_id_pos + 11   # Where in path does ID end
    return path[sub_id_pos:id_end_pos]


def correlate_pairwise(all_IDs, matrix_paths, cli_args):

    # Local vars
    mx_pair = dict()  # Each pair of matrices for each subject
    num_subjs = len(matrix_paths[1])  # Total number of subjects
    correls = dict() # [None] * num_subjs  # Return value: List of correlations
    corr_fn = (spearman_rho if cli_args.spearman_rho else
               lambda x, y: np.corrcoef(x, y)[0, 1])  # Correlation function
    progress = track_progress(1, [1] * num_subjs)  # Estimate time this'll take
    started_at = now()  # Timestamp where matrix correlation started

    # Correlate every matrix with its corresponding matrix
    for i in range(num_subjs):
        for j in (matrix_paths.keys()):
            mx_pair[j] = load_matrix_from(matrix_paths[j][i])
        correls[all_IDs[i]] = corr_fn(mx_pair[1], mx_pair[2])
        if not i % count_digits_of(i):
            progress = update_progress(
                progress, "calculating pairwise correlations", 1, started_at
            )
    return correls
        

if __name__ == "__main__":
    main()