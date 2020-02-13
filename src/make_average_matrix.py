#! /usr/bin/env python3

"""
Average matrix creation for automated subset analysis script
Greg Conan: conan@ohsu.edu
Created 2019-11-20
Last Updated 2019-02-13
"""

##################################
#
# Script to make average matrices of two groups of data, to use as inputs for
# automated_subset_analysis.py
#
##################################

# Imports
import argparse
from conan_tools import *
import nibabel
import numpy as np
import os
import pandas as pd
import socket
import sys

# Constants
PWD = get_pwd()
DEFAULT_DEM_VAR_PCONNS = "pconn10min"


def main():

    # Store and print the date and time when this script started running
    starting_timestamp = get_and_print_timestamp_when(sys.argv[0], "started")

    # Get and validate all command-line arguments from user
    cli_args = get_cli_args(
        "Script to get the average matrices of two subject sets",
        ["example_file", "group_1_avg_file", "group_2_avg_file", 
         "inverse_fisher_z", "matrices_conc_1", "matrices_conc_2", "output"],
        PWD, validate_cli_args
    )
    
    # Make and save average matrix for each group
    # The strings passed to build_avg_matrix below are the names of the group
    # files' columns: One for matrix file path and one for subject ID
    paths_col = "scalar"
    for group_num in ("1", "2"):      
        print("Making average matrix for group {}".format(group_num))
        save_to_cifti2(get_average_matrix(get_total_matrix_from_conc(
            cli_args, group_num, paths_col, "id_redcap"
        ), paths_col, cli_args), cli_args, group_num)

    # Print when this script started and finished
    print(starting_timestamp)
    get_and_print_timestamp_when(sys.argv[0], "finished")


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all given command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    # Matrix path .conc files are required
    if not (getattr(cli_args, "matrices_conc_1", None)
            and getattr(cli_args, "matrices_conc_2", None)):
        parser.error("Please include the --matrices-conc-1 and "
                     "--matrices-conc-2 arguments. Even though they are shown"
                     "as optional, this script requires them.")

    # Get an arbitrary cifti2 file to use as template
    if not cli_args.example_file:
        with open(cli_args.matrices_conc_1) as infile:
            try:
                example = rename_exacloud_path(infile.readline().strip())
                cli_args.example_file = valid_readable_file(example)
            except argparse.ArgumentTypeError:
                parser.error("Tried and failed to open {} from the first line "
                             "of {} as a file path. Please use the "
                             "--example-file argument.".format(example, 
                             cli_args.matrices_conc_1))
    if os.path.splitext(cli_args.example_file)[1] != ".nii":
        parser.error("All input data files must be cifti2 files (*.nii)")
    return cli_args


def get_total_matrix_from_conc(cli_args, gp_num, paths_col, id_col):
    """
    Get a DataFrame with all input data to be averaged
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the group_1_conc and group_2_conc arguments.
    :param gp_num: String that is just the group number of the matrix to get
    :param col_name: String naming the column with all of the matrix file paths
    :return: pandas.DataFrame with 2 columns (subjectkey and scalar file path),
             and a String with the 2-part file extension of the scalar files
    """
    # Import the file with a list of matrix/scalar file paths into a DataFrame
    group_scalars_df = pd.read_csv(
        getattr(cli_args, "matrices_conc_{}".format(gp_num)),
        sep="\n", header=None
    ).rename(columns=lambda x: paths_col)

    # Rename Exacloud paths
    if "exa" not in socket.gethostname():
        group_scalars_df[paths_col] = (group_scalars_df[paths_col]
                                       .apply(rename_exacloud_path))

    # Return the DataFrame with both columns, and also the file extension
    group_scalars_df[id_col] = (group_scalars_df[paths_col]
                                .apply(extract_subject_id_from))
    return group_scalars_df


def save_to_cifti2(avg_matrix, cli_args, gp_num):
    """
    Save the average matrix as a cifti2 file (*.nii) by importing an arbitrary
    cifti2 matrix and saving a copy of it with its data replaced by the data of
    the average matrix
    :param avg_matrix: numpy.ndarray with data to save into a cifti2 file
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the group_1_conc, group_2_conc,
                     --example_file, and --output arguments.
    :param gp_num: String that is just the group number of the matrix to save
    :return: N/A
    """
    group_file = getattr(cli_args, "matrices_conc_{}".format(gp_num))

    # Import an arbitrary input .nii file
    nii_matrix = nibabel.cifti2.load(cli_args.example_file)
    
    # Replace the arbitrary file's data with the average matrix's data
    avg_nii_matrix = nibabel.cifti2.cifti2.Cifti2Image(
        dataobj = avg_matrix,
        header = nii_matrix.header,
        nifti_header = nii_matrix.nifti_header,
        file_map = nii_matrix.file_map
    )

    # Get the file path to save the average matrix into
    cli_args = add_default_avg_matr_path_to(cli_args, gp_num,
                                            argparse.ArgumentParser(),
                                            DEFAULT_DEM_VAR_PCONNS)
    avg_matr = getattr(cli_args, "group_{}_avg_file".format(gp_num))

    # Save the average matrix file
    print("Saving group {} average matrix to {}".format(gp_num, avg_matr))
    nibabel.cifti2.save(avg_nii_matrix, avg_matr)


if __name__ == "__main__":
    main()
