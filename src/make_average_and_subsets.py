#! /usr/bin/env python3

"""
Greg Conan: conan@ohsu.edu
Created 2019-12-12
Last Updated 2019-12-18
"""

##################################
#
# Run make_average_matrix.py to get an average matrix of each entire group, and
# then run automated_subset_analysis.py using those average matrices
#
##################################

# Imports
from conan_tools import *
import shutil
import subprocess
import os

# Get the path to the Automated Subset Analysis (ASA) directory
PWD = get_pwd()
ASA_DIR = (os.path.abspath(os.path.dirname(PWD))
           if os.path.split(PWD)[1] == "src" else PWD)


def main():
    try:
        get_and_print_timestamp_when(sys.argv[0], "started")

        # Get all command-line arguments from user
        arg_names = remove_if_in(get_ASA_arg_names(), "continuous_variables")
        cli_args =  make_avg_paths_absolute(get_cli_args(
            "Make average matrices of groups and use them in subset analyses",
            arg_names, PWD
        ))

        # Get average matrix paths
        for gp_num in (1, 2):
            cli_args = add_default_avg_matr_path_to(cli_args, gp_num)

        # Make average matrices for each entire group
        run_ASA_or_MAM_script(
            os.path.join(PWD, "make_average_matrix.py"),
            ("matrices_conc_1", "matrices_conc_2", "group_1_avg_file",
             "group_2_avg_file", "output", "example_file", "inverse_fisher_z"),
             cli_args
        )
        print("Average matrices created: {}".format(*[
            getattr(cli_args, "group_{}_avg_file".format(x)) for x in (1, 2)
        ]))

        # Use those average matrices to run automated_subset_analysis.py
        asa_args = remove_if_in(arg_names, "example_file")
        run_ASA_or_MAM_script("automated_subset_analysis.py",
                              asa_args, make_avg_paths_absolute(cli_args))

        get_and_print_timestamp_when(sys.argv[0], "finished")

    except Exception as e:
        get_and_print_timestamp_when(sys.argv[0], "crashed")
        raise e


def make_avg_paths_absolute(cli_args):
    # Validate that average matrix file arguments are absolute paths
    for gp_num in (1, 2):
        avg_matr = "group_{}_avg_file".format(gp_num)
        avg_matr_path = getattr(cli_args, avg_matr, None)
        if avg_matr_path:
            setattr(cli_args, avg_matr, os.path.join(cli_args.output,
                                                     avg_matr_path))
    return cli_args
    

def empty_output_dir(cli_args):
    """
    Ensure that output dir is empty. If it isn't, then ask the user whether to 
    delete all of its contents or quit the script.
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses --output    
    :return: N/A
    """
    # If the output is empty, do nothing; otherwise ask user what to do
    files_in_output = os.listdir(cli_args.output)
    if files_in_output:
        print("Please use an empty directory as the --output argument.")
        clear_dir = "foo"
        while clear_dir.lower() != "y" and clear_dir.lower() != "n":
            clear_dir = input("Do you want to delete everything in the {} "
                              "directory? Y or N: ".format(cli_args.output))

            # Clear output dir if user says to
            if clear_dir.lower() == "y":
                for filename in files_in_output:
                    filepath = os.path.join(cli_args.output, filename)
                    try:
                        shutil.rmtree(filepath)
                    except OSError:
                        os.remove(filepath)
    
            # Otherwise quit
            elif clear_dir.lower() == "n":
                sys.exit(1)


def run_ASA_or_MAM_script(to_run, to_include, cli_args):
    """
    Run a script in the automated_subset_analysis directory
    :param to_run: String which is the base name of the python3 script to run
    :param to_include: List of strings which name the args for to_run to use
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    cmd = ["python3", os.path.join(ASA_DIR, to_run)]
    for arg in to_include:
        val = getattr(cli_args, arg, None)
        if val:
            if "demo_file" not in arg:
                cmd.append("--{}".format(arg).replace("_", "-"))
            if val is not True:
                if isinstance(val, list):
                    val = " ".join((str(v) for v in val))
                cmd.append(str(val))
    os.system(" ".join(cmd))


def remove_if_in(remove_from, to_remove):
    """
    :param remove_from: List of elements to remove something from
    :param to_remove: Object to remove from subtract_from, if it is in the list
    :return: remove_from, but without the object to_remove
    """
    if to_remove in remove_from:
        remove_from.remove(to_remove)
    return remove_from


if __name__ == "__main__":
    main()