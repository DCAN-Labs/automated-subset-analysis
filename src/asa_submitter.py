#! /usr/bin/env python3

"""
SBATCH job submitter for automated subset analysis
Greg Conan: conan@ohsu.edu
Created 2020-01-03
Updated 2020-02-11
"""

##################################
#
# Script to submit a SBATCH job, running many instances of
# automated_subset_analysis.py, for parallel processing on Exacloud
#
##################################

# Imports
import argparse
from automated_subset_analysis import validate_cli_args
import os
from src.conan_tools import *
import subprocess
import time

# Constant: Directory holding automated_subset_analysis.py file
PWD = get_pwd()


def main():

    # Store and print the date and time when this script started running
    starting_timestamp = get_and_print_timestamp_when(sys.argv[0], "started")

    # Get and validate all command-line arguments from user
    cli_args = get_submitter_cli_args(
        ("Script to run many instances of automated_subset_analysis.py in "
         "parallel."), get_ASA_arg_names(), PWD
    )

    try:
        submit_batch_job(cli_args)
        
    except Exception as e:
        get_and_print_timestamp_when(sys.argv[0], "crashed")
        raise e

    # Print the date and time when this script started and finished running
    print(starting_timestamp)
    get_and_print_timestamp_when(sys.argv[0], "finished")
    
    
def get_submitter_cli_args(script_description, arg_names, pwd, validate=None):
    """
    Get and validate all args from command line using argparse.
    :param script_description: String describing the basic purpose of a script,
                               shown when the user runs the script with --help
    :param arg_names: List of strings, each of which names a flag which the
                      user can call the script with
    :param pwd: String which is a valid path to the parent directory of the
                script currently being run
    :param validate: Function to pass output namespace and its parser into to
                     validate all user inputs
    :return: Namespace containing all validated inputted command line arguments
    """
    # Create arg parser, and fill it with parameters shared with other scripts
    parser = initialize_subset_analysis_parser(argparse.ArgumentParser(
        description=script_description
    ), pwd, arg_names)
    
    # This block differs from conan_tools.get_cli_args by adding a new 
    # argument and converting cli_args into a dictionary
    default_jobs = 100
    default_seconds = 10
    parser.add_argument(
        "-sbatch",
        "--sbatch-string",
        required=True,
        help=("String with all of the required parameters to run an SBATCH "
              "command. All other input arguments will be appended to this "
              "string, and then the result will be executed. This flag must "
              "be included with a command to run.")
    )
    parser.add_argument(
        "-seconds",
        "--seconds-to-wait",
        type=valid_whole_number,
        help=("Waiting period between batch job submissions by this script. "
              "Enter the number of seconds to wait, after submitting one "
              "job, to submit the next job. By default, this script will "
              "wait {} seconds.".format(default_seconds))
    )
    parser.add_argument(
        "-print-cmd",
        "--print-command",
        type="store_true",
        help=("Include this flag to print every command that is run to submit "
              "an automated_subset_analysis.py batch job.")
    )
    return vars(validate(parser.parse_args(), parser)
                if validate else parser.parse_args())


def get_asa_options(cli_args):
    """
    :param cli_args: argparse namespace with all validated command-line
                     arguments, all of which are used by this function
    :return: List of some cli_args optional arguments and their values
    """
    asa_optional_args = []
    for arg in get_ASA_arg_names():
        if arg not in (gp_demo.format(1), gp_demo.format(2),
                       "n_analyses", "output", "subset_size"):
            if cli_args[arg]:
                asa_optional_args.append(as_cli_arg(arg))
                if isinstance(cli_args[arg], list):
                    for el in cli_args[arg]:
                        asa_optional_args.append(str(el))
                elif not isinstance(cli_args[arg], bool):
                    asa_optional_args.append(str(cli_args[arg]))
    return asa_optional_args


def submit_batch_job(cli_args):
    """
    Submit batch job to run automated_subset_analysis in parallel
    :param cli_args: argparse namespace with all validated command-line
                     arguments, all of which are used by this function
    :return: N/A
    """
    gp_demo = "group_{}_demo_file"                # ASA required args
    asa_optional_args = get_asa_options(cli_args) # Some ASA optional args

    # Track how long submission took so far to estimate how long it will take
    progress = track_progress(cli_args.n_analyses, len(cli_args.subset_size))
    start_time = now()

    # Call batch command with all needed parameters for SBATCH command and 
    # automated_subset_analysis call
    for i in range(1, cli_args["n_analyses"] + 1):
        for subset_size in cli_args["subset_size"]:

            # Command to run automated_subset_analysis batch job
            cmd = (cli_args["sbatch_string"].split(" ") + [
                cli_args[gp_demo.format(1)],
                cli_args[gp_demo.format(2)],
                "--output",
                os.path.join(cli_args["output"], "output{}".format(i)),
                "--n-analyses", "1",
                "--parallel", PWD,
                "--subset-size", str(subset_size)
            ] + asa_optional_args)
            
            # Print information to inform user about process progress
            update_progress(progress,
                            "submitting automated_subset_analysis batch jobs",
                            1, start_time)
            if cli_args.print_command:
                print("Running {}".format(cmd))

            # Call command and wait before printing the next one
            subprocess.check_call(cmd)
            time.sleep(cli_args.seconds_to_wait)
        
if __name__ == "__main__":
    main()
