#! /usr/bin/env python3

"""
SBATCH job submitter for automated subset analysis
Greg Conan: conan@ohsu.edu
Created 2020-01-03
Updated 2020-02-12
"""

##################################
#
# Script to submit many automated_subset_analysis.py SBATCH jobs for parallel
# processing on the Exacloud server
#
##################################

# Imports
import argparse
from automated_subset_analysis import validate_cli_args
import os
from src.conan_tools import *
import subprocess
import time

# Constants: Demographics and job argument names, automated_subset_analysis dir
GP_DEMO_FILE = "group_{}_demo_file"
JOB_SHORTNAME = "automate"
PWD = get_pwd()

def main():

    # Store and print the date and time when this script started running
    starting_timestamp = get_and_print_timestamp_when(sys.argv[0], "started")

    # Get and validate all command-line arguments from user
    cli_args = get_submitter_cli_args(
        ("Script to run many instances of automated_subset_analysis.py in "
         "parallel."), get_ASA_arg_names(), PWD
    )

    cli_args["sbatch"] = ["sbatch", "--time={}".format(cli_args["time"]), 
                          "--mem=1gb", "-c", "1", "-A", "fnl_lab",
                          os.path.join(PWD, "automated_subset_analysis.py")]

    try:
        submit_batch_jobs(cli_args)
        
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
    default_sleep = 60
    default_time_limit = "04:00:00"
    parser.add_argument(
        "-print-cmd",
        "--print-command",
        action="store_true",
        help=("Include this flag to print every command that is run to submit "
              "an automated_subset_analysis.py batch job.")
    )
    parser.add_argument(
        "-q", "-queue",
        "--queue-max-size",
        type=valid_whole_number,
        default=default_jobs,
        help=("The maximum number of jobs to run simultaneously. By default, "
              "a maximum of {} jobs will run at once.".format(default_jobs))
    )
    parser.add_argument(
        "-sleep",
        "--seconds-between-jobs",
        dest="sleep",
        type=valid_whole_number,
        default=default_sleep,
        help=("Number of seconds to wait between batch job submissions. The "
              "default number is {}.".format(default_sleep))
    )
    parser.add_argument(
        "-time",
        "--job-time-limit",
        dest="time",
        type=valid_time_str,
        default=default_time_limit,
        help=("Time limit for each automated_subset_analysis batch job. The "
              "time limit must be formatted specifically as HH:MM:SS where HH "
              "is hours, MM is minutes, and SS is seconds. {} is the default "
              "time limit.".format(default_time_limit))
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
        if arg not in (GP_DEMO_FILE.format(1), GP_DEMO_FILE.format(2),
                       "n_analyses", "output", "subset_size"):
            if cli_args[arg]:
                asa_optional_args.append(as_cli_arg(arg))
                if isinstance(cli_args[arg], list):
                    for el in cli_args[arg]:
                        asa_optional_args.append(str(el))
                elif not isinstance(cli_args[arg], bool):
                    asa_optional_args.append(str(cli_args[arg]))
    return asa_optional_args


def valid_time_str(in_arg):
    """
    :param in_arg: Object to check if it's a time string in the right format
    :return: True if in_arg is a string representing a time limit in the format
             HH:MM:SS; otherwise False
    """
    try:
        split = in_arg.split(":")
        assert len(split) == 3
        for each_num in split:
            assert each_num.isdigit()
            assert int(each_num) >= 0
        return in_arg
    except (TypeError, AssertionError, ValueError):
        raise argparse.ArgumentTypeError("Invalid time string.")


def count_jobs_running():
    """
    :return: Integer counting how many ASA batch jobs are running right now
    """
    return subprocess.check_output("squeue", universal_newlines=True
                                   ).count(JOB_SHORTNAME)


def get_batch_command(cli_args, out_num, subset_size):
    """
    Get command to run automated_subset_analysis batch job
    :param cli_args: argparse namespace with all validated command-line
                     arguments, all of which are used by this function
    :param out_num: Integer from 1 to cli_args["n_analyses"] representing which
                    analysis this batch command is
    :param subset_size: Integer which is an element of cli_args["subset_size"]
    :return: List of strings which can be called as a command to run an
             automated_subset_analysis batch job
    """
    return (cli_args["sbatch"] + [
        cli_args[GP_DEMO_FILE.format(1)],
        cli_args[GP_DEMO_FILE.format(2)],
        "--output",
        os.path.join(cli_args["output"],"output{}".format(out_num)),
        "--n-analyses", "1",
        "--parallel", PWD,
        "--subset-size", str(subset_size)
    ] + get_asa_options(cli_args))


def submit_batch_jobs(cli_args):
    """
    Submit automated_subset_analysis batch jobs to run in parallel
    :param cli_args: argparse namespace with all validated command-line
                     arguments, all of which are used by this function
    :return: N/A
    """   
    all_jobs_subset_sizes = cli_args["subset_size"] * cli_args["n_analyses"]
    keep_adding_jobs = True
    out_num = 0
    while len(all_jobs_subset_sizes) > 0:
        with open("./submissions.txt", "a+") as infile:  # TODO remove this print?
            infile.write("keep_adding_jobs: {}, len(all_jobs_subset_sizes): "
                         "{}, out_num: {}, running: {}, queue_max: {}\n"
                         .format(keep_adding_jobs, len(all_jobs_subset_sizes),
                                 out_num, count_jobs_running(),
                                 cli_args["queue_max_size"]))
        if keep_adding_jobs:
            if all_jobs_subset_sizes[-1] == cli_args["subset_size"][-1]:
                out_num += 1
            subprocess.check_call(get_batch_command(
                cli_args, out_num, all_jobs_subset_sizes.pop()
            ))
        time.sleep(cli_args["sleep"])
        keep_adding_jobs = count_jobs_running() < cli_args["queue_max_size"]


if __name__ == "__main__":
    main()

