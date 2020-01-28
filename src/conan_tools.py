#! /usr/bin/env python3

"""
Conan Tools
Greg Conan: conan@ohsu.edu
Created 2019-11-26
Updated 2020-01-28
"""

##################################
#
# Common source for utility functions used by multiple scripts in the
# automated_subset_analysis folder/module
#
##################################

# Imports
import argparse
import datetime
import math
import nibabel
import numpy as np
import os
import pandas as pd
import random
from scipy import stats
from scipy.spatial import distance
import socket
import sys
from timeit import Timer

# Constants: Names of common automated_subset_analysis argparse parameters
GP_AV_FILE = "group_{}_avg_file"
GP_DEMO_FILE = "group_{}_demo_file"
GP_MTR_FILE = "matrices_conc_{}"
EXAMPLE_FILE = "example_file"


def add_default_avg_matr_path_to(cli_args, gp_num, parser, default):
    """
    Get paths to average matrices for each group, if not given
    :param cli_args: argparse namespace with most needed command-line
                     arguments. This function uses the --matrices-conc-{} and
                     --example-file arguments, but only the former is required
    :param gp_num: Integer which is the group number
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :param default: String naming the demographic data spreadsheet's column
                    with all of the paths to matrix files
    :return: cli_args, but with the group_{}_avg_file argument for group gp_num
    """
    gp_avg_arg = GP_AV_FILE.format(gp_num)
    if not getattr(cli_args, gp_avg_arg, None):
        matr_conc = getattr(cli_args, GP_MTR_FILE.format(gp_num))
        if not matr_conc:
            parser.error("Please provide the {} argument or the {} argument. "
                         .format(as_cli_arg(GP_AV_FILE.format(gp_num)),
                                 as_cli_arg(GP_MTR_FILE.format(gp_num))))
        example = getattr(cli_args, EXAMPLE_FILE, None)
        if not example:
            with open(matr_conc) as matr_conc_file_obj:
                example = matr_conc_file_obj.readline().strip()
        else:
            parser.error("Please provide the {} argument or the "
                         "--example-file argument.".format(as_cli_arg(
                            GP_MTR_FILE.format(gp_num)
                         )))
        setattr(cli_args, gp_avg_arg, os.path.join(cli_args.output, "".join((
            os.path.splitext(os.path.basename(matr_conc))[0],
            "_AVG", get_2_exts_of(example)
        ))))
    return cli_args


def as_cli_arg(arg_str, gp_num=None):
    """
    :param arg_str: String naming a stored argument taken from the command line
    :param gp_num: Integer to be interpolated into arg_str, to differentiate
                   the group that arg_str represents data about
    :return: String which is the command-line argument form of arg_str
    """
    if gp_num is not None:
        arg_str = arg_str.format(gp_num)
    return "--" + arg_str.replace("_", "-")


def chdir_to(folder):
    """
    Change the current working directory to folder unless already there
    :param folder: String that is a valid path to the directory to chdir to
    :return: N/A
    """
    abs_dir = os.path.abspath(folder)
    if os.path.abspath(os.getcwd()) != abs_dir:
        os.chdir(abs_dir)


def count_digits_of(a_num):
    """
    :a_num: Numeric value
    :return: Integer which is 10 to the power of the number of digits in a_num
    """
    return 10**int(math.log10(a_num))


def default_vis_titles():
    """
    :return: Dictionary mapping the keyword identifying each visualization to
             the text to use as that visualization's default title
    """
    return {"sub1_sub2": "Correlations Between Average Subsets",
            "sub1_all2": "Group 1 Subset to Group 2 Correlation",
            "sub2_all1": "Group 1 to Group 2 Subset Correlation",
            None: "Correlation Between Unknown Groups"}


def drop_nan_rows_from(group, group_num, nan_threshold, parser):
    """
    Check how many rows have NaN values. If it's under the threshold,
    then drop all of them with NaNs. Otherwise, raise a parser error
    :param group: pandas.DataFrame with all data from group demographics file
    :param nan_threshold: Float between 0 and 1 which represents the percentage
                          of cells with NaNs below which a row will be rejected
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: group without any of its rows that contained NaNs, if group has
             few enough NaNs to pass below the nan_threshold
    """
    nans = sum(group.apply(lambda x: 1 if x.isnull().any() else 0))
    nans_percent = (nans / len(group.index))
    print("Percentage of rows in group {} demographics file with NaN "
          "values: {:.1%}".format(group_num, nans_percent))
    if nans_percent < nan_threshold:
        group = group.dropna(how="any", axis=0)
        print("Dropping all {} rows with NaN values.".format(nans))
    else:
        parser.error("Too many rows with NaN values at threshold of "
                     "{:.1%}".format(nan_threshold))
    return group


def extract_subject_id_from(path):
    """
    Accepts the path to a subject's matrix file, and gets that subject's ID
    :param path: String which is a valid path to a CIFTI2 matrix file
    :return: String which is the subject ID which was in path
    """
    # Where in path does subject ID start, where in path would the "_" between 
    # "NDAR" and "INV" be, and where in path does the subject ID end
    sub_id_pos = path.find("NDAR")  
    id_mid_pos = sub_id_pos + 4
    id_end_pos = sub_id_pos + 15

    return path[sub_id_pos:id_end_pos] if path[id_mid_pos] == "_" else (
        "_".join((path[sub_id_pos:id_mid_pos], path[id_mid_pos:id_end_pos]))
    )


def fit_strings_to_width(elements, widths):
    """
    Given a list of elements to print as a row, and an equally sized list of
    integers where each element is the width of the column of the element at
    the same index, return one string which is a row where each element has its
    correct width.
    :param elements: List of objects that can all be turned into a strings
    :param widths: List of integers, the sum of which is the total output width
    :return: String which joins everything in elements together in order
    """
    result = []
    num_els = len(elements)
    try:
        assert num_els == len(widths)
        for i in range(num_els):
            string = str(elements[i])
            result.append(string + " " * (widths[i] - len(string)))
        return "".join(result)
    except AssertionError:
        print("Error: Inconsistent lengths in fit_strings_to_widths. "
              "\nElements: {}\nWidths:{}".format(elements, widths))


def get_2_exts_of(path):
    """
    :param path: String representing a file path with two extensions
    :return: String with just those two extensions (like ".dscalar.nii")
    """
    splitext = os.path.splitext(path)
    return os.path.splitext(splitext[0])[-1] + splitext[-1]


def get_and_print_timestamp_when(script, completion):
    """
    Print and return a string showing the exact date and time when the current 
    running script reached a certain part of its process
    :param script: String which is name of script that started/finished
    :param completion: String which is a past tense verb describing what script 
                       did at the timestamp, like "started" or "finished"
    :return: String with an easily human-readable message showing when a script
             either started or finished
    """
    timestamp = "\n{} {} at {}.".format(
        script, completion, now().strftime("%H:%M:%S on %b %d, %Y")
    )
    print(timestamp)
    return timestamp


def get_ASA_arg_names():
    """
    :return: List of strings where each names an argument accepted by 
             automated_subset_analysis.py as a cli_args argparse parameter
    """
    return [GP_DEMO_FILE.format(1), GP_DEMO_FILE.format(2), "axis_font_size", 
            "columns", "euclidean", "fill", GP_AV_FILE.format(1),
            GP_AV_FILE.format(2), GP_MTR_FILE.format(1), GP_MTR_FILE.format(2),
            "hide_legend", "marker_size", "n_analyses", "nan_threshold",
            "no_matching", "only_make_graphs", "output", "parallel", 
            "skip_subset_generation", "subset_size", "graph_title", 
            "title_font_size", "y_range", "inverse_fisher_z"]


def get_average_matrix(subset, paths_col, fisherz=None):
    """
    Build and return a matrix averaging over every matrix in the subset by
    adding every matrix to a running total (to avoid loading every matrix in 
    the subset into memory at once)
    :param subset: pandas.DataFrame with a column of paths to matrix files
    :param paths_col: String naming the column of paths to matrix files
    :param fisherz: Object which is truthy to do an inverse Fisher-Z 
                    transformation on the matrix data, or falsy not to
    :return: numpy.ndarray averaging all of the subset's subjects' matrices
    """
    subset_size = len(subset.index)

    # Get one matrix file to initialize the running total matrix
    subject_matrix_paths = subset[paths_col].iteritems()
    running_total = load_matrix_from(next(subject_matrix_paths)[1])

    # Iteratively add every matrix to the running total
    counter = 1  # (Manual counter because iteritems has no len)
    just_printed = 0
    for subj in subject_matrix_paths:
        counter += 1
        running_total = np.add(running_total, load_matrix_from(subj[1])) 

        # Display progress to user <= 100 times
        percent_done = counter / subset_size
        if (math.floor(100 * percent_done)
                != math.floor((100 * just_printed - 1) / subset_size)):
            print("Remaining matrices: {}. Progress: {:.1%} done."
                  .format(subset_size - counter, percent_done))
            just_printed = counter + 1
 
    # Divide running total matrix by number of matrices, then return it
    divisor_matrix = np.ndarray(running_total.shape)
    divisor_matrix.fill(subset_size)
    return np.divide(running_total, divisor_matrix)


def get_cli_args(script_description, arg_names, pwd, validate=None):
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
    return (validate(parser.parse_args(), parser)
            if validate else parser.parse_args())


def get_confidence_interval(series, confidence=0.95):
    """
    :param series: pandas.Series filled with numeric data
    :param confidence: Float that represents the percentage for the confidence
                       interval; 0.95 (95%) by default
    :return: Tuple of series's confidence interval: (lower bound, upper bound)
    """
    # Confidence interval = mean +/- (standard deviation * Z-score / sqrt(n))
    distance_from_mean = (series.std() * stats.norm.ppf(1-(1-confidence)/2)
                          / np.sqrt(len(series.index)))
    series_mean = series.mean()
    return series_mean - distance_from_mean, series_mean + distance_from_mean



def get_family_info(subject, family_vars):
    """
    Get the demographic information identifying a subject wrt their family
    :param subject: pandas.Series representing one subject who has a value
                    for each variable name in family_vars
    :param family_vars: List of strings where each names a demographic
                        variable identifying a subject wrt their family
    :return: Dictionary mapping each string in family_vars to the value of
             the variable by that name in subject
    """
    return {prop: int(subject[prop]) for prop in family_vars}


def get_group_demographics(cli_args, gp_num, gp_demo_str, parser):
    """
    :param cli_args: argparse namespace with most needed command-line arguments
    :param gp_num: Integer which is the group number
    :param gp_demo_str: String naming the .csv file path cli_args attribute
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: pandas.DataFrame with all of the group's demographic data read in
             from the .csv file, dropping all NaN values
    """
    return drop_nan_rows_from(pd.read_csv(
        getattr(cli_args, gp_demo_str.format(gp_num)),
        na_values=(" ", 777, 999)
    ), gp_num, cli_args.nan_threshold, parser)


def get_ID_string(subject_series, id_name):
    """
    :param sub: pandas.Series representing a subject
    :param id_name: String naming the subject ID demographic variable
    :return: String which is sub's subject ID
    """
    id_str = subject_series.get(id_name)
    return (id_str if isinstance(id_str, str)
            else id_str.to_string(index=False))


def get_pwd():
    """
    :return: String which is a valid path to the directory containing the file
             which is currently being run
    """
    try:
        pwd = os.path.dirname(os.path.abspath(sys.argv[0]))
    except OSError:
        pwd = os.getcwd()
    print("Running {} while in directory {}".format(sys.argv[0], pwd))
    return pwd


def get_subset_of(group, subset_size):
    """
    Randomly select, validate, and return a subset of a given size from group
    :param group: pandas.DataFrame++ which is an entire group of subjects
    :param subset_size: Integer which is the amount of subjects to randomly 
                        select from group to put into a subset
    :return: pandas.DataFrame which is a valid subset of group
    """
    # Constants: names of demographic variables for IDs of subject and family
    ID = "id_redcap"
    FAMILY = {"FAM": "rel_family_id", "REL": "rel_relationship",
              "GRP": "rel_group_id"}

    # Variables: Randomly selected subset from the overall group, and its 
    # members who have siblings not in the subset
    subset = group.sample(n=subset_size)
    subs_missing_sibs = set()

    # Local function to add one subject who has missing siblings to a set
    def get_invalid_subset_member(sub, subs_missing_sibs):
        """
        If a subject has a sibling not in the subset, add that subject to a set
        :param sub: pandas.Series representing a subject, 1 row of a DataFrame
        :param subs_missing_sibs: Set of strings; all are invalid subjects' IDs
        :return N/A:
        """
        if has_missing_siblings(sub, subset, FAMILY):
            subs_missing_sibs.add(sub[ID])

    # Local function to run get_invalid_subset_member with a default set
    # of subjects missing siblings
    def collect_invalid_members_of(subset_df, missing=subs_missing_sibs):
        """
        Put all subjects who have missing siblings into one set (missing)
        :param subset_df: pandas.DataFrame with all subjects in the subset
        :param missing: Set of strings which are all the subject IDs of the 
                        subjects in subset_df with siblings not in subset_df
        :return N/A:
        """
        subset_df.apply(lambda sub: get_invalid_subset_member(sub, missing),
                        axis="columns")

    # Any subset with one member of a specific family but not others is invalid
    collect_invalid_members_of(subset)
    return make_subset_valid(subs_missing_sibs, collect_invalid_members_of,
                             subset, group, FAMILY["REL"], ID, subset_size)


def get_which_str_in_filename(filename, possible_names):
    """
    :param filename: String that is the name of a readable file
    :param possible_names: List of names which might be in filename
    :return: Element of possible_names which is in the filename string
    """
    found = None
    while not found and len(possible_names) > 0:
        next_name = possible_names.pop()
        if next_name in filename:
            found = next_name
    return found


def has_missing_siblings(subj_row, subset, ids):
    """
    Check whether a subject in a subset has a sibling not in the subset.
    :param subj_row: pandas.Series representing one subject
    :param subset: pandas.DataFrame with one subject per row
    :param ids: Dictionary with the labels of relevant columns in subset:
                REL (relationship), FAM (family ID), and GRP (sibling group ID)
    :return: True if subj_row represents a subject who has siblings not in
             subset; otherwise False
    """
    # If subject has no family in the group, they have no missing siblings
    mis_sibs = False    
    subj_rel = int(subj_row[ids["REL"]])
    if subj_rel > 0:

        # Sibling missing if subject has family in the group but not the subset
        family = subset[subset[ids["FAM"]] == subj_row[ids["FAM"]]]
        if len(family[ids["GRP"]]) < 2:
            mis_sibs = True

        # Sibling missing if a twin has no family in the subset
        elif subj_rel == 2:
            twins = int(family[
                family[ids["GRP"]] == subj_row[ids["GRP"]]
            ][ids["GRP"]].tolist()[0])
            if subj_rel != twins:
                mis_sibs = True

        # Sibling missing if a triplet has <3 family members in subset
        elif len(family[family[ids["REL"]] == subj_row[ids["REL"]]]) <3:
            mis_sibs = True
    return mis_sibs


class HiddenPrints:
    """
    Object used to temporarily suppress printing text. This code is copied from
    https://stackoverflow.com/a/45669280
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def initialize_subset_analysis_parser(parser, pwd, to_add):
    """
    :param parser: argparse.ArgumentParser without any parameters
    :param pwd: String which is a path to the present working directory
    :param to_add: List of strings; each names a parameter to add to parser
    :return: parser, but with all of the parameters used by a script
    """
    # Default values of input parameters, to use if not specified by user
    default_cols = ['demo_comb_income_v2b', 'demo_ed_v2', 'demo_prnt_ed_v2b', 
                    'demo_sex_v2b', 'ehi_y_ss_scoreb', 'interview_age',
                    'medhx_9a', 'race_ethnicity', 'rel_relationship',
                    'site_id_l']
    choices_fill = ["all", "confidence_interval"]
    default_continuous_vars = ['demo_prnt_ed_v2b', 'interview_age', 
                               'rel_group_id', 'rel_relationship']
    default_dims = (1, 2)
    default_euclid_vals = [-0.44897407617376806, 3.7026679337563486]
    default_marker_size = 5
    default_n_analyses = 1
    default_nan_threshold = 0.1
    default_out_dir = os.path.join(pwd, "data")
    default_subset_size = [50, 100, 200, 300, 400, 500]
    default_text_size_axis = 30
    default_text_size_title = 40

    # Help messages used by multiple input parameters
    help_demo_file = ("Path to a .csv file containing all demographic "
                      "information about the subjects in group {}.")
    help_font_size = (
        "Font size of {0} text in visualization. Enter a positive integer for "
        "this argument. If it is excluded, then the default {0} font size "
        "will be {1}."
    )
    help_group_avg_file = (
        "Path to a .nii file containing the average matrix for group {0}. "
        "By default, this path will be to group{0}_10min_mean.pconn.nii "
        "file in this script's parent folder or the output folder."
    )
    help_matrices_conc = (
        "Path to a .conc file containing only a list of valid paths to group "
        "{0} matrix files. This flag is only needed if your group {0} "
        "demographics .csv file either does not have a column labeled "
        "'pconn10min' with paths to matrix files, or if it does include that "
        "column but you want to use different paths."
    )

    # Required: Full list of group 1 demographic data
    def group_1_demo_file():
        parser.add_argument(
            GP_DEMO_FILE.format(1),
            type=valid_readable_file,
            help=help_demo_file.format(1)
        )

    # Required: Full list of group 2 demographic data
    def group_2_demo_file():
        parser.add_argument(
            GP_DEMO_FILE.format(2),
            type=valid_readable_file,
            help=help_demo_file.format(2)
        )

    # Optional: Font size of axis text in visualization
    def axis_font_size():
        parser.add_argument(
            "-axis",
            "--axis-font-size",
            default=default_text_size_axis,
            type=valid_whole_number,
            help=help_font_size.format("axis", default_text_size_axis)
        )

    # Optional: Specify which columns to match on
    def columns():
        parser.add_argument(
            "-c",
            "--columns",
            nargs="*",
            default=default_cols,
            help=("Demographic variables to match subsets on. By default, the "
                  "subsets will be matched on these variables: {}"
                  .format(default_cols))
        )

    # Optional: Names of columns which are continuous, not categorical, vars
    def continuous_variables():
        parser.add_argument(
            "-con-vars",
            "--continuous-variables",
            dest="con_cols",
            default=default_continuous_vars,
            help=("All names of columns in the demographics .csv file which "
                  "have continuous instead of categorical data.")
        )

    # Optional: Custom logarithmic function for Euclidean distance threshold
    def euclidean():
        parser.add_argument(
            "-e",
            "--euclidean",
            nargs=2,
            type=float,
            metavar=("LOG_TERM", "CONSTANT_TERM"),
            default=default_euclid_vals,
            help=("Enter custom constant values for the logarithmic equation "
                  "used to estimate the Euclidean distance threshold used to "
                  "determine whether a randomly selected subset of a group is "
                  "significantly different from the other group. This must be "
                  "2 numbers, such that the first is the constant for the "
                  "logarithmic term and the second is the constant added to "
                  "it. For example, using the argument --euclidean -5 10 will "
                  "produce the equation y = -5*ln(x) + 10. Default values: {}"
                  .format(default_euclid_vals))
        )

    # Optional: One matrix file to get extension from & to save out avg matrix
    def example_file():
        parser.add_argument(
            "-ex",
            as_cli_arg(EXAMPLE_FILE),
            type=valid_readable_file,
            help="Path to one of the matrix files inputted to the average."
        )

    # Optional: Choose what to fill in the visualization
    def fill():
        parser.add_argument(
            "-f",
            "--fill",
            choices=choices_fill,
            default=choices_fill[1],
            help=("Choose which data to shade in the visualization. Choose {} "
                  "to shade in the area within the minimum and maximum "
                  "correlations in the dataset. Choose {} to only shade in the "
                  "95 percent confidence interval of the data. By default, "
                  "--fill will be {}.".format(*choices_fill, choices_fill[1]))
        )

    def graph_title():
        parser.add_argument(
            "-title",
            "--graph-title",
            help=("Include this argument with a custom title to show it at "
                  "the top of the output visualizations. Otherwise, each "
                  "visualization will have one of these default titles: {}"
                  .format(", ".join(default_vis_titles().values())))
        )

    # Optional: Path to average matrix .pconn file for group 1
    def group_1_avg_file():
        parser.add_argument(
            "-avg1",
            as_cli_arg(GP_AV_FILE, 1),
            help=help_group_avg_file.format(1)
        )

    # Optional: Path to average matrix .pconn file for group 2
    def group_2_avg_file():
        parser.add_argument(
            "-avg2",
            as_cli_arg(GP_AV_FILE, 2),
            help=help_group_avg_file.format(2)
        )

    def hide_legend():
        parser.add_argument(
            "-hide",
            "--hide-legend",
            action="store_true",
            help=("Include this flag to prevent the visualization from "
                  "displaying a legend in the corner.")
        )

    # Optional: Do inverse Fisher-Z transformation (hyperbolic tangent) on data
    def inverse_fisher_z():
        parser.add_argument(
            "-z",
            "--inverse-fisher-z",
            action="store_true",
            help=("Do inverse Fisher-Z transformation on the matrices "
                  "imported from the data matrix files before getting "
                  "correlations.")
        )

    def marker_size():
        parser.add_argument(
            "-marker",
            "--marker-size",
            type=valid_whole_number,
            default=default_marker_size,
            help=("Positive integer, the size of each data point in the "
                  "visualization that this script will create. The default "
                  "marker size is {}.".format(default_marker_size))
        )

    # Optional: .conc file with paths to group 1 matrix files
    def matrices_conc_1():
        parser.add_argument(
            "-conc1",
            as_cli_arg(GP_MTR_FILE, 1),
            type=valid_conc_file,
            help=help_matrices_conc.format(1)
        )

    # Optional: .conc file with paths to group 2 matrix files
    def matrices_conc_2():
        parser.add_argument(
            "-conc2",
            as_cli_arg(GP_MTR_FILE, 2),
            type=valid_conc_file,
            help=help_matrices_conc.format(2)
        )

    # Optional: Number of subsets
    def n_analyses():
            parser.add_argument(
            "-n",
            "--n-analyses",
            default=default_n_analyses,
            type=valid_whole_number,
            help="Number of times to generate and analyze a pair of subsets."
        )

    # Optional: NaN threshold
    def nan_threshold():

        def float_between_0_and_1(val):
            """
            :param val: Object to check, then throw an error if it is invalid
            :return: val if it is a float between 0 and 1 (otherwise invalid)
            """
            float_val = float(val)
            if 0 < float_val < 1:
                return float_val
            else:
                parser.error("NaN threshold must be a number between 0 and 1.")

        parser.add_argument(
            "-nan",
            "--nan-threshold",
            type=float_between_0_and_1,
            default=default_nan_threshold,
            help=("Enter a number between 0 and 1. If the percentage of rows "
                  "with NaN values in the data for either demographic file is "
                  "greater than the --nan-threshold, then the script will raise "
                  "an error. Otherwise, the script will drop every row "
                  "containing a NaN. The default NaN threshold is {}."
                  .format(default_nan_threshold))
        )

    # Optional: Don't demographically match subsets
    def no_matching():
        parser.add_argument(
            "--no-matching",
            action="store_true",
            help=("By default, each subset of one group is matched to the "
                  "other group so their demographics show no statistically "
                  "significant difference. Include this flag to skip matching "
                  "and include subsets regardless of demographic differences.")
        )

    # Optional: Import data from .csv files already created by this script
    # instead of creating new ones
    def only_make_graphs():
        parser.add_argument(
            "-graphs",
            "--only-make-graphs",
            nargs="+",
            type=valid_readable_file,
            metavar="SUBSET-CORRELATIONS-CSV",
            help=("Include this flag to import data from a .csv file already "
                  "made by this script instead of making a new one, just to "
                  "make a graph visualization of already-existing data. If "
                  "this flag is included, it must be a path to at least one "
                  "readable .csv file(s) with 2 columns: 'Subjects' (the "
                  "number of subjects in each subset) and 'Correlation' (the "
                  "correlation between each randomly generated subset in that "
                  "pair).")
        )

    # Optional: Output folder
    def output():
        parser.add_argument(
            "-out",
            "--output",
            type=valid_output_dir,
            default=default_out_dir,
            help=("Path to folder to save subset data in. The default folder "
                  "is " + default_out_dir)
        )
        
    # Optional: Parallel processing boolean flag
    def parallel():
        parser.add_argument(
            "--parallel",
            type=valid_readable_file,
            help=("Include this argument if you are running the "
                  "automated_subset_analysis script many times in parallel. "
                  "It should be a valid path to the directory containing the "
                  "automated_subset_analysis script. This argument will let "
                  "multiple running instances of this script append "
                  "correlation values to the same output .csv files.")
        )

    # Optional: Number of subjects in each subset pair
    def subset_size():
        parser.add_argument(
            "-subjects",
            "--subset-size",
            nargs="+",
            default=default_subset_size,
            type=valid_whole_number,
            help=("Number of subjects to include in subsets. Include a list "
                  "of whole numbers to generate subsets pairs of different "
                  "sizes. By default, the subset sizes will be {}."
                  .format(default_subset_size))
        )

    # Optional: Font size of title text in visualization
    def title_font_size():
        parser.add_argument(
            "-tfont",
            "--title-font-size",
            default=default_text_size_title,
            type=valid_whole_number,
            help=help_font_size.format("title", default_text_size_title)
        )

    # Optional: Make averages/correlations from existing subsets
    def skip_subset_generation():
        parser.add_argument(
            "-skip",
            "--skip-subset-generation",
            nargs="?",
            const="output",
            metavar="SUBSETS-DIR-PATH",
            help=("Include this flag to calculate correlations and create the "
                  "visualization using existing subsets instead of making new "
                  "ones. By default, the subsets to use will be assumed to "
                  "exist in the --output folder. If you want to load subsets "
                  "from a different folder, add the path to this flag as a "
                  "parameter.")
        )

    # Optional: Custom data range for visualization
    def y_range():
        parser.add_argument(
            "-y",
            "--y-range",
            nargs=2,
            type=float,
            metavar=("Y-MIN", "Y-MAX"),
            help=("Range of y_axis in visualization. By default, this script "
                  "will automatically set the y-axis boundaries to show all of "
                  "the correlation values and nothing else. If this argument "
                  "is used, then it should be two floating-point numbers: the "
                  "minimum and maximum values to be shown on the "
                  "visualizations' y-axes.")
        )

    # Add all requested parameters to parser
    for cli_arg in to_add:
        locals()[cli_arg]()
    return parser


def load_matrix_from(matrix_path):
    """
    :param matrix_path: String which is the absolute path to a matrix file
    :return: numpy.ndarray matrix of data from the matrix file
    """
    return np.array(nibabel.cifti2.load(matrix_path).get_data().tolist())


def look_for_file(orig_path, cli_args, arg_name, pwd, parser):
    """
    Given a filepath which might not be valid, look for files with its same
    base name in other likely locations (output dir and PWD).
    :param orig_path: String representing a file path
    :param cli_args: argparse namespace with all given command-line arguments
    :param arg_name: String which is the name of the cli_args attribute which
                     will store the found_path value
    :param pwd: String which is a path to the present working directory
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: String which is a valid path to the file being looked for
    """
    found_path = None
    if os.access(orig_path, os.R_OK):
        found_path = orig_path

    # Try to find file in other likely locations (output and PWD)
    else:
        filename = os.path.basename(orig_path) 
        for possible_dir in (cli_args.output, pwd):
            found_path = os.path.join(possible_dir, filename)
            if os.access(found_path, os.R_OK):
                break
            else:
                found_path = None

        # Raise parser error if file is not found
        if not found_path:
            parser.error("{} is not a valid path for the --{} argument."
                         .format(orig_path, arg_name))
    return found_path


def make_subset_valid(subs_missing_sibs, collect, subset, group, rel, id_var,
                      subset_size):
    """
    Given a subset, keep checking and swapping its members around until no 
    member of the subset has siblings outside of the subset
    :param subs_missing_sibs: Set of strings; each is an invalid subject's ID
    :param collect: Function to get all of the subjects in subset who have
                    missing siblings
    :param subset: pandas.DataFrame with a randomly selected subset of group
    :param group: pandas.DataFrame which is the whole group of which subset_IDs 
                  represents a subset 
    :param rel: String naming the demographic variable describing each
                subject's relationship to their sibling(s)
    :param id_var: String naming the subject ID demographic variable
    :return: subset, but without any subjects who have missing siblings
    """
    # Variables: All subject IDs in subset, subject IDs of all subjects with
    # siblings in group, counter flag to prevent while loop getting stuck
    subset_IDs = set(subset[id_var].unique())
    all_sibling_IDs = set(group[group[rel] > 0][id_var].tolist())
    stuck_if_at_10 = 0

    # Local function to shuffle a subset if the swapping process gets stuck
    def shuffle_subsets(subset_IDs, shuffle_out, other_sibs, stuck_if):
        """
        :param subset_IDs: Set of strings (all of the subset's subject IDs)
        :param shuffle_out: Set of strings (all subject IDs to remove)
        :param other_sibs: Set of strings (all subject IDs not in subset)
        :param stuck_if: Integer used as counter to trigger shuffling
        :return: stuck_if, but reset to 0 if a subset with fewer invalid 
                 subjects was created by shuffling some subjects around
        """
        subset_IDs = shuffle_out_subset_of(subset_IDs, shuffle_out, other_sibs)
        shuffled_subset = group[group[id_var].isin(subset_IDs)]
        missing_post_shuffle = set()
        collect(subset, missing=missing_post_shuffle)
        return (0 if len(missing_post_shuffle) < len(subs_missing_sibs)
                else stuck_if)

    # Keep swapping missing siblings in/out of the subset until the subset has 
    # no subjects with siblings not in the subset, then return the valid subset
    while len(subs_missing_sibs) > 0:
        mis_sibs_before = len(subs_missing_sibs)
        subset_IDs = shuffle_out_subset_of(
            subset_IDs, subs_missing_sibs,
            all_sibling_IDs.difference(subset_IDs)
        )
        subset = group[group[id_var].isin(subset_IDs)]
        subs_missing_sibs = set()
        collect(subset)
        if len(subs_missing_sibs) == mis_sibs_before:
            stuck_if_at_10 += 1
        shuffle_out = set(random.sample(subset_IDs, 10))
        stuck_if_at_10 = shuffle_subsets(
            subset_IDs, shuffle_out,
            all_sibling_IDs.difference(subset_IDs, shuffle_out), stuck_if_at_10
        )

    # If an extra subject was added to the subset, then remove them
    if len(subset) > subset_size:
        to_remove = subset.sample(n=1)
        while to_remove[rel].iloc[0] > 0:
            to_remove = subset.sample(n=1)
        subset = subset.drop(to_remove.index)

    return subset 


def natural_log(x_val, coefs):
    """
    :param x_val: Float, the natural log of which will be used in the equation
    :param coefs: Tuple of 2 floats: 1) log coefficient, 2) constant term
    :return: Float, the output of plugging inputs into a logarithmic equation
    """
    return coefs[0] * np.log(x_val) + coefs[1]


def now():
    """
    :return: Current date and time as a datetime.datetime object
    """
    return datetime.datetime.now()


def print_col_headers_and_get_widths(headers):
    """
    Prints the headers of every column in the chart of generated subsets' chi-
    squared values and other details to the command line
    :param headers: List of strings which are the labels for each column
    :return: List of integers which are all column widths, in order
    """
    col_widths = [len(hdr) + 2 for hdr in headers]
    print("Subsets randomly generated:\n"
          + fit_strings_to_width(headers, col_widths))
    return col_widths


def randomly_select_subset(group, group_n, sub_n, diff_group,
                           cli_args, loop_check_fn, eu_threshold=None):
    """
    Randomly select subsets of size sub_n from group until finding one which
    doesn't significantly differ demographically from diff_group. Then return
    either the Euclidean distance between them, or the subset.
    :param group: pandas.DataFrame with all demographic data from one group
    :param group_n: Integer which is the group number of the group df
    :param sub_n: Integer which is how many subjects to put in the subset
    :param diff_group: pandas.DataFrame with all demographic data from the
                       other group
    :param cli_args: argparse namespace with all given command-line arguments
    :param loop_check_fn: Function to determine when a subset is found
    :param eu_threshold: Float which is the Euclidean distance threshold
    :return: Either pandas.DataFrame or float, depending on eu_threshold:
    - If eu_threshold provided, then return a randomly generated valid subset
    - Otherwise return the Euclidean distance between the subset and diff_group
    """
    # Get information about group total dataframe for chi square comparison,
    # excluding any column where all subjects have the same value (GROUP) or
    # nearly all have different values (rel_family_id)
    columns = group.select_dtypes(include=["number"]).drop([
        "GROUP", "rel_family_id"
    ], axis="columns")
    group_avgs = [diff_group[col].mean(skipna=True)
                  for col in columns.columns.tolist()]

    # Randomly pick subsets until finding 1 that demographically represents    
    # the overall group    
    loops = 0
    keep_looping = True
    col_widths = None
    while keep_looping:
        loops += 1

        # Generate subset of group and get its columns' averages as a list
        subset = (group.sample(n=sub_n) if cli_args.no_matching
                  else get_subset_of(group, sub_n))
        sub_avgs = [subset[col].mean(skipna=True)
                    for col in columns.columns.tolist()]

        # If any column in subset averages has mean of N/A, throw it out
        if any(np.isnan(el) for el in sub_avgs):
            print("Too much missing data, skipping to next.")
        else:

            # Get Euclidean distance
            eu_dist = distance.euclidean(group_avgs, sub_avgs)

            # Check significance of Euclidean distance, and show user
            extra_args = [eu_threshold] if eu_threshold else [
                subset, diff_group, columns.columns.tolist(),
                cli_args.con_cols, col_widths
            ]
            keep_looping = loop_check_fn(loops, eu_dist, group_n,
                                         *extra_args)
            if loops == 2 and keep_looping and not eu_threshold:
                col_widths = print_col_headers_and_get_widths((
                    "Subset", "Group", "Statistic", "P-Value",
                    "Significant Difference", "Euclidean Distance"
                ))

    return subset if eu_threshold else eu_dist


def read_file_into_list(filepath):
    """
    Given the path to a file, read all of the file's lines as a list
    :param filepath: Path to file to read
    :return: List of all lines in the file at filepath
    """
    with open(filepath, "r") as infile:
        return [line.strip() for line in infile] 
        

def rename_exacloud_path(path):
    """
    Convert valid Exacloud path to valid Rushmore path
    :param path: String representing a valid file path on Exacloud server
    :return:     String representing a valid file path on Rushmore server
    """
    return (path.replace("mnt/rose/shared", "home/exacloud/lustre1/fnl_lab")
            if "exa" in socket.gethostname() else 
            path.replace("home/exacloud/lustre1/fnl_lab", "mnt/rose/shared"))
    

def shuffle_out_subset_of(subset, to_shuffle_out, exclusive_pool):
    """
    Given a subset, remove all members of to_shuffle_out and add an equal
    number of randomly selected members from exclusive_pool
    :param subset: Set of strings
    :param to_shuffle_out: Set of strings in subset to remove from it
    :param exclusive_pool: Set of strings which are not in subset
    :return: subset, but without any member of to_shuffle_out and with some
             random new members of exclusive_pool; subset's len doesn't change
    """
    if len(exclusive_pool) < len(to_shuffle_out):
        sys.exit("Not enough subjects in population to randomly select a "
                 "sample with {} subjects, because {} subjects cannot be "
                 "randomly swapped out from a pool of {} subjects".format(
                     len(subset), len(to_shuffle_out), len(exclusive_pool)
                 ))
    return subset.difference(to_shuffle_out).union(set(
        random.sample(exclusive_pool, len(to_shuffle_out))
    ))


def stringify_timedelta(timestamp):
    """
    Turn a timedelta object, or a number representing one, into a string
    :param timestamp: datetime.timedelta object, or float representing a number
                      of seconds
    :return: string showing a timedelta in H:MM:SS.SSSSSS format
    """
    return (str(timestamp) if isinstance(timestamp, datetime.timedelta)
            else str(datetime.timedelta(seconds=timestamp)))
    

def time_since(start_time):
    """
    :param start_time: datetime.datetime representing when something started
    :return: Amount of time between start_time and now, as a datetime.timedelta
    """
    return now() - start_time


def timeit_comparison(fns_dict, runs):
    """
    Given a dictionary of functions and their arguments, and a number of times 
    to run them, print and return how long it takes to run them that many times
    :param fns_dict: Dictionary mapping a function object to a list of values
                     accepted by that function as arguments/parameters
    :param runs: Integer which is the number of times to run each function
    :return: Dictionary mapping the name of each function in fns_dict to how
             long it took to run that function {runs} times
    """
    times = {}
    for fn, fn_args in fns_dict.items():
        times[fn.__name__] = timeit_fn(fn, runs, fn_args)
    for fn_name, fn_time in times.items():
        print("Time taken by function {}: {}".format(fn_name, fn_time))
    return times


def timeit_fn(fn, runs, fn_args):
    """
    Return how long it takes to run a function a given number of times
    :param fn: Callable function object which accepts fn_args as parameters
    :param runs: Integer which is the number of times to run fn
    :param fn_args: List of arguments passed to fn
    :return: Float which is how long it took to run fn(*fn_args) the number of 
             times specified in the runs parameter
    """
    return Timer(lambda: fn(*fn_args)).timeit(runs)


def track_progress(cli_args):
    """
    Initialize dictionary to track how long making average matrices will take
    :return: Dictionary tracking progress so far, with 3 string:integer pairs
        {"matrices_left":  Matrices left to make at some point in the process
         "total_matrices": Matrices to make, total
         "seconds_taken": Seconds the script spent making matrices so far}
    """
    matrices_to_do = cli_args.n_analyses * sum(cli_args.subset_size)
    return {"total_matrices": matrices_to_do, "matrices_left": matrices_to_do,
            "seconds_taken": 0}


def update_progress(progress, doing_what, subset_size, start_time):
    """
    Gets and prints how far the script has gotten processing multiple subsets
    :param progress: Dictionary keeping track of progress so far
    :param doing_what: String describing what the script did in the time
                       represented by each cell of seconds_taken
    :param subset_size: Integer which is how many subjects are in the subset
    :param start_time: datetime.datetime when this subset started processing
    :return: progress, after updating to account for the last step finishing
    """
    progress["matrices_left"] -= subset_size
    progress["seconds_taken"] += time_since(start_time).total_seconds()
    seconds_left = (progress["matrices_left"] * (progress["seconds_taken"] / (
                    progress["total_matrices"] - progress["matrices_left"])))
    print("       Time spent so far {0} (H:MM:SS.SSSSSS) is {1}\n"
          "Estimated time remaining {0} (H:MM:SS.SSSSSS) is {2}\n"
          .format(doing_what, stringify_timedelta(progress["seconds_taken"]),
                  stringify_timedelta(seconds_left)))
    return progress


def valid_conc_file(path):
    """
    Throw argparse exception unless parameter is a valid .conc file path
    :param path: String to check if it represents a valid filename
    :return: String representing a valid path to a readable .conc file
    """
    try:
        assert os.path.splitext(path)[1] == ".conc"
        return valid_readable_file(path)
    except (AssertionError, OSError, TypeError):
        raise argparse.ArgumentTypeError("{} is not a .conc file".format(path))


def valid_output_dir(path):
    """
    Try to create a directory to write files into at the given path, and throw 
    argparse exception if that fails
    :param path: String which is a valid (not necessarily real) folder path
    :return: String which is a validated absolute path to real writeable folder
    """
    try:
        out_dir = os.path.abspath(path)
        os.makedirs(out_dir, exist_ok=True)
        assert os.access(out_dir, os.W_OK)
        return out_dir
    except (AssertionError, OSError, TypeError):
        raise argparse.ArgumentTypeError("Cannot write to directory at {}"
                                         .format(path))


def valid_readable_file(path):
    """
    Throw exception unless parameter is a valid readable filename string. This
    is used instead of argparse.FileType("r") because the latter leaves an open
    file handle, which has caused problems.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename
    """
    try:
        assert os.access(path, os.R_OK)
        return os.path.abspath(path) 
    except (AssertionError, OSError, TypeError):
        raise argparse.ArgumentTypeError("Cannot read file at {}".format(path))


def valid_whole_number(to_validate):
    """
    Throw argparse exception unless to_validate is an integer greater than 0
    :param to_validate: Object to test whether it is an integer greater than 0
    :return: to_validate if it is an integer greater than 0
    """
    try:
        to_validate = int(to_validate)
        assert to_validate > 0
        return to_validate
    except (AssertionError, TypeError):
        raise argparse.ArgumentTypeError("{} is not a positive integer."
                                         .format(to_validate))
