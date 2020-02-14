#! /usr/bin/env python3

"""
Automated subset selection and analysis for ABCD resource paper
Greg Conan: conan@ohsu.edu
Created 2019-09-17
Updated 2020-02-14
"""

##################################
#
# Script to randomly select pairs of subsets of data of given sizes, find the
# correlation between the average matrices of each subset, and plot all of
# those correlations for each subset size
#
##################################

# Standard Python imports
import argparse
import numpy as np
import os
import pandas as pd
import plotly
import pprint
import re
import sys

# Ensure that this script can find its local imports if parallel processing
if "--parallel" in sys.argv:
    parallel_flag_pos = 0
    while sys.argv[parallel_flag_pos] != "--parallel":
        parallel_flag_pos += 1
    sys.path.append(os.path.abspath(sys.argv[parallel_flag_pos + 1]))

# Local custom imports
from src.conan_tools import *

# Constants: Default demographic variable names and PWD
DEFAULT_DEM_VAR_MATR = "matrix_file"
DEFAULT_DEM_VAR_PCONNS = "pconn10min"
DEFAULT_DEM_VAR_SUBJID = "id_redcap"
GP_DEMO_STR = "group_{}_demo"
PWD = get_pwd()


def main():

    # Store and print the date and time when this script started running
    starting_timestamp = get_and_print_timestamp_when(sys.argv[0], "started")

    # Get and validate all command-line arguments from user
    cli_args = get_cli_args(
        ("Script to randomly select pairs of subsets of data of given sizes, "
         "find the correlation between the average matrices of each subset, "
         "and plot all of those correlations for each subset size."),
        get_ASA_arg_names(), PWD, validate_cli_args
    )

    # If user said to skip the subset generation and use pre-existing subset
    # correlations, then make pd.DataFrame to visualize those
    try:
        if cli_args.only_make_graphs:
            only_make_graphs(cli_args)
        else:

            # Make and save all subsets and their correlations, unless said to
            # get pre-existing subsets' correlations instead
            get_subs = (skip_subset_generation
                        if getattr(cli_args, "skip_subset_generation", None)
                        else save_and_get_all_subsets)
            all_subsets = get_subs(cli_args, "subset_{}_with_{}_subjects.csv")

            # Go to output dir to make visualizations of subset correlations
            if not cli_args.parallel:
                chdir_to(cli_args.output)
            for correls_df_name, correls_df in get_correl_dataframes(
                    all_subsets, cli_args).items():
                if not cli_args.parallel:
                    make_visualization(correls_df, cli_args, correls_df_name)

    except Exception as e:
        get_and_print_timestamp_when(sys.argv[0], "crashed")
        raise e

    # Print the date and time when this script started and finished running
    print(starting_timestamp)
    get_and_print_timestamp_when(sys.argv[0], "finished")


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all given command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    try:
        # If user said to only make graphs, validate that the path they gave is
        # a file with correlations instead of a path to a directory
        if cli_args.only_make_graphs:
            for correl_file in cli_args.only_make_graphs:
                if os.path.isdir(correl_file):
                    parser.error(correl_file + " is a directory, not a file. "
                                 "Please enter the name of a readable file as "
                                 "the --only-make-graphs argument.")
        else:
            # If user said to get correlations from existing subsets, get the
            # path to the directory to save correlations in
            path_skip_sub = getattr(cli_args, "skip_subset_generation", None)
            if path_skip_sub == "output":
                cli_args.skip_subset_generation = cli_args.output
            elif path_skip_sub:
                valid_readable_file(path_skip_sub)  # Raise error unless valid

            # For each group, get the path to the directory with its .pconn
            # files, and the path to the file containing its demographic data
            cli_args = add_pconn_paths_to(cli_args, [1, 2], parser)

        return cli_args
    except (OSError, argparse.ArgumentTypeError) as e:
        parser.error(str(e))


def add_pconn_paths_to(cli_args, group_nums, parser):
    """
    Get necessary paths to .pconn.nii files and add them to the cli_args object
    :param cli_args: argparse namespace with most needed command-line arguments
    :param group_nums: List of whole numbers; each element is a group number
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: cli_args, but with all needed arguments
    """
    for gp_num in group_nums:

        # Import all of this group's demographics; drop empty columns
        group_demo_str = GP_DEMO_STR.format(gp_num)
        demographics = get_group_demographics(cli_args, gp_num,
                                              group_demo_str + "_file", parser)

        # Add paths to matrices if the user gave paths in separate .conc file
        matr_conc = getattr(cli_args, "matrices_conc_{}".format(gp_num), None)
        setattr(cli_args, group_demo_str, demographics if not matr_conc else
                replace_paths_column(demographics, matr_conc))

        # Get average matrix for each group
        group_avg_file_str = "group_{}_avg_file".format(gp_num)
        cli_args = add_default_avg_matr_path_to(cli_args, gp_num, parser,
                                                DEFAULT_DEM_VAR_PCONNS)
        try:
            valid_readable_file(getattr(cli_args, group_avg_file_str))
        except argparse.ArgumentTypeError as e:
            parser.error(str(e))
        setattr(cli_args, "group_{}_avg".format(gp_num),
                load_matrix_from(getattr(cli_args, group_avg_file_str)))
    return cli_args


def replace_paths_column(demographics, matr_conc):
    """
    :param demographics: pandas.DataFrame with demographic information,
                         including 0 or 1 column(s) with paths to matrix files
    :param matr_conc: String which is a valid path to a .conc file with paths
                      to matrix files
    :return: demographics, but replacing its previous paths column with the
             contents of the matr_conc file
    """
    matrix_paths = pd.read_csv(matr_conc, converters={
        0: rename_exacloud_path
    }, sep="\n", header=None).rename(columns=lambda x: DEFAULT_DEM_VAR_MATR)
    matrix_paths[DEFAULT_DEM_VAR_SUBJID] = matrix_paths.apply(lambda x: (
        extract_subject_id_from(x.loc[DEFAULT_DEM_VAR_MATR])
    ), axis="columns")
    if DEFAULT_DEM_VAR_PCONNS in demographics:
        demographics = demographics.drop(DEFAULT_DEM_VAR_PCONNS, axis=1)
    return demographics.merge(matrix_paths, on=DEFAULT_DEM_VAR_SUBJID)
    

def only_make_graphs(cli_args):
    """
    If user said to, skip subset generation and skip correlation calculation,
    then make visualizations from already-existing correlation files
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the --output and --only_make_graphs 
                     arguments, and passes cli_args to make_visualization.
    :return: N/A
    """
    chdir_to(cli_args.output)
    titles = [k for k in default_vis_titles().keys() if k]
    for correls_csv in cli_args.only_make_graphs:
        make_visualization(pd.read_csv(correls_csv), cli_args,
                           get_which_str_in_filename(correls_csv, titles))


def skip_subset_generation(cli_args, subsets_file_name):
    """
    Get all subsets that have already been generated if user skipped subset
    generation on this run. If no subset files exist, then terminate the script
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the --group_1_demo, --group_2_demo, and
                     --skip_subset_generation arguments.
    :param subsets_file_name: String with the format of subset file names
    :return: List of dictionaries where each has two elements mapping a group
             number to a subset of that group, and also has an element mapping 
             the string "subset_size" to the number of subjects in both subsets
    """
    all_subsets = []  # Return value: List of subsets

    # Local function to get a group's subset's demographics data
    def get_demographics_of_subset_of_gp(gp_num):
        """
        :param gp_num: String with the group number, either "1" or "2"       
        :return: pandas.Series with all of a subset's data from group gp_num
        """
        gp_demo = getattr(cli_args, GP_DEMO_STR.format(gp_num))
        return gp_demo[gp_demo[DEFAULT_DEM_VAR_SUBJID
                               ].isin(subset_df.pop(gp_num).tolist())]

    # Get every subset based on cli_args, excluding other files in the same dir
    subset_name_parts = subsets_file_name.split("{}")
    for file_name in os.listdir(cli_args.skip_subset_generation):
        subset_csv = os.path.join(cli_args.skip_subset_generation, file_name)
        if is_subset_csv(subset_csv, subset_name_parts, cli_args.n_analyses):

            # Read subset from file into pandas.DataFrame 
            subset_df = pd.read_csv(subset_csv)
            if ("1" in subset_df and "2" in subset_df 
                    and len(subset_df.index) in cli_args.subset_size):
                all_subsets.append({1: get_demographics_of_subset_of_gp("1"),
                                    2: get_demographics_of_subset_of_gp("2"),
                                    "subset_size": len(subset_df.index)})
    if len(all_subsets) == 0:
        raise FileNotFoundError("No subsets found at {}"
                                .format(cli_args.skip_subset_generation))
    return all_subsets


def is_subset_csv(path, subset_filename_parts, n_analyses):
    """
    Check if a path points to a subset .csv file made by this script within the
    --n-analyses variable
    :param path: String which should be a valid path to a readable file
    :param subsets_filename_parts: List of strings which each have part of the
                                   format of subset file names
    :param n_analyses: Integer which is the --n-analyses argument value
    :return: True if path is to a readable .csv file following this script's 
             subset file naming conventions, where the analysis number <= 
             n_analyses, with two columns labeled '1' and '2'; otherwise False
    """
    try:
        path = valid_readable_file(path)
        assert os.path.splitext(path)[1] == ".csv"
        with open(path, "r") as infile:
            row_1 = infile.readline().strip().split(",")
        name = os.path.basename(path)
        match = re.search(r"(\d+)", name)
        result = (len(row_1) == 2 and row_1[0] == "1" and row_1[1] == "2"
                  and all(part in name for part in subset_filename_parts)
                  and match and (int(match.group()) <= n_analyses))
    except (OSError, argparse.ArgumentTypeError, AssertionError):
        result = False
    return result


def save_and_get_all_subsets(cli_args, subsets_file_name):
    """
    Randomly generate all pairs of subsets, save each of them to a .csv file,
    and then return dictionaries with pairs of subsets
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the --n-analyses, --output, & --subset-size 
                     arguments; and passes cli_args to randomly_select_subset.
    :param subsets_file_name: String with the format of subset file names
    :return: List of dictionaries, each of which maps a subset's group number
             to the subset for one pair of subsets
    """
    all_subsets = []  # Return value: List of subsets
    progress = track_progress(cli_args.n_analyses, cli_args.subset_size)

    # Get average correlation from user-defined number of pairs of average
    # matrices of randomly generated subsets
    for i in range(cli_args.n_analyses):
        for subset_size in cli_args.subset_size:
            start_time = now()
            print("Making randomly selected subset pair {} out of {} with {} "
                  "subjects.".format(i + 1, cli_args.n_analyses, subset_size))

            # Iterate over both groups to get a subset of each
            subsets = {1: None, 2: None}
            print("Estimated Euclidean distance threshold for statistical "
                  "significance for subset with {} subjects: {}".format(
                      subset_size, natural_log(subset_size, cli_args.euclidean)
                  ))

            # Select subsets and make average matrices from .pconns for each
            for group_n in subsets.keys():
                subsets[group_n] = randomly_select_subset(
                    getattr(cli_args, GP_DEMO_STR.format(group_n)),
                    group_n, subset_size,
                    getattr(cli_args,
                            GP_DEMO_STR.format(1 if group_n == 2 else 2)),
                    cli_args, check_keep_looping,
                    natural_log(subset_size, cli_args.euclidean)
                )

            # Save randomly generated subsets
            save_subsets(subsets, cli_args.output, i + 1,
                         subsets_file_name, DEFAULT_DEM_VAR_SUBJID)
            subsets["subset_size"] = subset_size
            all_subsets.append(subsets)

            # Print how long this has taken, and roughly how much time is left
            progress = update_progress(progress, "making subsets", subset_size,
                                       start_time)
    return all_subsets


def check_keep_looping(loops, eu_dist, gp_n, eu_threshold):
    """ 
    If a Euclidean distance threshold was given, use it to determine whether to
    stop randomly generating subsets
    :param loops: Integer which is how many random subsets have been generated
    :param eu_dist: Float which is the Euclidean distance between the latest 
                    randomly generated subset and the total group it is from
    :param gp_n: Integer which is the total group's group number
    :param eu_threshold: Float above which a Euclidean distance between two 
                         groups is declared to show a significant difference
    :return: True if another subset should be randomly made; otherwise False
    """
    keep_looping = eu_dist >= eu_threshold
    if not keep_looping or not loops % count_digits_of(loops):  
        print("{} subsets of group {} randomly generated.".format(loops, gp_n))
    return keep_looping


def save_subsets(subs, output_dir, sub_num, subsets_file_name, sub_id_col):
    """
    Given a pair of subsets, save them to a .csv file
    :param subs: Dictionary mapping each subset's group number to that subset
    :param output_dir: String path to directory to save subset .csv files into
    :param sub_num: Integer that is arbitrary unique number so that multiple
                    to prevent a group's subsets' filenames from conflicting
    :param subsets_file_name: String with the format of subset file names
    :param sub_id_col: String naming the subject ID demographic variable/column
    :return: N/A
    """
    # Local function to get all subject ID strings in a group's subset
    def extract_subj_ids(gp_num):
        """
        :param gp_num: Integer which is the number of the group to get IDs from
        :return: pandas.Series with all of group gp_num's subset's subject IDs
        """
        return subs[gp_num][sub_id_col].reset_index(drop=True)

    # Save the pair of subsets to .csv
    to_save = pd.DataFrame({1: extract_subj_ids(1), 2: extract_subj_ids(2)})
    to_save.to_csv(os.path.join(
        output_dir, subsets_file_name.format(sub_num, len(to_save.index))
    ), index=False)


def get_correl_dataframes(all_subsets, cli_args):
    """
    :param all_subsets: List of dictionaries, each of which maps both group
                        numbers to a randomly selected subset of that group
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the --group_1_avg, --group_2_avg, and 
                     --output arguments. It also passes cli_args to 
                     get_avg_matrices_of_subsets and get_correls_between.
    :return: Dictionary of 3 string:pandas.DataFrame pairs such that each
             DataFrame has a column of subset sizes and one of correlations:
       {sub1_sub2: Correlations between the average matrices of both subsets
        sub1_all2: Correls between group 1 subset avg matrix and group 2 total
        sub2_all1: Correls between group 2 subset avg matrix and group 1 total}
    """
    # Return value: Dict of correlation lists to become pandas.DataFrames
    correl_lists = {subset_id: [] for subset_id in
                    default_vis_titles().keys() if subset_id}

    # Keep track of how long this function takes, to show the user during loop
    progress = track_progress(cli_args.n_analyses, cli_args.subset_size)

    # Get each pair of average matrices, their correlation with each other, and
    # each one's correlation with the other group's average matrix
    for s in range(len(all_subsets)):
        start_time = now()
        sub_pair = all_subsets[s]
        sub1_avg, sub2_avg = get_avg_matrices_of_subsets(sub_pair.copy(),
                                                         cli_args).values()

        # If data is 2-dimensional, flatten it to make it 1-dimensional
        subsets = {"sub1": sub1_avg, "all1": cli_args.group_1_avg,
                   "sub2": sub2_avg, "all2": cli_args.group_2_avg}
        for set_name, set_avg in subsets.items():
            subsets[set_name] = set_avg.flatten()

        # Get and show all subset correlations; put them in the dict to return
        subset_size = sub_pair.pop("subset_size")
        correl_lists = get_sub_pair_correls(subset_size, correl_lists, subsets)

        # Print how long this has taken, and roughly how much time is left
        progress = update_progress(progress, "making average matrices",
                                   subset_size, start_time)

    return {name: save_correlations_and_get_df(
                correls, "correlations_{}.csv".format(name), cli_args
            ) for name, correls in correl_lists.items()}


def get_avg_matrices_of_subsets(subsets_dict, cli_args):
    """
    Make and return an average matrix for each subset in subsets_dict
    :param subsets_dict: Dict mapping the group number to its subset DataFrame
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses --inverse_fisher_z
    :return: Dictionary matching each group's number to its subset's average
             matrix from .pconn files of all its subjects
    """
    avg_matrices = {}  # Return value: Average matrices of both groups
    sub_size = subsets_dict.pop("subset_size")  # Number of subjects per group

    # Get all data from .pconn files of every subject in the subset
    for sub_num, subset in subsets_dict.items():
        col = (DEFAULT_DEM_VAR_MATR if getattr(
                   cli_args, "matrices_conc_{}".format(sub_num), None
               ) else DEFAULT_DEM_VAR_PCONNS) 
        print("Making average matrix for group {} subset.".format(sub_num))
        avg_matrices[sub_num] = get_average_matrix(subset, col, cli_args)
        print("Group {} subset's average matrix: \n{}"
              .format(sub_num, avg_matrices[sub_num]))
    return avg_matrices


def get_sub_pair_correls(subset_size, correl_lists, subsets):
    """
    Get, print, and collect the correlations between a subset pair
    :param subset_size: Integer which is how many subjects are in each subset
    :param correl_lists: List of dictionaries holding all subset correlations
    :param subsets: Dictionary mapping label string to average subsets
    :return: correl_lists, plus the correlations between the subset pair
    """
    for df_name in correl_lists.keys():
        sub_keys = df_name.split("_")
        params = [subsets[sub] for sub in sub_keys]
        params.append(subset_size)
        correl_lists[df_name].append(get_correls_between(*params))
        title = ("group {0}'s subset and group {1}".format(
                     sub_keys[0][-1], sub_keys[1][-1]
                 ) if "all" in df_name else "both subsets")
        print("Correlations between average matrices of {}:\n{}"
              .format(title, pprint.pformat(correl_lists[df_name])))
    return correl_lists


def get_correls_between(arr1, arr2, num_subjects):
    """
    :param arr1: np.ndarray with only numeric values
    :param arr2: np.ndarray with only numeric values
    :param num_subjects: Integer, number of subjects in each group
    :return: Dictionary mapping "Subjects" to num_subjects and mapping
             "Correlation" to the correlation between arr1 and arr2 
    """
    return {"Subjects": num_subjects,
            "Correlation": np.corrcoef(arr1, arr2)[0, 1]}


def save_correlations_and_get_df(correls_list, correl_file_name, cli_args):
    """
    Save correlations to .csv file
    :param correls_list: List of dictionaries where each dictionary contains
    "Subjects" and "Correlation" as a key with a numerical value
    :param correl_file_name: Path to the file to save correlations into
    :return: pandas.DataFrame with a header row ("Subjects", "Correlation"),
    subset sizes in one column, and correlations in its second column
    """
    out_dir = (os.path.dirname(cli_args.output)
               if getattr(cli_args, "parallel", False) else cli_args.output)
    append_rows_to_file(correls_list, os.path.join(out_dir, correl_file_name))
    return pd.DataFrame(correls_list, columns=["Subjects", "Correlation"])
    
    
def append_rows_to_file(correls_list, filename):
    """
    Append a list of pairs of subset sizes and correlations to a file
    :param correls_list: List of dictionaries where each dictionary contains
                         "Subjects" and "Correlation" as a key with a number
    :param filename: String that is a valid path to an existing .csv file
    :return: N/A
    """
    mode = "a" if os.access(filename, os.W_OK) else "w+"
    with open(filename, mode) as out:
        if mode == "w+":
            out.write("Subjects,Correlation\n")
        for row in correls_list:
            out.write("{},{}\n".format(row["Subjects"], row["Correlation"]))


def make_visualization(correls_df, cli_args, corr_df_name):
    """
    Create a graph visualization from correlation data
    :param correls_df: pandas.DataFrame with one column titled "Subjects" and
                       another titled "Correlations"; both have numeric values
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the --axis_font_size, --fill, --output, 
                     --title_font_size, and --y_range arguments.
    :param corr_df_name: String identifying the visualization to create
    :return: N/A
    """
    # Colors as RGBA strings, for the visualizations' lines and shading
    clear = "rgba(0,0,0,0)"
    def red(opacity):
        return "rgba(255,0,0,{})".format(opacity)

    # Make scatter plot mapping subset size to pairwise correlations
    scatter_plot = plotly.graph_objs.Scatter(
        x=correls_df["Subjects"], y=correls_df["Correlation"], mode="markers",
        name="All correlations", line_color=red(1),
        marker={"size": cli_args.marker_size}
    )

    # Add average lines to plot using averages of each subset size
    avgs = correls_df.groupby(["Subjects"]).agg(lambda x: 
                                                x.unique().sum() / x.nunique())
    last_avg = float(avgs.tail(1).values)
    avgs_plot = plotly.graph_objs.Scatter(
        x=avgs.index.values, y=avgs["Correlation"], mode="lines",
        name="Average correlations"
    )

    # Get and display the visualization title, and averages
    vis_title = (cli_args.graph_title if cli_args.graph_title
                 else default_vis_titles()[corr_df_name])
    vis_file = "".join(vis_title.replace("<br>", "_").split()) + ".html"
    i = 1
    while os.access(vis_file, os.R_OK):
        i += 1
        vis_file = "{}_{}.html".format("".join(vis_title.split()), i)
    print("{}:\n{}".format(vis_title, avgs))

    # Add upper & lower bounds (all data or confidence intervals) of shaded
    # area to plot as lines
    bounds = get_shaded_area_bounds(correls_df, cli_args.fill)
    bounds_params = ({"showlegend": False} if cli_args.fill == "all" else
                     {"name": "95% confidence interval", "showlegend": True})
    lower_plot = plotly.graph_objs.Scatter(
        x=avgs.index.values, y=bounds[0], fill="tonexty", line_color=clear, 
        fillcolor=red(0.2), **bounds_params
    )
    upper_plot = plotly.graph_objs.Scatter(x=avgs.index.values, y=bounds[1],
                                           line_color=clear, showlegend=False)

    # Show image & export it as a .png file, but suppress the massive block of
    # text that plotly.offline.plot() would normally print to the command line
    with HiddenPrints():
        plotly.offline.init_notebook_mode()
        plotly.offline.plot({
            "data": [scatter_plot, avgs_plot, upper_plot, lower_plot],
            "layout": get_plot_layout(get_layout_args(cli_args, correls_df,
                                                      vis_title, last_avg))
        }, image="png", filename=vis_file)
    print("Saving .png image of {} offline using browser.".format(vis_file))


def get_layout_args(cli_args, correls_df, title, last_avg):
    """
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the --axis_font_size, --title_font_size, and 
                     --y_range arguments.
    :param correls_df: pandas.DataFrame with one column titled "Subjects" and
                       another titled "Correlations"; both have numeric values
    :param title: String to put at the top of the visualization
    :return: Tuple of all input arguments for make_visualization function
    """
    # Use default y-axis range or one chosen by user
    y_axis_min, y_axis_max = cli_args.y_range if cli_args.y_range else (
        correls_df["Correlation"].min(), correls_df["Correlation"].max()
    )
    return (y_axis_min, y_axis_max, title, last_avg, cli_args.title_font_size,
            cli_args.axis_font_size, correls_df["Subjects"].mean(),
            not cli_args.hide_legend)


def get_shaded_area_bounds(all_data_df, to_fill):
    """
    :param all_data_df: pandas.DataFrame with two columns of numeric values, 
                        labeled "Subjects" and "Correlation"
    :param to_fill: String that is either "all" to shade the area between the
                    min and max correlations for each number of subjects in
                    all_data_df or "confidence_interval" to shade the area 
                    within a 95% confidence interval of the correlations for
                    each number of subjects in all_data_df
    :return: Tuple of 2 pandas.Series objects where the first is the lower
             boundary of the shaded area and the second is the upper boundary
    """
    # Local function to aggregate correlation values
    def aggregate_data(quantile):
        """
        :param quantile: Integer representing the quantile of data to aggregate
        :return: pandas.Series with all aggregated correlation values
        """
        return all_data_df.groupby(["Subjects"]).agg(
            lambda x: x.quantile(quantile)
        )["Correlation"]

    if to_fill == "confidence_interval":
        intervals = all_data_df.groupby(["Subjects"]).agg(
            lambda x: get_confidence_interval(x)
        )
        intervals = pd.DataFrame(intervals["Correlation"].tolist(),
                                 index=intervals.index)
        result = (intervals[0], intervals[1])
    elif to_fill == "all":
        result = (aggregate_data(0), aggregate_data(1))
    else:
        raise ValueError("Invalid value for --fill parameter.")
    return result


def get_plot_layout(all_args):
    """
    Return all format settings for creating a pretty plot visualization.
    :param layout_args: Tuple of all arguments needed for visualization format
    :return: Nested dictionary containing all needed plot attributes
    """
    # Local variables from unpacked list of args to use in layout: 
    #   Lowest and highest y-values to show, graph title, last y-value to show,
    #   graph title and axis title font sizes, average subset size
    y_min, y_max, title, last_y, title_size, axis_font, x_avg, show = all_args

    # Others: RGBA colors as well as space buffer above y_max and below y_min
    black = "rgb(0, 0, 0)"
    white = "rgb(255, 255, 255)"
    y_range_step = (y_max - y_min) / 10

    def get_axis_layout(title_txt, **kwargs):
        """
        Get all of the parameters that "xaxis" and "yaxis" have in common for 
        get_plot_layout in order to avoid redundant code.
        :param title: String of text to be displayed at the top of the graph
        :param kwargs: Dictionary containing all of the parameters that "xaxis"
                       and "yaxis" do not have in common
        :return: Dictionary with all parameters needed for formatting an axis
        """
        result = {"title": {"font": {"size": title_size}, "text": title_txt},
                  "tickcolor": black, "ticklen": 15, "ticks": "outside",
                  "tickfont": {"size": axis_font}, "tickwidth": 2,
                  "showline": True, "linecolor": black, "linewidth": 2}
        result.update(kwargs)
        return result

    return {"title": {"text": title, "x": 0.5, "xanchor": "center",
                      "font": {"size": title_size}},
            "paper_bgcolor": white, "plot_bgcolor": white, "showlegend": show,
            "legend": {"font": {"size": axis_font}, "x": 0.5,
                       # Place the legend in white space away from the graph
                       "y": 0.9 if last_y < ((y_max + y_min) / 2) else 0.1},
            "xaxis": get_axis_layout(
                title_txt="Sample Size (n)", tick0=0, tickmode="linear",
                dtick=count_digits_of(x_avg)
            ),
            "yaxis": get_axis_layout(
                title_txt="Correlation (r)", tickmode="auto", nticks=5,
                range=(y_min - y_range_step, y_max + y_range_step) 
            )}


if __name__ == '__main__':
    main()

