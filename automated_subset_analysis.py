#! /usr/bin/env python

"""
Automated subset selection and analysis for ABCD resource paper
Greg Conan: conan@ohsu.edu
Created 2019-09-17
Updated 2020-11-25
"""

##################################
#
# Script to randomly select pairs of subsets of data of given sizes, find the
# correlation between the average matrices of each subset, and plot all of
# those correlations for each subset size
#
##################################
import argparse
import numpy as np
import os
import pandas as pd
import plotly
import pprint
import subprocess
import sys

# Ensure that this script can find its local imports if parallel processing
if "--parallel" in sys.argv:
    parallel_flag_pos = 0
    while sys.argv[parallel_flag_pos] != "--parallel":
        parallel_flag_pos += 1
    sys.path.append(os.path.abspath(sys.argv[parallel_flag_pos + 1]))

# Local custom imports
from src.conan_tools import *

# Constants: Default demographic variable names, PWD, and MATLAB plot code info
DEFAULT_CF = "correlations_{}.csv"
DEFAULT_DEM_VAR_MATR = "matrix_file"
DEFAULT_DEM_VAR_PCONNS = "pconn10min"
DEFAULT_DEM_VAR_SUBJID = "id_redcap"
DEFAULT_MAT_CSV = "matlab_parameters_{}.csv"
DEFAULT_OPACITY = 0.3
GP_DEMO_STR = "group_{}_demo"
PWD = get_pwd()
MATLAB_PLOTTER = os.path.join(PWD, "src", "run_MultiShadedBars.sh")


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
        
        # Get and save correlations or effect sizes, then make visualizations
        correl_fn = (save_subset_effect_size_matrices if cli_args.calculate ==
                     "effect-size" else get_correl_dataframes)
        for corr_df_name, corr_df in correl_fn(all_subsets, cli_args).items():
            if not cli_args.parallel:  # must call correl_fn even if --parallel
                plot_args = (corr_df, cli_args, corr_df_name)
                make_visualization(get_traces_to_plot(*plot_args), *plot_args)

    print(starting_timestamp)  # Print when script started and finished running
    get_and_print_timestamp_when(sys.argv[0], "finished")


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all given command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    try:  # Validate that --y-range has exactly 2 numbers, min and max
        if (cli_args.y_range and isinstance(cli_args.y_range[0], float)
                and len(cli_args.y_range) != 2):
            parser.error("--y-range must be only two numbers.")
        
        # Validate RGBA & threshold values for MATLAB plotting code
        if getattr(cli_args, "matlab_rgba", None):
            if not (3 <= len(cli_args.matlab_rgba) <= 5):
                parser.error("--matlab-rgba must only have 3-5 parameters.")
            elif len(cli_args.matlab_rgba) == 3:
                cli_args.matlab_rgba += [DEFAULT_OPACITY]

        if cli_args.only_make_graphs:  # Check that each trace has a title
            if not cli_args.trace_titles:
                cli_args.trace_titles = ["correlations"] 
            if len(cli_args.trace_titles) != len(cli_args.only_make_graphs):
                parser.error("--trace-titles must have exactly as many input "
                             "arguments as --only-make-graphs.") 
            for correl_file in cli_args.only_make_graphs:
                if os.path.isdir(correl_file):  # Check correlation file path
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

    except (OSError, argparse.ArgumentTypeError) as e:
        parser.error(str(e))
    return cli_args


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
        matr_conc = getattr(cli_args, GP_MTR_FILE.format(gp_num), None)
        setattr(cli_args, group_demo_str, demographics if not matr_conc else
                replace_paths_column(demographics, matr_conc))

        # Get average matrix for each group
        gp_av_arg = GP_AV_FILE.format(gp_num)
        cli_args = add_and_validate_gp_file(cli_args, gp_num, parser,
                                            DEFAULT_DEM_VAR_PCONNS, gp_av_arg)
        setattr(cli_args, "group_{}_avg".format(gp_num),
                load_matrix_from(getattr(cli_args, gp_av_arg)))

        # Validate group variance matrix file paths
        if cli_args.calculate in ("variance", "effect-size"):
            gp_v_f = GP_VAR_FILE.format(gp_num)
            cli_args = add_and_validate_gp_file(cli_args, gp_num, parser,
                                                DEFAULT_DEM_VAR_PCONNS, gp_v_f)
            fname = ("group_{}_variance_matrix{}"
                     .format(gp_num, get_2_exts_of(GP_AV_FILE.format(gp_num))))
            if not getattr(cli_args, gp_v_f, None):
                setattr(cli_args, gp_v_f, os.path.join(cli_args.output, fname))
    return cli_args


def replace_paths_column(demographics, matr_conc):
    """
    :param demographics: pandas.DataFrame with demographic information,
                         including 0 or 1 column(s) with paths to matrix files
    :param matr_conc: String, valid path to a .conc file with matrix file paths
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
    :param cli_args: argparse namespace with --output and --only-make-graphs 
                     command-line arguments to pass to make_visualization.
    :return: N/A
    """
    chdir_to(cli_args.output)
    titles = [k for k in default_vis_titles().keys() if k]
    to_plot = []
    for i in range(len(cli_args.only_make_graphs)):
        correls_csv = cli_args.only_make_graphs[i]
        plot_ttl = cli_args.trace_titles[i]
        corr_df = pd.read_csv(correls_csv)
        corr_df = corr_df.loc[corr_df["Subjects"].isin(cli_args.subset_size)]
        to_plot += get_traces_to_plot(corr_df, cli_args, plot_ttl, i)
    make_visualization(to_plot, corr_df, cli_args,
                       get_which_str_in_filename(correls_csv, titles))


def skip_subset_generation(cli_args, subsets_file_name):
    """
    Get all subsets that have already been generated if user skipped subset
    generation on this run. If no subset files exist, then terminate the script
    :param cli_args: argparse namespace with all command-line arguments
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
    # List of subsets to return, and progress tracker to estimate time left
    all_subsets = []
    progress = track_progress(cli_args.n_analyses, cli_args.subset_size)

    # Get average correlation from user-defined number of pairs of average
    # matrices of randomly generated subsets
    for i in range(cli_args.n_analyses):
        for sub_n in cli_args.subset_size:
            start_time = now()
            print("Making randomly selected subset pair {} out of {} with {} "
                  "subjects.".format(i + 1, cli_args.n_analyses, sub_n))

            # Select subsets from each group and get each's average matrix
            subsets = {1: None, 2: None}
            print("Estimated Euclidean distance threshold for statistical "
                  "significance for subset with {} subjects: {}"
                  .format(sub_n, natural_log(sub_n, cli_args.euclidean)))
            for group_n in subsets.keys():
                subsets[group_n] = randomly_select_subset(
                    getattr(cli_args, GP_DEMO_STR.format(group_n)),
                    group_n, sub_n, getattr(
                        cli_args, GP_DEMO_STR.format(other_group_n(group_n))
                    ), cli_args, check_keep_looping,
                    natural_log(sub_n, cli_args.euclidean)
                )

            # Save randomly generated subsets
            save_subsets(subsets, cli_args.output, i + 1,
                         subsets_file_name, DEFAULT_DEM_VAR_SUBJID)
            subsets["subset_size"] = sub_n
            all_subsets.append(subsets)
            progress = update_progress(progress, "making subsets", sub_n,
                                       start_time)
    return all_subsets


def check_keep_looping(loops, eu_dist, gp_n, eu_threshold):
    """ 
    Use Euclidean threshold to check whether to stop randomly making subsets
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


def save_subset_effect_size_matrices(all_subsets, cli_args):
    """
    For every subset, calculate its effect size matrix and save it to a .csv
    :param all_subsets: List including every pair of subsets
    :param cli_args: argparse namespace with all command-line arguments
    :return: Dictionary mapping name of subset pair to its effect sizes
    """
    # Get both groups' average matrices, group sizes, and variance matrices
    group_averages, gp_vars, gp_sizes = get_groups_avg_var_and_size(
        cli_args, get_matr_file_col_name(cli_args)
    )

    # Get effect size matrices for all subset pairs
    progress = track_progress(cli_args.n_analyses, cli_args.subset_size)
    effect_sizes = dict()
    VAR = "variance"
    AVG = "average"
    for s in range(len(all_subsets)):
        start_time = now()
        sub_pair = all_subsets[s]
        subset_size = sub_pair.pop("subset_size")
        subset_data = dict()

        # Effect size matrices for each subset to the other group
        for gp_num, subset in sub_pair.items():

            # Get group's and subset's variance matrices
            other_gp = other_group_n(gp_num)
            subset_data[gp_num] = get_average_matrix(
                subset, get_matr_file_col_name(cli_args), cli_args
            )
            subset_data[gp_num][VAR] = pd.DataFrame(subset_data[gp_num][VAR])

            # Make subset-to-group effect sizes matrix
            gp_ids = "sub{}_all{}".format(gp_num, other_gp)
            effect_sizes = touch_dict(effect_sizes, gp_ids, [])
            effect_sizes[gp_ids] = make_effect_size_matrix(
                subset_data[gp_num][VAR], pd.DataFrame(gp_vars[other_gp]),
                subset_size, gp_sizes[other_gp], subset_data[gp_num][AVG],
                group_averages[other_gp], effect_sizes[gp_ids]
            )

        # Make subset-to-subset effect sizes matrix
        gp_ids = "sub1_sub2"
        effect_sizes = touch_dict(effect_sizes, gp_ids, [])
        effect_sizes[gp_ids] = make_effect_size_matrix(
            subset_data[1][VAR], subset_data[2][VAR], subset_size, subset_size,
            subset_data[1][AVG], subset_data[2][AVG], effect_sizes[gp_ids]
        )
        progress = update_progress(progress, "getting effect sizes",
                                   subset_size, start_time)
    
    # Save subset-size-to-effect-size lists as .csv files then return them
    for gp_id, dicts_list in effect_sizes.items():
        effect_sizes[gp_id] = save_correlations_and_get_df(
            cli_args, dicts_list, "effect_sizes_{}.csv".format(gp_id)
        )
    return effect_sizes


def make_effect_size_matrix(sub_vars, var2, sub_size, size2, sub_avg, avg2,
                            effect_sizes):
    """
    Create a matrix where each cell is the effect size comparing that cell of
    a subset to that cell of another group, then return those effect sizes
    :param sub_vars: pandas.DataFrame with the subset's variance for each cell
    :param var2: pandas.DataFrame with another group's variance for each cell 
    :param subset_size: Integer which is how many subjects are in the subset
    :param size2: Integer which is how many subjects are in the other group
    :param sub_avg: pandas.DataFrame with average values of each subset cell
    :param avg2: pandas.DataFrame with the other group's cells' average values
    :param effect_sizes: List of dictionaries, each of which has "Subjects" and
                         "Correlation" (the latter actually means effect size)
    :return: effect_sizes, with the effect sizes just calculated added in
    """
    # Make pooled standard deviation matrix
    pool_stdev = sub_vars.apply(lambda x: dual_df_apply(
        x, var2, get_pooled_stdev, sub_size, size2
    ))

    # Ignore invalid division warnings, then replace invalid quotients with 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        effect_size_matrix =  np.divide(np.subtract(sub_avg, avg2), pool_stdev
                                        ).replace([np.inf, np.nan, -np.inf], 0)

    # Collect all effect sizes mapped to their subject size, for visualization
    effect_size_matrix.applymap(lambda x: effect_sizes.append({
        "Subjects": sub_size, "Correlation": x
    }))
    return effect_sizes


def get_correl_dataframes(all_subsets, cli_args):
    """
    :param all_subsets: List of dictionaries, each of which maps both group
                        numbers to a randomly selected subset of that group
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses the --group-1-avg, --group-2-avg, and 
                     --output arguments. It also passes cli_args to 
                     get_avg_matrices_of_subsets and get_correls_between.
    :return: Dictionary of 3 string:pandas.DataFrame pairs such that each
             DataFrame has a column of subset sizes and one of correlations:
       {sub1_sub2: Correlations between the average matrices of both subsets
        sub1_all2: Correls between group 1 subset avg matrix and group 2 total
        sub2_all1: Correls between group 2 subset avg matrix and group 1 total}
    """
    # Dict of correlation lists to return as DataFrames, and progress tracker
    correl_lists = {subset_id: [] for subset_id in
                    default_vis_titles().keys() if subset_id}
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
        correl_lists = get_sub_pair_correls(subset_size, correl_lists, subsets,
                                            cli_args.spearman_rho)
        progress = update_progress(progress, "making average matrices",
                                   subset_size, start_time)
    return {name: save_correlations_and_get_df(
                cli_args, correls, DEFAULT_CF.format(name)
            ) for name, correls in correl_lists.items()}


def get_avg_matrices_of_subsets(subsets_dict, cli_args):
    """
    Make and return an average matrix for each subset in subsets_dict
    :param subsets_dict: Dict mapping the group number to its subset DataFrame
    :param cli_args: argparse namespace with all command-line arguments.
    :return: Dictionary matching each group's number to its subset's average
             matrix from .pconn files of all its subjects
    """
    avg_matrices = dict()  # Return value: Average matrices of both groups
    subsets_dict.pop("subset_size") 
    
    # Get all data from .pconn files of every subject in the subset
    for sub_num, subset in subsets_dict.items():
        col = get_matr_file_col_name(cli_args) 
        print("Making average matrix for group {} subset.".format(sub_num))
        avg_matrices[sub_num] = get_average_matrix(subset, col, cli_args)
    return avg_matrices


def get_matr_file_col_name(cli_args):
    """
    :param cli_args: argparse namespace with --matrices-conc-1
    :return: String naming demographics .csv column with .nii matrix
    """
    return (DEFAULT_DEM_VAR_MATR if getattr(cli_args, "matrices_conc_1", None)
            else DEFAULT_DEM_VAR_PCONNS) 


def get_sub_pair_correls(subset_size, correl_lists, subsets, rho=None):
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
        if rho:
            params.append(rho)
        correl_lists[df_name].append(get_correls_between(*params))
        title = ("group {0}'s subset and group {1}".format(
                     sub_keys[0][-1], sub_keys[1][-1]
                 ) if "all" in df_name else "both subsets")
        print("Correlations between average matrices of {}:\n{}"
              .format(title, pprint.pformat(correl_lists[df_name])))
    return correl_lists


def save_correlations_and_get_df(cli_args, correls_list, correl_file_name):
    """
    Save correlations to .csv file
    :param cli_args: argparse namespace with all command-line arguments
    :param correls_list: List of dictionaries, each has "Subjects" and
                         "Correlation" as a key with a numerical value
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


def get_traces_to_plot(corr_df, cli_args, ttl, dataset_num=0):
    """
    Get averages, standard devs, confidence intervals, &/or scatter plot traces
    :param corr_df: pandas.DataFrame with one column titled "Subjects" and 
                    another titled "Correlations"; both have numeric values
    :param cli_args: argparse namespace with all command-line arguments
    :param ttl: String naming each dataset to plot on the figure
    :return: List of plotly.graph_objs.Scatter traces to add to plotly figure
    """
    COLORS = plotly.colors.DEFAULT_PLOTLY_COLORS
    nxt_clr = COLORS[dataset_num % len(COLORS)]
    avg_plot, avgs = get_stdev_trace_and_avg(cli_args, corr_df, nxt_clr, ttl)
    to_plot = [avg_plot]
    if "scatter" in cli_args.plot:
        to_plot.append(get_scatter_plot_trace(cli_args, corr_df, nxt_clr, ttl))
    if cli_args.fill:
        to_plot += get_confidence_interval_traces(cli_args, avgs, corr_df)
    return to_plot


def get_stdev_trace_and_avg(cli_args, correls_df, nxt_clr, title):
    """
    :param cli_args: argparse namespace with all command-line arguments
    :param correls_df: pandas.DataFrame with one column titled "Subjects" and
                       another titled "Correlations"; both have numeric values
    :param nxt_clr: String with the RGBA code of an opaque color
    :param title: String naming the dataset to plot on the figure
    :return: plotly.graph_objs.Scatter object made from correls_df data
    """
    avgs = correls_df.groupby(["Subjects"]).agg(lambda x: 
                                                x.unique().sum() / x.nunique())
    stdev = dict()  
    if "stdev" in cli_args.plot:
        avgs["StDev"] = get_shaded_area_bounds(correls_df, "stdev")
        stdev = {"type": "data", "array": avgs["StDev"], "visible": True}
    print(avgs)
    return plotly.graph_objs.Scatter(
        x=avgs.index.values, y=avgs["Correlation"], mode="lines", error_y=stdev,
        name="Average {}".format(title), line_color=nxt_clr
    ), avgs


def get_scatter_plot_trace(cli_args, correls_df, nxt_clr, title):
    """
    :param cli_args: argparse namespace with all command-line arguments
    :param correls_df: pandas.DataFrame with one column titled "Subjects" and
                       another titled "Correlations"; both have numeric values
    :param nxt_clr: String with the RGBA code of an opaque color
    :param title: String naming the dataset to plot on the figure
    :return: plotly.graph_objs.Scatter object made from correls_df data
    """
    if cli_args.rounded_scatter: # Round to reduce # of points
        digits = correls_df["Correlation"].apply(lambda x: len(str(x))-2)
        correls_df = correls_df.round(decimals=int(digits.max()**(1/4))
                                      ).drop_duplicates()
    return plotly.graph_objs.Scatter(
        x=correls_df["Subjects"], y=correls_df["Correlation"],
        name="All {}".format(title), line_color=nxt_clr,
        marker={"size": cli_args.marker_size}, mode="markers"
    )


def get_confidence_interval_traces(cli_args, avgs, correls_df):
    """
    Add upper & lower bounds (all data or CI) of shaded area to plot as lines
    :param cli_args: argparse namespace with all command-line arguments
    :param avgs: pandas.Series aggregating averages across correls_df sub_sizes 
    :param correls_df: pandas.DataFrame with one column titled "Subjects" and
                       another titled "Correlations"; both have numeric values
    :return: List of 2 plotly.graph_objs.Scatter objects: (lower, upper)
    """
    bounds = get_shaded_area_bounds(correls_df, cli_args.fill)
    bounds_params = ({"showlegend": False} if cli_args.fill == "all" else
                     {"name": "95 percent confidence interval",
                      "showlegend": True})
    return [plotly.graph_objs.Scatter(
        x=avgs.index.values, y=bounds[0], line_color=rgba("red", 0), # clear
        fillcolor=rgba("red", 0.2), fill="tonexty", **bounds_params
    ), plotly.graph_objs.Scatter(x=avgs.index.values, y=bounds[1],
                                 line_color=rgba("red", 0), showlegend=False)]


def make_visualization(traces, corr_df, cli_args, corr_df_name):
    """
    Show image & export it as a .png file 
    :param traces: List of plotly.graph_objs.Scatter objects to plot
    :param corr_df: pandas.DataFrame with one column titled "Subjects" and
                       another titled "Correlations"; both have numeric values
    :param cli_args: argparse namespace with all command-line arguments
    :param corr_df_name: String identifying the visualization to create
    """
    vis_ttl, vis_file = get_vis_title_and_fname(cli_args, corr_df_name)
    if cli_args.plot_with_matlab:  # Make visualization w/ compiled MATLAB code
        corrfile_dir = (os.path.dirname(cli_args.only_make_graphs[0]) if
                        cli_args.only_make_graphs else cli_args.output)
        corrfile = os.path.join(corrfile_dir, DEFAULT_CF.format(corr_df_name))
        nowID = ''.join(c for c in str(now()) if c.isalnum())
        paramcsv = os.path.join(cli_args.output, DEFAULT_MAT_CSV.format(nowID))
        with open(paramcsv, "w+") as outf:  # Write param file for MATLAB code
            outf.write(",".join([corrfile,
                                *[str(x) for x in cli_args.matlab_rgba]])+"\n")
        os.chmod(paramcsv, 775)
        mat_params = []
        pm_IDs = {"LowerBound": getattr(cli_args, "matlab_lower_bound", None),
                  "NoEdge": getattr(cli_args, "matlab_no_edge", None),
                  "ShowThreshold": getattr(cli_args, "matlab_show", None),
                  "UpperBound": getattr(cli_args, "matlab_upper_bound", None)} 
        for pm_name, pm_val in pm_IDs.items():
            if pm_val:
                mat_params += [pm_name, str(pm_val)]
        subprocess.check_call((
            MATLAB_PLOTTER, cli_args.plot_with_matlab, paramcsv, *mat_params,
            "Outputfile",
            os.path.join(cli_args.output,  vis_file.split(".")[0] + ".tif"),
        )) 
    else:  # Make visualization w/ python plotly package
        layout = get_plot_layout(get_layout_args(cli_args, corr_df, vis_ttl))
        with HiddenPrints():  # suppress plotly.offline.plot()'s huge textwall 
            plotly.offline.init_notebook_mode()
            plotly.offline.plot({"data": traces, "layout": layout},
                                image="png", filename=vis_file)


def get_vis_title_and_fname(cli_args, corr_df_name):
    """
    :param cli_args: argparse namespace with --graph-title argument
    :param corr_df_name: String identifying the visualization to create
    :return: Tuple of 2 strings, the visualization's title and filename
    """
    vis_title = (cli_args.graph_title if cli_args.graph_title
                 else default_vis_titles()[corr_df_name])
    vis_file = graph_title_to_vis_fname(vis_title) + ".html"
    i = 1
    while os.access(vis_file, os.R_OK):
        i += 1
        vis_file = "{}_{}.html".format(graph_title_to_vis_fname(vis_title), i)
    return vis_title, vis_file


def graph_title_to_vis_fname(title):
    """
    :param title: String naming the title of the visualization
    :return: title, but converted into a valid filename
    """
    return "".join(title.replace("<br>", "_").replace("/", "_").split())


def get_layout_args(cli_args, correls_df, title):
    """
    :param cli_args: argparse namespace with all command-line arguments. This
                     function uses axis and title font sizes, plus y-range.
    :param correls_df: pandas.DataFrame with one column titled "Subjects" and
                       another titled "Correlations"; both have numeric values
    :param title: String to put at the top of the visualization
    :return: Tuple of all input arguments for get_plot_layout function
    """
    # Use default y-axis range or one chosen by user
    y_axis_min, y_axis_max = cli_args.y_range if cli_args.y_range else (
        correls_df["Correlation"].min(), correls_df["Correlation"].max()
    )
    return (y_axis_min, y_axis_max, title, cli_args.place_legend,
            cli_args.title_font_size, cli_args.axis_font_size,
            correls_df["Subjects"].mean(), not cli_args.hide_legend,
            ("Effect Size (d)" if cli_args.calculate == "effect-size"
             else "Correlation (r)"))


def get_plot_layout(args):
    """
    Return all format settings for creating a pretty plot visualization.
    :param args: Tuple of all arguments needed for visualization format
    :return: Nested dictionary containing all needed plot layout attributes
    """
    # Local variables from unpacked list of args to use in layout: 
    # Lowest and highest y-values to show, graph title, y-position of legend,
    # graph title and axis title font sizes, average subset size, y-axis title
    y_min, y_max, title, lgnd_y, ttl_size, axis_font, x_avg, show, y_ttl = args
    y_range_step = (y_max - y_min) / 10 # Space buffer above y_max and below y_min

    def get_axis_layout(title_txt, **kwargs):
        """
        Get all of the parameters that "xaxis" and "yaxis" have in common for 
        get_plot_layout in order to avoid redundant code.
        :param title: String of text to be displayed at the top of the graph
        :param kwargs: Dictionary containing all of the parameters that "xaxis"
                       and "yaxis" do not have in common
        :return: Dictionary with all parameters needed for formatting an axis
        """
        result = {"title": {"font": {"size": ttl_size}, "text": title_txt},
                  "tickcolor": rgba("black"), "ticklen": 15, "ticks": "outside",
                  "tickfont": {"size": axis_font}, "tickwidth": 2,
                  "showline": True, "linecolor": rgba("black"), "linewidth": 2}
        result.update(kwargs)
        return result

    return {"title": {"text": title, "x": 0.5, "xanchor": "center",
                      "font": {"size": ttl_size}},
            "paper_bgcolor": rgba("white"), "plot_bgcolor": rgba("white"),
            "legend": {"font": {"size": axis_font}, "x": 0.5, "y": lgnd_y},
            "xaxis": get_axis_layout(
                title_txt="Sample Size (n)", tick0=0, tickmode="linear",
                dtick=count_digits_of(x_avg)
            ), "showlegend": show,
            "yaxis": get_axis_layout(
                title_txt=y_ttl, tickmode="auto", nticks=6,
                range=(y_min - y_range_step, y_max + y_range_step) 
            )}


if __name__ == '__main__':
    main()
