#! /usr/bin/env python3

"""
Automated subset selection and analysis for ABCD resource paper
Greg Conan: conan@ohsu.edu
Created 2019-09-17
Last Updated 2019-10-18
"""

##################################
#
# Script to randomly select pairs of subsets of data of given sizes, find the
# correlation between the average matrices of each subset, and plot all of
# those correlations for each subset size
#
##################################

import argparse
import nibabel
import numpy as np
import os
import pandas as pd
import plotly
from pprint import pprint
import psutil  # TODO: Test full script after removing this to see if it's unneeded
from scipy.spatial.distance import euclidean
from scipy.stats import sem, t
import subprocess

# Constants
DEFAULT_AVG_GROUP_1 = ("/mnt/rose/shared/projects/ABCD/avg_pconn_maker/"
                       "group1_10min_mean.pconn.nii")
DEFAULT_AVG_GROUP_2 = ("/mnt/rose/shared/projects/ABCD/avg_pconn_maker/"
                       "group2_10min_mean.pconn.nii")
DEFAULT_DEMO_TOTAL_FILE = "total_group_matching_stats.csv"
DEFAULT_FOLDER = "subset-analysis-output"
DEFAULT_N_ANALYSES = 1
DEFAULT_SUBSET_SIZE = [50, 100, 200, 300, 400, 500]
VISUALIZATION_TITLES = {
    "sub1_sub2": "Correlations Between Average Subsets",
    "sub1_all2": "Group 1 Subset to Group 2 Correlation",
    "sub2_all1": "Group 1 to Group 2 Subset Correlation"
}
PWD = os.getcwd()


def main():

    # Get and validate all command-line arguments from user
    cli_args = get_cli_args()

    started_at = subprocess.check_output("date").decode("utf-8")
    print("Started at " + started_at)

    # If user said to skip the subset generation and use pre-existing subset
    # correlations, then make pd.DataFrame to visualize those
    if cli_args.only_make_graphs:
        for correls_csv in cli_args.only_make_graphs:
            correls_df = pd.read_csv(correls_csv)
            csv_sets_name = get_vis_name_from(correls_csv, list(
                VISUALIZATION_TITLES.keys()))
            csv_sets_name = (
                VISUALIZATION_TITLES[csv_sets_name] if csv_sets_name in
                VISUALIZATION_TITLES else "Correlations Between Unknown Groups"
            )
            make_visualizations(correls_df, cli_args.fill, cli_args.y_range,
                                csv_sets_name)
    else:

        # If user said to skip subset generation but get correlations for pre-
        # existing subsets, then do that
        if getattr(cli_args, "skip_subset_generation", None) is not None:
            all_subsets = []
            for subset_csv in os.listdir(cli_args.skip_subset_generation):
                subset_df = pd.read_csv(os.path.join(
                    cli_args.skip_subset_generation, subset_csv))
                if "1" in subset_df and "2" in subset_df:
                    all_subsets.append({
                        1: cli_args.group_1_demo[cli_args.group_1_demo[
                            "subjectkey"].isin(subset_df.pop("1").tolist())],
                        2: cli_args.group_2_demo[cli_args.group_2_demo[
                            "subjectkey"].isin(subset_df.pop("2").tolist())],
                        "subset_size": len(subset_df.index)
                    })

        # Otherwise (normal use case), generate and save all subsets and the
        # correlations between them
        else:
            all_subsets = save_and_get_all_subsets(cli_args)
        for correls_df_name, correls_df in get_correl_dataframes(
                all_subsets, cli_args).items():
            make_visualizations(correls_df, cli_args.fill, cli_args.y_range,
                                VISUALIZATION_TITLES[correls_df_name])

    print("Started at " + started_at)
    print("Finished at " + subprocess.check_output("date").decode("utf-8"))


def get_vis_name_from(file_name, vis_names):
    """
    :param file_name:
    :param vis_names:
    :return: Element of vis_names which is in the file_name string
    """
    found = None
    while not found and len(vis_names) > 0:
        next_vis_name = vis_names.pop()
        if next_vis_name in file_name:
            found = next_vis_name
    return found


def get_cli_args():
    """
    Get and validate all args from command line using argparse.
    :return: Namespace containing all validated inputted command line arguments
    """
    # Create arg parser
    parser = argparse.ArgumentParser(
        description=("Script to randomly select pairs of subsets of data of "
                     "given sizes, find the correlation between the average "
                     "matrices of each subset, and plot all of those "
                     "correlations for each subset size.")
    )

    # Required: Full list of group 1 demographic data
    parser.add_argument(
        "group_1_demo_file",
        type=validate_readable_file,
        help=("Path to a .csv file containing all demographic information "
              "about the subjects in group 1.")
    )

    # Required: Full list of group 2 demographic data
    parser.add_argument(
        "group_2_demo_file",
        type=validate_readable_file,
        help=("Path to a .csv file containing all demographic information "
              "about the subjects in group 2.")
    )

    # Optional: Choose what to fill in the visualization
    parser.add_argument(
        "-f",
        "--fill",
        choices=["confidence_interval", "all"],
        default="confidence_interval",
        help=("Choose which data to shade in the visualization. Choose 'all' "
              "to shade in the area within the minimum and maximum "
              "correlations in the dataset. Choose 'confidence_interval' to "
              "only shade in the 95 percent confidence interval of the data. "
              "By default, this argument will be 'confidence_interval'.")
    )

    # Optional: Import data from .csv files already created by this script
    # instead of creating new ones
    parser.add_argument(
        "-g",
        "--only_make_graphs",
        nargs="+",
        type=validate_readable_file,
        help=("Include this flag to import data from a .csv file already made "
              "by this script instead of making a new one, just to make a "
              "graph visualization of already-existing data. If this flag "
              "is included, it must be a path to a readable .csv file with 2 "
              "columns: 'Subjects' (the number of subjects in each subset) "
              "and 'Correlation' (the correlation between each randomly "
              "generated subset in that pair).")
    )

    # Optional: Path to average matrix .pconn file for group 1
    parser.add_argument(
        "--group_1_avg_file",
        type=validate_readable_file,
        default=DEFAULT_AVG_GROUP_1,
        help=("Path to a .pconn file containing the average matrix for group "
              "1. By default, this path will be " + DEFAULT_AVG_GROUP_1)
    )

    # Optional: Path to average matrix .pconn file for group 2
    parser.add_argument(
        "--group_2_avg_file",
        type=validate_readable_file,
        default=DEFAULT_AVG_GROUP_2,
        help=("Path to a .pconn file containing the average matrix for group "
              "2. By default, this path will be " + DEFAULT_AVG_GROUP_2)
    )

    # Optional: Number of subsets
    parser.add_argument(
        "-n",
        "--n_analyses",
        default=DEFAULT_N_ANALYSES,
        type=validate_whole_number,
        help="Number of times to generate a pair of subsets and analyze them."
    )

    # Optional: Output folder
    parser.add_argument(
        "-o",
        "--output",
        default=os.path.join(PWD, DEFAULT_FOLDER),
        help=("Path to folder to save subset data in. The default folder is "
              + DEFAULT_FOLDER)
    )

    # Optional: Path to parent directory of folders with .pconn files
    parser.add_argument(
        "-p"
        "--parent_path",
        type=validate_readable_file,
        help=("Path to parent directory containing .pconn files. The "
              "group_1_demo_file and the group_2_demo_file should both have "
              "a column named 'pconn_10min' where each row has a path to the "
              ".pconn file for the subject in that row. That path should be "
              "a relative path from the parent_path directory, the value "
              "given as this argument. So, parent_path must be a valid path "
              "to a readable directory. By default, if this flag is excluded, "
              "then parent_path will be the parent directory of "
              "group_1_demo_file and group_2_demo_file.")
    )

    # Optional: Number of subjects in each subset pair
    parser.add_argument(
        "-s",
        "--subset_size",
        nargs="+",
        default=DEFAULT_SUBSET_SIZE,
        type=validate_whole_number,
        help=("Number of subjects to include in subsets. Include a list of "
              "whole numbers to generate subsets pairs of different sizes. By "
              "default, the subset sizes will be " + str(DEFAULT_SUBSET_SIZE))
    )

    # Optional: Make averages/correlations from existing subsets
    parser.add_argument(
        "-u",
        "--skip_subset_generation",
        nargs="?",
        const="output",
        help=("Include this flag to calculate correlations and create the "
              "visualization using existing subsets instead of making new "
              "ones. By default, the subsets to use will be assumed to exist "
              "in the --output folder. If you want to load subsets from "
              "a different folder, add the path to this flag as a parameter.")
    )

    # Optional: Custom data range for visualization
    parser.add_argument(
        "-y",
        "--y_range",
        nargs=2,
        type=float,
        help=("Range of y_axis in visualization. By default, this script will "
              "automatically set the y-axis to show all of the correlation "
              "values and nothing else. If this argument is used, then it "
              "should be two floating-point numbers: the minimum and maximum "
              "values to be shown on the visualizations' y-axes.")
    )

    # Optional: Do inverse Fisher-Z transformation (hyperbolic tangent) on data
    parser.add_argument(
        "-z",
        "--inverse_fisher_z",
        action="store_true",
        help=("Do inverse Fisher-Z transform on the matrices imported from "
              "the .pconn files of the data before getting correlations.")
    )

    return validate_cli_args(parser.parse_args(), parser)


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    try:
        if not cli_args.only_make_graphs:

            if (hasattr(cli_args, "skip_subset_generation") and
                    cli_args.skip_subset_generation is not None):
                if cli_args.skip_subset_generation == "output":
                    cli_args.skip_subset_generation = cli_args.output
                else:
                    validate_readable_file(cli_args.skip_subset_generation)

            # Make directory to save matrices in
            os.makedirs(os.path.abspath(cli_args.output), exist_ok=True)

            # Get path to .pconn parent directory
            if not hasattr(cli_args, "parent_path"):
                cli_args.parent_path = os.path.commonpath((
                    cli_args.group_1_demo_file, cli_args.group_2_demo_file))

            # For each group, get the path to the directory with its .pconn
            # files, and the path to the file containing its demographic data
            for group_num in ("group_1_", "group_2_"):

                # Import all of this group's demographics; drop empty columns
                setattr(cli_args, group_num + "demo", pd.read_csv(
                    getattr(cli_args, group_num + "demo_file"), na_values=" "
                ).dropna(how="all", axis=1))

                # Get average matrix for each group
                setattr(cli_args, group_num + "avg", load_matrix_from_pconn(
                    getattr(cli_args, group_num + "avg_file"), cli_args))

        return cli_args

    except OSError as e:
        parser.error(str(e))


def validate_readable_file(path):
    """
    Throw exception unless parameter is a valid readable filename string. This
    is used instead of argparse.FileType("r") because the latter leaves an open
    file handle, which has caused problems.
    :param path: Parameter to check if it represents a valid filename
    :return: A valid filename as a string
    """
    try:
        assert os.access(path, os.R_OK)
        return path if os.path.isabs(path) else os.path.abspath(path)
    except (AssertionError, OSError, TypeError):
        argparse.ArgumentTypeError("Cannot read file at " + path)


def validate_whole_number(to_validate):
    """
    Throw argparse exception unless to_validate is a positive integer
    :param to_validate: Object to test whether it is a positive integer
    :return: to_validate if it is a positive integer
    """
    try:
        to_validate = int(to_validate)
        assert to_validate > 0
        return to_validate
    except (AssertionError, TypeError):
        raise argparse.ArgumentTypeError(
            str(to_validate) + " is not a positive integer.")


def save_and_get_all_subsets(cli_args):
    """
    Randomly generate all pairs of subsets, save each of them to a .csv file,
    and then return a pandas.DataFrame with the correlation between each pair
    :param cli_args: argparse namespace with all command-line arguments
    :return: List of dictionaries, each of which maps a subset's group number
    to the subset for one pair of subsets
    """
    # Return value: List of subsets
    all_subsets = []

    for subset_size in cli_args.subset_size:
        print("Generating {} randomly selected group(s) of {} subjects."
              .format(cli_args.n_analyses, subset_size))

        # Get average correlation from user-defined number of pairs of average
        # matrices of randomly generated subsets
        for i in range(cli_args.n_analyses):

            # Select subsets and make average matrices from .pconns for each
            subsets = randomly_select_subsets(cli_args, subset_size)

            # Save randomly generated subsets to .csv files
            save_subsets(subsets, cli_args.output, i)

            all_subsets.append(subsets)
    return all_subsets


def get_correls_between(arr1, arr2, num_subjects):
    """
    :param arr1: np.ndarray with only numeric values
    :param arr2: np.ndarray with only numeric values
    :param num_subjects: Number of subjects to randomly select from each group
    :return: Dictionary mapping "Subjects" to the number of rows in set1
    and mapping "Correlation" to the correlation between both sets flattened
    """
    return {"Subjects": num_subjects,
            "Correlation": np.corrcoef(arr1.flatten(), arr2.flatten())[0, 1]}


def save_subsets(subs, output_dir, sub_num):
    """
    Given a pair of subsets, save them to a .csv file
    :param subs: Dictionary mapping each subset's group number to that subset
    :param output_dir: Path to directory to save subset .csv files into
    :param sub_num: Arbitrary number so that multiple subsets of the same group
    can be saved without conflicting filenames
    :return: N/A
    """
    to_save = pd.DataFrame({
        1: subs[1]["subjectkey"].reset_index(drop=True),
        2: subs[2]["subjectkey"].reset_index(drop=True)
    })
    to_save.to_csv(os.path.join(
        output_dir, "_".join(("subset", str(sub_num + 1), "with", str(
            len(to_save.index)), "subjects.csv"))), index=False)


def read_file_into_list(filepath):
    """
    Given the path to a file, read all of the file's lines as a list
    :param filepath: Path to file to read
    :return: List of all lines in the file at filepath
    """
    with open(filepath, "r") as infile:
        return [line.strip() for line in infile]


def randomly_select_subsets(cli_args, subset_size):
    """
    :param cli_args: argparse namespace with all command-line arguments
    :param subset_size: Number of subjects to randomly select from each group
    :return: A dictionary mapping group numbers to randomly selected subsets
    of size subset_size from those groups, where each subset is demographically
    representative of the overall group as shown by chi square. The dictionary
    also includes the subset size.
    """
    # Iterate over both groups to get a subset of each
    subsets = {1: None, 2: None, "subset_size": subset_size}
    euclidean_threshold = get_est_euclid_chi_sq_threshold(subset_size)
    print("Estimated Euclidean distance threshold for chi-squared "
          "significance: " + str(euclidean_threshold))
    for group_num in subsets.keys():
        group_all = getattr(cli_args, str(group_num).join(("group_", "_demo")))
        other_group = getattr(cli_args, "group_{}_demo".format(
                              1 if group_num == 2 else 2))

        # Get information about group total dataframe for comparison
        columns = group_all.select_dtypes(include=["number"]).columns.values

        # Compare subset of one group to the other's total
        group_avgs = [other_group[col].mean(skipna=True) for col in columns]

        # Randomly pick subsets until finding 1 that demographically represents
        # the overall group
        loops = 0
        subset_matches_total = False
        while not subset_matches_total:
            loops += 1

            # Generate subset of group
            subsets[group_num] = group_all.sample(n=subset_size)

            # Get subset demographics to compare to total demo of other group
            sub_avgs = [subsets[group_num][col].mean(skipna=True)
                        for col in columns]

            # If any column in subset averages has mean of N/A, throw it out
            if not any(np.isnan(el) for el in sub_avgs):
                eu_dist = euclidean(group_avgs, sub_avgs)
                subset_matches_total = (eu_dist < euclidean_threshold)
            if subset_matches_total or loops % 100 == 0:
                print("{} subsets of group {} randomly generated.".format(
                      loops, group_num))
            if subset_matches_total:
                print("Euclidean distance: " + str(eu_dist))

    return subsets


def get_est_euclid_chi_sq_threshold(subjects):
    """
    Estimated Euclidean distance threshold for chi-squared significance
    :param subjects: Number of subjects in a subset (int)
    :return: The estimated threshold of a Euclidean distance between a subset
    and a total set such that the subset is not statistically different from
    the total set according to a chi-squared test.
    """
    return 19.0225676337-(4.6882289084*np.log(subjects/100))


def get_average_matrices_of_subsets(subsets, cli_args):
    """
    :param subsets: Dictionary mapping the group number to its pandas.DataFrame
    with a subjectkey column listing subject IDs
    :param cli_args: argparse namespace with all command-line arguments
    :return: Dictionary matching each group's number to its subset's average
    matrix from .pconn files of all its subjects
    """
    # Return value: Average matrices of group 1 and group 2
    average_matrices = {}

    sub_size = subsets.pop("subset_size")

    # Local function to append a subject's 2D matrix to the subset's 3D matrix
    def add_pconn_to_3d_matrix(subject_pconn, full_matrix):
        """
        :param subject_pconn: Path to a .pconn with a subject's 2D data matrix
        :param full_matrix: 3D matrix to append the imported data matrix to
        :return: full_matrix, with subject's matrix added on the 3rd dimension
        """
        matrix_to_add = load_matrix_from_pconn(subject_pconn, cli_args)
        if cli_args.inverse_fisher_z:
            matrix_to_add = np.tanh(matrix_to_add)
        return np.dstack((full_matrix, matrix_to_add))

    # Get all data from .pconn files of every subject in the subset
    print("Importing .pconn files from " + cli_args.parent_path)
    for subset_num, subset in subsets.items():

        # Get one .pconn file to initialize 3D matrix
        subject_pconns = subset["pconn_10min"].iteritems()
        first_pconn = load_matrix_from_pconn(next(subject_pconns)[1], cli_args)

        # Initialize the 3D matrix from that .pconn file
        all_matrices = first_pconn.reshape((
            first_pconn.shape[0], first_pconn.shape[1], 1
        ))

        # Build the rest of the 3D matrix
        for subject in subject_pconns:
            all_matrices = add_pconn_to_3d_matrix(subject[1], all_matrices)
            print("Matrix size: {}. Progress: {:.1%} done.".format(
                all_matrices.shape, all_matrices.shape[2] / sub_size))

        # Compute average matrix
        average_matrix = all_matrices.mean(axis=2)
        average_matrices[subset_num] = average_matrix
        print("Average matrix:")
        print(average_matrix)
    return average_matrices


def load_matrix_from_pconn(pconn, cli_args):
    """
    Import a subject's 2D data matrix from a .pconn file
    :param pconn: Path to .pconn file of subject
    :param cli_args: argparse namespace with all command-line arguments
    :return: 2D matrix of data from the subject's .pconn
    """
    return np.array(nibabel.cifti2.load(os.path.join(
                    cli_args.parent_path, pconn)).get_data().tolist())


def get_correl_dataframes(all_subsets, cli_args):
    """
    :param all_subsets: List of dictionaries, each of which maps both group
    numbers to a randomly selected subset of that group
    :param cli_args: argparse namespace with all command-line arguments
    :return: Dictionary of 3 string:pandas.DataFrame pairs such that each
    DataFrame has a column of subset sizes and a column of correlations:
       {sub1_sub2: Correlations between the average matrices of both subsets
        sub1_all2: Correls between group 1 subset avg matrix and group 2 total
        sub2_all1: Correls between group 2 subset avg matrix and group 1 total}
    """
    # Return value: Dict of correlation lists to become pandas.DataFrames
    result = {"sub1_sub2": [], "sub1_all2": [], "sub2_all1": []}

    # Get each pair of average matrices, their correlation with each other, and
    # each one's correlation with the other group's average matrix
    for sub_pair in all_subsets:
        sub1_avg, sub2_avg = get_average_matrices_of_subsets(
            sub_pair.copy(), cli_args).values()
        subset_size = sub_pair.pop("subset_size")
        print("subset_size: " + str(subset_size))

        result["sub1_sub2"].append(get_correls_between(sub1_avg, sub2_avg,
                                                       subset_size))
        result["sub1_all2"].append(get_correls_between(cli_args.group_2_avg,
                                                       sub1_avg, subset_size))
        result["sub2_all1"].append(get_correls_between(cli_args.group_1_avg,
                                                       sub2_avg, subset_size))
        pprint("\nCorrelations between average matrices of ".join((
            "", "both subsets:\n" + str(result["sub1_sub2"]),
            "subset 1 and group 2:\n" + str(result["sub1_all2"]),
            "subset 2 and group 1:\n" + str(result["sub2_all1"])
        )))

    return {name: save_correlations_and_get_df(
            correls, cli_args, "correlations_{}.csv".format(name))
            for name, correls in result.items()}


def save_correlations_and_get_df(correls_list, cli_args, correl_file):
    """
    Save correlations to .csv file
    :param correls_list: List of dictionaries where each dictionary contains
    "Subjects" and "Correlation" as a key with a numerical value
    :param cli_args: argparse namespace with all command-line arguments
    :param correl_file: Name of the file to save correlations into
    :return: pandas.DataFrame with a header row ("Subjects", "Correlation"),
    subset sizes in one column, and correlations in its second column
    """
    n_vs_correls_df = pd.DataFrame(correls_list, columns=[
        "Subjects", "Correlation"])
    n_vs_correls_df.to_csv(os.path.join(cli_args.output, correl_file),
                           index=False)
    return n_vs_correls_df


def make_visualizations(correls_df, fill_area, y_range, vis_title):
    """
    Create a graph visualization from correlation data
    :param correls_df: pandas.DataFrame with one column titled "Subjects" and
    another titled "Correlations" where both have numeric values
    :param fill_area: A string that is either "all" or "confidence_interval"
    :param y_range: None if the y-axis values range will be automatically
    calculated based on the data, or a tuple of 2 floats if the user entered a
    custom y-axis range.
    :param vis_title: Title of visualization to create
    :return: N/A
    """
    # Make scatter plot mapping subset size to pairwise correlations
    scatter_plot = plotly.graph_objs.Scatter(
        x=correls_df["Subjects"], y=correls_df["Correlation"], mode="markers",
        name="All correlations", line_color="rgba(255,0,0,1)",
        marker={"size": 8})

    # Add average lines to plot using averages of each subset size
    averages = correls_df.groupby(["Subjects"]).agg(
        lambda x: x.unique().sum()/x.nunique())
    print("Average correlations for {}:\n{}".format(vis_title, averages))
    avgs_plot = plotly.graph_objs.Scatter(
        x=averages.index.values, y=averages["Correlation"], mode="lines",
        name="Average correlations")

    # Add upper & lower bounds (all data or confidence intervals) of shaded
    # area to plot as lines
    bounds = get_shaded_area_bounds(correls_df, fill_area)
    bounds_params = ({"showlegend": False} if fill_area == "all" else
                     {"name": "95% confidence interval", "showlegend": True})
    lower_plot = plotly.graph_objs.Scatter(
        x=averages.index.values, y=bounds[0], fill="tonexty",
        fillcolor="rgba(255,0,0,0.2)", line_color="rgba(0,0,0,0)",
        **bounds_params
    )
    upper_plot = plotly.graph_objs.Scatter(
        x=averages.index.values, y=bounds[1], line_color="rgba(0,0,0,0)",
        showlegend=False
    )

    # Show plots, either using custom y-axis range or calculating a range
    if y_range:
        y_axis_min, y_axis_max = y_range
    else:
        y_axis_min = correls_df["Correlation"].min()
        y_axis_max = correls_df["Correlation"].max()
    plotly.io.show({
        "data": [scatter_plot, avgs_plot, upper_plot, lower_plot],
        "layout": get_plot_layout(y_axis_min, y_axis_max, vis_title)
    })


def get_shaded_area_bounds(all_data_df, to_fill):
    """
    :param all_data_df: pandas.DataFrame with a "Subjects" column and a
    "Correlation" column, where both only have numeric values
    :param to_fill: A string that is either "all" to shade the area between the
    minimum and maximum correlations for each number of subjects in all_data_df
    or "confidence_interval" to shade the area within a 95% confidence interval
    of the correlations for each number of subjects in all_data_df
    :return: A tuple of 2 pandas.Series objects where the first is the lower
    boundary of the shaded area and the second is the upper boundary
    """
    if to_fill == "confidence_interval":
        intervals = all_data_df.groupby(["Subjects"]).agg(
            lambda x: get_confidence_interval(x))
        intervals = pd.DataFrame(intervals["Correlation"].tolist(),
                                 index=intervals.index)
        result = (intervals[0], intervals[1])
    elif to_fill == "all":
        result = (all_data_df.groupby(["Subjects"]).agg(
                      lambda x: x.quantile(0))["Correlation"],
                  all_data_df.groupby(["Subjects"]).agg(
                      lambda x: x.quantile(1))["Correlation"])
    else:
        raise ValueError("Invalid value for --fill parameter.")
    return result


def get_confidence_interval(series, confidence=0.95):
    """
    :param series: pandas.Series filled with numeric data
    :param confidence: Percentage for confidence interval; default 95%
    :return: Tuple representing series's confidence interval: (start, end)
    """
    mean = series.mean()
    h = sem(series) * t.ppf((1 + confidence) / 2, len(series) - 1)
    return mean - h, mean + h


def get_plot_layout(y_min, y_max, graph_title):
    """
    Return all format settings for creating a pretty plot visualization. This
    function needs its parameters to determine the range of the y-axis.
    :param y_min: Lowest y-value to be displayed on the graph
    :param y_max: Highest y-value to be displayed on the graph
    :param graph_title: String of text to be displayed at the top of the graph
    :return: Nested dictionary containing all needed plot attributes
    """
    title_font_size = 35
    font_size = 25
    black = "rgb(0, 0, 0)"
    white = "rgb(255, 255, 255)"
    y_range_step = (y_max - y_min)/10  # Space buffer above y_max & below y_min

    def get_axis_layout(title, **kwargs):
        result = {
            "title": {"font": {"size": title_font_size}, "text": title},
            "tickcolor": black, "ticklen": 15, "tickfont": {"size": font_size},
            "ticks": "outside", "tickwidth": 2, "showline": True,
            "linecolor": black, "linewidth": 2
        }
        result.update(kwargs)
        return result

    return {
        "title": {"text": graph_title, "x": 0.5, "xanchor": "center",
                  "font": {"size": title_font_size}},
        "paper_bgcolor": white,
        "plot_bgcolor": white,
        "legend": {"font": {"size": font_size}, "y": 0.4, "x": 0.5},
        "xaxis": get_axis_layout(
            title="Sample Size (n)", tick0=0, dtick=100, tickmode="linear"
        ),
        "yaxis": get_axis_layout(
            title="Correlation (r)", nticks=5, tickmode="auto",
            range=(y_min - y_range_step, y_max + y_range_step),
        )
    }


if __name__ == '__main__':
    main()