#! /usr/bin/env python3

"""
Euclidean distance threshold finder for subset analysis
Greg Conan: conan@ohsu.edu
Created 2019-09-17
Last Updated 2019-10-14
"""

##################################
#
# Script to randomly select pairs of subsets of data of a given size and get
# the Euclidean distance threshold below which the subsets pass chi-square test
# showing that each is not significantly different from the other's total
# population, then get an equation to estimate that threshold given group size
#
##################################

import argparse
import numpy as np
import os
import pandas as pd
from scipy import stats
from scipy.spatial.distance import euclidean as euclidean_distances
import subprocess

# Constants
DEFAULT_DEMO_TOTAL_FILE = "total_group_matching_stats.csv"
DEFAULT_MINUTES_LIMIT = 10
DEFAULT_N_ANALYSES = 1
DEFAULT_SUBSET_SIZE = [50, 100, 200, 300, 400, 500]
PWD = os.getcwd()
RAW_DATA_FOLDER = os.path.join(PWD, "abcd-raw-data-for-subset-analysis")


def main():

    # Get and validate all command-line arguments from user
    cli_args = get_cli_args()

    started_at = subprocess.check_output("date").decode("utf-8")
    print("Started at " + started_at)

    frequency_info = []
    freq_maximums = []
    for subset_size in cli_args.subset_size:
        print(" ".join((
              "Running", str(cli_args.n_analyses), "analyses of randomly "
              "selected groups of", str(subset_size), "subjects.")))

        for i in range(cli_args.n_analyses):

            # Select subsets and get the Euclidean distances for each subset
            freq_1, freq_2 = randomly_select_subsets(cli_args, subset_size)
            frequency_info.append(freq_1)
            frequency_info.append(freq_2)
            print("All frequency info: " + str(frequency_info))

            # Find the maximum Euclidean distance for each number of subjects
            for freq in (freq_1, freq_2):
                if (not freq_maximums or freq["Subjects"] != freq_maximums[-1][
                        "Subjects"]):
                    print("Adding a maximum for subjects="
                          + str(freq["Subjects"]))
                    freq_maximums.append(freq)
                elif (freq["Euclidean Distance"] >
                      freq_maximums[-1]["Euclidean Distance"]):
                    print(str(freq) + " > " + str(freq_maximums[-1]))
                    freq_maximums[-1] = freq
                print("All maximums: " + str(freq_maximums))

    # Collect and print all
    all_freq_info = pd.DataFrame(frequency_info, columns=[
        "Group", "Subjects", "Euclidean Distance"])
    print("All frequency info:\n" + str(all_freq_info))

    # Get logarithmic best-fit curve of maximum Euclidean distances for each
    # number of subjects
    maximums_df = pd.DataFrame(freq_maximums)
    print("All maximums:\n" + str(maximums_df))
    print("Logarithmic best-fit equation: " + get_logarithmic_fit(
        maximums_df["Subjects"], maximums_df["Euclidean Distance"]
    ))

    print("Started at " + started_at)
    print("Finished at " + subprocess.check_output("date").decode("utf-8"))


def get_cli_args():
    """
    Get and validate all args from command line using argparse.
    :return: Namespace containing all validated inputted command line arguments
    """
    # Create arg parser
    parser = argparse.ArgumentParser(
        description=("Script to randomly select pairs of subsets of data of a "
                     "given size and get the Euclidean distance threshold "
                     "below which the subsets pass chi-square test.")
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

    # Optional: Threshold to distinguish categorical and continuous demo vars
    parser.add_argument(
        "-e",
        "--exclude",
        type=validate_whole_number,
        default=None,
        help=("Threshold to distinguish categorical and continuous "
              "demographic variables. If this flag is included, then any "
              "demographic variable with at least this many unique values "
              "will be excluded from chi-squared comparison on the basis of "
              "being a continuous instead of a categorical variable.")
    )

    # Optional: Number of times to generate new subsets until being done
    parser.add_argument(
        "-l",
        "--loops",
        type=validate_whole_number,
        default=None,
        help="How many times to generate another subset before giving up and "
             "printing the frequency distribution"
    )

    # Optional: Max number of rows to display
    parser.add_argument(
        "-m",
        "--max_rows",
        type=validate_whole_number,
        default=None,
        help="Maximum number of rows to display when printing a Pandas object."
    )

    # Optional: Number of subsets
    parser.add_argument(
        "-n",
        "--n_analyses",
        type=validate_whole_number,
        default=DEFAULT_N_ANALYSES,
        help="Number of times to generate a pair of subsets and analyze them."
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

    # Optional: Subset size
    parser.add_argument(
        "-s",
        "--subset_size",
        nargs="+",
        default=DEFAULT_SUBSET_SIZE,
        type=validate_whole_number,
        help="Number of subjects to include in subsets."
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
        # Read in all lines from both files listing all subject paths
        setattr(cli_args, "demo_total", pd.read_csv(cli_args.demo_total_file))

        # For each group, get the path to the directory containing its .pconn
        # files, and the path to the file containing its demographic data
        for group_num in ("group_1_", "group_2_"):

            # Import all of this group's demographics; drop empty columns
            setattr(cli_args, group_num + "demo", pd.read_csv(
                getattr(cli_args, group_num + "demo_file"), na_values=" "
            ).dropna(how="all", axis=1))

        # Get path to .pconn parent directory
        if not hasattr(cli_args, "parent_path"):
            cli_args.parent_path = os.path.commonpath((
                cli_args.group_1_demo_file, cli_args.group_2_demo_file))

        # Ensure that all subjects in a subset can be displayed
        if cli_args.max_rows:
            pd.set_option("display.max_rows", cli_args.max_rows)
        else:
            pd.set_option("display.max_rows", max(cli_args.subset_size))

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
            str(to_validate) + " is not a positive integer."
        )


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
    :return: A dictionary mapping group numbers to randomly selected subsets
    of size subset_size from those groups, where each subset is
    demographically representative of the overall group as shown by chi square
    """
    # Iterate over both groups to get a subset of each
    subsets = {1: None, 2: None}

    # Return value: List of dicts where each one will be a row in a pandas
    # DataFrame. Each row will have subset_size, Euclidean distance, and group
    freq_info = []

    for group_num in subsets.keys():
        group_all = getattr(cli_args, str(group_num).join(("group_", "_demo")))
        other_group = getattr(cli_args, {1: "2", 2: "1"}[group_num].join(
            ("group_", "_demo")))

        # Euclidean distances between each subset and the other group's total
        eu_distances = []

        # Get information about group total dataframe for chi square comparison
        categorical_columns = group_all.select_dtypes(include=["number"])
        if cli_args.exclude:
            continuous_columns = []
            uniques = group_all.nunique()
            for col in categorical_columns:
                if uniques.loc[col] >= cli_args.exclude:
                    continuous_columns.append(categorical_columns.pop(col))
            print("Excluded " + str(len(continuous_columns)) + " demographic "
                  "variables.")
        columns = categorical_columns

        # Compare subset of one group to the other's total
        group_avgs = [other_group[col].mean(skipna=True) for col in columns]

        # Randomly pick subsets until finding 1 that demographically represents
        # the overall group
        loops = 0
        loop_cond = True
        while loop_cond:
            print(str(loops+1) + " subsets of group " + str(group_num)
                  + " randomly generated.", end=" ")

            # Generate subset of group
            subsets[group_num] = group_all.sample(n=subset_size)

            # Compare subset demographics to total demo of other group
            sub_avgs = [subsets[group_num][col].mean(skipna=True)
                        for col in columns]

            # If any column in subset averages has mean of N/A, throw it out
            if any(np.isnan(el) for el in sub_avgs):
                print("Too much missing data, skipping to next.")
            else:

                # Get chi-squared value and Euclidean distance
                loops += 1
                loop_cond = chi_square_difference(
                    other_group, subsets[group_num], columns
                )
                eu_dist = euclidean_distances(group_avgs, sub_avgs)
                eu_distances.append(eu_dist)
                if loop_cond:
                    print("Euclidean distance=" + str(round(eu_dist)))
                else:
                    print("Match found. Subset:")
                    print(subsets[group_num]["subjectkey"])
                    print("Match's Euclidean distance=" + str(round(eu_dist)))

            if cli_args.loops:
                loop_cond = bool(loops < cli_args.loops)

        print("Group " + str(group_num) + " frequency distribution info:")
        info = {"Group": group_num, "Subjects": subset_size,
                "Euclidean Distance": eu_distances[-1]}
        freq_info.append(info)
        print(info, end="\n\n")

    return freq_info[0], freq_info[1]


def chi_square_difference(total_df, subset_df, columns):
    """
    :param total_df: pandas.DataFrame
    :param subset_df: pandas.DataFrame
    :param columns: List of columns in both DataFrame objects to check
    :return: None if every numeric column shows no significant difference
    based on a chi-squared comparison of that column in total_df and subset_df;
    otherwise the numeric column showing a significant difference
    """
    difference = None
    subset_size = len(subset_df.index)
    for column in columns:
        try:
            subset_counts, total_counts = fill_missing_indices_of(
                get_count_distribution_in_column(subset_df, column, subset_size),
                get_count_distribution_in_column(total_df, column, subset_size)
            )
        except ValueError:
            print("\n\nSKIPPING COLUMN " + column + " OF CHI SQ TEST\n\n")
            columns = [col for col in columns if col != column]

        [chi], [p] = get_chi_squared(
            subset_counts, total_counts)  # counts.sort_values(subset_counts.columns[0])

        # Ignore columns which cannot be compared using chi-square test
        if pd.isna(p):
            # print("Ignoring column " + column + " for chi-square test")
            pass
        elif p <= 0.05:
            difference = column
            print("chi = {0:.3g} and p = {1:.3g} for ".format(chi, p) + column,
                  end=". ")
            break
    return difference


def get_count_distribution_in_column(df, column, subset_size):
    """
    :param df: pandas.DataFrame
    :param column: The column of df to count instances of each index value in
    :param subset_size: Number of subjects in subset_df
    :return: A new pandas.DataFrame with one column, mapping each value from
    column in df to the number of occurrences of that value in column in df,
    having removed the space character (" ") as an index
    """
    result = pd.crosstab(index=df[column], columns="count",
                         normalize="columns") * subset_size
    while " " in result.index:
        result = result.drop([" "])
    return result


def get_chi_squared(observed, expected):
    """
    :param observed: pandas.DataFrame
    :param expected: pandas.DataFrame
    :return: Chi-squared comparison between observed and expected DataFrames
    """
    # If both DataFrames have the same indices, get their chi-square comparison
    try:
        chisquare = stats.chisquare(observed, expected)

    # Otherwise, for each index that one DataFrame has and the other does not,
    # add a zero at the missing index
    except ValueError:
        print("Filled missing values with 0")
        observed, expected = fill_missing_indices_of(observed, expected)
        chisquare = stats.chisquare(observed, expected)
    return chisquare


def fill_missing_indices_of(df1, df2):
    """
    :param df1: pandas.DataFrame object with one column
    :param df2: pandas.DataFrame object with one column
    :return: Both dfs, but with a zero value added at each missing index
    """
    for ix in df1.index.difference(df2.index):
        df2.loc[ix] = [0]
    for ix in df2.index.difference(df1.index):
        df1.loc[ix] = [0]
    return df1, df2


def get_logarithmic_fit(x, y):
    """
    Given two arrays of numbers which have equal length, return a logarithmic
    equation to predict y-values given x-values
    :param x: numpy array of numbers
    :param y: numpy array of numbers
    :return: The logarithmic best-fit curve equation for x and y, as a string
    """
    coefs = np.polyfit(np.log(x), y, 1)
    return " + ".join((str(coefs[1]), str(coefs[0]) + " * ln(x)"))


if __name__ == '__main__':
    main()