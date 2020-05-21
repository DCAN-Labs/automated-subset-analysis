#! /usr/bin/env python

"""
Euclidean distance threshold finder for automated_subset_analysis.py
Greg Conan: conan@ohsu.edu
Created 2019-09-17
Updated 2020-05-05
"""

##################################
#
# Script to randomly select pairs of subsets of data of a given size and get
# the Euclidean distance threshold below which the subsets are not significantly
# different from the other's total population, then get an equation to estimate
# that threshold given group size
#
##################################

# Imports
import argparse
from conan_tools import *
import numpy as np
import os
import pandas as pd
import scipy

# Constants
ARG_NAMES = ["group_1_demo_file", "group_2_demo_file", "columns", "n_analyses", 
             "nan_threshold", "output", "subset_size", "calculate",
             "continuous_variables", "no_matching"]
FREQ_COL_NAMES = ["Group", "Subjects", "Euclidean Distance"]
PWD = get_pwd()


def main():
    # Store and print the date and time when this script started running
    starting_timestamp = get_and_print_timestamp_when(sys.argv[0], "started")

    # Get and validate all command-line arguments from user
    cli_args = get_cli_args(
        ("Script to randomly select pairs of subsets of data of a given size "
         "and get the Euclidean distance threshold below which the subsets "
         "pass a statistical test."), ARG_NAMES, PWD, validate_cli_args
    )

    # Get all subset frequencies and maximums
    frequency_info, freq_maximums = do_all_freq_analyses(cli_args,
                                                         FREQ_COL_NAMES)
    all_freq_info = pd.DataFrame(frequency_info, columns=FREQ_COL_NAMES)
    maximums_df = pd.DataFrame(freq_maximums, columns=FREQ_COL_NAMES)

    # Print a line, then print all subset frequencies and maximums to show user
    print(sum([len(f)+2 for f in FREQ_COL_NAMES]) * "_")
    print("\nComplete frequency distribution:\n{}\n"
          "\nAll maximum Euclidean distances:\n{}\n"
          .format(all_freq_info.to_string(index=False),
                  maximums_df.to_string(index=False)))

    # Calculate and print the logarithmic best-fit curve of maximum Euclidean 
    # distances for each number of subjects
    log_fn, log_fn_str = get_logarithmic_fit(maximums_df["Subjects"],
                                             maximums_df["Euclidean Distance"])
    print("Logarithmic best-fit equation: {}\n"
          "Correlation between the expected values from the best-fit equation "
          "and the subsets' max Euclidean distances: {}".format(log_fn_str, 
              np.corrcoef(maximums_df["Euclidean Distance"], 
                          [log_fn(y) for y in maximums_df["Subjects"]])[0, 1]))

    # Print when this script started and finished
    print(starting_timestamp)
    get_and_print_timestamp_when(sys.argv[0], "finished")


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    try:
        # Import each group's demographic data, minus any rows with NaNs
        for gp_num in (1, 2):
            gp_str = "group_{}_demo".format(gp_num)
            setattr(cli_args, gp_str, get_group_demographics(
                cli_args, gp_num, gp_str + "_file", parser
            ))
        return cli_args
    except OSError as e:
        parser.error(str(e))


def do_all_freq_analyses(cli_args, freq_cols):
    """
    Compute, collect, and return all Euclidean distances
    :param cli_args: argparse namespace with all command-line arguments
    :param freq_cols: List of strings which are the names of the keys/columns
                      in frequency info data structures
    :return: Tuple of 2 lists of dictionaries: The first list has one dict for
             every analysis, and the second list has only the dicts which
             have the maximum Euclidean distances for their subset size
    """
    # Return values: Lists of one dict per analysis where each dict includes
    # a Euclidean distance
    frequency_info = []  # Distances for all analyses
    freq_maximums = []   # One distance per subset size; only maximum Euclid dist

    def print_frequencies(label, freqs):
        print("{}:\n{}".format(label, pd.DataFrame(freqs, columns=freq_cols
                                                   ).to_string(index=False)))

    # Run the number of analyses in --n_analyses for each subset size given
    progress = track_progress(cli_args.n_analyses, cli_args.subset_size)
    for subset_size in cli_args.subset_size:
        print("Running {} analyses of randomly selected groups of {} subjects."
              .format(cli_args.n_analyses, subset_size))
        for i in range(cli_args.n_analyses):

            # Select subsets and get the Euclidean distances for each subset
            start_time = now()
            freq_1, freq_2 = get_random_subsets(cli_args, subset_size,
                                                loop_checker)
            frequency_info.append(freq_1)
            frequency_info.append(freq_2)

            # Find the maximum Euclidean distance for each number of subjects
            for freq in (freq_1, freq_2):
                if (not freq_maximums
                        or freq["Subjects"] != freq_maximums[-1]["Subjects"]):
                    print("Adding a maximum for {}-subject subset."
                          .format(freq["Subjects"]))
                    freq_maximums.append(freq)
                elif (freq["Euclidean Distance"] >
                      freq_maximums[-1]["Euclidean Distance"]):
                    freq_maximums[-1] = freq
            progress = update_progress(
                progress, "calculating maximum Euclidean distances",
                subset_size, start_time
            )
            print_frequencies("All maximums so far", freq_maximums)

    return frequency_info, freq_maximums


def loop_checker(loops, eu_dist, group_number, subset, diff_group,
                 column_names, continuous_column_names, col_widths):
    """
    Function to pass as an object parameter to randomly_generate_subsets. This
    checks whether a given subset demographically matches the other group using
    either a chi-squared test or a t-test. Print info to user every 10 subjects
    :param loops: Integer which is how many subsets were randomly generated
    :param eu_dist: Float which is the Euclidean distance between the average  
                    of the subset's columns and that of the other group's
    :param group_number: Integer which is the group number (1 or 2)
    :param subset: pandas.DataFrame which has the subset
    :param diff_group: pandas.DataFrame which has the other group
    :param column_names: List of Strings where each names the column in the
                         same position in subset.columns and diff_group.columns
    :param continuous_column_names: List of Strings where each names a column
                                    with continuous data in the DataFrames 
    :param col_widths: List of integers where each is the maximum length of its
                       corresponding column (max for all subsets generated)
    :return: True if subset and diff_group significantly differ; else False
    """ 
    diff_col, stat, p = get_statistical_difference(
        diff_group, subset, column_names, continuous_column_names
    )

    keep_looping = bool(diff_col) # If no column significantly differs, stop
    FLOAT_FORMAT = "{0:.3g}"      # String format to print a decimal

    # Tell user if a match was found 
    if not keep_looping:
        print("Match found. Its Euclidean distance="
              + FLOAT_FORMAT.format(eu_dist))

    # If no match was found, print the next row of the table to the
    # command line for the user to track the script's progress
    elif not loops % count_digits_of(loops):
        if col_widths:
            print(fit_strings_to_width((
                loops, group_number, FLOAT_FORMAT.format(stat),
                FLOAT_FORMAT.format(p), diff_col, FLOAT_FORMAT.format(eu_dist)
            ), col_widths))

    return keep_looping


def get_random_subsets(cli_args, subset_size, loop_check_fn, eu_thresh=None):
    """
    :param cli_args: argparse namespace with all command-line arguments
    :param subset_size: Number of subjects per subset
    :return: Tuple of 2 dictionaries mapping group numbers to randomly selected 
             subsets of size subset_size from those groups, where each subset
             demographically represents the overall group, shown by chi square
    """
    # Iterate over both groups to get a subset of each
    subsets = {1: None, 2: None}

    # Return value: List of dicts where each one will be a row in a pandas
    # DataFrame. Each row will have subset_size, Euclidean distance, and group
    freq_info = []

    # Get a random subset of each group which does not significantly differ
    # demographically from the other group
    for group_num in subsets.keys():
        group_all = getattr(cli_args, str(group_num).join(("group_", "_demo")))
        other_group = getattr(cli_args, {1: "2", 2: "1"}[group_num].join(
            ("group_", "_demo")))

        eu_dist = randomly_select_subset(group_all, group_num, subset_size,
                                         other_group, cli_args, loop_check_fn)
        info = {"Group": group_num, "Subjects": subset_size,
                "Euclidean Distance": eu_dist} 
        freq_info.append(info)

    # Print and return the current subsets' Euclidean distances & etc
    print("Frequency distribution so far:\n{}".format(pd.DataFrame(
            freq_info, columns=FREQ_COL_NAMES
    ).to_string(index=False)))
    return freq_info[0], freq_info[1]


def get_statistical_difference(total_df, subset_df, columns, continuous_cols):
    """
    :param total_df: pandas.DataFrame
    :param subset_df: pandas.DataFrame
    :param columns: List of columns in both DataFrame objects to check
    :param continuous_cols: List of names of columns containing continuous data
    :return: Tuple with 3 elements:
        1) None if every numeric column shows no significant difference in a
           statistical comparison of that column in total_df and subset_df;
           otherwise, the numeric column showing a significant difference
        2) chi-squared or t statistic
        3) p-value from statistical test
    """
    difference = None  # Return value
    subset_size = len(subset_df.index)
    for col in columns:

        # Independent samples t-test if column is continuous
        if col in continuous_cols:
            stat, p = scipy.stats.ttest_ind(
                total_df[col], subset_df[col],
                nan_policy="raise", equal_var=False
            )

        # Chi-squared comparison if column is categorical
        else:
            # Get count distribution of each value in each column
            try:
                subset_counts, total_counts = fill_missing_indices_of(
                    get_count_distribution_in_col(subset_df, col, subset_size),
                    get_count_distribution_in_col(total_df, col, subset_size)
                )
            except ValueError:
                print("\n\nSKIPPING COLUMN {} OF STAT TEST\n\n".format(col))
                columns = [column for column in columns if col != column]

            # Get chi-squared statistic
            [stat], [p] = get_chi_squared(subset_counts, total_counts)

        if pd.isna(p) or p <= 0.05:
            difference = col
            break
    return difference, stat, p


def get_count_distribution_in_col(df, column, subset_size):
    """
    :param df: pandas.DataFrame
    :param column: The column of df to count instances of each index value in
    :param subset_size: Number of subjects in subset_df
    :return: A new pandas.DataFrame with one column, mapping each value from
             column in df to the number of occurrences of that value in column 
             in df, having removed the space character (" ") as an index
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
        chisquare = scipy.stats.chisquare(observed, expected)

    # Otherwise, for each index that one DataFrame has and the other does not,
    # add a zero at the missing index
    except ValueError:
        print("Filled missing values with 0")
        observed, expected = fill_missing_indices_of(observed, expected)
        chisquare = scipy.stats.chisquare(observed, expected)
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
    equation to predict y-values given x-values - as a function and as a string
    :param x: numpy array of numbers
    :param y: numpy array of numbers
    :return: A tuple of the logarithmic best-fit curve equation for x and y in
             two forms: as a function and as a string
    """
    coefs = np.polyfit(np.log(x), y, 1)
    return (lambda x: natural_log(x, coefs),
            "f(x) = {} * ln(x) + {}".format(coefs[0], coefs[1]))


if __name__ == '__main__':
    main()
