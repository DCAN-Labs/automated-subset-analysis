#! /usr/bin/env python3

"""
Pairwise Correlator
Greg Conan: gconan@umn.edu
Created 2020-08-03
Updated 2020-08-14
"""

##################################
#
# Compare 2 groups of pconns subject-by-subject using .conc files
#
##################################

# Imports
import argparse
from conan_tools import (
    count_digits_of, get_cli_args, get_group_demographics, get_pwd,
    GP_DEMO_FILE, GP_MTR_FILE, load_matrix_from, now, PATH_EXA, PATH_RUSH,
    read_file_into_list, rename_exacloud_path, spearman_rho, track_progress,
    update_progress, valid_conc_file, valid_output_dir, valid_readable_file
)
import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go
import subprocess
import time
import warnings

# Constants: Meanings of demographic data values, this script's parent folder,
# and the name of the demographic .csv column with subject IDs
DATA_LABELS = {
    "site_id_l": "Site ID", "race_ethnicity": "Race/Ethnicity",
    "interview_age": "Age in Months at Interview", "demo_sex_v2b": "Sex", 
    "demo_prnt_ed_v2b": "Highest Grade Level Passed by Parent"
}
DATA_VALUES = {
    "site_id_l": {x: x for x in range(1, 22)},
    "interview_age": {x: "{} months".format(x) for x in range(104, 134)},
    "race_ethnicity": {1: "White", 2: "Black", 3: "Hispanic",
                       4: "Asian", 5: "Other"},
    "demo_sex_v2b": {1: "Male", 2: "Female"},
    "demo_prnt_ed_v2b": {x: x for x in range(1, 22)}
}
PWD = get_pwd()
SUBJ_ID_COL = "id_redcap"


def main():
    # Get command line arguments from user
    cli_args = get_cli_args("Pairwise connectivity matrix correlator", (
        "axis_font_size", GP_DEMO_FILE.format(1), GP_MTR_FILE.format(1),
        GP_MTR_FILE.format(2), "graph_title", "hide_legend", "marker_size",
        "nan_threshold", "only_make_graphs", "output", "spearman_rho",
        "title_font_size", "y_range"
    ), PWD, validate_cli_args)

    # Import each subject's matrices from both .conc files,
    # compare them, and save correlations into .csv files
    if not cli_args.only_make_graphs:
        all_IDs, matrix_paths = import_concs(cli_args)
        correlate_pairwise(all_IDs, matrix_paths, cli_args)

    # Skip correlation and only make graphs, if user said to
    else:
        prep_data_and_make_visualizations(cli_args)


def validate_cli_args(cli_args, parser):
    """
    Check that all command line arguments will allow this script to work.
    :param cli_args: argparse namespace with all given command-line arguments
    :param parser: argparse ArgumentParser to raise error if anything's invalid
    :return: Validated command-line arguments argparse namespace
    """
    # Matrix path .conc files are required
    if not (cli_args.only_make_graphs or (
            getattr(cli_args, "matrices_conc_1", None)
            and getattr(cli_args, "matrices_conc_2", None))):
        parser.error("Please include the --matrices-conc-1 and "
                     "--matrices-conc-2 arguments. Even though they are shown "
                     "as optional, this script requires them.")
    
    # Import all of this group's demographics; drop empty columns
    setattr(cli_args, GP_DEMO_FILE.format(1)[:-5], get_group_demographics(
        cli_args, 1, GP_DEMO_FILE.format(1), parser
    ))

    return cli_args


def import_concs(cli_args):
    """
    :param cli_args: argparse namespace with all given command-line arguments
    :return: Tuple of 2 objects, a list of subject ID strings in order and a
             a dict mapping the group number to all of that group's file paths
    """
    # Local variables: String for error to show if .conc files' subjects do not
    # match exactly, and also collections of subject IDs/paths to return
    MATCH_ERR = ("Line {0} of {1} has subject {2}, but line {0} of {3} has "
                 "subject {4}. Both lines should have the same subject.")
    CONC_1 = os.path.basename(getattr(cli_args, GP_MTR_FILE.format(1)))
    CONC_2 = os.path.basename(getattr(cli_args, GP_MTR_FILE.format(2)))
    all_IDs = []
    matrix_paths = dict()

    # Import file paths from .conc files
    for i in (1, 2):
        matrix_paths[i] = read_file_into_list(getattr(cli_args,
                                                      GP_MTR_FILE.format(i)))

    # Ensure that both .conc files have the same number of matrix file paths
    if len(matrix_paths[1]) != len(matrix_paths[2]):
        raise ValueError("Both matrix .conc files must have exactly the same "
                         "number of file paths, 1 per line.")  
    
    # Validate all paths in .conc files
    for i in range(len(matrix_paths[1])):
    
        # Ensure that both matrix file path lists have the same subjects
        id1 = get_subject_id_from(matrix_paths[1][i])
        id2 = get_subject_id_from(matrix_paths[2][i])
        if id1 == id2:
            all_IDs.append(id1)
        else:
            raise ValueError(MATCH_ERR.format(i, CONC_1, id1, CONC_2, id2))
            
        # Ensure that both matrix file path lists point to real files
        for j in (1, 2):
            if not os.path.exists(matrix_paths[j][i]):
                matrix_paths[j][i] = rename_exacloud_path(os.path.join(
                    PATH_RUSH, matrix_paths[j][i]
                ))
            assert os.path.exists(matrix_paths[j][i])
                    
    return all_IDs, matrix_paths


def get_subject_id_from(path):
    """
    Accepts the path to a subject's matrix file, and gets that subject's ID
    :param path: String which is a valid path to a CIFTI2 matrix file
    :return: String which is the subject ID which was in path
    """
    sub_id_pos = path.find("INV")              # Where in path ID starts
    return path[sub_id_pos:(sub_id_pos + 11)]  # through where in path ID ends


def correlate_pairwise(all_IDs, matrix_paths, cli_args):
    """
    Calculate and save the correlation between each subject's conn matrix pair
    :param all_IDs: List of strings which are all the subject IDs in order
    :param matrix_paths: Dictionary mapping each group's number to all of its
                         connectivity matrix file paths
    :param cli_args: argparse namespace with all given command-line arguments
    :return: N/A                    
    """
    # Local variables
    mx_pair = dict()  # Each pair of matrices for each subject
    num_subjs = len(matrix_paths[1])  # Total number of subjects
    corr_fn = (spearman_rho if cli_args.spearman_rho else
               lambda x, y: np.corrcoef(x, y)[0, 1])  # Correlation function
    progress = track_progress(1, [1] * num_subjs)  # Estimate time this'll take
    started_at = now()  # Timestamp where matrix correlation started

    # Correlate every matrix with its corresponding matrix and write to file
    filepath = cli_args.only_make_graphs[0] if cli_args.only_make_graphs else (
        os.path.join(cli_args.output, "correlations_out_{}.csv".format(now()))
    ).replace(" ", "-").replace("\\", "-").replace(":", "_").replace(".", "-")
    if not os.path.exists(filepath):  # Make a new file if one doesn't exist
        with open(filepath, 'w+') as outfile: 
            outfile.write("Subject ID,Correlation\n")
    for i in range(num_subjs):
        for j in (1, 2):
            mx_pair[j] = load_matrix_from(matrix_paths[j][i]).flatten()
        correl = corr_fn(mx_pair[1], mx_pair[2])
        print("Subject {}'s correlation is {}"
              .format(all_IDs[i], correl))
        with open(filepath, "a+") as outf:
            outf.write("{},{}\n".format(all_IDs[i], correl))
        progress = update_progress(
            progress, "calculating pairwise correlations", 1, started_at
        )


def prep_data_and_make_visualizations(cli_args):
    """
    :param cli_args: argparse namespace with all given command-line arguments
    :return: N/A
    """
    for correl_csv in cli_args.only_make_graphs:
    
        # Add correlations column to demographics .csv
        all_correls = pd.read_csv(correl_csv) # getattr(cli_args, GP_DEMO_FILE.format(1))
        all_correls["Subject ID"] = (all_correls["Subject ID"]
                                     .apply(lambda x: "NDAR_" + x))
        full_df = all_correls.merge(getattr(
            cli_args, GP_DEMO_FILE.format(1)[:-5]
        ), left_on="Subject ID", right_on=SUBJ_ID_COL, how="inner").dropna()

        # Split .csv into 3 groups by correl: high and low outliers vs. med
        outliers_lo = full_df.loc[full_df["Correlation"] <
                                  np.percentile(full_df["Correlation"], 5)]
        outliers_hi = full_df.loc[full_df["Correlation"] > 
                                  np.percentile(full_df["Correlation"], 95)]
        med_range = full_df.drop(outliers_hi.index).drop(outliers_lo.index)
        corr_gps = (outliers_lo, med_range, outliers_hi)
        bar_x_vals = ["Low-Correlation Outliers", "Middle Range",
                      "High-Correlation Outliers"]

        make_visualizations(cli_args, full_df, bar_x_vals, corr_gps)

    
def count_val_in_col(the_df, col_name, val):
    """
    :param the_df: pandas.DataFrame with a col_name column
    :param col_name: String naming a column in the_df
    :param val: Element to count the occurrences of in the_df in col_name
    :return: Float, the proportion of the_df rows which have val in the column
    """
    return len(the_df.loc[the_df[col_name] == val])/len(the_df)


def get_y_range(cli_args, full_df):
    """
    Get y-axis range to display on figure
    :param cli_args: argparse namespace with all given command-line arguments
    :param full_df: pandas.DataFrame with a "Correlation" column of floats
                    ranging from 0 to 1
    :return: Tuple of 2 floats, min and max y-values to show on figure
    """
    y_min = min(full_df["Correlation"])
    y_max = max(full_df["Correlation"])
    if y_max > 0.98:
        y_max = min(full_df.nlargest(5, "Correlation")["Correlation"])
    y_step = (y_max - y_min)/10
    return (cli_args.y_range if cli_args.y_range
            else y_min - y_step, y_max + y_step)


def make_visualizations(cli_args, full_df, bar_x_vals, corr_groups):
    """
    :param cli_args: argparse namespace with all given command-line arguments
    :param full_df: pandas.DataFrame with a "Correlation" column of floats
                    ranging from 0 to 1
    :param bar_x_vals: List of strings where each labels a bar plot column
    :param corr_groups: List of 3 pandas.DataFrame objects:
                        1) outliers at or under the 5th percentile
                        2) data between the 5th and 95th percentiles
                        3) outliers at or over the 95th percentile
    :return: N/A
    """
    outliers_lo, med_range, outliers_hi = corr_groups
    y_range = get_y_range(cli_args, full_df)

    # Make 1 visualization per column/variable in demographics .csv
    for column in DATA_LABELS.keys():
        box_plot = format_plot_y_axis(go.Figure())
        bars = []

        # Group and sort the correlation dataframes for display
        correls_grouped = full_df.groupby([column]).groups
        sorted_keys = list(correls_grouped.keys())
        sorted_keys.sort()
        for uniq in sorted_keys:
            tracename = DATA_VALUES[column][int(uniq)]

            # Add each trace to the box plot 
            box_plot.add_trace(go.Box(
                jitter=0.5, name=tracename,
                marker_size=2, boxpoints="all", pointpos=-2,
                y=(full_df.loc[full_df[column] == uniq].reset_index(
                   drop=True)["Correlation"].tolist())
            )) 

            # Add each trace to the bar plot
            bar_y_vals = [count_val_in_col(outliers_lo, column, uniq),
                          count_val_in_col(med_range, column, uniq),
                          count_val_in_col(outliers_hi, column, uniq)]
            bars.append(go.Bar(
                name=tracename, x=bar_x_vals, y=bar_y_vals,
                textposition="auto",
                text=["{}: {:.2%}".format(DATA_VALUES[column][uniq], x)
                        for x in bar_y_vals]
            ))

        # Format and save each box plot figure
        box_plot.update_layout(
            title=cli_args.graph_title + ": " + column, showlegend=False,
            yaxis=dict(title="Correlation", range=y_range, autorange=False),
            xaxis=dict(title=DATA_LABELS[column], dtick=1) #, tickangle=45)
        )
        save_plot_to_html(box_plot, cli_args.output, column + "_by_demo")

        # Format and save each bar plot figure
        bar_plot = format_plot_y_axis(go.Figure(data=bars))
        bar_plot.update_layout(
            title=cli_args.graph_title + " Outlier Demographics: " + column,
            showlegend=True, barmode="stack",
            yaxis=dict(autorange=True, # range=[0,1],
                       title=DATA_LABELS[column] + " by Outlier Grouping"),
            xaxis=dict(title="Outlier Groupings (2 std. dev.)", dtick=1)
        )
        save_plot_to_html(bar_plot, cli_args.output, column + "_by_outlier")


def format_plot_y_axis(a_plot):
    """
    :param a_plot: plotly.graph_objects.Figure 
    :return: a_plot, but with the y-axis formatted in its layout
    """
    a_plot.update_layout(yaxis=dict(
        showgrid=True, gridwidth=1, dtick=0.1, showticklabels=True,
        tickfont=dict(size=18, color="black", family="Courier New, monospace")
    ))
    return a_plot


def save_plot_to_html(plot_to_save, out_dir, fname):
    """
    :param plot_to_save: plotly.graph_objects.Figure
    :param out_dir: String, a valid path pointing to a real directory
    :param fname: String, the base name of the file to save (without extension)
    :return: N/A
    """
    plot_to_save.write_html(os.path.join(out_dir, fname + ".html"))


if __name__ == "__main__":
    main()
