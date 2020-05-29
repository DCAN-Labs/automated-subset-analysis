#! /usr/bin/env python3

"""
Multiple visualization creator
Greg Conan: conan@ohsu.edu
Created 2020-05-11
Updated 2020-05-11
"""

##################################
#
# Script to create multiple visualizations using automated subset analysis
#
##################################

# Imports
import argparse
from conan_tools import GP_DEMO_FILE, get_cli_args, get_pwd, as_cli_arg
import os
import subprocess
import time


PWD = get_pwd()
ASA_DIR = (os.path.abspath(os.path.dirname(PWD))
           if os.path.split(PWD)[1] == "src" else PWD)
ASA_SCRIPT = os.path.join(ASA_DIR, "automated_subset_analysis.py")

def main():

    # Get arguments from user via command line
    args_to_get = [
        GP_DEMO_FILE.format(1), GP_DEMO_FILE.format(2), "axis_font_size", 
        "fill", "hide_legend", "marker_size", "n_analyses", 
        "output", "subset_size", "title_font_size", "y_range"
    ]
    cli_args = get_cli_args(
        ("Script to create multiple subset reliability visualizations. Note "
         "that --output must be the same folder that the .csv files are in."),
        args_to_get, PWD
    ) 

    # Make a visualization of every .csv file with correlations
    for csv_file in os.scandir(cli_args.output):
        if csv_file.is_file() and os.path.splitext(csv_file.path)[1] == ".csv":
            print(csv_file.path)
            graph_title = get_graph_title(csv_file)
            run_ASA_script(args_to_get, csv_file.path, cli_args, graph_title)


def get_graph_title(csv_file):
    """
    :param csv_file: os.DirEntry of a correlations .csv file 
    :return: Title of visualization generated from csv_file based on filename
    """
    csv_filename = os.path.splitext(csv_file.name)[0]
    return "<br>".join((
        {"pconn": "Parcel. Connectivity Matrix",
         "curv_map": "Curvature Maps",
         "cortical": "Cortical Thickness",
         "myelin_cortical": "Myelin Maps to Cortical Thickness",
         "myelin_myelin": "Myelin Maps",
         "sulcus": "Sulcus Depth"}[csv_filename[:-10]],
        {"sub2_all1": "Group 1 to Subset 2",
         "sub1_all2": "Group 2 to Subset 1",
         "sub1_sub2": "Both Subsets"}[csv_filename[-9:]]
    ))


def run_ASA_script(to_include, csv_file, cli_args, graph_title):
    """
    Run automated_subset_analysis script with --only-make-graphs flag
    :param to_include: List of strings which name the args for to_run to use
    :param csv_file: String with the path to a correlations .csv file
    :param cli_args: argparse namespace with all command-line arguments
    :return: N/A
    """
    cmd = ["python3", ASA_SCRIPT]
    for arg in to_include:
        val = getattr(cli_args, arg, None)
        if val:
            if "demo_file" not in arg:
                cmd.append(as_cli_arg(arg))
            if val is not True:
                if isinstance(val, list):
                    val = " ".join((str(v) for v in val))
                cmd.append(str(val))
    print("Graph title: {}".format(graph_title))
    print("Filename: {}".format(csv_file))
    cmd += ["--only-make-graphs", csv_file, "--graph-title", graph_title.join(("'", "'"))]
    os.system(" ".join(cmd))


if __name__ == '__main__':
    main()
