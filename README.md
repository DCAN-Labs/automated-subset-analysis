# Automated Subset Analysis

This script was written for the ABCD resource paper. It uses Python's `argparse` package so that it can be run from the BASH command line and accept command-line arguments.

## Purpose

Randomly selects and compares pairs of subsets of data from two groups. `automated_subset_analysis.py` makes and saves the subsets, the correlation values between them and the groups, and graph visualizations of those correlations.

### Brief Explanation of Steps

1. `automated_subset_analysis.py` accepts demographic data for 2 datasets. For each dataset, it randomly selects a subset which is not significantly different from the other dataset's total demographics.
1. Once this script has made subsets of both groups, it finds correlations between the average of each subset and the other group. It also finds the correlation between average of both subsets. This process repeats a certain number of times for subsets which include a certain number of subjects. After finding the correlations, the script saves them to `.csv` files.
1. All of the correlation values are then plotted as data points on graph visualizations. Those graphs are saved when `automated_subset_analysis` finishes executing.

For a more detailed explanation, see this document's `Explanation of Process` section.

## Dependencies

1. [Python 3.5.2](https://www.python.org/downloads/release/python-352) or greater

### Python Packages

1. `nibabel`
1. `numpy`
1. `pandas`
1. `plotly`
1. `pprint`
1. `scipy`

## Usage

### Required Arguments (2)

1. `group_1_demo_file` is a path to a text file which contains only a list of paths (1 per line) to the `.nii` files of all subjects in the first group to analyze a subset of.

1. `group_2_demo_file` is a path to a text file which contains only a list of paths (1 per line) to the `.nii` files of all subjects in the second group to analyze a subset of.

Example of a basic call to this script:
```
csv1=/home/user/conan/data/group1_pconn.csv
csv2=/home/user/conan/data/group2_pconn.csv
python3 automated_subset_analysis.py ${csv1} ${csv2}
```

### Optional Arguments (21)

#### File Paths with Default Values (5)

1. `--group-1-avg-file` takes one valid path to a readable `.nii` file containing the average matrix for the entire group 1. By default, this path will be to the `group1_10min_mean.pconn.nii` file in this script's parent folder or in the `--output` folder.

1. `--group-2-avg-file` also takes one valid path to a readable `.nii` file, just like `--group-1-avg-file`. By default, it will point to the `group2_10min_mean.pconn.nii` file in one of the same places.

1. `--matrices-conc-1` takes one path to a readable `.conc` file containing only a list of valid paths to group 1 matrix files. This flag is only needed if your group 1 demographics `.csv` file either does not have a column labeled `'pconn10min'` with paths to matrix files, or if it does include that column but you want to use different paths.

1. `--matrices-conc-2` also takes one path to a readable `.conc` file. It is just like `--matrices-conc-1`, but for group 2.

1. `--output` takes one file path to a directory where all files produced by this script will be saved. If the directory already exists, then this script will add files to it and only overwrite files with conflicting filenames. If not, then this script will create a directory at the `--output` path. If this flag is excluded, then the script will save files to a new subdirectory of the present working directory, `./data/`.

#### Numerical Values for Data Processing (3)

1. `--n-analyses` takes one positive integer, the number of times to generate a pair of subsets and analyze them. For every integer given in the `--subset-size` list, this script will randomly generate `--n-analyses` subsets, creating `--subset-size * --n-analyses` total `.csv` files.

1. `--nan-threshold` takes one floating-point number between 0 and 1. If the percentage of rows with Not-a-Number (NaN) values in the data for either group's demographic file is greater than the `--nan-threshold`, then the script will raise an error. Otherwise, the script will drop every row containing a NaN. The default NaN threshold is `0.1`, meaning that if over 10% of rows have a NaN value, the script will crash.

1. `--subset-size` takes one or more positive integers, the number of subjects to include in subsets. Include a list of whole numbers to generate subsets pairs of different sizes. By default, the subset sizes will be `[50, 100, 200, 300]`. An example of entering a different list of sizes is `--subset-size 100 300 500 1000`.

#### Flags to Skip Steps of Process (2)

1. `--only-make-graphs` takes one or more paths to readable `.csv` files as a parameter. Include this flag to import average correlations data from `.csv` files instead of making any new ones. Given this flag, `automated_subset_analysis.py` only makes graph visualizations of already-existing data.

    - If this flag is included, it must include paths to readable `.csv` files with 2 columns: `Subjects` (the number of subjects in each subset) and `Correlation` (the correlation between each randomly generated subset in that pair).

1. `--skip-subset-generation` takes either no parameters or one path to a readable directory as a parameter. Include this flag to calculate correlations and create the visualization using existing subsets instead of randomly generating new ones. By default, the subsets to use for calculating the correlations between average matrices and producing a visualization will be assumed to exist in the `--output` folder. To load subsets from a different folder, add the path to this flag as a parameter.

#### Notes for Step-Skipping Flags

- If `--skip-subset-generation` or `--only-make-graphs` is included, then `--subset-size` and `--n-analyses` will do nothing.
- If `--only-make-graphs` is included, then `--skip-subset-generation` will do nothing.
- Unless the `--only-make-graphs` flag is used, the `.csv` file(s) with subsets' average correlations will/must be called `correlations_sub1_sub2.csv`, `correlations_sub1_all2.csv`, and/or `correlations_sub2_all1.csv`.

#### Visualization Formatting Arguments (7)

1. `--axis-font-size` takes one positive integer, the font size of the text on both axes of the visualizations that this script will create. If this argument is excluded, then by default, the font size will be `30`.

1. `--fill` takes one parameter, a string that is either `all` or `confidence-interval`. Include this flag to choose which data to shade in the visualization. Choose `all` to shade in the area within the minimum and maximum correlations in the dataset. Choose `confidence-interval` to only shade in the 95% confidence interval of the data. By default, this argument will be `confidence-interval`.

1. `--graph-title` takes one string, the title at the top of all output visualizations and the name of the output `.html` visualization files. Unless this flag is included, each visualization will have one of these default titles:
    - "Correlations Between Average Subsets"
    - "Group 1 Subset to Group 2 Correlation"
    - "Group 1 to Group 2 Subset Correlation"
    - "Correlation Between Unknown Groups"

1. `--hide-legend` takes no parameters. Unless this flag is included, the output visualization(s) will display a legend in the top- or bottom-right corner showing the name of each thing plotted on the graph: data points, average trendline, confidence interval, and/or entire data range.

1. `--marker-size` takes one positive integer to determine the size (in pixels) of each data point in the output visualization. The default size is 5.

1. `--title-font-size` takes one positive integer. It is just like `--axis-font-size`, except for the title text in the visualizations. This flag determines the size of the title text above the graph as well as both axis labels. If this argument is excluded, then by default, the font size will be `40`.

1. `--y-range` takes two floating-point numbers, the minimum and maximum values to be displayed on the y-axis of the graph visualizations that this script will create. By default, this script will automatically set the y-axis boundaries to show all of the correlation values and nothing else.

#### Other Flags (5)

1. `--columns` takes one or more strings. Each should be the name of a column in the demographics `.csv` which contains numerical data to include in the subset correlations analysis. By default, the script will assume that both input demographics `.csv` files have columns of numerical data with these names:
    ```
    demo_comb_income_v2b, demo_ed_v2, demo_prnt_ed_v2b, demo_sex_v2b, ehi_y_ss_scoreb interview_age, medhx_9a, race_ethnicity, rel_relationship, site_id_l
    ```

1. `--correlate-variances` takes no parameters. By default, subset analysis will calculate correlations between subsets'/groups' average values. Include this flag to correlate the subsets' variances instead.

1. `--inverse-fisher-z` takes no parameters. Include this flag to do an inverse Fisher-Z transformation on the matrices imported from the `.pconn` files of the data before getting correlations.

1. `--no-matching` takes no parameters. Include this flag to match subsets on every demographic variable *except* family relationships. Otherwise, subsets will be matched on all demographic variables.

    - By default, `automated_subset_analysis.py`  checks that every subset of one group has the same proportion of twins, triplets, or other siblings as the other group. It also checks that no one in the subset has family members outside the subset. `--no-matching` will skip both checks. 

    - Use this flag if `--subset-size` includes any number under 25, because family relationship matching takes a very long time for small subsets.

1. `--parallel` takes one valid path, the directory containing `automated_subset_analysis.py`. It should be included to simultaneously run multiple different instances of `automated_subset_analysis.py` as a batch command executed by `asa_submitter.py`. Otherwise, this flag is not needed.

For more information, including the shorthand flags for each option, run this script with the `--help` command: `python3 automated_subset_analysis.py --help`

### Advanced Usage Examples

Generate 50 subsets and save a `.csv` file of each, including 10 subsets each of sizes 50, 100, 300, 500, and 1000:

```
python3 automated_subset_analysis.py ${csv1} ${csv2} --subset-size 50 100 300 500 1000 --n-analyses 10
```

Calculate the correlations between average matrices of already-generated subsets in the `./subsets/` folder, then save the correlations and a visualization of them to the `./correls/` folder:

```
python3 automated_subset_analysis.py ${csv1} ${csv2} --skip-subset-generation ./subsets/ --output ./correls/
```

## Explanation of Process

Two `.csv` files, each with demographic data about subjects from a group, are given by the user. One subset is randomly selected from each group repeatedly. The amount of subjects in each subset depends on `--subset-size`, and the number of times that amount is selected depends on `--n-analyses`.

Once a pair has been selected, the script calculates the Euclidean distance between the average demographics of each subset and the average demographics of the whole other group. If the Euclidean distance between the subset and the total is higher than the estimated maximum value for significance (given in the equation calculated by `./src/euclidean_threshold_estimator.py`<sup> 1 </sup>), then another subset is randomly generated and tested for significance. Otherwise, the subset pair is deemed valid and saved to a `.csv` file. The `.csv` has one subset per column and one subject ID per row, excluding the header row which only contains the group number of each subset.

After finding a valid pair of subsets, the script calculates the correlation between the subset of group 1 and the subset of group 2. This correlation value is stored with the number of subjects in both subsets described by the correlation. Once the correlation values are all calculated, each correlation value is saved out to a `.csv` file with a name starting with `correlations`. That `.csv` has two columns. It has one row per subset pair, excluding the header row which contains only the names of the columns: `Subjects` and `Correlation`.

Once the correlation `.csv` files are made, the script will make a graph visualization for each one. That graph will plot how the number of subjects in a subset pair (x-axis) relates to the correlation between the subsets in that pair (y-axis).<sup> 2 </sup>

### Notes

<sup>1</sup> The equation currently used in `automated_subset_analysis.py` to predict significant Euclidean distance threshold using subset size was found using this function call:
```
python3 ./src/euclidean_threshold_estimator.py ./raw/ABCD_2.0_group1_data_10minpconns.csv ./raw/ABCD_2.0_group2_data_10minpconns.csv -con-vars ./automated_subset_analysis_files/continuous_variables.csv --subset-size 1300 1200 1100 1000 900 800 700 600 500 400 300 200 100 90 80 70 60 50 --n-analyses 10
```
The data used to calculate that equation can be found in `./src/euclidean_threshold_estimate_data/est-eu-thresh-2019-12-12`.

<sup>2</sup> The output visualization will include:
1. Each correlation value as 1 data point,
1. A trendline using the average correlation values of each subset size,
1. A shaded region showing a certain data range (confidence interval or all data), and
1. A legend to identify all of these parts (unless `--hide-visualization` is used).

## Metadata

Information about this `README` file:

- Created by Greg Conan, 2019-10-03
- Updated by Greg Conan, 2020-02-06
