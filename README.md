# Automated Subset Analysis

These scripts were written for the ABCD resource paper. They use Python's `argparse` package so that they can be run from the BASH command line and accept command-line arguments. They should both be run from within the `automated_subset_analysis` directory.

## Purpose of Scripts

- `automated_subset_analysis.py`

Randomly selects and compares pairs of subsets of data. It outputs the subsets, the correlation values between them, and a graph visualization of those correlations.
 
This script accepts demographic data for 2 datasets. For each dataset, it randomly selects a subset which is not significantly different from the other dataset's total demographics. Once this script has the  find the correlation between the average matrices of each subset, and plot all of those correlations for each subset size

- `euclidean_threshold_estimator.py`

Uses the same code as `automated_subset_analysis.py` to randomly generate subsets, but instead of analyzing them, it estimates the maximum Euclidean distance required to create a subset which is not significantly different from a total set. Prints out a list of these estimated Euclidean distances, including a logarithmic regression equation to estimate the Euclidean distance threshold for a specific number of subjects. (*Warning: The calculated regression equation may not be accurate yet*)

## Dependencies

1. [Python 3.5.2](https://www.python.org/downloads/release/python-352) or greater

### Python Packages

1. `nibabel`
1. `numpy`
1. `pandas`
1. `plotly`
1. `scipy`

## Usage: `automated_subset_analysis.py`

### Required Arguments

1. `group_1_demo_file` is a path to a text file which contains only a list of paths (1 per line) to the .pconn files of all subjects in the first group to analyze a subset of.

1. `group_2_demo_file` is a path to a text file which contains only a list of paths (1 per line) to the .pconn files of all subjects in the second group to analyze a subset of.

Example of a basic call to this script:
```
python3 automated_subset_analysis.py /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp1_10min_pconn.csv /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp2_10min_pconn.csv
```

### Optional Arguments: File Paths with Default Values

1. `--output` takes one file path to a directory where all files produced by this script will be saved. If the directory already exists, then this script will add files to it and only overwrite files with conflicting filenames. If not, then this script will create a directory at the `--output` path. If this flag is excluded, then the script will save files to a subdirectory called `subset-analysis-output` of the present working directory.

1. `--parent_path` takes one file path to the parent directory of a directory containing .pconn files. The `group_1_demo_file` and the `group_2_demo_file` should both have a column named `pconn_10min` where each row has a path to the `.pconn` file for the subject in that row. That path should be a relative path from the `--parent_path` directory, the value given as this argument. So, `--parent_path` must be a valid path to a readable directory. By default, if this flag is excluded, then `--parent_path` will be the longest common path directory of `group_1_demo_file` and `group_2_demo_file`.

### Optional Arguments: Numerical Values

1. `--n_analyses` takes one positive integer, the number of times to generate a pair of subsets and analyze them. For every integer given in the `--subset_size` list, this script will randomly generate `--n_analyses` subsets, creating `--subset_size` * `--n_analyses` subset `.csv` files in total.

1. `--subset_size` takes one positive integer or a list of positive integers, the number of subjects to include in subsets. Include a list of whole numbers to generate subsets pairs of different sizes. By default, the subset sizes will be `[50, 100, 200, 300]`.

### Optional Arguments: Runtime Options

1. `--skip_subset_generation` takes either no parameters or one path to a readable directory as a parameter. Include this flag to calculate correlations and create the visualization using existing subsets instead of randomly generating new ones. By default, the subsets to use for calculating the correlations between average matrices and producing a visualization will be assumed to exist in the `--output` folder. To load subsets from a different folder, add the path to this flag as a parameter.

1. `--only_make_graphs` takes either no parameters or one path to a readable `.csv` file as a parameter. Include this flag to import data from a `.csv` file instead of making a new one, just to make a graph visualization of already-existing data. If this flag is included, it must be a path to a readable `.csv` file with 2 columns: `Subjects` (the number of subjects in each subset) and `Correlation` (the correlation between each randomly generated subset in that pair). If the `.csv` was already made by this script, it will be called `correlations.csv`. 

1. `--inverse_fisher_z` takes no parameters. Include this flag to do an inverse Fisher-Z transformation on the matrices imported from the `.pconn` files of the data before getting correlations.

If `--skip_subset_generation` or `--only_make_graphs` is included, then `--subset_size` and `--n_analyses` will do nothing. If `--only_make_graphs` is included, then `--skip_subset_generation` will do nothing.

For more information, including the shorthand flags for each option, run this script with the `--help` command: `python3 automated_subset_analysis.py --help`

### Advanced Usage Examples

Generate 50 subsets and make a `.csv` file of each, including 10 subsets each of sizes 50, 100, 300, 500, and 1000:

```
python3 automated_subset_analysis.py /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp1_10min_pconn.csv /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp2_10min_pconn.csv --subset_size 50 100 300 500 1000 --n_analyses 10 
```

Calculate the correlations between average matrices of already-generated subsets in the `./subsets` folder, then save the correlations and a visualization of them to the `./correls` folder:

```
python3 automated_subset_analysis.py /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp1_10min_pconn.csv /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp2_10min_pconn.csv --skip_subset_generation subsets --output correls
```

## Usage: `euclidean_threshold_estimator.py`

This script uses many of the same arguments as `automated_subset_analysis.py`, which work exactly the same way:

- `--group_1_demo_file` 
- `--group_2_demo_file`
- `--parent_path`
- `--n_analyses`
- `--subset_size`

However, some of the optional arguments from `automated_subset_analysis.py` are not used at all by `euclidean_threshold_estimator.py`:

- `--output`
- `--only_make_graphs`
- `--skip_subset_generation`
- `--inverse_fisher_z`

And conversely, `euclidean_threshold_estimator.py` takes some optional arguments that `automated_subset_analysis.py` does not:

- `--exclude` takes one positive integer representing the threshold at which to distinguish categorical and continuous demographic variables. If this flag is included, then any demographic variable with at least this many unique values will be excluded from chi-squared comparison on the basis of being a continuous instead of a categorical variable. Otherwise, every variable will be treated as categorical and included in the chi-squared comparison.

- `--max_rows` takes one positive integer, the maximum number of rows to display in the terminal when printing the contents of a Pandas object.

- `--loops` takes one positive integer representing how many times to generate another subset before giving up on finding a subset which passes the chi-squared test and printing the frequency distribution of already-generated subsets. If this flag is excluded, then the script will keep generating subsets until one of them passes.

For more information, including the shorthand flags for each option, run this script with the `--help` command: `python3 euclidean_threshold_estimator.py --help`

## Explanation of Process

- ### `automated_subset_analysis.py`

Two `.csv` files, each with demographic data about subjects from a group, are given by the user. One subset is randomly selected from each group repeatedly; the amount of subjects in each group depends on `--subset_size` and the number of times that amount is selected depends on `--n_analyses`.

Once a pair has been selected, the script calculates the Euclidean distance between the average demographics of each subset and the average demographics of the whole other group. If the Euclidean distance between the subset and the total is higher than the estimated maximum value for significance (given in the equation calculated by `euclidean_threshold_estimator.py`*), then another subset is randomly generated and tested for significance. Otherwise, the subset pair is deemed valid and saved to a `.csv` file. The `.csv` has one subset per column and one subject ID per row, excluding the header row which contains only the group number of each subset.

After finding a valid pair of subsets, the script calculates the correlation between the subset of group 1 and the subset of group 2. This correlation value, and the number of subjects in both subsets described by the correlation, are stored. Once they are all calculated, they are all saved out to a file called `correlations.csv`. That `.csv` has two columns, and one row per subset pair, excluding the header row which contains only the names of the columns: `Subjects` and `Correlation`.

Once `correlations.csv` is made, the script will then make a visualization plotting how the number of subjects in a subset pair (x-axis) relates to the correlation between the subsets in that pair (y-axis). 

- ### `euclidean_threshold_estimator.py`

Two `.csv` files, each with demographic data about subjects from a group, are given by the user. One subset is randomly selected from each group repeatedly; the amount of subjects in each group depends on `--subset_size` and the number of times that amount is selected depends on `--n_analyses`.

Once a pair has been selected, a chi-squared test of significance is run comparing the average demographics of each subset to the average demographics of the whole other group. If the chi-squared test shows a significant difference, then another subset is randomly generated and tested for significance. Otherwise, the subset pair is deemed valid. If the `--exclude` flag is used, then only those demographic variables with fewer unique values than the value of `--exclude` are included in the chi-squared test. 

After finding a valid pair of subsets, the script calculates the Euclidean distance between each subset's average demographics and the whole other group's average demographics. It collects these Euclidean distances for however many subset sizes and number of times that the user specified.

Once all of the Euclidean distances are calculated, the script takes the maximum Euclidean distance for each subset size and calculates a logarithmic regression equation to predict the maximum Euclidean distance given only a subset size. This equation is used by `automated_subset_analysis.py` to determine whether a subset's average demographics is significantly different from a total's, because calculating Euclidean distance is so much faster than trying to run a chi-squared test for every necessary subset comparison.

### Notes

*The equation currently used in `automated_subset_analysis.py` to predict significant Euclidean distance threshold using subset size was found using the equivalent of this function call:
```
python3 euclidean_threshold_estimator.py --subset_size 100 200 300 400 500 1000 1500 2000 --exclude 10 --n_analyses 3
```
The data and methods used to calculate that equation can be found in `euclidean_threshold_raw_data.txt` and `euclidean_threshold_estimation.ods`, which are both in the `automated_subset_analysis_files` subdirectory of this folder.  

## Meta

This `README` was created on 2019-10-03 by Greg Conan and last updated on 2019-10-08 by Greg Conan.