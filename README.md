# Automated Subset Analysis

This script was written for the ABCD resource paper. It uses Python's `argparse` package so that it can be run from the BASH command line and accept command-line arguments.

## Purpose

Randomly selects and compares pairs of subsets of data. It outputs the subsets, the correlation values between them, and a graph visualization of those correlations.

`automated_subset_analysis.py` accepts demographic data for 2 datasets. For each dataset, it randomly selects a subset which is not significantly different from the other dataset's total demographics. Once this script has the  find the correlation between the average matrices of each subset, and plot all of those correlations for each subset size

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

1. `group_1_demo_file` is a path to a text file which contains only a list of paths (1 per line) to the `.pconn` files of all subjects in the first group to analyze a subset of.

1. `group_2_demo_file` is a path to a text file which contains only a list of paths (1 per line) to the `.pconn` files of all subjects in the second group to analyze a subset of.

Example of a basic call to this script:
```
python3 automated_subset_analysis.py /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp1_10min_pconn.csv /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp2_10min_pconn.csv
```

### Optional Arguments (18)

#### File Paths with Default Values

1. `--group-1-avg-file` takes one valid path to a readable `.nii` file containing the average matrix for the entire group 1. By default, this path will be to the `group1_10min_mean.pconn.nii` file in this script's parent folder or in the `--output` folder.

1. `--group-2-avg-file` also takes one valid path to a readable `.nii` file, just like `--group-1-avg-file`. By default, it will point to the `group1_10min_mean.pconn.nii` file in one of the same places.

1. `--matrices-conc-1` takes one path to a readable `.conc` file containing only a list of valid paths to group 1 matrix files. This flag is only needed if your group 1 demographics `.csv` file either does not have a column labeled `'pconn10min'` with paths to matrix files, or if it does include that column but you want to use different paths.

1. `--matrices-conc-2` also takes one path to a readable `.conc` file. It is just like `--matrices-conc-1`, but for group 2.

1. `--output` takes one file path to a directory where all files produced by this script will be saved. If the directory already exists, then this script will add files to it and only overwrite files with conflicting filenames. If not, then this script will create a directory at the `--output` path. If this flag is excluded, then the script will save files to a new subdirectory of the present working directory, `./data/`.

#### Numerical Values for Data Processing

1. `--dimensions` takes an integer, either 1 or 2, representing how many dimensions are in the input matrices. `.*conn.nii` files are 2-dimensional, whereas `.*scalar.nii` files are 1-dimensional. This script can process either kind of file, but by default will assume that a given file is 2-dimensional.

1. `--n-analyses` takes one positive integer, the number of times to generate a pair of subsets and analyze them. For every integer given in the `--subset-size` list, this script will randomly generate `--n-analyses` subsets, creating `--subset-size * --n-analyses` total `.csv` files.

1. `--nan-threshold` takes one floating-point number between 0 and 1. If the percentage of rows with Not-a-Number (NaN) values in the data for either group's demographic file is greater than the `--nan-threshold`, then the script will raise an error. Otherwise, the script will drop every row containing a NaN. The default NaN threshold is `0.1`, meaning that if over 10% of rows have a NaN value, the script will crash.

1. `--subset-size` takes one or more positive integers, the number of subjects to include in subsets. Include a list of whole numbers to generate subsets pairs of different sizes. By default, the subset sizes will be `[50, 100, 200, 300]`. An example of entering a different list of sizes is `--subset-size 100 300 500 1000`.

#### Flags to Skip Steps of Process

1. `--only-make-graphs` takes one or more paths to readable `.csv` files as a parameter. Include this flag to import average correlations data from `.csv` files instead of making any new ones. Given this flag, `automated_subset_analysis.py` only makes graph visualizations of already-existing data. If this flag is included, it must include paths to readable `.csv` files with 2 columns: `Subjects` (the number of subjects in each subset) and `Correlation` (the correlation between each randomly generated subset in that pair). The `.csv` files must be called `correlations_sub1_sub2.csv`, `correlations_sub1_all2.csv`, and/or `correlations_sub2_all1.csv`.

1. `--skip-subset-generation` takes either no parameters or one path to a readable directory as a parameter. Include this flag to calculate correlations and create the visualization using existing subsets instead of randomly generating new ones. By default, the subsets to use for calculating the correlations between average matrices and producing a visualization will be assumed to exist in the `--output` folder. To load subsets from a different folder, add the path to this flag as a parameter.

If `--skip-subset-generation` or `--only-make-graphs` is included, then `--subset-size` and `--n-analyses` will do nothing. If `--only-make-graphs` is included, then `--skip-subset-generation` will do nothing.

#### Numerical Values for Visualization Format

1. `--axis-font-size` takes one positive integer, the font size of the text on both axes of the visualizations that this script will create. If this argument is excluded, then by default, the font size will be `30`.

1. `--title-font-size` takes one positive integer. It is just like `--axis-font-size`, except for the title text in the visualizations. Title text includes the title at the top as well as both axis labels. If this argument is excluded, then by default, the font size will be `40`.

1. `--y-range` takes two floating-point numbers, the minimum and maximum values to be displayed on the y-axis of the graph visualizations that this script will create. By default, this script will automatically set the y-axis boundaries to show all of the correlation values and nothing else.

#### Other Flags

1. `--columns` takes one or more strings. Each should be the name of a column in the demographics `.csv` which contains numerical data to include in the subset correlations analysis. By default, the script will assume that both input demographics `.csv` files have columns of numerical data with these names:
    ```
    demo_comb_income_v2b, demo_ed_v2, demo_prnt_ed_v2b, demo_sex_v2b, ehi_y_ss_scoreb interview_age, medhx_9a, race_ethnicity, rel_relationship, site_id_l
    ```

1. `--fill` takes one parameter, a string that is either `all` or `confidence-interval`. Include this flag to choose which data to shade in the visualization. Choose `all` to shade in the area within the minimum and maximum correlations in the dataset. Choose `confidence-interval` to only shade in the 95% confidence interval of the data. By default, this argument will be `confidence-interval`.

1. `--inverse-fisher-z` takes no parameters. Include this flag to do an inverse Fisher-Z transformation on the matrices imported from the `.pconn` files of the data before getting correlations.

1. `--parallel` takes one valid path, the directory containing `automated_subset_analysis.py`. It should be included to simultaneously run multiple different instances of `automated_subset_analysis.py` as a batch command executed by `asa_submitter.py`, but otherwise it is not needed.

For more information, including the shorthand flags for each option, run this script with the `--help` command: `python3 automated_subset_analysis.py --help`

### Advanced Usage Examples

Generate 50 subsets and make a `.csv` file of each, including 10 subsets each of sizes 50, 100, 300, 500, and 1000:

```
python3 automated_subset_analysis.py /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp1_10min_pconn.csv /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp2_10min_pconn.csv --subset-size 50 100 300 500 1000 --n-analyses 10
```

Calculate the correlations between average matrices of already-generated subsets in the `./subsets/` folder, then save the correlations and a visualization of them to the `./correls/` folder:

```
python3 automated_subset_analysis.py /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp1_10min_pconn.csv /mnt/rose/shared/projects/ABCD/avg_pconn_maker/gp2_10min_pconn.csv --skip-subset-generation ./subsets/ --output ./correls/
```

## Explanation of Process

Two `.csv` files, each with demographic data about subjects from a group, are given by the user. One subset is randomly selected from each group repeatedly. The amount of subjects in each subset depends on `--subset-size`, and the number of times that amount is selected depends on `--n-analyses`.

Once a pair has been selected, the script calculates the Euclidean distance between the average demographics of each subset and the average demographics of the whole other group. If the Euclidean distance between the subset and the total is higher than the estimated maximum value for significance (given in the equation calculated by `./src/euclidean_threshold_estimator.py`<sup> 1 </sup>), then another subset is randomly generated and tested for significance. Otherwise, the subset pair is deemed valid and saved to a `.csv` file. The `.csv` has one subset per column and one subject ID per row, excluding the header row which contains only the group number of each subset.

After finding a valid pair of subsets, the script calculates the correlation between the subset of group 1 and the subset of group 2. This correlation value is stored with the number of subjects in both subsets described by the correlation. Once the correlation values are all calculated, they are saved out to a file called `correlations.csv`. That `.csv` has two columns, and one row per subset pair, excluding the header row which contains only the names of the columns: `Subjects` and `Correlation`.

Once `correlations.csv` is made, the script will then make a visualization plotting how the number of subjects in a subset pair (x-axis) relates to the correlation between the subsets in that pair (y-axis).

### Note

<sup>1</sup> The equation currently used in `automated_subset_analysis.py` to predict significant Euclidean distance threshold using subset size was found using this function call:
```
python3 ./src/euclidean_threshold_estimator.py ./raw/ABCD_2.0_group1_data_10minpconns.csv ./raw/ABCD_2.0_group2_data_10minpconns.csv -con-vars ./automated_subset_analysis_files/continuous_variables.csv --subset-size 1300 1200 1100 1000 900 800 700 600 500 400 300 200 100 90 80 70 60 50 --n-analyses 10
```
The data used to calculate that equation can be found in `./src/euclidean_threshold_estimate_data/est-eu-thresh-2019-12-12`.

## Metadata

Information about this `README` file:

- Created by Greg Conan, 2019-10-03
- Updated by Greg Conan, 2020-01-08
