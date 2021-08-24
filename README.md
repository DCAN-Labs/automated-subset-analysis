# Automated Subset Analysis

These automated split-half subset reliability analysis scripts were written for the ABCD resource paper. They use Python's `argparse` package so that it can be run from the BASH command line and accept command-line arguments.

*The rest of this README focuses on `automated_subset_analysis.py`. If you are looking for a description of `asa_submitter.py`, then see `src/README.md`.*

## Purpose

Randomly selects and compares pairs of subsets of data from two groups. `automated_subset_analysis.py` makes and saves the subsets, the correlation values between them and the groups, and graph visualizations of those correlations.

### Brief Explanation of Steps

1. `automated_subset_analysis.py` accepts demographic data for 2 datasets. For each dataset, it randomly selects many subsets which are not significantly different from the other dataset's total demographics.
1. Once this script has made subsets of both groups, it finds correlations between the average of each subset and the other group. It also finds the correlation between the averages of both subsets. This process repeats a specified number of times for subsets which include specified numbers of subjects. After finding the correlations, the script saves them to `.csv` files.
1. The correlation values are then plotted as data points on graph visualizations. Those graphs are saved when `automated_subset_analysis` finishes executing.

For a more detailed explanation, see this document's `Explanation of Process` section.

## Installation

Installation should be simple:

1. Clone this repository to a location on your local filesystem.
1. Run `pip install -r requirements.txt` from within the new `automated_subset_analysis` directory.
1. To verify that the code is set up, run `python3 automated_subset_analysis.py --help` within the same directory.

## Requirements

### Dependencies

- The only dependency is [Python 3.6.8](https://www.python.org/downloads/release/python-368) or greater.

### Python Packages

- All of the Python packages required by this script which are not default can be found in this directory's `requirements.txt` file.

### Limitations

- This repository's Python code should theoretically be able to run on any operating system. However, so far it has only been tested on `*nix` systems. It may be unable to run on Windows or Mac. 

## Usage

### Required Arguments (2)

1. `group_1_demo_file` is a path to a `.csv` file which contains demographic data about all subjects in group 1. By default, the script will assume that the group 1 input demographics `.csv` file has a column of numerical data under each of these names:

    ```
    demo_comb_income_v2b, demo_ed_v2, demo_prnt_ed_v2b, demo_sex_v2b, ehi_y_ss_scoreb interview_age, medhx_9a, race_ethnicity, rel_relationship, site_id_l
    ```


    The last column of the demographics `.csv` file should be a list of paths (1 per line) to the `.nii` files of all subjects in group 1.

1. `group_2_demo_file` is a path to a `.csv` file just like `group_1_demo_file`, but containing demographic data about all subjects in group 2.

Example of a basic call to this script:

```sh
demo1=/home/user/conan/data/group1_pconn.csv
demo2=/home/user/conan/data/group2_pconn.csv
python3 automated_subset_analysis.py ${demo1} ${demo2}
```

### Optional Arguments (34)

#### File Paths with Default Values (7)

1. `--group-1-avg-file` takes one valid path to a readable `.nii` file containing the average matrix for the entire group 1. By default, this path will be to the `group1_10min_mean.pconn.nii` file in this script's parent folder or in the `--output` folder.

1. `--group-2-avg-file` also takes one valid path to a readable `.nii` file, just like `--group-1-avg-file`. By default, it will point to the `group2_10min_mean.pconn.nii` file in one of the same places.

1. `--group-1-var-file` also takes a valid `.nii` file path pointing to group 1's total variance matrix. By default, it will point to the `group1_variance_matrix.*.nii` file in one of the same places.

1. `--group-2-var-file` also takes a valid `.nii` file path pointing to group 2's total variance matrix. By default, it will point to the `group2_variance_matrix.*.nii` file in one of the same places.

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

    - Giving this flag multiple `.csv` files will put all of their correlations onto one visualization.

1. `--skip-subset-generation` takes either no parameters or one path to a readable directory as a parameter. Include this flag to calculate correlations and create the visualization using existing subsets instead of randomly generating new ones. By default, the subsets to use for calculating the correlations between average matrices and producing a visualization will be assumed to exist in the `--output` folder. To load subsets from a different folder, add the path to this flag as a parameter.

##### Notes for Step-Skipping Flags

- If `--skip-subset-generation` or `--only-make-graphs` is included, then `--subset-size` and `--n-analyses` will do nothing.
- If `--only-make-graphs` is included, then `--skip-subset-generation` will do nothing.
- Unless the `--only-make-graphs` flag is used, the `.csv` file(s) with subsets' average correlations will/must be called `correlations_sub1_sub2.csv`, `correlations_sub1_all2.csv`, and `correlations_sub2_all1.csv`.

#### Optional Plotly Visualization Elements (4)

1. `--fill` takes one parameter, a string that is either `all` or `confidence-interval`. Include this flag to choose which data to shade in the visualization. Choose `all` to shade in the area within the minimum and maximum correlations in the dataset. Choose `confidence-interval` to only shade in the 95% confidence interval of the data. By default, neither will be shaded. This argument cannot be used if `--only-make-graphs` has multiple parameters.

1. `--hide-legend` takes no parameters. Unless this flag is included, the output visualization(s) will display a legend in the top- or bottom-right corner showing the name of each thing plotted on the graph: data points, average trendline, confidence interval, and/or entire data range.

1. `--plot` takes one or more strings: `scatter` and/or `stdev`. By default, a visualization will be made with only the average value for each subset size. Include this flag with the parameter `scatter` to also plot all data points as a scatter plot, and/or with the parameter `stdev` to also plot standard deviation bars for each subset size.

1. `--rounded-scatter` takes no parameters. Include this flag to reduce the total number of data points plotted on any scatter-plot visualization by only including points at rounded intervals. This flag does nothing unless `--plot` includes `scatter`. 

#### Plotly Visualization Formatting Arguments (7)

1. `--axis-font-size` takes one positive integer, the font size of the text on both axes of the visualizations that this script will create. If this argument is excluded, then by default, the font size will be `30`.

1. `--graph-title` takes one string, the title at the top of all output visualizations and the name of the output `.html` visualization files. To break the title into two lines, include `<br>` in the `--graph-title` string. Unless this flag is included, each visualization will have one of these default titles:
    - "Correlations Between Average Subsets"
    - "Group 1 Subset to Group 2 Correlation"
    - "Group 1 to Group 2 Subset Correlation"
    - "Correlation Between Unknown Groups"

1. `--marker-size` takes one positive integer to determine the size (in pixels) of each data point in the output visualization. The default size is 5.

1. `--place-legend` takes one number between 0 and 1, the location of the legend on the y-axis in the output visualization. 0 is the very bottom of the visualization and 1 is the very top. By default, this value will be 0.05.

1. `--title-font-size` takes one positive integer. It is just like `--axis-font-size`, except for the title text in the visualizations. This flag determines the size of the title text above the graph as well as both axis labels. If this argument is excluded, then by default, the font size will be `40`.

1. `--trace-titles` takes one or more strings. Each will label one dataset in the output visualization. Each should be the title of one of the `.csv` files given to `--only-make-graphs`. Include exactly as many titles as there are `--only-make-graphs` parameters, in exactly the same order as those parameters, to match titles to datasets correctly. This argument only does anything when running the script in `--only-make-graphs` mode. 

1. `--y-range` takes two floating-point numbers, the minimum and maximum values to be displayed on the y-axis of the graph visualizations that this script will create. By default, this script will automatically set the y-axis boundaries to show all of the correlation values and nothing else.

#### MATLAB Visualization Arguments (6)

The following arguments only apply when making a visualization using the compiled MATLAB code instead of the Python Plotly code. So, they do nothing unless the `--plot-with-matlab` argument is included.

1. `--plot-with-matlab` takes one string, a valid path to an existing directory for the MATLAB Runtime Environment v9.4. Include this flag to create the output visualization using compiled MATLAB "MultiShadedBars" code (see `src`). Otherwise, none of the `matlab` flags will do anything and the subset analysis code will produce an output visualization using `plotly`.

1. `--matlab-lower-bound` takes one decimal number between 0 and 1, the lower bound of data to display on the MATLAB output visualization.

1. `--matlab-no-edge` takes no parameters. By default, the output visualization will display an edge. Include this flag to hide that edge.

1. `--matlab-show` takes no parameters. Include this flag to display the threshold as a line on the output visualization. Otherwise, the line will not be shown.

1. `--matlab-upper-bound` takes one decimal number between 0 and 1, the upper bound of data to display on the MATLAB output visualization.

1. `--matlab-rgba` takes 3 to 5 3 to 5 numbers between 0 and 1, the RGBA values and line threshold for producing the visualization. Respectively those numbers are the red value, green value, blue value, (optional) alpha opacity value, and (optional) threshold to include a line at on the visualization.

#### Other Flags (5)

1. `--columns` takes one or more strings. Each should be the name of a column in the demographics `.csv` which contains numerical data to include in the subset correlations analysis. By default, the script will assume that both input demographics `.csv` files have columns of numerical data with these names:

    ```
    demo_comb_income_v2b, demo_ed_v2, demo_prnt_ed_v2b, demo_sex_v2b, ehi_y_ss_scoreb interview_age, medhx_9a, race_ethnicity, rel_relationship, site_id_l
    ```

1. `--calculate` takes one string to define the output metric. With its default value of `mean`, the subset analysis will calculate correlations between subsets'/groups' average values. Use `--calculate variance` to correlate the subsets' variances instead, or `--calculate effect-size` to measure the effect size of the difference between each subset and the total group. 

1. `--inverse-fisher-z` takes no parameters. Include this flag to do an inverse Fisher-Z transformation on the matrices imported from the `.pconn` files of the data before getting correlations.

1. `--no-matching` takes no parameters. Include this flag to match subsets on every demographic variable *except* family relationships. Otherwise, subsets will be matched on all demographic variables.

    - By default, `automated_subset_analysis.py`  checks that every subset of one group has the same proportion of twins, triplets, or other siblings as the other group. It also checks that no one in the subset has family members outside the subset. `--no-matching` will skip both checks. 

    - Use this flag if `--subset-size` includes any number under 25, because family relationship matching takes a *very* long time for small subsets.

1. `--parallel` takes one valid path, the directory containing `automated_subset_analysis.py`. It will automatically be included by `asa_submitter.py` to simultaneously run multiple different instances of `automated_subset_analysis.py` as a batch command. Otherwise, this flag is not needed. Do not use this flag, because it will be included automatically if needed.

For more information, including the shorthand flags for each option, run this script with the `--help` command: `python3 automated_subset_analysis.py --help`

### Advanced Usage Examples

Generate 50 subsets and save a `.csv` file of each, including 10 subsets each of sizes 50, 100, 300, 500, and 1000:

```sh
python3 automated_subset_analysis.py ${demo1} ${demo2} --subset-size 50 100 300 500 1000 --n-analyses 10
```

Calculate the correlations between average matrices of already-generated subsets in the `./subsets/` folder, then save the correlations and a visualization of them to the `./correls/` folder:

```sh
python3 automated_subset_analysis.py ${demo1} ${demo2} --skip-subset-generation ./subsets/ --output ./correls/
```

## Explanation of Process

Two `.csv` files, each with demographic data about subjects from a group, are given by the user. One subset is randomly selected from each group repeatedly. The amount of subjects in each subset depends on `--subset-size`, and the number of times that amount is selected depends on `--n-analyses`.

Once a pair has been selected, the script calculates the Euclidean distance between the average demographics of each subset and the average demographics of the whole other group. If the Euclidean distance between the subset and the total is higher than the estimated maximum value for significance (given in the equation calculated by `./src/euclidean_threshold_estimator.py`<sup> 1 </sup>), then another subset is randomly generated and tested for significance. Otherwise, the subset pair is deemed valid and saved to a `.csv` file. The `.csv` has one subset per column and one subject ID per row, excluding the header row which only contains the group number of each subset.

After finding a valid pair of subsets, the script calculates the correlation between the subset of group 1 and the subset of group 2. This correlation value is stored with the number of subjects in both subsets described by the correlation. Once the correlation values are all calculated, each correlation value is saved out to a `.csv` file with a name starting with `correlations`. That `.csv` has two columns. It has one row per subset pair, excluding the header row which contains only the names of the columns: `Subjects` and `Correlation`.

Once the correlation `.csv` files are made, the script will make a graph visualization for each one. That graph will plot how the number of subjects in a subset pair (x-axis) relates to the correlation between the subsets in that pair (y-axis).<sup> 2 </sup>

### Notes

<sup>1</sup> The equation currently used in `automated_subset_analysis.py` to predict significant Euclidean distance threshold using subset size was found using this Bash code:

```sh
python3 ./src/euclidean_threshold_estimator.py ./raw/ABCD_2.0_group1_data_10minpconns.csv ./raw/ABCD_2.0_group2_data_10minpconns.csv -con-vars ./automated_subset_analysis_files/continuous_variables.csv --subset-size 1300 1200 1100 1000 900 800 700 600 500 400 300 200 100 90 80 70 60 50 --n-analyses 10
```

The data used to calculate that equation can be found in `./src/euclidean_threshold_estimate_data/est-eu-thresh-2019-12-12`.

<sup>2</sup> If `--plot-with-matlab` is not used, the output visualization will include:

1. One trendline using the average correlation values of each subset size (or more if `--only-make-graphs` includes multiple parameters),
1. A shaded region showing a data range, either the confidence interval or all data (if `--fill` is used),
1. Each correlation value as 1 data point (if `--plot` includes `scatter`; if `--rounded-scatter` is used, only correlation values at specific intervals will be plotted),
1. Standard deviation bars above and below each data point (if `--plot` includes `stdev`), and
1. A legend to identify all of these parts (unless `--hide-legend` is used).

## Metadata

Information about this `README` file:

- Created by Greg Conan, 2019-10-03
- Updated by Greg Conan, 2021-08-24
 
