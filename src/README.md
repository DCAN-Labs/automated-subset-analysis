# Automated Subset Analysis: Other Scripts

These are a collection of scripts which supplement `automated_subset_analysis.py`, using much of the same functionality:
- `euclidean_threshold_estimator.py` calculates the equation used by `automated_subset_analysis.py` to validate random subsets, using the latter script's `--euclidean` argument.
- `make_average_matrix.py` takes groups of matrices and makes one average matrix for each group. Those average matrices can then be used for the `--group-1-avg-file` and `--group-2-avg-file` arguments of `automated_subset_analysis.py`.
- `make_average_and_subsets.py` runs `make_average_matrix.py` to get average matrices, and then immediately afterwards runs `automated_subset_analysis.py` using those average matrices.
- `asa_submitter.py` runs many instances of `automated_subset_analysis.py` in parallel to speed up its data processing.
- `conan_tools.py` does nothing, but it contains all of the functionality used by multiple scripts in this repository.

## Purpose

- ### `euclidean_threshold_estimator.py`

Randomly generates subsets in the same way as `automated_subset_analysis.py`, but instead of analyzing them and producing visualizations, it estimates the maximum Euclidean distance required to create a subset which is not significantly different from a total set. It then prints out a list of these estimated Euclidean distances, including a logarithmic regression equation to estimate the Euclidean distance threshold for a specific number of subjects.

- ### `make_average_matrix.py`

Accepts paths to `.nii` data files, then creates new `.nii` files with the average values of the originals. Each group's average matrix is saved into its own new `.nii` file. The output matrices will have the same dimensions as the original `.nii` files, but each point in the average matrix will have the average of all of the values at the same point in the input matrices from a group. Using this script on every matrix in group 1 provides the `--group-1-avg-file` file to use in `automated_subset_analysis.py`. The same applies to group 2 and `--group-2-avg-file`.

- ### `make_average_and_subsets.py`

This script is simply a wrapper to run `make_average_matrix.py` on two groups, and then run `automated_subset_analysis.py` using the new average matrices. Immediately after `make_average_matrix.py` finishes running, `automated_subset_analysis.py` will correlate the newly created average matrices with subsets of the two groups, and then produce exactly the same outputs that you would expect from running `automated_subset_analysis.py` normally.

- ### `asa_submitter.py`

This script is another wrapper, but for a different reason: Instead of running two scripts in a row, this runs many instances of `automated_subset_analysis.py` at once for improved speed. It was written to run a batch job using the [slurm job scheduler](https://slurm.schedmd.com) on the Exacloud server, so it may not work on using other servers or job schedulers without minor changes to the code.

- ### `conan_tools.py`

Since the scripts above share so much functionality with `automated_subset_analysis.py`, all of the scripts in this repo import their shared functionality from the common source script `conan_tools.py`. It has no `main` function, so is not meant to be run on its own.

## Usage

Just like `automated_subset_analysis.py`, all of the scripts in this directory (except `conan_tools.py`) are run by calling them from the command line using the `python3` command.

### Examples

- #### `euclidean_threshold_estimator.py`

```
python3 ./euclidean_threshold_estimator.py ../raw/group_1_demographics.csv ../group_2_demographics.csv --n-analyses 4 --subset-size 50 75 100 --inverse-fisher-z
```

- #### `make_average_matrix.py`

```
python3 ./make_average_matrix.py --matrices-conc-1 ../raw/group_1_matrix_paths.conc --matrices-conc-2 ../raw/group_2_matrix_paths.conc --output ../data/
```

- #### `make_average_and_subsets.py`

```
python3 ./make_average_and_subsets.py ../raw/group_1_demographics.csv ../raw/group_2_demographics.csv --matrices-conc-1 ../raw/group_1_matrix_paths.conc --matrices-conc-2 ../raw/group_2_matrix_paths.conc --group-1-avg-file ../raw/gp1_avg_matrix.dscalar.nii --group-2-avg-file ../raw/gp2_avg_matrix.dscalar.nii --output ../data/ --n-analyses 10 --subset-size 100
```

- ### `asa_submitter.py`

```
python3 ./asa_submitter.py ../raw/group_1_demographics.csv ../raw/group_2_demographics.csv --sbatch-string "sbatch --time=1:00:00 --cpus 1 -A my_group ./automated_subset_analysis.py" --output ./data/ --n-analyses 5 --subset-size 50 100 200 
```


## Command-Line Arguments

Most of the command-line arguments used by the scripts in this directory are the same as those used by `automated_subset_analysis.py`, and work in the same way. However, each script accepts a different combination of those arguments, shown in the chart below. There are also two flags accepted by scripts in this directory but not by `automated_subset_analysis.py`, which are described below the combinations chart.

### Chart: Which Scripts Accept Which Arguments

| Argument Name | automated... | euclidean... | make...matrix | make...subsets | asa... |
|----------------------------|:-:|:-:|:-:|:-:|:-:|
| `group_1_demo_file`        | *Required* | *Required* | | *Required* | *Required* |
| `group_2_demo_file`        | *Required* | *Required* | | *Required* | *Required* |
| `--columns`                | Optional | Optional | | Optional | Optional |
| `--n-analyses`             | Optional | Optional | | Optional | Optional |
| `--nan-threshold`          | Optional | Optional | | Optional | Optional |
| `--subset-size`            | Optional | Optional | | Optional | Optional |
| `--continuous-variables`   | | Optional | | | |
| `--example-file`           | | | Optional | Optional | | 
| `--inverse-fisher-z`       | Optional | | Optional | Optional | Optional
| `--matrices-conc-1`        | Optional | | *Required* | *Required* | Optional |
| `--matrices-conc-2`        | Optional | | *Required* | *Required* | Optional |
| `--output`                 | Optional | | Optional | Optional | Optional |
| `--group-1-avg-file`       | Optional | | Optional | Optional | Optional |
| `--group-2-avg-file`       | Optional | | Optional | Optional | Optional |
| `--axis-font-size`         | Optional | | | Optional | Optional |
| `--dimensions`             | Optional | | | Optional | Optional |
| `--euclidean`              | Optional | | | Optional | Optional |
| `--fill`                   | Optional | | | Optional | Optional |
| `--title-font-size`        | Optional | | | Optional | Optional |
| `--y-range`                | Optional | | | Optional | Optional |
| `--only-make-graphs`       | Optional | | | | Optional |
| `--skip-subset-generation` | Optional | | | | Optional |
| `--parallel`               | Optional<sup>1</sup> | | | | |
| `--sbatch-string`          | | | | | *Required* |

<sup>1</sup>The `--parallel` argument does not need to be used for `asa_submitter.py`, because that script will add it automatically when running `automated_subset_analysis.py`.

### Arguments for Other Scripts

The flags listed above which are used by the main `automated_subset_analysis.py` script take exactly the same values for it as for the other scripts in this repo. However, 3 of those flags are not used by the main script:

#### `euclidean_threshold_estimator.py`

- `--continuous-variables` takes one or more strings which name columns in the demographics `.csv` files with continuous instead of categorical data. The data from any variable from this list will be tested for statistical significance using a T-test. Any other data will be tested for significance using a chi-squared test. By default, if this flag is excluded, these four demographic variables will be considered continuous: `demo_prnt_ed_v2b`, `interview_age`, `rel_group_id`, and `rel_relationship`.

#### `make_average_matrix.py`

- `--example-file` takes one valid path to a matrix file. It should be the same type of matrix file as the ones listed in `--matrices-conc-1` and `--matrices-conc-2`, so that when the entire group's average matrix is created, the example file can be used as a template to save out the average matrix to a `*.nii` file. By default, if this flag is excluded, the first path listed in `--matrices-conc-1` will be used as the template.

#### `asa_submitter.py`

- `--sbatch-string` is a string containing all of the required parameters to run multiple instances of `automated_subset_analysis.py` in parallel as an SBATCH command. All other input arguments will be appended to this string, and then the result will be executed. This flag must be included with a command string in order to run `asa_submitter.py`.

## Explanation of Process

- ### `euclidean_threshold_estimator.py`

Two `.csv` files, each with demographic data about subjects from a group, are given by the user. A subset is randomly selected from each group repeatedly; the amount of subjects in each group depends on `--subset-size` and the number of times that amount is selected depends on `--n-analyses`.

Once a pair has been selected, a statistical (chi-squared- or T-) test of significance is run comparing the average demographics of each subset to the average demographics of the whole other group. If the statistical test shows a significant difference, then another subset is randomly generated and tested for significance. Otherwise, the subset pair is deemed valid.

After finding a valid pair of subsets, the script calculates the Euclidean distance between each subset's average demographics and the whole other group's average demographics. It collects these Euclidean distances for however many subset sizes and number of times that the user specified.

Once all of the Euclidean distances are calculated, the script takes the maximum Euclidean distance for each subset size and calculates a logarithmic regression equation to predict the maximum Euclidean distance given only a subset size. This equation is used by `automated_subset_analysis.py` to determine whether a subset's average demographics is significantly different from a total's, because calculating Euclidean distance is so much faster than trying to run a chi-squared test or t-test for every necessary subset comparison.

- ### `make_average_matrix.py`

Two `.conc` files, one for each group, are given with lists of paths to matrix files. The script will take these and a path to an `--output` directory. For each `.conc` file, the script will save a matrix averaging over every matrix in the `.conc` file into the `--output` directory.

- ### `make_average_and_subsets.py`

Normally, `automated_subset_analysis.py` accepts already-existing average matrices of both groups as its `--group-1-avg-file` and `--group-2-avg-file` arguments. However, `make_average_and_subsets.py` runs `make_average_matrix.py` to create new ones and then passes them to `automated_subset_analysis.py`.

- ### `asa_submitter.py`

The `--sbatch-string` is executed together with all of the other input parameters to run a batch job of many (specifically, `--n-analyses` * the number of parameters in `--subset-size`) instances of `automated_subset_analysis.py` at once. Each instance will save its output subsets into a subdirectory of the `--output` directory. The number of subdirectories created will be equl to `--n-analyses`, numbered from 1 to `--n-analyses` like so:
  - `./data/output1`
  - `./data/output2`
  - ... 
  - `./data/output20` (if `--n-analyses` is 20)

Each instance will append its output correlations to 3 `.csv` files in the `--output` directory: 1 for correlating the subsets to each other, and 1 each correlating a group to the subsets of the other group. If there are already correlation `.csv` files in the `--output` directory from previous runs, the output correlation values of any subsequent runs will be appended to the already-existing files. Otherwise, new `.csv` files will be created. 

## Metadata

Information about this `README` file:

- Created by Greg Conan, 2019-12-16
- Updated by Greg Conan, 2020-01-08
