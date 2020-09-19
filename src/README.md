# Automated Subset Analysis: Other Scripts

These are a collection of scripts which supplement `automated_subset_analysis.py`, using much of the same functionality:

- `euclidean_threshold_estimator.py` calculates the equation used by `automated_subset_analysis.py` to validate random subsets, using the latter script's `--euclidean` argument.
- `make_average_matrix.py` takes groups of matrices and makes one average matrix for each group. Those average matrices can then be used for the `--group-1-avg-file` and `--group-2-avg-file` arguments of `automated_subset_analysis.py`. It can also make variance matrices with the `--calculate variance` argument.
- `make_average_and_subsets.py` runs `make_average_matrix.py` to get average matrices, and then immediately afterwards runs `automated_subset_analysis.py` using those average matrices.
- `asa_submitter.py` runs many instances of `automated_subset_analysis.py` in parallel to speed up its data processing.
- `pairwise_correlator.py` compares 2 groups of `.nii` files, where both groups have the exact same subjects, and outputs the within-subject correlations broken down by each demographic variable in a bar plot and in a box plot.
- `conan_tools.py` does nothing, but it contains all of the functionality used by multiple scripts in this repository.

## Purpose

- ### `euclidean_threshold_estimator.py`

Randomly generates subsets in the same way as `automated_subset_analysis.py`, but instead of analyzing them and producing visualizations, it estimates the maximum Euclidean distance required to create a subset which is not significantly different from a total set. It then prints out a list of these estimated Euclidean distances, including a logarithmic regression equation to estimate the Euclidean distance threshold for a specific number of subjects.

- ### `make_average_matrix.py`

Accepts paths to `.nii` data files, then creates new `.nii` files with the average values of the originals. Each group's average matrix is saved into its own new `.nii` file. The output matrices will have the same dimensions as the original `.nii` files, but each point in the average matrix will have the average of all of the values at the same point in the input matrices from a group. Using this script on every matrix in group 1 provides the `--group-1-avg-file` file to use in `automated_subset_analysis.py`. The same applies to group 2 and `--group-2-avg-file`.

- ### `make_average_and_subsets.py`

This script is simply a wrapper to run `make_average_matrix.py` on two groups, and then run `automated_subset_analysis.py` using the new average matrices. Immediately after `make_average_matrix.py` finishes running, `automated_subset_analysis.py` will correlate the newly created average matrices with subsets of the two groups, and then produce exactly the same outputs that you would expect from running `automated_subset_analysis.py` normally.

- ### `asa_submitter.py`

*NOTE: `asa_submitter.py` must be moved into the parent directory before being run, since it uses functions from `automated_subset_analysis.py`*

This script is another wrapper, but for a different reason: Instead of running two scripts in a row, this runs many instances of `automated_subset_analysis.py` at once for improved speed. It was written to run batch jobs on OHSU's Exacloud server. So, it may not work on using other servers or job schedulers without changing the code.

- ### `conan_tools.py`

Since the scripts above share so much functionality with `automated_subset_analysis.py`, all of the scripts in this repo import their shared functionality from the common source script `conan_tools.py`. It has no `main` function, so is not meant to be run on its own.

## Usage

Just like `automated_subset_analysis.py`, all of the scripts in this directory (except `conan_tools.py`) are run by calling them from the command line using the `python3` command.

### Examples

- #### `euclidean_threshold_estimator.py` examples

```sh
python3 ./euclidean_threshold_estimator.py ../raw/group_1_demographics.csv ../group_2_demographics.csv --n-analyses 4 --subset-size 50 75 100 --inverse-fisher-z
```

- #### `make_average_matrix.py` examples

```sh
python3 ./make_average_matrix.py --matrices-conc-1 ../raw/group_1_matrix_paths.conc --matrices-conc-2 ../raw/group_2_matrix_paths.conc --output ../data/
```

- #### `make_average_and_subsets.py` examples

```sh
python3 ./make_average_and_subsets.py ../raw/group_1_demographics.csv ../raw/group_2_demographics.csv --matrices-conc-1 ../raw/group_1_matrix_paths.conc --matrices-conc-2 ../raw/group_2_matrix_paths.conc --group-1-avg-file ../raw/gp1_avg_matrix.dscalar.nii --group-2-avg-file ../raw/gp2_avg_matrix.dscalar.nii --output ../data/ --n-analyses 10 --subset-size 100
```

- ### `asa_submitter.py` examples

```sh
python3 ./asa_submitter.py ../raw/group_1_demographics.csv ../raw/group_2_demographics.csv --sbatch-string "sbatch --time=1:00:00 --cpus 1 -A my_group ./automated_subset_analysis.py" --output ./data/ --n-analyses 5 --subset-size 50 100 200 
```

## Command-Line Arguments

Most of the command-line arguments used by the scripts in this directory are the same as those used by `automated_subset_analysis.py`, and work in the same way. However, each script accepts a different combination of those arguments, shown in the chart below. There are also two flags accepted by scripts in this directory but not by `automated_subset_analysis.py`, which are described below the combinations chart.

### Chart: Which Scripts Accept Which Arguments

        "nan_threshold", "only_make_graphs", "output", "spearman_rho",


| Argument Name | automated... | euclidean... | make...matrix | make...subsets | asa... | pairwise...
|----------------------------|:-:|:-:|:-:|:-:|:-:|:-:|
| `group_1_demo_file`        | **Required** | **Required** | | **Required** | **Required** | **Required** |
| `group_2_demo_file`        | **Required** | **Required** | | **Required** | **Required** | |
| `--columns`                | Optional | Optional | | Optional | Optional |
| `--n-analyses`             | Optional | Optional | | Optional | Optional |
| `--nan-threshold`          | Optional | Optional | | Optional | Optional | Optional |
| `--no-matching`            | Optional | Optional | | Optional | Optional |
| `--subset-size`            | Optional | Optional | | Optional | Optional |
| `--continuous-variables`   | | Optional | | | | |
| `--example-file`           | | | Optional | Optional | |
| `--inverse-fisher-z`       | Optional | | Optional | Optional | Optional
| `--matrices-conc-1`        | Optional | | **Required** | **Required** | Optional | **Required** |
| `--matrices-conc-2`        | Optional | | **Required** | **Required** | Optional | **Required** |
| `--calculate`              | Optional | | Optional | Optional | Optional |
| `--group-1-avg-file`       | Optional | | Optional | Optional | Optional |
| `--group-2-avg-file`       | Optional | | Optional | Optional | Optional |
| `--group-1-var-file`       | Optional | | Optional | Optional | Optional |
| `--group-2-var-file`       | Optional | | Optional | Optional | Optional |
| `--output`                 | Optional | | Optional | Optional | Optional | Optional |
| `--spearman-rho`           | Optional | | | Optional | Optional | Optional |
| `--axis-font-size`         | Optional | | | Optional | | Optional |
| `--euclidean`              | Optional | | | Optional | | |
| `--fill`                   | Optional | | | Optional | | |
| `--graph-title`            | Optional | | | Optional | | Optional |
| `--hide-legend`            | Optional | | | Optional | | Optional |
| `--marker-size`            | Optional | | | Optional | | Optional |
| `--plot`                   | Optional | | | Optional | | |
| `--rounded-scatter`        | Optional | | | Optional | | |
| `--title-font-size`        | Optional | | | Optional | | Optional |
| `--trace-titles`           | Optional | | | Optional | | |
| `--y-range`                | Optional | | | Optional | | Optional |
| `--parallel`               | Optional<sup>1</sup> | | | | |
| `--only-make-graphs`       | Optional | | | | Optional | Optional | Optional |
| `--skip-subset-generation` | Optional | | | | Optional | 
| `--job-time-limit`         | | | | | Optional |
| `--print-command`          | | | | | Optional |
| `--queue-max-size`         | | | | | Optional |
| `--seconds-between-jobs`   | | | | | Optional |

<sup> 1 </sup>The `--parallel` flag probably never needs to be added by the user. `asa_submitter.py` will add it automatically when running `automated_subset_analysis.py` batch jobs, and that is the flag's only intended use.

### Arguments for Other Scripts

The flags listed above which are used by the main `automated_subset_analysis.py` script take exactly the same values for it as for the other scripts in this repo. However, several of those flags are not used by the main script:

#### `euclidean_threshold_estimator.py` arguments

- `--continuous-variables` takes one or more strings which name columns in the demographics `.csv` files with continuous instead of categorical data. The data from any variable from this list will be tested for statistical significance using a T-test. Any other data will be tested for significance using a chi-squared test. By default, if this flag is excluded, these four demographic variables will be considered continuous: `demo_prnt_ed_v2b`, `interview_age`, `rel_group_id`, and `rel_relationship`.

#### `make_average_matrix.py` arguments

- `--example-file` takes one valid path to a matrix file. It should be the same type of matrix file as the ones listed in `--matrices-conc-1` and `--matrices-conc-2`, so that when the entire group's average matrix is created, the example file can be used as a template to save out the average matrix to a `*.nii` file. By default, if this flag is excluded, the first path listed in `--matrices-conc-1` will be used as the template.

#### `asa_submitter.py` arguments

- `--job-time-limit` takes a time string formatted as `HH:MM:SS`, e.g. `07:31:52` for 7 hours, 31 minutes, and 52 seconds. The time string will determine how long to run each `automated_subset_analysis.py` batch job submitted by `asa_submitter.py`. The default time limit is `04:00:00`.

- `--print-command` takes no arguments. Include this flag to print every Bash command that is run to submit an `automated_subset_analysis.py` batch job.

- `--queue-max-size` takes one positive integer, the maximum number of  `automated_subset_analysis.py` jobs to run in parallel as an SBATCH command. By default, a maximum of 100 jobs will run at once.

- `--seconds-between-jobs` takes one positive integer, the number of seconds to wait after submitting one `automated_subset_analysis.py` batch job before submitting the next. By default, this script will wait 60 seconds between job submissions. However, if the number of jobs running is `--queue-max-size` or greater, then this script will not submit another job until queue space clears up.

## Explanation of Process

- ### `euclidean_threshold_estimator.py` explanation

Two `.csv` files, each with demographic data about subjects from a group, are given by the user. A subset is randomly selected from each group repeatedly; the amount of subjects in each group depends on `--subset-size` and the number of times that amount is selected depends on `--n-analyses`.

Once a pair has been selected, a statistical (chi-squared- or T-) test of significance is run comparing the average demographics of each subset to the average demographics of the whole other group. If the statistical test shows a significant difference, then another subset is randomly generated and tested for significance. Otherwise, the subset pair is deemed valid.

After finding a valid pair of subsets, the script calculates the Euclidean distance between each subset's average demographics and the whole other group's average demographics. It collects these Euclidean distances for however many subset sizes and number of times that the user specified.

Once all of the Euclidean distances are calculated, the script takes the maximum Euclidean distance for each subset size and calculates a logarithmic regression equation to predict the maximum Euclidean distance given only a subset size. This equation is used by `automated_subset_analysis.py` to determine whether a subset's average demographics is significantly different from a total's, because calculating Euclidean distance is so much faster than trying to run a chi-squared test or t-test for every necessary subset comparison.

- ### `make_average_matrix.py` explanation

Two `.conc` files, one for each group, are given with lists of paths to matrix files. The script will take these and a path to an `--output` directory. For each `.conc` file, the script will save a new matrix into the `--output` directory with either the mean or the variance of every cell position in the matrices. If `--calculate` is `mean`, then the script will average out every matrix in the `.conc` file. If `--calculate` is `variance`, then the script will save the variance between all of the matrices in one `.conc` file into one new matrix.

- ### `make_average_and_subsets.py` explanation

Normally, `automated_subset_analysis.py` accepts already-existing average matrices of both groups as its `--group-1-avg-file` and `--group-2-avg-file` arguments. However, `make_average_and_subsets.py` runs `make_average_matrix.py` to create new ones and then passes them to `automated_subset_analysis.py`.

- ### `pairwise_correlator.py` explanation

`pairwise_correlator.py` accepts 2 `.conc` files in its `--matrices-conc-*` arguments. Other scripts use those arguments for comparing `.nii` files from 2 lists of totally different subjects. However, `pairwise_correlator.py` accepts 2 `.conc` files which have lists of different `.nii` files for the *same* subjects. Each line in each `.conc` file must have a path to a `.nii` file for the same subject as that same line in the other `.conc` file. 

For the same reason, `pairwise_correlator.py` only uses the `group_1_demo_file` argument and not `group_2_demo_file`. Because each subject only belongs to one group, and both `.conc` files include exactly the same subjects, this script can only be run on one group at a time. Ff the subjects are from group 2, the user should pass the path to group 2's demographics `.csv` file as the `group_1_demo_file` parameter.

The script can be run in one of two stages: correlating the matrices pairwise and then making visualizations. For each subject, the first stage finds the correlation between its matrix in `--matrices-conc-1` and its matrix in `--matrices-conc-2`. That correlation is saved out to a `.csv` file with all of the others. By default, that file will be called `correlations_out_{CURRENT-DATE-AND-TIME}.csv`.

Stage two makes visualizations of the correlation values and saves them into `.html` files. Each demographic variable to analyze gets two visualizations: a box-and-whisker plot showing how the correlations differ between demographic groups, and a bar plot showing how the correlations differ between each demographic value's median, low outliers, and high outliers. Each box-and-whisker plot is saved as `{demographic variable}_by_demo.html`. Each bar plot is saved as `{demographic variable}_by_outlier.html`. 

- ### `asa_submitter.py` explanation

This script was designed for parallel processing using the [SLURM job scheduling cluster](https://slurm.schedmd.com/overview.html) on the OHSU Exacloud server. `asa_submitter.py` runs many (specifically, `--n-analyses` * the number of parameters in `--subset-size`) batch jobs of `automated_subset_analysis.py` at once. Each job will save its output subsets into a subdirectory of the `--output` directory. The number of subdirectories created will be equl to `--n-analyses`, numbered from 1 to `--n-analyses` like so:

- `./data/output1`
- `./data/output2`
- ...
- `./data/output20` (if `--n-analyses` is 20)

Each job will append its output correlations to 3 `.csv` files in the `--output` directory: 1 for correlating the subsets to each other, and 1 each correlating a group to the subsets of the other group. If there are already correlation `.csv` files in the `--output` directory from previous runs, the output correlation values of any subsequent runs will be appended to the already-existing files. Otherwise, new `.csv` files will be created.

## Metadata

Information about this `README` file:

- Created by Greg Conan, 2019-12-16
- Updated by Greg Conan, 2020-09-18
