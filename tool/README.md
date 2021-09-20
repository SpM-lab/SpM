# Wrapper＆score

## About This Directory

The configuration of this directory is shown below.

```
.
├── wrapper_cv_Gtau.py  # Wrapper for generating dataset for CV
└── calc_cv_score.py    # Calculating score
```

## Usage

### Making datasets for Cross-Validation

`wrapper_cv_Gtau.py` will generates input files for CV.

To use this script, first, the data set of imaginary-time Green's functions must be prepared.
The name of files stored imaginary-time Green's functions must be Gtau_xx.in (xx is the integer).
Then, `wrapper_cv_Gtau.py` generates training and validation files for CV by averaging imaginary-time Green's functions.
You can also do SpM calculation if you set `[job]` section in the input file of `wrapper_cv_Gtau.py`.
You can select two methods, Leave one out and Kfold, for CV.

### Input file

Parameters of the input file for `wrapper_cv_Gtau.py` is explained below.

- \[cond\] section
  - cv_type
    - type: str
    - Description: Set the Cross validation method from LOO or Kfold.
  - nsamples
    - type: int
    - Description: The total number of samples.
- \[param\] section
  - k_fold
    - type: int
    - Description: The number of folds used in Kfold method.
- \[file\] section
  - work_dir
    - type: str
    - Description: The name of work directory
- \[job\] section
  - cmd
    - type: list
    - Description: Command of executing SpM in working directory.  
      ex.) \["./SpM.out"\], \["sh", "./spm.sh"\], \["qsub", "./spm.sh"\]
  - path_to_job
    - type: str
      - Description: Path to base script or execution file.
        This file will be copied to each working directory.

If `cmd` in \[job\] section is not defined, only input files are generated.

### Example

Sample files (`*.toml` and `input/Gtau_xx`) are available at `sample/tool/cv`.

1. Leave one out (LOO)

    In LOO method, one sample is chosen from all samples to make a data set in all combinations.
    Type as follows if the numbar of all samples is 20.

        $ python wrapper_cv_Gtau.py -i input_LOO.toml

    If succeeded, data sets including `G_test.in` and `G_train.in`  are created in `input/20samples_leave_one_out/` directory.

2. Kfold

    In Kfold method, all samples are divided into K parts in order from the head, and each of them is considered as a unit to create a learning data set of K combinations.
    Type as follows if the number of all samples is 20 and K = 10.

        $ python wrapper_cv_Gtau.py -i input_10fold.toml

    If succeeded, data sets are created in `input/20samples_10fold/` directory.

## score

### Calculating score
To calculate the scores, `calc_cv_score.py` is prepared.
The input file which is same as the above wrapper for generating the CV set is needed.
For example, for LOO, you can get score by typing

    $ python calc_cv_score.py -i input_LOO.toml

If succeeded, `input/20samples_10fold/score.dat` is created.
