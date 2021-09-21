import os
import subprocess
import argparse
import numpy as np
import shutil
import glob

def read_Gtau(path_to_file, n_seed):
    """
    Get Header, tau, green's function from Gtau_*.in.

    Parameters
    ----------
    path_to_file : str
        Path to the parent directory of the samples
    n_seed : int
        The number of samples

    Returns
    -------
    header : list
        List of Setting-informations
    tau : ndarray
        Array of imaginary time (real number)
    g_tau : ndarray
        Array of imaginary time Green's function (complex number)
    """
    file_name = os.path.join(path_to_file, "Gtau_{}.in".format(n_seed))
    with open(file_name, "r") as fr:
        lines = fr.readlines()
        header = lines[:5]
        num = len(lines[5:])
        g_tau = np.zeros(num, dtype=np.complex128)
        tau = np.zeros(num)
        for idx, line in enumerate(lines[5:]):
            values = line.split()
            tau[idx] = values[0]
            g_tau[idx] = float(values[1]) + 1j * float(values[2])
    return header, tau, g_tau


def write_Gtau(path_to_file, file_name, header, tau, g_tau):
    """
    Write Header, tau, green's function in the file

    Parameters
    ----------
    path_to_file : str
        Path to the directory where the file will be written
    file_name : str
        Name of the file in which to write the green's function
    header : list
        List of Setting-informations
    tau : ndarray
        Array of imaginary time (real number)
    g_tau : ndarray
        Array of imaginary time Green's function (complex number)
    """
    path_to_file = os.path.join(path_to_file, file_name)
    with open(path_to_file, "w") as fw:
        for line in header:
            fw.write(line)
        for idx, value in enumerate(tau):
            fw.write("{} {} {}\n".format(value, g_tau[idx].real, g_tau[idx].imag))


def make_Gtau_LOO(path_to_file, nsamples):
    """
    Create directories of datasets for Leave-one-out Cross Validation.

    Parameters
    ----------
    path_to_file : str
        path to the parent directory of the samples
    nsamples : int
        the number of samples used to create datasets
    """
    parent_dir_fold = os.path.join(
        path_to_file, "{}samples_leave_one_out".format(nsamples)
    )
    os.makedirs(parent_dir_fold, exist_ok=True)
    g_tau_list = []
    for idx in range(nsamples):
        header, tau, g_tau = read_Gtau(path_to_file, idx)
        g_tau_list.append(g_tau)

    for idx in range(nsamples):
        g_train = np.zeros_like(g_tau_list[0])
        g_test = np.zeros_like(g_tau_list[0])
        os.makedirs(os.path.join(parent_dir_fold, str(idx)), exist_ok=True)
        path_to_fold = os.path.join(parent_dir_fold, str(idx))
        for tmp_idx in range(nsamples):
            if idx == tmp_idx:
                g_test = g_tau_list[tmp_idx]
            else:
                g_train += g_tau_list[tmp_idx]
        g_train /= nsamples - 1
        write_Gtau(path_to_fold, "G_train.in", header, tau, g_train)
        write_Gtau(path_to_fold, "G_test.in", header, tau, g_test)


def make_Gtau_Kfold(path_to_file, nsamples, kfold):
    """
    Create directories of datasets for K-Fold Cross Validation.

    Parameters
    ----------
    path_to_file : str
        path to the parent directory of the samples
    nsamples : int
        the number of samples used to create datasets
    kfold : int
        the number of folds
    """
    parent_dir_fold = os.path.join(
        path_to_file, "{}samples_{}fold".format(nsamples, kfold)
    )
    os.makedirs(parent_dir_fold, exist_ok=True)
    ksample = int(nsamples / kfold)
    g_tau_list = []
    for idx in range(nsamples):
        header, tau, g_tau = read_Gtau(path_to_file, idx)
        g_tau_list.append(g_tau)

    for idx in range(kfold):
        g_train = np.zeros_like(g_tau_list[0])
        g_test = np.zeros_like(g_tau_list[0])
        os.makedirs(os.path.join(parent_dir_fold, str(idx)), exist_ok=True)
        path_to_fold = os.path.join(parent_dir_fold, str(idx))
        for tmp_idx in range(nsamples - (nsamples % ksample)):
            if idx == int(tmp_idx / ksample):
                g_test += g_tau_list[tmp_idx]
            else:
                g_train += g_tau_list[tmp_idx]
        g_test /= int(nsamples / kfold)
        g_train /= nsamples - (nsamples % ksample) - int(nsamples / kfold)
        write_Gtau(path_to_fold, "G_train.in", header, tau, g_train)
        write_Gtau(path_to_fold, "G_test.in", header, tau, g_test)


def make_Gtau_CV(path_to_file, cv_type, nsamples, kfold):
    """
    Create directories of datasets for Cross Validation.

    Parameters
    ----------
    path_to_file : str
        path to the parent directory of the samples
    cv_type : str
        the type of Crass Validation (LOO, Kfold or MC, LPO)
    nsamples : int
        the number of samples used to create datasets
    kfold : int
        the number of folds (for Kfold)
    """
    if cv_type == "Kfold":
        if kfold < 2:
            print("Error: kfold = {} must be greater than 2.".format(kfold))
            exit(1)
        make_Gtau_Kfold(path_to_file, nsamples, kfold)
    elif cv_type == "LOO":
        make_Gtau_LOO(path_to_file, nsamples)
    else:
        print("Error: cv_type {} is incorrect.".format(cv_type))
        exit(1)


def throw_job_CV(path_to_file, cv_type, nsamples, kfold, path_to_script, cmd):
    """
    Perform SpM calculations on a data set

    Parameters
    ----------
    path_to_file : str
        Path to the parent directory of the samples
    CV_type : str
    nsamples : int
        The number of samples
    kfold : int
    psamples : int
    iter_rand : int
    path_to_script : str
        Path to the script file
    """
    if cv_type == "Kfold":
        fold_dir = "{}samples_{}fold".format(nsamples, kfold)
        num_idx = kfold
    elif cv_type == "LOO":
        fold_dir = "{}samples_leave_one_out".format(nsamples)
        num_idx = nsamples
    else:
        print("Error: cv_type {} is incorrect.".format(cv_type))
        exit(1)

    parent_dir_fold = os.path.join(os.getcwd(), path_to_file, fold_dir)
    files = glob.glob(os.path.join(path_to_script,"*"))
    for idx in range(num_idx):
        os.makedirs(os.path.join(parent_dir_fold, str(idx)), exist_ok=True)
        path_to_fold = os.path.join(parent_dir_fold, str(idx))
        os.chdir(path_to_fold)
        for src_file in files:
            shutil.copy(src_file, os.path.join(path_to_fold))
        subprocess.run(cmd)

parser = argparse.ArgumentParser(
    description="To generate Gtau for CV with white noise from a spectrum",
    add_help=True,
)

parser.add_argument(
    "-i", dest="input_file_toml", default="input.toml", type=str, help="toml file",
)


import toml

args = parser.parse_args()
dict = toml.load(args.input_file_toml)

input_parent_dir = os.path.join(os.getcwd(), dict["file"]["work_dir"])
current_dir = os.getcwd()
path_to_fold = input_parent_dir
make_Gtau_CV(
    path_to_fold,
    dict["cond"]["cv_type"],
    dict["cond"]["nsamples"],
    dict["param"].get("k_fold", 1),
)

if "job" in dict:
    if dict["job"].get("cmd") is not None:
        path_to_job = dict["file"].get("path_to_job", "./spm.sh")
        fullpath_to_job = os.path.join(os.getcwd(), path_to_job)
        throw_job_CV(
            path_to_fold,
            dict["cond"]["cv_type"],
            dict["cond"]["nsamples"],
            dict["param"].get("k_fold", 1),
            fullpath_to_job,
            dict["job"]["cmd"],
        )
    else:
        print("Error: cmd does not exist in job section.")
