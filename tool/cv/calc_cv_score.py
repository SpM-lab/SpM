import os
import argparse
import numpy as np
import glob
from sklearn.metrics import mean_squared_error


def read_Gtau_test(path_to_file):
    file_name = os.path.join(path_to_file)
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


def read_ytw(path_to_file):
    with open(os.path.join(path_to_file, "y_tw.dat"), "r") as fr:
        lines = fr.readlines()
        num = len(lines)
        g_tau = np.zeros(num, dtype=np.complex128)
        for idx, line in enumerate(lines):
            values = line.split()
            g_tau[idx] = float(values[2]) + 1j * float(values[3])
    return g_tau


def get_cv_type_name(cv_type, nsamples, kfold):
    if cv_type == "Kfold":
        cv_type_name = "{}samples_{}fold".format(nsamples, kfold)
        num_idx = kfold
    elif cv_type == "LOO":
        cv_type_name = "{}samples_leave_one_out".format(nsamples)
        num_idx = nsamples
    else:
        print("Error: cv_type {} is incorrect.".format(cv_type))
        exit(1)
    return cv_type_name, num_idx


parser = argparse.ArgumentParser(
    description="Calculating score functions", add_help=True
)

parser.add_argument(
    "-i", dest="input_file_toml", default="input.toml", type=str, help="toml file",
)


import toml

args = parser.parse_args()
dict = toml.load(args.input_file_toml)

input_parent_dir = os.path.join(os.getcwd(), dict["file"]["work_dir"])
current_dir = os.getcwd()
cv_type_name, num_idx = get_cv_type_name(
    dict["cond"]["cv_type"], dict["cond"]["nsamples"], dict["param"].get("k_fold", 1)
)

path_to_fold = os.path.join(input_parent_dir, cv_type_name)
rmse = []
lambda_list = [
    os.path.basename(p)[7:]
    for p in glob.glob(os.path.join(path_to_fold, str(0), "output", "lambda", "*"))
]
if lambda_list == []:
    print(
        "Error: output directory does not exist in {}.".format(
            os.path.join(path_to_fold, "0")
        )
    )
    exit(1)
score_list = np.zeros((len(lambda_list), 3))
green_list = []

with open(os.path.join(path_to_fold, "score.dat"), "w") as fw:
    for l_idx, dlambda in enumerate(lambda_list):
        print("     lambda: {}/{}".format(l_idx + 1, len(lambda_list)))
        lambda_dir = "lambda_{}".format(dlambda)
        score_list_lambda = np.zeros(num_idx)
        for idx in range(num_idx):
            print("       fold: {}/{}".format(idx + 1, num_idx))
            path_to_train = os.path.join(path_to_fold, str(idx), "G_test.in")
            header, tau, g_tau_test = read_Gtau_test(path_to_train)
            path_to_train = os.path.join(path_to_fold, str(idx), "G_train.in")
            header, tau, g_tau_train = read_Gtau_test(path_to_train)
            path_to_lambda = os.path.join(
                path_to_fold, str(idx), "output", "lambda", lambda_dir
            )
            g_train_opt = -read_ytw(path_to_lambda)
            score_list_lambda[idx] = 1.0 - mean_squared_error(
                g_tau_test.real, g_train_opt.real
            ) / mean_squared_error(
                g_tau_test.real, np.full(g_tau_test.real.shape, g_tau_test.mean().real)
            )
        score_list[l_idx][0] = dlambda
        score_list[l_idx][1] = np.mean(score_list_lambda)
        score_list[l_idx][2] = np.std(score_list_lambda)
        fw.write(
            "{} {} {}\n".format(
                score_list[l_idx][0], score_list[l_idx][1], score_list[l_idx][2]
            )
        )
        fw.flush()
