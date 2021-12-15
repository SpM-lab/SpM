from typing import Iterable, Union

import os
import os.path
import copy
import subprocess
import shutil
import numpy as np

default_spm_input_params = {
    # INPUT/OUTPUT
    # "statistics": "fermion",
    # "beta": 4,
    "filein_G": "gtau.in",
    "column": 1,
    "fileout_spec": "spectrum.dat",
    # OMEG,
    # "Nomega": 1001,
    # "omegamin": -10,
    # "omegamax": 10,
    # ADMM
    "tolerance": 1e-10,
    "maxiteration": 1000,
}


# If None, parameter should be set by user
valid_kwargs = {
    # MODEL
    "statistics": "fermion",
    "beta": None,
    # OMEGA
    "Nomega": 1001,
    "omegamin": None,
    "omegamax": None,
    # ADMM
    "SVmin": 1.0e-10,
    "lambdalogbegin": 2,
    "lambdalogend": -8,
    "lambdalogmesh": 0.2,
    # SpM-Pade
    "PadeEta": 0.0,
    "filein_G_sigma": "gtau.in",
    "column_sigma": 1,
    "g_sigma": "inf",
}


class SpMSolver:
    def __init__(
        self, exec_path: str = "SpM.out", work_dir: str = "work_SpM", **kwargs
    ):
        if os.path.dirname(exec_path) != "":
            self.exec_path = os.path.abspath(exec_path)
            if not os.access(self.exec_path, mode=os.X_OK):
                raise RuntimeError(f"{exec_path} is not found")
        else:
            # solve PATH
            found = False
            for P in os.environ["PATH"].split(":"):
                self.exec_path = os.path.join(P, exec_path)
                if os.access(self.exec_path, mode=os.X_OK):
                    found = True
                    break
            if not found:
                raise RuntimeError(f"{exec_path} is not found")
        self.work_dir = work_dir
        self.stdout_err = None

        additional_params = copy.deepcopy(valid_kwargs)
        for k, v in kwargs.items():
            if k not in valid_kwargs.keys():
                raise ValueError(f"Invalid key {k}!")
            additional_params[k] = v

        self.spm_input_params = copy.deepcopy(default_spm_input_params)
        for k, v in additional_params.items():
            if v is None:
                raise RuntimeError(f"kwarg {k} must be set!")
            self.spm_input_params[k] = v

    def solve(self, gtau: Union[Iterable, np.ndarray]) -> np.ndarray:
        """
        Perform analytic continuation

        gtau: 1D array of float
            Green's function to be continued
        return: 1D array of float
            Spectral function
        """
        gtau = np.asarray(gtau)
        if np.iscomplexobj(gtau):
            raise RuntimeError("gtau is complex array")
        if gtau.ndim != 1:
            N = len(gtau)
            for n in gtau.shape:
                if n != N and n != 1:
                    raise RuntimeError("gtau is not 1D array")
            gtau = gtau.reshape(N)

        if os.path.exists(self.work_dir):
            if os.path.isdir(self.work_dir):
                shutil.rmtree(self.work_dir)
            else:
                raise RuntimeError(f"Path {self.work_dir} is not a directory!")

        dir_org = os.getcwd()
        try:
            os.mkdir(self.work_dir)
            os.chdir(self.work_dir)

            with open("param.in", "w") as f:
                for k, v in self.spm_input_params.items():
                    print(f"{k} = {v}", file=f)

            xs = np.linspace(0, 1, gtau.size)
            np.savetxt("gtau.in", np.vstack((xs, gtau)).T)

            commands = [self.exec_path, "param.in"]
            with open("output.txt", "w") as f:
                return_code = subprocess.call(
                    commands, stdout=f, stderr=subprocess.STDOUT
                )
            with open("output.txt", "r") as f:
                self.stdout_err = f.readlines()

            if return_code:
                raise RuntimeError("Subprocess of SpM returned error!")

            res = np.loadtxt("output/spectrum.dat")

            return res[:, 1]
        finally:
            os.chdir(dir_org)
