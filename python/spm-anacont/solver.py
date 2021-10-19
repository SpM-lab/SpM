import numpy as np
import os
import copy
import subprocess


default_spm_input_params = {
   # INPUT/OUTPUT
   #"statistics": "fermion",
   #"beta": 4,
   "filein_G": "gtau.txt",
   "column": 1,
   "fileout_spec": "spectrum.dat",
   # OMEG,
   #"Nomega": 1001,
   #"omegamin": -10,
   #"omegamax": 10,
   # ADMM
   "lambdalogbegin": 2,
   "lambdalogend": -8,
   "tolerance": 1e-10,
   "maxiteration": 1000,
}


valid_kwargs = {
    'statistics': 'fermion',
    'beta': None,
    'Nomega': None,
    'omegamin': None,
    'omegamax':  None,
}


class SpMSolver:
    def __init__(self, exec_path='SpM.out', work_dir='work_SpM', **kwargs):
        self.exec_path = exec_path
        self.work_dir = work_dir

        additional_params = copy.deepcopy(valid_kwargs)
        for k, v in kwargs:
            if k not in valid_kwargs.keys():
                raise ValueError(f"Invalid key {k}!")
            additional_params[k] = v

        self.spm_input_params = copy.deepcopy(default_spm_input_params)
        for k, v in additional_params:
            if v is None:
                raise RuntimeError("kward {k} must be set!")
            self.spm_input_params[k] = v


    def solve(self):
        dir_org = os.getcwd()
        if os.path.exists(self.work_dir):
            if os.path.isdir(self.work_dir):
                os.removedirs(self.work_dir)
            else:
                raise RuntimeError("Path {self.work_dir} exits but it's not a directory!")

        try:
            os.chdir(self.work_dir)

            with open('param.in', 'w') as f:
                for k, v in self.spm_input_params:
                    print(f'{k} = {v}', file=f)

            commands = [self.exec_path, 'param.in']
            return_code = subprocess.call(commands, stdout='stdout.txt', stderr=subprocess.STDOUT)

            if return_code:
                raise RuntimeError("Subprocess of SpM returned error!")
            
            res = np.loadtxt('spectrum.dat')
            # Load additional results

            os.chdir(dir_org)

            return res[:,1] + 1j * res[:,2]

        except Exception as e:
            os.chdir(dir_org)
            return e