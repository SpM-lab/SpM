import numpy as np
import os
import copy
import subprocess


default_spm_input_params = {
   # INPUT/OUTPUT
   #"statistics": "fermion",
   #"beta": 4,
   "filein_G": "gtau.in",
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
        self.stdout_err = None

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


    def solve(self, gtau):
        """
        Perform analytic continuation

        gtau: 1D array of float
            Green's function to be continued
        return: 1D array of float
            Spectral function
        """
        gtau = np.asarray(gtau)
        assert not np.iscomplexobj(gtau)
        assert gtau.ndim == 1

        dir_org = os.getcwd()
        if os.path.exists(self.work_dir):
            if os.path.isdir(self.work_dir):
                os.removedirs(self.work_dir)
            else:
                raise RuntimeError("Path {self.work_dir} is not a directory!")

        try:
            os.chdir(self.work_dir)

            with open('param.in', 'w') as f:
                for k, v in self.spm_input_params:
                    print(f'{k} = {v}', file=f)

            xs = np.linspace(0, 1, gtau.size)
            np.savetxt('gtau.in', np.vstack((xs, gtau)).T)

            commands = [self.exec_path, 'param.in']
            return_code = subprocess.call(commands, stdout='output.txt', stderr=subprocess.STDOUT)
            with open('output.txt') as f:
                self.stdout_err = f.readlines()

            if return_code:
                raise RuntimeError("Subprocess of SpM returned error!")
            
            res = np.loadtxt('spectrum.dat')
            # TODO: Load additional results

            os.chdir(dir_org)

            return res[:,1]

        except Exception as e:
            os.chdir(dir_org)
            return e