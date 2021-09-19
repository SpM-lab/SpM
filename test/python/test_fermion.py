import unittest
import diff
import subprocess
import os

class TestFermion(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestFermion, self).__init__(*args, **kwargs)
        self.file_exe = "../../../../c++/src/SpM.out"
        self.file_path = "./data/fermion/"
        self.ref_path = "./ref/"
        # {filename:header_number}
        self.data_list = {"find_lambda_opt.dat":1,
                          "SV.dat":0,
                          "lambda_dep.dat":0,
                          "pade.dat":0,
                          "spectrum.dat":1}
        
    def test_fermion(self):
        os.chdir(self.file_path)
        subprocess.call(self.file_exe+ " param.in", shell = True)
        #for _key in self.data_list.keys():
            #diff_max=diff.get_diff_max(self.ref_path+_key, "./output/"+_key, header_number = self.data_list[_key])
            #self.assertLess(diff_max, 1e-8)

if __name__ == '__main__':
    unittest.main()
