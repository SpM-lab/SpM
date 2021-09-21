# -*- coding: utf-8 -*-
import sys
import numpy as np

def _read_data(_filename, header_number=0):
    values = []
    with open(_filename, "r") as f:
        for line in f.readlines()[header_number:]:
            values.append([float(s) for s in line.split()])
    return np.array(values)

def _get_diff_max(ref_data, org_data):
    # ref_data and org_data must be ndarray
    diff_abs = np.absolute(ref_data - org_data)
    return np.amax(diff_abs)

def get_diff_max(_filename_ref, _filename_org, header_number=0):
    ref_data = _read_data(_filename_ref, header_number)
    org_data = _read_data(_filename_org, header_number)
    return( _get_diff_max(ref_data, org_data))
    
    
#if __name__ == "__main__":
    #print(get_diff_max("find_lambda_opt.dat", "find_lambda_opt.dat", header_number = 1))
    
