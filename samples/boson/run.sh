#!/bin/sh

file_exe="../SpM.out"
dir_plt="../plt"

# INPUT/OUTPUT
statistics="boson"  # "fermion" or "boson"
beta=100
filein_G="Gtau.in"
column=1  # start from 0
fileout_spec="spectrum.dat"

# OMEGA
N_omega=1001
omega_min=-4
omega_max=4

# SVD
SV_min=1e-10  # 0 for no truncation

# ADMM
N_lambda=51
lambda_begin=1e+2
lambda_end=1e-8
penalty=1.  # if negative, penalty is optimized during the iteration starting with its absolute value
tolerance=1e-10
max_iter=1000
print_level=1  # 0: minimum, 1: moderate, 2: verbose

# CROSS VALIDATION
cross_validation="n"  # y/n

# =========================== RUN
$file_exe $statistics $beta $filein_G $column $fileout_spec $N_omega $omega_min $omega_max $SV_min $N_lambda $lambda_begin $lambda_end $penalty $tolerance $max_iter $print_level $cross_validation
# ===========================

# PLOT
echo "plotting..."
cd output
gnuplot ../$dir_plt/*.plt
cd lambda_opt
gnuplot ../../$dir_plt/lambda_fix/*.plt

# execute below if you want pdf for all data
# cd ../lambda
# for dir in `ls -F | grep /`; do
# 	cd $dir
# 	gnuplot ../../../$dir_plt/lambda_fix/*.plt
# 	cd ..
# done
