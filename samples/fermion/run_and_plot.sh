#!/bin/sh

# =========================
file_exe="../SpM.out"
dir_plt="../plt"
# =========================

# CHECK FILE AND DIRECTORY
if [ ! -e $file_exe ]; then
    echo "Error: file '$file_exe' not found"
    exit
fi

abs_plt=$(cd $dir_plt; pwd)  # get full path

# RUN
echo "### Running..."
$file_exe -i param.in

# PLOT
echo "### Plotting..."
cd output
gnuplot $abs_plt/*.plt

cd lambda_opt
gnuplot $abs_plt/lambda_fix/*.plt

# Comment-out below if you want PDFs for all data
# cd ../lambda
# for dir in `ls -F | grep /`; do
# 	cd $dir
# 	gnuplot $abs_plt/lambda_fix/*.plt
# 	cd ..
# done
