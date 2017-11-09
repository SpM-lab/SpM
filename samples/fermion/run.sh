#!/bin/sh

# =========================
file_exe="../SpM.out"
dir_plt="../plt"
output_pdf='n'  # 'y' to get pdf instead of eps
plot_all='n'  # 'y' if you want figures for all data
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
[ $output_pdf = 'y' ] && opt="-e flag_pdf=1"  # set flag to generate pdf
echo "### Plotting..."
cd output
gnuplot $opt $abs_plt/*.plt

# plot results for the optimal value of lambda
cd lambda_opt
gnuplot $opt $abs_plt/lambda_fix/*.plt

# plot results for all values of lambda
if [ $plot_all = 'y' ]; then
	cd ../lambda
	for dir in `ls -F | grep /`; do
		cd $dir
		gnuplot $opt $abs_plt/lambda_fix/*.plt
		cd ..
	done
fi
