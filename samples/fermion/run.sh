#!/bin/sh

# =========================
file_exe="../SpM.out"
dir_plt="../plt"
plot_level=1  # 0: no plot, 1: plot major data, 2: plot all data
eps_or_pdf='eps'  # 'eps' or 'pdf' (epstopdf is required)
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
[ $eps_or_pdf = 'pdf' ] && opt="-e flag_pdf=1"  # set flag to generate pdf
if [ $plot_level -ge 1 ]; then
    echo "### Plotting..."
    cd output
    gnuplot $opt $abs_plt/*.plt

    # plot results for the optimal value of lambda
    cd lambda_opt
    gnuplot $opt $abs_plt/lambda_fix/*.plt
fi
if [ $plot_level -ge 2 ]; then
    # plot results for all values of lambda
    if [ $plot_all = 'y' ]; then
    cd ../lambda
    for dir in `ls -F | grep /`; do
        cd $dir
        gnuplot $opt $abs_plt/lambda_fix/*.plt
        cd ..
    done
    fi
fi
