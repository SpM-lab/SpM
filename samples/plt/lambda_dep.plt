# To get PDF instead of EPS, execute this script with the option -e "flag_pdf=1"

reset
set key

set style line 1 lt 1 lw 2 lc 1 pt 6 ps 2.0
set style line 2 lt 1 lw 2 lc 4 pt 2 ps 2.0
set style line 3 lt 1 lw 2 lc 3 pt 4 ps 2.0



set terminal postscript eps enhanced color font "Times-Roman, 22" size 5, 3
# set output "| epstopdf -f -o=lambda_dep.pdf"
set output (exists("flag_pdf")) ? "| epstopdf -f -o=lambda_dep.pdf" : "lambda_dep.eps"
set xlabel "{/Symbol-Oblique l}"
set ylabel ""
set xrange[*:*]
set yrange[*:*]
set logscale
set format "10^{%L}"
set grid

plot "lambda_dep.dat"u 1:2 title"RMSE (from x')" w lp ls 1,\
 ""u 1:3 title"RMSE (from x)" w lp ls 2,\
 ""u 1:4 title"L1" w lp ls 3,\

set output
