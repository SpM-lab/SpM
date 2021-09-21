# To get PDF instead of EPS, execute this script with the option -e "flag_pdf=1"

reset
set key

set style line 1 lt 1 lw 2 lc 1 pt 6 ps 2.0
set style line 2 lt 1 lw 2 lc 4 pt 2 ps 2.0
set style line 3 lt 1 lw 2 lc 3 pt 4 ps 2.0



set terminal postscript eps enhanced color font "Times-Roman, 22" size 5, 3
# set output "| epstopdf -f -o=find_lambda_opt.pdf"
set output (exists("flag_pdf")) ? "| epstopdf -f -o=find_lambda_opt.pdf" : "find_lambda_opt.eps"
set xlabel "log_{10} {/Symbol-Oblique l}"
set ylabel "Reduction  log_{10} ( {/Times-Italic f} / {/Symbol-Oblique c}^2 )"
set xrange[*:*]
set yrange[*:*]
# set logscale
# set format "10^{%L}"
set grid

plot "find_lambda_opt.dat"u 1:2 title"" w lp ls 1

set output
