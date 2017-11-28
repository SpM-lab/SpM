# To get PDF instead of EPS, execute this script with the option -e "flag_pdf=1"

reset
set key

set style line 1 lt 1 lc 1 pt 2
set style line 2 lt 1 lc 3 pt 2

set style line 99 lt 1 lc 4 pt 7 ps 2.0



set terminal postscript eps enhanced color font "Times-Roman, 22" size 5, 3
# set output "| epstopdf -f -o=SV.pdf"
set output (exists("flag_pdf")) ? "| epstopdf -f -o=SV.pdf" : "SV.eps"
set xlabel "{/Times-Italic l}"
set ylabel "{/Times-Italic s_l}"
# set xlabel "i"
# set ylabel "S_i"
set xrange[0:*]
set yrange[0:*]
set grid

plot "SV.dat"u 1:2 title"" w p ls 99

set output



# set output "| epstopdf -f -o=SV_log.pdf"
set output (exists("flag_pdf")) ? "| epstopdf -f -o=SV_log.pdf" : "SV_log.eps"

set yrange[*:*]
set logscale y
set format y "10^{%L}"

replot

set output
