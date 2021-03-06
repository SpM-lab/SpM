reset
set key
set grid

set style line  1 lt 1 pt 1 lw 3 ps 2.2 lc rgb "dark-cyan"  # input (with noise)
set style line  2 lt 1 pt 6 lw 3 ps 2.4 lc 1  # result
set style line  3 lt 1 pt 2 lw 3 ps 1.8 lc 3  # exact
set style line  4 lt 1 pt 8 lw 3 ps 2.4 lc 2  # result

set style line 99 lt 1 lc 4 pt 7 ps 1.0



set terminal postscript eps enhanced color font "Times-Roman, 22" size 5, 3
# set output "| epstopdf -f -o=x_sv.pdf"
set output (exists("flag_pdf")) ? "| epstopdf -f -o=x_sv.pdf" : "x_sv.eps"
set xlabel "{/Times-Italic l}"
# set ylabel "| {/Times-Italic x'_l} |"
# set ylabel "| {/Times-Italic x'_l} | / {/Symbol Dw}^{1/2}"
set ylabel "| {/Symbol-Oblique r}@{/Symbol \\242}_{/Times-Italic l} |"
set xrange[0:60]
set yrange[0:*]

plot "../../Gtau.in.sv_basis"u 1:(abs($2)) title"exact" w p ls 3,\
 "x_sv.dat"u 1:(abs($3)) title"z'" w p ls 2
 
#  "x_sv.dat"u 1:(abs($2)) title"x'" w p ls 4,\

set output



# set output "| epstopdf -f -o=x_sv-log.pdf"
set output (exists("flag_pdf")) ? "| epstopdf -f -o=x_sv-log.pdf" : "x_sv-log.eps"

set yrange[1e-10:*]
set logscale y
set format y "10^{%L}"

replot
# replot "../SV.dat"u 1:2 title"SV" w p ls 99

set output
