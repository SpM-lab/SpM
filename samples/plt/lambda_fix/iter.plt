reset
set key

set style line 1 lt 1 lw 2 lc 1 pt 2 ps 2.0
set style line 2 lt 1 lw 2 lc 3 pt 2 ps 2.0
set style line 3 lt 1 lw 2 lc 2 pt 2 ps 2.0
set style line 4 lt 1 lw 2 lc 4 pt 2 ps 2.0
set style line 5 lt 1 lw 2 lc 5 pt 2 ps 2.0



set terminal postscript eps enhanced color font "Times-Roman, 22" size 5, 3
set output "| epstopdf -f -o=iter-1.pdf"
set xlabel "iteration"
set ylabel ""
set xrange[0:*]
set yrange[*:*]
set logscale y
set format y "10^{%L}"
set grid

plot "iter.dat"u 1:2 title"diff(x, x_{old})" w lp ls 1,\
 ""u 1:3 title"(for L1 regular.) primary residual" w lp ls 2,\
 ""u 1:4 title"dual residual" w lp ls 3,\
 ""u 1:5 title"(for non-negative) primary residual" w lp ls 4,\
 ""u 1:6 title"dual residual" w lp ls 5,\

set output




set output "| epstopdf -f -o=iter-2.pdf"

plot "iter.dat"u 1:7 title"RMSE" w lp ls 1,\
 ""u 1:8 title"|| x' ||_1" w lp ls 2,\
 ""u 1:9 title"sum(x)" w lp ls 3,\
 ""u 1:10 title"negative weight" w lp ls 4,\

set output
