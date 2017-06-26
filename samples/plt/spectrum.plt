reset
set key
set zeroaxis lw 1

set style line 1 lt 1 lc 1 lw 3 pt 1
set style line 2 lt 1 lc 3 lw 3



set terminal postscript eps enhanced color font "Times-Roman, 22" size 5, 3
set output "| epstopdf -f -o=spectrum.pdf"
set xlabel "{/Symbol-Oblique w}"
set ylabel "{/Symbol-Oblique r}({/Symbol-Oblique w})"
set xrange[*:*]
#set yrange[0:*]

plot "spectrum.dat"u 1:2 title"" w l ls 1,\
 "../Gtau.in.dos"u 1:2 title"exact" w l ls 2,

set output
