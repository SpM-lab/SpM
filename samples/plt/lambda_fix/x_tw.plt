reset
set key
set zeroaxis lw 1

set style line 1 lt 1 lc 1 lw 4 pt 6
set style line 2 lt 1 lc 2 lw 4 pt 1
set style line 3 lt 1 lc 4 lw 4 pt 1

set style line 99 lt 1 lc 3 lw 2



set terminal postscript eps enhanced color font "Times-Roman, 22" size 5, 3
# set output "| epstopdf -f -o=x_tw.pdf"
set output (exists("flag_pdf")) ? "| epstopdf -f -o=x_tw.pdf" : "x_tw.eps"
set xlabel "{/Symbol-Oblique w}"
# set ylabel "{/Times-Italic x}"
set ylabel "{/Symbol-Oblique r}({/Symbol-Oblique w})"
set xrange[*:*]
set yrange[*:*]

plot "x_tw.dat"u 1:2 title"V x'" w l ls 1,\
 ""u 1:3 title"V z'" w l ls 2,\
 ""u 1:4 title"z" w l ls 3,\
 "../../Gtau.in.dos"u 1:2 title"exact" w l ls 99

set output




# set output "| epstopdf -f -o=x_tw-zoom.pdf"
# set xrange[-0.4:0.4]
# set yrange[*:*]

# replot

# set output
