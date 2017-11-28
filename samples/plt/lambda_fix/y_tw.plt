reset
set key

set style line 1 lt 1 lc rgb "dark-cyan" lw 2 pt 1 ps 2.0  # input
# set style line 2 lt 3 lc 1 lw 6 dashtype (5,5)  # recovered
set style line 2 lt 3 lc 1 lw 6  # recovered
set style line 3 lt 3 lc 3 lw 7  # exact



set terminal postscript eps enhanced color font "Times-Roman, 22" size 5, 3
# set output "| epstopdf -f -o=y_tw.pdf"
set output (exists("flag_pdf")) ? "| epstopdf -f -o=y_tw.pdf" : "y_tw.eps"
set xlabel "{/Symbol-Oblique t} / {/Symbol-Oblique b}"
# set ylabel "{/Times-Italic y}"
set ylabel "{/Times-Italic G}({/Symbol-Oblique t})"
set xrange[0:1]
#set yrange[0:*]

plot "y_tw.dat"u 1:2 title"input" w p ls 1,\
 "../../Gtau.in"u 1:(abs($3)) title"exact" w l ls 3,\
 "y_tw.dat"u 1:3 title"recovered from x'" w l ls 2

set output




# set output "| epstopdf -f -o=y_tw-1.pdf"
# set xrange[0:0.1]
# set yrange[*:*]
# replot
# set output



# set output "| epstopdf -f -o=y_tw-2.pdf"
# set xrange[0.3:0.7]
# set yrange[*:*]
# replot
# set output



# set output "| epstopdf -f -o=y_tw-log.pdf"
# set xrange[0:1]
# set yrange[*:*]
# set logscale y
# replot
# set output
