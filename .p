set terminal postscript font "Helvetica,20" enhanced
set grid
show grid
set xlabel "Stromdichte (A/cm^2)" font "Helvetica,20"
set ylabel "Zellspannung (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "./.dat" using 1:2 with lines lt 1 lw 5 lc -1 title ""
set terminal postscript font "Helvetica,20" enhanced
set grid
show grid
set xlabel "Stromdichte (A/cm^2)" font "Helvetica,20"
set ylabel "Zellspannung (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "./.dat" using 1:2 with lines lt 1 lw 5 lc -1 title ""
