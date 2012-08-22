reset
set terminal postscript font "Helvetica,20" enhanced
set output "PEM_Elektrolyse_Pol.ps"
set grid
show grid
set xlabel "Stromdichte (A/cm^2)" font "Helvetica,20"
set ylabel "Zellspannung (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica, 20"
set key left top
plot "Data/onemm.dat" using 1:2 with lines lt 1 lw 5 lc -1 title "1 mm Channel", "Data/ponemm.dat" using 1:2 with lines lt 2 lw 5 lc 1 title "0.1 mm Channel"
