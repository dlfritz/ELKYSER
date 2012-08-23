set terminal postscript font "Helvetica,20" enhanced
set output "Figures/TMMEEComp.ps"
set grid
show grid
set xlabel "Stromdichte (A/cm^2)" font "Helvetica,20"
set ylabel "Zellspannung (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica, 20"
set key left top
#set xrange [0:1.8]
plot "Data/TMLW.dat" using 1:2 with lines lt 1 lw 5 lc -1 title "Liquid Water", "Data/TMMV.dat" using 1:2 with lines lt 1 lw 5 lc 1 title "5% Vapor"
#,"Data/CNOChannel.dat" using 1:2 with lines lt 1 lw 5 lc 2 title "No Channels on Cathode Side","Data/BNOChannel.dat" using 1:2 with lines lt 1 lw 5 lc 3 title "No Channels"
