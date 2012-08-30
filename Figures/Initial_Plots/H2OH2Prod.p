set terminal postscript font "Helvetica,20" enhanced
set output "Figures/Consumption.ps"
set grid
show grid
set xlabel "Stromdichte (A/cm^2)" font "Helvetica,20"
set ylabel "Zellspannung (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica, 20"
set key left top
plot "Data/Kickoff1.dat" using 1:3 with lines lt 1 lw 5 lc -1 title "Hydrogen Production Nl/m^3", "Data/Kickoff1.dat" using 1:4 with lines lt 1 lw 5 lc 1 title "Water Consumption" 
#,"Data/CNOChannel.dat" using 1:2 with lines lt 1 lw 5 lc 2 title "No Channels on Cathode Side","Data/BNOChannel.dat" using 1:2 with lines lt 1 lw 5 lc 3 title "No Channels"
