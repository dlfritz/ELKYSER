set terminal postscript font "Helvetica,20" enhanced
set output "Cathode_Catalyst_10_Pol.ps"
set grid
show grid
set xlabel "Current Density (A/cm^2)" font "Helvetica,20"
set ylabel "Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "../../Data/Cathode_Catalyst_10/Cathode_Catalyst_10.dat" using 1:2 with lines lt 1 lw 5 lc -1 title "Cathode_Catalyst_10"
reset
set terminal postscript font "Helvetica,20" enhanced
set output "Cathode_Catalyst_10_H2Prod.ps"
set grid
show grid
set xlabel "Current Density (A/cm^2)" font "Helvetica,20"
set ylabel "Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "../../Data/Cathode_Catalyst_10/Cathode_Catalyst_10.dat" using 1:3 with lines lt 1 lw 5 lc -1 title "Cathode_Catalyst_10"
reset
set terminal postscript font "Helvetica,20" enhanced
set output "Cathode_Catalyst_10_H2OCons.ps"
set grid
show grid
set xlabel "Current Density (A/cm^2)" font "Helvetica,20"
set ylabel "Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "../../Data/Cathode_Catalyst_10/Cathode_Catalyst_10.dat" using 1:4 with lines lt 1 lw 5 lc -1 title "Cathode_Catalyst_10"
