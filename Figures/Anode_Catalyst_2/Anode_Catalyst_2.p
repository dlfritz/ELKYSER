set terminal postscript font "Helvetica,20" enhanced
set output "Anode_Catalyst_2_Pol.ps"
set grid
show grid
set xlabel "Current Density (A/cm^2)" font "Helvetica,20"
set ylabel "Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "../../Data/Anode_Catalyst_2/Anode_Catalyst_2.dat" using 1:2 with lines lt 1 lw 5 lc -1 title "Anode_Catalyst_2"
reset
set terminal postscript font "Helvetica,20" enhanced
set output "Anode_Catalyst_2_H2Prod.ps"
set grid
show grid
set xlabel "Current Density (A/cm^2)" font "Helvetica,20"
set ylabel "Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "../../Data/Anode_Catalyst_2/Anode_Catalyst_2.dat" using 1:3 with lines lt 1 lw 5 lc -1 title "Anode_Catalyst_2"
reset
set terminal postscript font "Helvetica,20" enhanced
set output "Anode_Catalyst_2_H2OCons.ps"
set grid
show grid
set xlabel "Current Density (A/cm^2)" font "Helvetica,20"
set ylabel "Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "../../Data/Anode_Catalyst_2/Anode_Catalyst_2.dat" using 1:4 with lines lt 1 lw 5 lc -1 title "Anode_Catalyst_2"
