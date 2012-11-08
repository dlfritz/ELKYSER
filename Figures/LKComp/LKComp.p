set terminal postscript font "Helvetica,20" enhanced
set output "LKComp_Pol.ps"
set grid
show grid
set xlabel "Current Density (A/cm^2)" font "Helvetica,20"
set ylabel "Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set key left top
plot "../../Data/L1K1/L1K1.dat" using 1:2 with lines lt 1 lw 5 lc -1 title "Land/Channel 1:1", "../../Data/L2K1/L2K1.dat" using 1:2 with lines lt 1 lw 5 lc 1 title "Land/Channel 2:1", "../../Data/L3K1/L3K1.dat" using 1:2 with lines lt 1 lw 5 lc 2 title "Land/Channel 3:1"


