reset
set terminal postscript font "Helvetica,20" enhanced
set output "EKOLYSER.ps"
set grid
show grid
set xlabel "Current Density(A/cm^2)" font "Helvetica,20"
set ylabel "Cell Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica, 20"
set key left top
plot "../../Data/Initial_Data/ANOChannel.dat" using 1:2 with lines lt 1 lw 5 lc -1 title "No Flow-field Anode Side", "../../Data/Initial_Data/KT80.dat" using 1:2 with lines lt 2 lw 5 lc 1 title "Serpentine Flow-field"
