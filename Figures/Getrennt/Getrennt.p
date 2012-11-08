set terminal postscript font "Helvetica,20" enhanced dl 2 
set output "Getrennt_Pol.ps"
set grid
show grid

# Defining Colors
dblue="#006f94"
lblue="#558ed5"
dgrey="#646667"
lgrey="#c5c6cb" 
white="#ffffff"
orang="#e46c0a"
green="#a0ce63"
black="#000000"

# Axis Properties
set xlabel "Current Density (A/cm^2)" font "Helvetica,20"
set ylabel "Voltage (V)" font "Helvetica,20"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
set yrange [1:3]
set key off

# Set arrow style for the voltage ranges
set style arrow 1 head filled size graph 0.02,15,60 lt 1 lw 3 lc rgb black

# Define arrow for IR
set arrow from 0.75,1.548 to 0.75,1.233 as 1
set arrow from 0.75,1.69 to 0.75,1.955 as 1
 
# Define arrow for anodic losses
set arrow from 0.75,2.035 to 0.75,1.971 as 1
set arrow from 0.75,2.12 to 0.75,2.188 as 1

# Define arrow for cathodic losses
set arrow from 0.75,2.315 to 0.75,2.22 as 1
set arrow from 0.75,2.4 to 0.75,2.505 as 1

# Define arrow for the theoretical voltage
set arrow from 0.45,1.15 to 0.48,1.22 as 1

# Draw labels
set label "{/Times-Italic i}R" at 0.735,1.619 front
set label "{/Symbol h}_a" at 0.74,2.085 front
set label "{/Symbol h}_c" at 0.74,2.36 front
set label "V@_{r}" at 0.425,1.08 front

# Draw arrows
show arrow 

plot "../../Data/Getrennt/Getrennt.dat" using 1:8 with lines lt 1 lw 6 lc rgb black title "V_0","" using 1:5 with lines lt 3 lw 4 lc rgb dblue title "IR","" using 1:6 with lines lt 5 lw 4 lc rgb dgrey title "Cathode Act.","" using 1:7 with lines lt 1 lw 4 lc rgb black title "Anode Act."

