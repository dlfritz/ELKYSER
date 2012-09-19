#!/bin/bash


# Create empty file for the gnuplot script.
touch $1.p

# Add commands to the gnuplot script
echo 'set terminal postscript font "Helvetica,20" enhanced' >> $1.p
echo 'set output "'$1'_Pol.ps"' >> $1.p
echo 'set grid' >> $1.p
echo 'show grid' >> $1.p
echo 'set xlabel "Current Density (A/cm^2)" font "Helvetica,20"' >> $1.p
echo 'set ylabel "Voltage (V)" font "Helvetica,20"' >> $1.p
echo 'set xtics font "Helvetica,20"' >> $1.p
echo 'set ytics font "Helvetica,20"' >> $1.p
echo 'set key left top' >> $1.p
echo 'plot "./'$1'.dat" using 1:2 with lines lt 1 lw 5 lc -1 title "'$1'"' >> $1.p

echo 'reset' >> $1.p
echo 'set terminal postscript font "Helvetica,20" enhanced' >> $1.p
echo 'set output "'$1'_H2Prod.ps"' >> $1.p
echo 'set grid' >> $1.p
echo 'show grid' >> $1.p
echo 'set xlabel "Current Density (A/cm^2)" font "Helvetica,20"' >> $1.p
echo 'set ylabel "Voltage (V)" font "Helvetica,20"' >> $1.p
echo 'set xtics font "Helvetica,20"' >> $1.p
echo 'set ytics font "Helvetica,20"' >> $1.p
echo 'set key left top' >> $1.p
echo 'plot "./'$1'.dat" using 1:3 with lines lt 1 lw 5 lc -1 title "'$1'"' >> $1.p

echo 'reset' >> $1.p
echo 'set terminal postscript font "Helvetica,20" enhanced' >> $1.p
echo 'set output "'$1'_H2OCons.ps"' >> $1.p
echo 'set grid' >> $1.p
echo 'show grid' >> $1.p
echo 'set xlabel "Current Density (A/cm^2)" font "Helvetica,20"' >> $1.p
echo 'set ylabel "Voltage (V)" font "Helvetica,20"' >> $1.p
echo 'set xtics font "Helvetica,20"' >> $1.p
echo 'set ytics font "Helvetica,20"' >> $1.p
echo 'set key left top' >> $1.p
echo 'plot "./'$1'.dat" using 1:4 with lines lt 1 lw 5 lc -1 title "'$1'"' >> $1.p


#plot the pol curve
gnuplot	$1.p

# Edit gnuplot file to work from its new location.
rm $1.p
touch $1.p

# Add commands to the gnuplot script
echo 'set terminal postscript font "Helvetica,20" enhanced' >> $1.p
echo 'set output "'$1'_Pol.ps"' >> $1.p
echo 'set grid' >> $1.p
echo 'show grid' >> $1.p
echo 'set xlabel "Current Density (A/cm^2)" font "Helvetica,20"' >> $1.p
echo 'set ylabel "Voltage (V)" font "Helvetica,20"' >> $1.p
echo 'set xtics font "Helvetica,20"' >> $1.p
echo 'set ytics font "Helvetica,20"' >> $1.p
echo 'set key left top' >> $1.p
echo 'plot "../../Data/'$1'/'$1'.dat" using 1:2 with lines lt 1 lw 5 lc -1 title "'$1'"' >> $1.p

echo 'reset' >> $1.p
echo 'set terminal postscript font "Helvetica,20" enhanced' >> $1.p
echo 'set output "'$1'_H2Prod.ps"' >> $1.p
echo 'set grid' >> $1.p
echo 'show grid' >> $1.p
echo 'set xlabel "Current Density (A/cm^2)" font "Helvetica,20"' >> $1.p
echo 'set ylabel "Voltage (V)" font "Helvetica,20"' >> $1.p
echo 'set xtics font "Helvetica,20"' >> $1.p
echo 'set ytics font "Helvetica,20"' >> $1.p
echo 'set key left top' >> $1.p
echo 'plot "../../Data/'$1'/'$1'.dat" using 1:3 with lines lt 1 lw 5 lc -1 title "'$1'"' >> $1.p

echo 'reset' >> $1.p
echo 'set terminal postscript font "Helvetica,20" enhanced' >> $1.p
echo 'set output "'$1'_H2OCons.ps"' >> $1.p
echo 'set grid' >> $1.p
echo 'show grid' >> $1.p
echo 'set xlabel "Current Density (A/cm^2)" font "Helvetica,20"' >> $1.p
echo 'set ylabel "Voltage (V)" font "Helvetica,20"' >> $1.p
echo 'set xtics font "Helvetica,20"' >> $1.p
echo 'set ytics font "Helvetica,20"' >> $1.p
echo 'set key left top' >> $1.p
echo 'plot "../../Data/'$1'/'$1'.dat" using 1:4 with lines lt 1 lw 5 lc -1 title "'$1'"' >> $1.p


# Create directory for output storage
mkdir Figures/$1
mkdir Data/$1

# Move files to simulation directory
mv *.p Figures/$1
mv *.ps Figures/$1
mv *.dat Data/$1
cp PEM_einstellung.in Data/$1
