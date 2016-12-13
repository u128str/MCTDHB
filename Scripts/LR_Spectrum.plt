#!/usr/bin/gnuplot 
print "to use for different Morb you need to execute: gnuplot -e \" Morb=4; Emin=0.0001; Emax=5.0 \"  thisScript"
print " Full usage: gnuplot -e \"Morb=4\"   $mctdhb_dir/Scripts/LR_Spectrum.plt"
print "where: Morb= MCTDHB(Morb) levels"
print "     : Emin minimal value of excitation energy "
print "     : Emax maximal value of excitation energy "
print "     : $mctdhb_dir is the installation directory of the MCTDHB_V3* package"

dirIN="DATA/getLR/"
dirOUT="./media/"

filein=sprintf("%sMC_anlsplot.out",dirIN)
fileout=sprintf("%sLR_Spectrum.jpeg",dirOUT)

if (!exists("Morb")) Morb=4
if (!exists("Emin")) Emin=0.00001
if (!exists("Emax")) Emax=5


#set terminal jpeg enh size 1600,1200  font "arial,24"
set terminal jpeg enh size 800,600  font "arial,14"
lwv=1
set out fileout
#set terminal postscript landscape enhanced color solid defaultplex "Helvetica" 20 
#set terminal postscript landscape enhanced color dash defaultplex "Helvetica" 20 
#set output 'LR_Spectrum.ps'
#set key title "{/Symbol=20 r}(k,t); N=2"
set title 'Excitation spectrum obtained with LR-MCTDHB M='.Morb
unset title
set key at 3.0,10 box ls -1 lw lwv
set nokey 
set samples 1000,1000
set border 4095 lt -1 lw 1.00
set size 1.0, 1.0
###########################################
###########    LEFT
#########################################
set multiplot
set size 1.0,0.50
set origin 0.0, 0.0
unset label
set fit noerrorvariables
set xlabel "LR-excitation energy"
set ylabel "Intensity (resp. ampl. to <x>)" offset 0.8,+0
#set grid 
set ytics 0,1,2
set mytics 4
set label 1 "Response to ungerade <x> perturbation" at 0.1,1.1
set y2label "  " offset 0.8,0
set y2tics (" " 0)
plot [Emin:Emax][0:*] \
filein u ($2/sqrt(1)):($6) w p lt 1  pt 9 ps 1.7  ,\
filein u ($2/sqrt(1)):($6) w i lt 1 lw lwv+2

unset label 
set title 'Excitation spectrum obtained with LR-MCTDHB(M) M='.Morb
set xlabel " "
set ylabel " " offset 0.8,0
set ylabel "Intensity (resp. ampl. to <x^2>)" offset 0.8,+0
set size 1.0,0.55
set origin 0.0, 0.43
set mx2tics 
set mxtics 20 
set xtics scale 1,1 
set xtics -1,10,100
set ytics 0,1,2
set mytics 4
set label "Response to gerade <x^2> perturbation" at 0.1,1.1
set y2label "  " offset 0.8,0
set y2tics (" " 0)
plot [Emin:Emax] [0:*]\
filein u ($2/sqrt(1)):7 w i lt 2 lw lwv+2,\
filein u ($2/sqrt(1)):7 w p lt 2  pt 5 ps 1.2


print "LR spectrum vs Response amplitudes is plotted"
print "File name  in: ", filein 
print "File name out: ", fileout


###########################################
###########    UPPER INSET
#########################################

#set boxwidth 0.05 absolute
#set style fill solid 0.3 noborder
#
#unset label 
#unset xlabel 
#unset ylabel
#unset title
#unset arrow
#unset key
##set arrow from 1.0,0.0 to  1.0,25 nohead lt -1 lw 1.5 
##set arrow from 1.1067,0.0 to  1.1067,25 nohead lt -1 lw 1.5
#set origin 0.0,0.72
#set size 1.0,0.25
#set noxtics 
##set ytics (-1,1,0.2)
#set ytics ("" 0, "" 0.25,"" 0.5,"" 0.75,"  " 1)
#set y2tics ("0" 0, "" 0.25,"" 0.5,"" 0.75,"1" 1)
##set mytics 4
##set noytics
#
#set xlabel "" 
#set ylabel "||{CI}||" offset 0.8,0
#set ylabel "  " offset 0.8,0
#set y2label "Type" offset 0.8,0
##set noylabel 
##set y2label "Type" offset 0.8,0
#set mx2tics 
#set mxtics 20 
#set xtics scale 1,1 
#set xtics -1,10,100
#
#
#plot [0.0001:3.1][0:1.1] filein u ($2/sqrt(2)):(($2<0.1) ? 1/0 :$5) with boxes lt -1
#
###################### Trap and Density
#set size 0.3,0.2
#set origin 0.6,0.58
#set notics
#set noborder
#set nolabel
#set noy2label
#plot [-3:3][0:0.7] "./M2/20.0000000time.dat" u ($1):($5/4) with l lt -1, "" u ($1):($6) with l lt 1 lw 3
###################### Trap and Density
#set size 0.30,0.25
#set origin 0.6, 0.211
#set notics
#set noborder
#set nolabel
#set noy2label
#plot [-3:3][0:0.7] "./M2/20.0000000time.dat" u ($1):($5/4) with l lt -1, "" u ($1):($6) with l lt 1 lw 3
