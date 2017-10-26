#!/usr/bin/gnuplot 
print "to use for different Morb you need to execute: gnuplot -e \"Morb=2\"  thisScript"
print " Full usage: gnuplot -e \"Time_Bgn=0; Time_Fnl=15; dt=0.1; Morb=4\"  $mctdhb_dir/Scripts/Err_of_t.plt"
print "where: Time_Bgn=0 a time point to start from"
print "where: Time_Fnl=0 a time point where to finish"
print "where: dt=0  time step"
print "where: Morb=4 number of MCTDHB orbitals used"
print "where: $mctdhb_dir is the installation directory of the MCTDHB_V3* package"

dirIN="./"
dirOUT="./media/"

if (!exists("Morb")) Morb=1
print "Morb: ", Morb
if (!exists("Time_Bgn")) Time_Bgn='*'
if (!exists("Time_Fnl")) Time_Fnl='*'

Tmin= Time_Bgn
Tmax= Time_Fnl
#Xmin=$minX
#Xmax=$maxX
#Xmin=-0.00001
#Xmax=+0.00001
#Morb=2

#set function style lines
set ticslevel 0.5
set mxtics 10  
set mytics 10
set xtics border nomirror norotate 
set ytics border nomirror norotate 
set ytics border nomirror norotate 
set ztics border nomirror norotate 
set nox2tics 
set noy2tics
set zero 1e-010
set lmargin -1
set bmargin -1
set rmargin -1
set tmargin -1
set locale "C"


#set terminal jpeg medium size 1600x1200
# kommemnt oksa 
#set terminal jpeg  size 1200,600 
#set terminal jpeg enhanced 20 size 1200,600
set terminal jpeg enhanced 20 size 800,600
set terminal jpeg  enh size 1600,1200  font "arial,24"
set grid
set nobar 
set key  bot right
#set term po po enh col solid "Times-Roman" 18
#set term po land  enh col solid "Times-Roman" 18
Occmin=0.00001
Occmax=1.1
X1=Tmin
X2=Tmax
Y1=Occmin
Y2=Occmax
#unset key
zero=1.0e-16
#=============== Err =====================
set title 'Evolution of the Error,  MCTDHB M='.Morb 
set nologscale y 
set xlabel "t" 
set format y "%g"
set ylabel "Absolute Error " 
set nologscale y 

filein=sprintf("%sNO_PR.out",dirIN)
fileout1=sprintf("%sfig_Err_of_t_PLT.jpeg",dirOUT)
set out fileout1

set ylabel "Error " 
set key top left
set format y "%2.1e"
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(Morb+3)) t "Assumed Error" w l lt -1 lw 3,\
filein  using 1:(column(Morb+5)) t "Actual Error" w l lt 1  lw 3


set key bottom left

fileout2=sprintf("%sfig_Err_of_t_Log_Decomp_PLT.jpeg",dirOUT)
set out fileout2

set ylabel 'Error' 
set title "Partioning Error (Log scale), MCTDHB M=".Morb 
set logscale y 
plot  [X1:X2] [1.0e-18:*] \
filein  using 1:Morb+3 t "Assumed Err" w l lt -1 lw 3,\
filein  using 1:Morb+6 t "CI Error %"  w l lt  3 lw 3,\
filein  using 1:Morb+7 t "PSI Error %" w l lt  4 lw 3 

set key top left

fileout3=sprintf("%sfig_Err_of_t_Decomp_PLT.jpeg",dirOUT)
set out fileout3

set title "Partioning of Error to CI and ORB parts " 
set nologscale y 
plot  [X1:X2] [1.0e-18:*] \
filein  using 1:Morb+3 t "Assumed Err" w l lt -1 lw 3,\
filein  using 1:Morb+6 t "CI Error %" w l lt 3 lw 3 ,\
filein  using 1:Morb+7 t "PSI Error %" w l lt 4 lw 3 



print "Error analysis is  plotted"
print "File name  in: ", filein
print "File name out: ", fileout1
print "File name out: ", fileout2
print "File name out: ", fileout3

