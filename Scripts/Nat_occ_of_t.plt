#!/usr/bin/gnuplot 
print "to use for different Morb you need to execute: gnuplot -e \"Morb=2\"  thisScript"
print " Full usage: gnuplot -e \"Time_Bgn=0; Time_Fnl=15; dt=0.1; Morb=4\"  $mctdhb_dir/Scripts/Nat_occ_of_t.plt"
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
#Tmin=zero
#Tmax=0.05

#set function style lines
#set ticslevel 0.5
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

#set format y "%4.2e"

#set terminal jpeg medium size 1600x1200
#set terminal jpeg  enh size 1600,1200  font "arial,24"
set terminal jpeg  enh size 800,600  font "arial,14"
lwv=1
set grid
set nobar 
set key  bot right
#set term po po enh col solid "Times-Roman" 18
#set term po land  enh col solid "Times-Roman" 18

filein="NO_PR.out"
set title 'Evolution of natural occupations: MCTDHB(M) M='.Morb
set xlabel "t"
set ylabel "Log n_i/N" 
set logscale y 10

#zero=1.0e-16
#============== Nocc(i) ======================
# oksa add _PLT
fileout1="./media/fig_Nat_occ_of_t_log_PLT.jpeg"
set out fileout1
err=9.0e-10
Occmin=err
Occmax=1.1
X1=Tmin
X2=Tmax
Y1=Occmin
Y2=Occmax
plot  [X1:X2] [Y1:Y2] 0.5 t "50x50" w l lt -1,\
for[i=Morb:1:-1] "NO_PR.out" using 1:(column(i+1)) t 'n_'.(Morb-i+1) with lines linetype i lw lwv


unset logscale 
set ylabel "n_i/N" 
fileout2="./media/fig_Nat_occ_of_t_PLT.jpeg"
set out fileout2
set out "./media/fig_Nat_occ_of_t_PLT.jpeg"
plot  [X1:X2] [Y1:Y2] 0.5 t "50x50" w l lt -1,\
for[i=Morb:1:-1] "NO_PR.out" using 1:(column(i+1)) t 'n_'.(Morb-i+1) with lines linetype i lw lwv

print "Natural occupation n_i(t) are plotted"
print "File name  in: ", filein, X1,X2
print "File name out: ", fileout1
print "File name out: ", fileout2
