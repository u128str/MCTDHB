#!/usr/bin/gnuplot 
print "to use for different Morb you need to execute: gnuplot -e \"Morb=2\"  thisScript"
print " Full usage: gnuplot -e \"Time_Bgn=0; Time_Fnl=15; dt=0.1; Morb=4\"  $mctdhb_dir/Scripts/E_of_t.plt"
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
set terminal jpeg enhanced 20 size 800,600
set terminal jpeg  enh size 1600,1200  font "arial,24"
set grid
set nobar 
set key  bot right
#set term po po enh col solid "Times-Roman" 18
#set term po land  enh col solid "Times-Roman" 18
fileout1="./media/fig_E_of_t_PLT.jpeg"
set out fileout1
Occmin=0.00001
Occmax=1.1
X1=Tmin
X2=Tmax
Y1=Occmin
Y2=Occmax
unset key
#=============== E(t) =====================
#set out "fig_E_of_t.ps"
set title 'Evolution of the system, MCTDHB M='.Morb
set nologscale y 
set ylabel "E(t)/N"
set xlabel "time"
set format y "%10.6e"
plot  [X1:X2] [] \
"NO_PR.out"  using 1:Morb+2 t "E_{total}" w l lt 1 lw 3

set format y "%5.2f"
unset format 
#OP_PR.out:
#   1:time    2:User OP1=x   3:User OP2=r*r   4: Potential V   5:Kinetic T   6:Interaction W  7: Total E
#   1:time 2:<x>  3:<y>  4:<z>   5:<x*x>2body 6:<x*x>1body 7:Var(x) 8:<y*y>2body 9:<y*y>1body 10:Var(y) 
#   11:<z*z>2body 12:<z*z>1body 13:Var(z)  14: <V>    15:<T>     16: <W>     17: <E>           
filein=sprintf("%sOP_PR.out",dirIN)
fileout2=sprintf("%sfig_E_of_t_T_V_W_PLT.jpeg",dirOUT)
set out fileout2

set ylabel "E(t), T(t), V(t), W(t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(14)) t "V(t)" w l lt 1 lw 3,\
filein  using 1:(column(15)) t "T(t)" w l lt 2 lw 3,\
filein  using 1:(column(16)) t "W(t)" w l lt 3 lw 3,\
filein  using 1:(column(17)) t "E(t)" w l lt 4 lw 3

fileout3=sprintf("%sfig_OP_of_t_R_PLT.jpeg",dirOUT)
set out fileout3
set ylabel "Exp. val. of Operators [<X>, <Y>, <Z>](t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(2)) t "Sum <x_i>" w l lt 1 lw 3,\
filein  using 1:(column(3)) t "Sum <y_i>" w l lt 2 lw 3,\
filein  using 1:(column(4)) t "Sum <z_i>" w l lt 3 lw 3


fileout4=sprintf("%sfig_OP_of_t_RR_PLT.jpeg",dirOUT)
set out fileout4
set ylabel "Exp. val. of Operators [<X^2>, <Y^2>, <Z^2>](t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(5) +column(6))  t "<(Sum x_i)^2>" w l lt 1 lw 3,\
filein  using 1:(column(8) +column(9))  t "<(Sum y_i)^2>" w l lt 2 lw 3,\
filein  using 1:(column(11)+column(12)) t "<(Sum z_i)^2>" w l lt 3 lw 3




print "Energy, V, T, W, <x>,<x*x> are plotted"
print "File name  in: ", filein 
print "File name out: ", fileout1
print "File name out: ", fileout2
print "File name out: ", fileout3
print "File name out: ", fileout4
