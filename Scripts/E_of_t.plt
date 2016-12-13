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
#set terminal jpeg enhanced 20 size 800,600
#set terminal jpeg  enh size 1600,1200  font "arial,24"
set terminal jpeg  enh size 800,600  font "arial,14"
lwv=3
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
"NO_PR.out"  using 1:Morb+2 t "E_{total}" w l lt 1 lw lwv

set format y "%5.2f"
unset format 
#OP_PR.out:
#  1: time  
#  2: <E>/N     3:  <T>/N    4: <V>/N  5: <W>/N  
#  6: <x>/N     7:  <y>/N    8:   <z>/N  
#  9: <x*x>/N  10: <y*y>/N  11: <z*z>/N 
# 12: Var(x)/N 13: Var(y)/N 14: Var(z)/N 
# 15:  <px>/N  16: <py>/N   17:  <pz>/N  
# 18:<px*px>/N 19:<py*py>/N 20:<pz*pz>/N
# 21:Var(px)/N 22:Var(py)/N 23:Var(pz)/N
filein=sprintf("%sOP_PR.out",dirIN)
fileout2=sprintf("%sfig_E_of_t_T_V_W_PLT.jpeg",dirOUT)
set out fileout2

set ylabel "E(t), T(t), V(t), W(t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(4)) t "V(t)/N" w l lt 1 lw lwv,\
filein  using 1:(column(3)) t "T(t)/N" w l lt 2 lw lwv,\
filein  using 1:(column(5)) t "W(t)/N" w l lt 3 lw lwv,\
filein  using 1:(column(2)) t "E(t)/N" w l lt 4 lw lwv

fileout3=sprintf("%sfig_OP_of_t_R_PLT.jpeg",dirOUT)
set out fileout3
set ylabel "Exp. val. of position operators [<X>, <Y>, <Z>](t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(6)) t "Sum <x_i>/N" w l lt 1 lw lwv,\
filein  using 1:(column(7)) t "Sum <y_i>/N" w l lt 2 lw lwv,\
filein  using 1:(column(8)) t "Sum <z_i>/N" w l lt 3 lw lwv 


fileout4=sprintf("%sfig_OP_of_t_RR_PLT.jpeg",dirOUT)
set out fileout4
set ylabel "Exp. val. [<X^2>, <Y^2>, <Z^2>](t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(9))  t "<(Sum x_i)^2>/N" w l lt 1 lw lwv,\
filein  using 1:(column(10)) t "<(Sum y_i)^2>/N" w l lt 2 lw lwv,\
filein  using 1:(column(11)) t "<(Sum z_i)^2>/N" w l lt 3 lw lwv

fileout5=sprintf("%sfig_OP_of_t_VarRR_PLT.jpeg",dirOUT)
set out fileout5
set ylabel "Var. [<X^2>-<X>^2, <Y^2>-<Y>^2, <Z^2>-<Z>^2](t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(12)) t "<(Sum x_i)^2>/N-(Sum<x_i>)^2/N" w l lt 1 lw lwv,\
filein  using 1:(column(13)) t "<(Sum y_i)^2>/N-(Sum<y_i>)^2/N" w l lt 2 lw lwv,\
filein  using 1:(column(14)) t "<(Sum z_i)^2>/N-(Sum<z_i>)^2/N" w l lt 3 lw lwv

fileout6=sprintf("%sfig_OP_of_t_dP_PLT.jpeg",dirOUT)
set out fileout6
set ylabel "Exp. val. of momentum operators [<Px>, <Py>, <Pz>](t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(15)) t "Sum <Px_i>/N" w l lt 1 lw lwv,\
filein  using 1:(column(16)) t "Sum <Py_i>/N" w l lt 2 lw lwv,\
filein  using 1:(column(17)) t "Sum <Pz_i>/N" w l lt 3 lw lwv 

fileout7=sprintf("%sfig_OP_of_t_dPP_PLT.jpeg",dirOUT)
set out fileout7
set ylabel "Exp. val. [<Px*Px>, <Py*Py>, <Pz*Pz>](t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(18)) t "(<Sum Px_i)^2>/N" w l lt 1 lw lwv,\
filein  using 1:(column(19)) t "(<Sum Py_i)^2>/N" w l lt 2 lw lwv,\
filein  using 1:(column(20)) t "(<Sum Pz_i)^2>/N" w l lt 3 lw lwv

fileout8=sprintf("%sfig_OP_of_t_VarPP_PLT.jpeg",dirOUT)
set out fileout8
set ylabel "Var. [<Px^2>-<Px>^2, <Py^2>-<Py>^2, <Pz^2>-<Pz>^2](t)" 
set key top left spacing 1.5
err=9.0e-3
plot  [X1:X2] [] \
filein  using 1:(column(21)) t "<(Sum Px_i)^2>/N-(Sum<Px_i>)^2/N" w l lt 1 lw lwv,\
filein  using 1:(column(22)) t "<(Sum Py_i)^2>/N-(Sum<Py_i>)^2/N" w l lt 2 lw lwv,\
filein  using 1:(column(23)) t "<(Sum Pz_i)^2>/N-(Sum<Pz_i>)^2/N" w l lt 3 lw lwv


print "Energy, V, T, W, <x>, <x*x>, <p> <p*p> ...  are plotted"
print "File name  in: ", filein 
print "File name out: ", fileout1
print "File name out: ", fileout2
print "File name out: ", fileout3
print "File name out: ", fileout4
print "File name out: ", fileout5
print "File name out: ", fileout6
print "File name out: ", fileout7
print "File name out: ", fileout8
