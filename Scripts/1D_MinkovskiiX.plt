#!/usr/bin/gnuplot 
#    G N U P L O T
# thisScript   1D_MinkovskiiX.plt
print "usage: gnuplot   -e \"Time_Bgn=0; Time_Fnl=15; dt=0.1; Morb=4; Time_xint=-10; Time_xfn=10;\"  thisScript"
print "where: Time_Bgn=0  a time point to start from"
print "where: Time_Fnl=15 a time point where to finish"
print "where: dt=0  time step"
print "where: Time_xint=-5 a X-point to start from"
print "where: Time_xfnl=5 a X-point point where to finish"
print "where: Morb=4 number of MCTDHB orbitals used"



dirIN="DATA/orb_R/"
dirOUT="media/"
if (!exists("Morb")) Morb=1
print "Morb: ", Morb
if (!exists("Time_Bgn")) Time_Bgn=0.0  #Ymin
if (!exists("Time_Fnl")) Time_Fnl=10.0 #Ymax

if (!exists("Time_xint")) Time_xint=-8.0 #Xmin
if (!exists("Time_xfnl")) Time_xfnl=+8.0 #Xmax


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
set zero 1e-08
set lmargin -1
set bmargin -1
set rmargin -1
set tmargin -1
set locale "C"


#set terminal jpeg medium size 1600x1200
#set terminal jpeg  enh size 1200,1200  font "arial,24"
set terminal jpeg enh size 800,600  font "arial,14"
set grid
set nobar 
set key  top right spacing 1.5

set ylabel "Time (t)" 
set xlabel "Space (X)" 

zero=1.0e-16
unset colorbox
#============== PSI(i) ======================
filein =sprintf("%sDNS_Minkovskii.dat",dirIN)
fileout=sprintf("%sDNS_Minkovskii.jpeg",dirOUT)
print "File namein: ", filein, " out: ", fileout
#fileout=sprintf("%10.8ftimeWrkOrb.jpeg",Time)
#filein=sprintf("%s%10.8ftime.dat",dir,Time)
dnsttl=sprintf("DNS MCTDHB(%i)",Morb)
#orbttl(i)=sprintf("Working |PSI_%i|^2",i)

Fulltitle=sprintf(" MCTDHB(%i): rho(x,x\'|t) evolution in  Minkovskii Space-Time representation",Morb)
set title Fulltitle 
set out fileout 
set size 1.0,1.0;
set pm3d map;
set pm3d interpolate 5,5 flush begin ftriangles nohidden3d corners2color mean


minX= Time_xint
maxX= Time_xfnl
minY= Time_Bgn
maxY= Time_Fnl

set xrange [ minX : maxX ] noreverse nowriteback	
set yrange [ minY : maxY ] noreverse nowriteback	
splot [minX:maxX][minY:maxY][0:1.01] filein using 1:4:3 t ""


print "DNS(X,t) in Mincovskii Space-Time representation is plotted"
print "File name  in: ", filein
print "File name out: ", fileout
