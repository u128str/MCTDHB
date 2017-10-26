#    G N U P L O T
#to use for different Morb you need to execute: gnuplot -e "Morb=2" -e "Time=1.2" thisScript
if (!exists("Morb")) Morb=2
if (!exists("Time")) Time=0.0 
print "Morb: ", Morb, " at t= ",Time
#                            DATA/orb_R/1.00000000time.dat
dirIN="DATA/orb_R/"
dirOUT="media/"
#dirOUT="DATA/orb_R/"

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


#set terminal jpeg medium size 1600x1200
set terminal jpeg  size 1600,1200  font "arial,24"
set terminal jpeg  enh size 1600,1200  font "arial,24"
set grid
set nobar 
set key  top right spacing 1.5

set xlabel "Re" font "arial,30"
#set ylabel "V(x,t),rho, Working Orbitals" font "arial,30" 
set ylabel "x" font "arial,30" 
set zlabel "Im" font "arial,30"

set view 120,30,1,1
#set xrange [-1:1] 

zero=1.0e-16
#============== PSI(i) ======================
filein1=sprintf("%-30.10f",Time)
filein=sprintf("%s%.10stime.dat",dirIN,filein1)
fileout=sprintf("%s%.10stimeWrkOrbReIm_3Dview.jpeg",dirOUT,filein1)
print "File namein: ", filein, " out: ", fileout

#fileout=sprintf("%10.8ftimeWrkOrbReIm.jpeg",Time)
#filein=sprintf("%s%10.8ftime.dat",dir,Time)
dnsttl=sprintf("DNS MCTDHB(%i)",Morb)
orbttl(i)=sprintf("Working |PSI_%i|^2",i)

Fulltitle=sprintf("Working MCTDHB(%i) Orbitals in Re-Im space at T=%f",Morb,Time)
set title Fulltitle 
#fileout=sprintf("%s",dir)
print "File name : ",fileout
set out fileout 
X1=-10
X2= 10
Y1=-0.1
Y2=Morb*0.75+1
scl=0
splot [][X1:X2][-1:1] \
filein using (0.0):($1):(0.02*$5-0.01) t "V(x,t)" w lp lt -1 pt 7 ps 0.0 ,\
for[i=1:Morb] "" using (column(2*i+8)/column(4)-i*scl):1:(column(2*i+9)/column(4)) \
t "[^{WK}psi_".(i)."]" w lp  ,\
"" using (0.0):1:($8) t dnsttl w lp lt 4 lw 3 pt 6 ps 0.7

#filein using (0.02*$5):1:(0.0) t "V(x,t)" w lp  lc rgb "black"  pt 6,\
#for[i=1:Morb] "" using (column(2*i+8)/column(4)-i*scl):($1):(column(2*i+9)/column(4)) t orbttl(i) w lp lt i  pt 6 ps 1.2 lw 2 ,\
#    "" using 8:1:(0.0) t dnsttl w lp lc rgb "#00008B" lw 4 pt 6
