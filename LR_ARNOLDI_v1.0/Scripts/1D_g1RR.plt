#    G N U P L O T
#to use for different Morb you need to execute: gnuplot -e "Morb=2" -e "Time=1.2" thisScript
if (!exists("Morb")) Morb=2
if (!exists("Time")) Time=0.0 
print "Morb: ", Morb, " at t= ",Time
#                            DATA/orb_R/1.00000000time.dat
dirIN="DATA/g1_RR/"
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
set terminal jpeg  enh size 1200,1200  font "arial,24"
set grid
set nobar 
set key  top right spacing 1.5

set xlabel "x" font "arial,35"
set ylabel "V(x,t),rho, Working Orbitals" font "arial,30" 
set ylabel "x'" font "arial,35" 

zero=1.0e-16
#============== PSI(i) ======================
filein1=sprintf("%-30.10f",Time)
filein =sprintf("%s%.10sx-correlations.dat",dirIN,filein1)
fileout=sprintf("%s%.10sx-correlations.jpeg",dirOUT,filein1)
print "File namein: ", filein, " out: ", fileout
#fileout=sprintf("%10.8ftimeWrkOrb.jpeg",Time)
#filein=sprintf("%s%10.8ftime.dat",dir,Time)
dnsttl=sprintf("DNS MCTDHB(%i)",Morb)
orbttl(i)=sprintf("Working |PSI_%i|^2",i)

Fulltitle=sprintf(" MCTDHB(%i):  |g^{(1)}(x\',x|t)|^2  at T=%f",Morb,Time)
set title Fulltitle 
#fileout=sprintf("%s",dir)
#print "File name : ",fileout
set out fileout 
set size 1.0,1.0;
set pm3d map;
a= 9
minX=-a
maxX= a
minY=-a
maxY= a
cutOFF=0.0000000001

set xrange [ minX : maxX ] noreverse nowriteback	
set yrange [ minY : maxY ] noreverse nowriteback	
splot [minX:maxX][minY:maxY][0:1.01] filein using \
1:4:(($7**2)<cutOFF || ($10**2)<cutOFF   ? 1/0 :(($8**2+$9**2)/($7*$10)>1 ? 1 : ($8**2+$9**2)/($7*$10))) t ""

#filein using ($1):(0.02*$5-0.01) t "V(x,t)" w filledcurves y1=-0.01 lc rgb "black" fs transparent solid 0.5 noborder,
#"$i" using 1:4:(((\$8**2+\$9**2)/(\$7*\$10)>1) ? 1 : (\$8**2+\$9**2)/(\$7*\$10)) t "X-Space g^{(1)} MCTDHB(\$4)"
#(($7**2)<cutOFF || ($10**2)<cutOFF   ? 1/0 :(($8**2+$9**2)/($7*$10)>1 ? 1 : ($8**2+$9**2)/($7*$10))) t ""

