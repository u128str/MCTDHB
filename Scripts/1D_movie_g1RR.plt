#!/usr/bin/gnuplot 
#    G N U P L O T
# thisscript 1D_movie_g1RR.plt
print "GIF usage: gnuplot  -e \"Morb=2; Time_Bgn=0.0; Time_Fnl=11; dt=0.2 \"  thisscript"
print "OUT animated gif: ../media/1D_g1RR.gif"
print "where: Time_Bgn=0 a time point to start from"
print "where: Time_Fnl=0 a time point where to finish"
print "where: dt=0  time step"
print "where: Morb=4 number of MCTDHB orbitals used"
print "where: $mctdhb_dir is the installation directory of the MCTDHB_V3* package"
print "JPG usage: gnuplot  -e \"Morb=2; Time_Bgn=0.0 \"  $mctdhb_dir/Scripts/thiscript"
print "OUT single-shot jpeg: media/0.00000000x-correlations.jpeg "

if (!exists("Morb")) Morb=2
if (!exists("Time_Bgn")) Time_Bgn=0
if (!exists("Time_Fnl")) Time_Fnl=10
Time=Time_Bgn
print "Morb: ", Morb, " at t= ",Time
PIC=2
if (!exists("dt")) {dt=0;Time_Fnl=Time_Bgn;PIC=1;Time=Time_Bgn}
Tmin= Time_Bgn
Tmax= Time_Fnl
dt=dt
n=1
if (dt!=0) n=1+(Tmax-Tmin)/dt #n frames
print "PIC=",PIC," n=",n, " dt=", dt," Time_Bgn=",Time_Bgn," Time_Fnl=",Time_Fnl
#                            DATA/orb_R/1.00000000time.dat
dirIN="DATA/g1_RR/"
dirOUT="media/"

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
#set terminal jpeg  size 1600,1200  font "arial,24"
#set terminal jpeg  enh size 1600,1200  font "arial,24"
if(PIC==1) {
#set terminal jpeg enh  size 1600,1200  font "arial,24"
set terminal jpeg enh size 800,600  font "arial,14"
filein1=sprintf("%-30.10f",Time)
filein =sprintf("%s%.10sx-correlations.dat",dirIN,filein1)
fileout1=sprintf("%s%.10sx-correlations.jpeg",dirOUT,filein1)
print "I am in JPG PIC=",PIC
}

if(PIC==2){
#set term gif enh animate size 1600,1200 font "arial,30"
set term gif enh animate size 1600,1200 font "arial,14"
set term gif enh animate size 800,600 font "arial,14"
filein1=sprintf("%-30.10f",Tmin) 
filein =sprintf("%s%.10sx-correlations.dat",dirIN,filein1)
fileout1="./media/1D_g1RR.gif"
print "I am in GIF PIC=",PIC
}
set out fileout1
print "File namein: ", filein, " out: ", fileout1
if (!exists("filein")) print "Input File: ",filein," does not exist!!!" 
if (!exists("filein")) exit 

set grid
set nobar 
set key  top right spacing 1.5
dnsttl=sprintf("DNS MCTDHB(%i)",Morb)
orbttl(i)=sprintf("|^{WK}PSI_%i|^2",i)

set xlabel "x" 
set ylabel "x'" 
zero=1.0e-16

set size 1.0,1.0;
set pm3d map;
set pm3d interpolate 2,2 flush begin ftriangles nohidden3d corners2color mean
a= 9
minX=-a
maxX= a
minY=-a
maxY= a
cutOFF=0.0000000001

set xrange [ minX : maxX ] noreverse nowriteback
set yrange [ minY : maxY ] noreverse nowriteback
#set cbrange [ 0.0 : 1.01] noreverse nowriteback


do for [j=1:n]{
Fulltitle=sprintf(" MCTDHB(%i):  |g^{(1)}(x\',x|t)|^2  at T=%2.3f",Morb,Time+dt*(j-1))
set title Fulltitle 
filein1=sprintf("%-30.10f",Tmin+(j-1)*dt) 
mvfilein(j)=sprintf("%s%.10sx-correlations.dat",dirIN,filein1)
print "File namein: ", mvfilein(j)

splot [minX:maxX][minY:maxY][0.0:1.01] mvfilein(j) using \
1:4:(($7**2)<cutOFF || ($10**2)<cutOFF   ? 1/0 :(($8**2+$9**2)/($7*$10)>1 ? 1 : ($8**2+$9**2)/($7*$10))) t ""
}

print "File name  out: ", fileout1 


