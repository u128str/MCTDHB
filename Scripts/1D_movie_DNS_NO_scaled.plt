#!/usr/bin/gnuplot 
#    G N U P L O T
# thisscript 1D_movie_DNS_Re_Im_NO_scaled.plt
print "GIF usage: gnuplot  -e \"Morb=2; Time_Bgn=0.0; Time_Fnl=11; dt=0.2 \"  thisscript"
print "OUT animated gif: ./media/1D_DNS_NO_scaled.gif"
print "where: Time_Bgn=0 a time point to start from"
print "where: Time_Fnl=0 a time point where to finish"
print "where: dt=0  time step"
print "where: Morb=4 number of MCTDHB orbitals used"
print "where: $mctdhb_dir is the installation directory of the MCTDHB_V3* package"
print "JPG usage: gnuplot  -e \"Morb=2; Time_Bgn=10.0 \"  $mctdhb_dir/Scripts/thiscript"
print "OUT single-shot jpeg: ./media/10.0000000timeSclNatOrb.jpeg"

if (!exists("Morb")) Morb=2
if (!exists("Time")) Time=0.0 
print "Morb: ", Morb, " at t= ",Time
if (!exists("Time_Bgn")) Time_Bgn=0
if (!exists("Time_Fnl")) Time_Fnl=10
PIC=2
if (!exists("dt")) {dt=0;Time_Fnl=Time_Bgn;PIC=1;Time=Time_Bgn}
Tmin= Time_Bgn
Tmax= Time_Fnl
dt=dt
n=1
if (dt!=0) n=1+(Tmax-Tmin)/dt #n frames
print "PIC=",PIC," n=",n, " dt=", dt," Time_Bgn=",Time_Bgn," Time_Fnl=",Time_Fnl
#                            DATA/orb_R/1.00000000time.dat
dirIN="DATA/orb_R/"
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
set terminal jpeg enh  size 800,600  font "arial,14"
filein1=sprintf("%-30.10f",Time)
filein=sprintf("%s%.10stime.dat",dirIN,filein1)
fileout1=sprintf("%s%.10stimeSclNatOrb.jpeg",dirOUT,filein1)
print "I am in JPG PIC=",PIC
}

if(PIC==2){
#set term gif enh animate size 1600,1200 font "arial,30"
set term gif enh animate size 800,600 font "arial,14"
filein1=sprintf("%-30.10f",Tmin) 
filein=sprintf("%s%.10stime.dat",dirIN,filein1)
fileout1="./media/1D_DNS_NO_scaled.gif"
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
orbttl(i)=sprintf("n(%i)|^{NO}PSI_%i|^2",Morb-i+1,Morb-i+1)
nocc(i)=4*Morb+10+i-1
set xlabel "x" 
set ylabel "V(x,t),rho(x,t),Natural Orbitals" 
X1=-10
X2= 10
Y1=-0.1
Y2=Morb*0.75+1
zero=1.0e-16
set xrange [X1:X2]
set yrange [Y1:Y2]
do for [j=1:n]{
Fulltitle=sprintf("Scaled Natural MCTDHB(%i) Orbitals at T=%2.3f",Morb,Time+dt*(j-1))
set title Fulltitle 
filein1=sprintf("%-30.10f",Tmin+(j-1)*dt) 
mvfilein(j)=sprintf("%s%.10stime.dat",dirIN,filein1)
print "File namein: ", mvfilein(j)
plot [X1:X2] [Y1:Y2] \
mvfilein(j)  using ($1):($8) t sprintf("rho(t=%2.3f)",Tmin+dt*(j-1)) \
 w filledcurves y1=0.0 lc rgb "#00008B" fs transparent solid 0.4 noborder,\
"" using ($1):(0.02*$5-0.01) t "V(x,t)" w filledcurves y1=-0.01 lc rgb "black" fs transparent solid 0.5 noborder,\
for[i=1:Morb] "" using 1: (column(nocc(i))*(column(2*Morb+2*(i-1)+10)**2+column(2*Morb+2*(i-1)+11)**2)/(column(4)*column(4)) +(Morb-i+1)*0.75) \
t orbttl(i) w filledcurves y1=(Morb-i+1)*0.75 lt Morb-i+1  fs trans fs transparent solid 0.8 noborder
}

print "File name  out: ", fileout1 


