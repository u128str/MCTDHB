#    G N U P L O T
# thisscript 1D_movie_DNS_Re_Im_WK_3Dview.plt
print "GIF usage: gnuplot  -e \"Morb=2; Time_Bgn=0.0; Time_Fnl=11; dt=0.2 \"  thisscript"
print "OUT animated gif: ./media/1D_DNS_Re_Im_WK_3Dview.gif"
print "where: Time_Bgn=0 a time point to start from"
print "where: Time_Fnl=0 a time point where to finish"
print "where: dt=0  time step"
print "where: Morb=4 number of MCTDHB orbitals used"
print "where: $mctdhb_dir is the installation directory of the MCTDHB_V3* package"
print "JPG usage: gnuplot  -e \"Morb=2; Time_Bgn=0.0 \"  $mctdhb_dir/Scripts/thiscript"
print "OUT single-shot jpeg: ./media/0.00000000timeWrkOrbReIm_3Dview.jpeg"

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
set terminal jpeg enh  size 1600,1200  font "arial,24"
filein1=sprintf("%-30.10f",Time)
filein=sprintf("%s%.10stime.dat",dirIN,filein1)
fileout1=sprintf("%s%.10stimeWrkOrbReIm_3Dview.jpeg",dirOUT,filein1)
print "I am in JPG PIC=",PIC
}

if(PIC==2){
set term gif enh animate size 1600,1200 font "arial,30"
filein1=sprintf("%-30.10f",Tmin) 
filein=sprintf("%s%.10stime.dat",dirIN,filein1)
fileout1="./media/1D_DNS_Re_Im_WK_3Dview.gif"
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
X1=-10
X2= 10
Y1=-0.75
Y2=0.75
Z1=-0.75
Z2=0.75
set view 120,95,1.2,1.2
set xlabel "Im" offset 0, 0, -3 
set ylabel "X"  offset 8, -1, 0   
set zlabel "Re" offset 0,0.5, 0   
set view 130,20,1.1,1.2
set xlabel "Re" offset 0, 0, -3 
set ylabel "X"  offset 8, -1, 0   
set zlabel "Im" offset 0,0.5, 0   
scl=0.0
zero=1.0e-16
do for [j=1:n]{
Fulltitle=sprintf("Re and Im parts of working MCTDHB(%i) Orbitals at T=%2.3f",Morb,Tmin+dt*(j-1))
set title Fulltitle 
filein1=sprintf("%-30.10f",Tmin+(j-1)*dt) 
mvfilein(j)=sprintf("%s%.10stime.dat",dirIN,filein1)
print "File namein: ", mvfilein(j)
splot  [Z1:Z2] [X1:X2] [Y1:Y2] \
mvfilein(j) using (0.02*$5-0.01):1:(0.0) t "V(x,t)" w lp lt -1 pt 7 ps 0.0 lw 4 ,\
        ""  using ($8):1:(0.0)  t dnsttl w lp lt 4 lw 3 pt 6 ps 0.7,\
for[i=1:Morb] \
"" using (column(2*i+8)/column(4)-i*scl):1:(column(2*i+9)/column(4)) \
t "[^{WK}psi_".(i)."]" w lp 
}

print "File name  out: ", fileout1 


