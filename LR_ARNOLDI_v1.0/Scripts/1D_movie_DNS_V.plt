#    G N U P L O T
#to use for different Morb you need to execute: gnuplot -e "Morb=2" -e "Time=1.2" thisScript
print "to use for different Morb you need to execute: gnuplot -e \"Morb=2\"  thisScript"
print " Full usage: gnuplot -e \"Time_Bgn=0; Time_Fnl=15; dt=0.1; Morb=4\"  $mctdhb_dir/Scripts/E_of_t.plt"
print "where: Time_Bgn=0 a time point to start from"
print "where: Time_Fnl=0 a time point where to finish"
print "where: dt=0  time step"
print "where: Morb=4 number of MCTDHB orbitals used"
print "where: $mctdhb_dir is the installation directory of the MCTDHB_V3* package"

if (!exists("Morb")) Morb=2
if (!exists("Time")) Time=0.0 
print "Morb: ", Morb, " at t= ",Time
if (!exists("Time_Bgn")) Time_Bgn=0
if (!exists("Time_Fnl")) Time_Fnl=10
if (!exists("dt")) dt=0.1
Tmin= Time_Bgn
Tmax= Time_Fnl
dt=dt
n=1+(Tmax-Tmin)/dt #n frames
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
#set terminal jpeg  size 1600,1200  font "arial,24"
#set terminal jpeg  enh size 1600,1200  font "arial,24"
#set terminal jpeg enh  size 1600,1200  font "arial,24"
set grid
set nobar 
set key  top right spacing 1.5

set xlabel "x" font "arial,30"
set ylabel "V(x,t),rho(x,t), Natural Orbitals" font "arial,30" 

zero=1.0e-16
#============== PSI(i) ======================
#filein1=sprintf("%-30.10f",Time)
#filein=sprintf("%s%.10stime.dat",dirIN,filein1)
#fileout=sprintf("%s%.10stimeNatOrbReIm.jpeg",dirOUT,filein1)
#
#print "File namein: ", filein, " out: ", fileout
#if (!exists("filein")) print "Input File: ",filein," does not exist!!!" 
#if (!exists("filein")) exit 

dnsttl=sprintf("DNS MCTDHB(%i)",Morb)
orbttl(i)=sprintf("Natural |PSI_%i|^2",Morb-i+1)
orbttl(i)=sprintf("Natural |PSI_%i|^2",Morb-i+1)
#41 echo "\"$i\" using (\$1):("$x1"+    ((\$"$i1"))/(\$4)) t \"NO n("$i2")|Re (PSI_"$i2")\" w lp lt "$i2" pt 6 ps 0.1 ,\\">>$i.gnu
#142 echo "\"$i\" using (\$1):("$x1+0.1"+((\$"$i0"))/(\$4)) t \"NO n("$i2")|Im (PSI_"$i2")\" w lp lt "$i2" pt 5 ps 0.5 ,\\">>$i.gnu


#Fulltitle=sprintf("Re and Im parts of natural MCTDHB(%i) Orbitals at T=%f",Morb,Time)
#set title Fulltitle 
#fileout=sprintf("%s",dir)
#print "File name : ",fileout
#set out fileout .
X1=-10
X2= 10
Y1=-0.1
Y2=Morb*0.75+1
#plot  [X1:X2] [Y1:Y2] \
#filein using ($1):(0.02*$5-0.01) t "V(x,t)" w filledcurves y1=-0.01 lc rgb "black" fs transparent solid 0.5 noborder,\
#for[i=1:Morb] \
#"" using 1: (column(2*Morb+2*(i-1)+11)/column(4) +(Morb-i+1)*0.75+0.1) \
#t "Im[^{NO}psi_".(Morb-i+1)."]" w filledcurves y1=(Morb-i+1)*0.75+0.1 lt 2*Morb-i+1+1  fs trans fs transparent solid 0.5 noborder,\
#for[i=1:Morb] \
#"" using 1: (column(2*Morb+2*(i-1)+10)/column(4) +(Morb-i+1)*0.75) \
##t "Re[^{NO}psi_".(Morb-i+1)."]" w filledcurves y1=(Morb-i+1)*0.75 lt Morb-i+1  fs trans fs transparent solid 0.5 noborder,\
#"" using ($1):($8) t dnsttl w filledcurves y1=0.0 lc rgb "#00008B" fs transparent solid 0.4 noborder 

#n=50 #n frames
#do for[i=1:n]{
#filein1=sprintf("%-30.10f",Time+i*dt) 
#mvfilein(i)=sprintf("%s%.10stime.dat",dirIN,filein1)
#print "File namein: ", mvfilein(i)
#}

#reset
set term gif enh animate size 1600,1200 font "arial,30"
set output "animate.gif"
set xrange [X1:X2]
set yrange [Y1:Y2]
do for [i=1:n]{
Fulltitle=sprintf("Re and Im parts of natural MCTDHB(%i) Orbitals at T=%2.3f",Morb,Tmin+dt*(i-1))
set title Fulltitle 
filein1=sprintf("%-30.10f",Tmin+(i-1)*dt) 
mvfilein(i)=sprintf("%s%.10stime.dat",dirIN,filein1)
print "File namein: ", mvfilein(i)
plot mvfilein(i)  using ($1):($8) t sprintf("rho(t=%2.3f)",Tmin+dt*(i-1)) w filledcurves y1=0.0 lc rgb "#00008B" fs transparent solid 0.4 noborder,\
"" using ($1):(0.02*$5-0.01) t "V(x,t)" w filledcurves y1=-0.01 lc rgb "black" fs transparent solid 0.5 noborder
}
