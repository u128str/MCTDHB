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
set grid
set nobar 
set key  top right spacing 1.5

set xlabel "x" font "arial,30"
set xlabel "y" font "arial,30"
set zlabel "V(x,t),rho" font "arial,30" 

zero=1.0e-16
#============== PSI(i) ======================
#fileout=sprintf("%10.8ftimeSclNatOrb.jpeg",Time)
#filein=sprintf("%s%10.8ftime.dat",dir,Time)
filein1=sprintf("%-30.10f",Time)
filein=sprintf("%s%.10stime.dat",dirIN,filein1)
fileout=sprintf("%s%.10stimeDNS2D.jpeg",dirOUT,filein1)
print "File namein: ", filein, " out: ", fileout



dnsttl=sprintf("DNS MCTDHB(%i)",Morb)
orbttl(i)=sprintf("Natural n(%i)|PSI_%i|^2",Morb-i+1,Morb-i+1)

nocc(i)=4*Morb+10+i-1

Fulltitle=sprintf("MCTDHB(%i) Density in 2D at T=%f",Morb,Time)
set title Fulltitle 
#fileout=sprintf("%s",dir)
print "File name : ",fileout
set out fileout 


set pm3d implicit at b
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
unset colorbox 
set loadpath 
set fontpath 
set fit noerrorvariables



X1=-10
X2= 10
Y1=-10
Y2= 10
Z1=-0.1
Z2=Morb*0.75+1
splot  [X1:X2] [Y1:Y2] \
filein using ($1):2:((0.02*$5-0.01)>0.5  ? 1/0 :(0.02*$5-0.01)) t "V(x,t)" w l ,\
filein using ($1):2:(($6<0.0000001 ? 1/0 :$6)) t dnsttl   w l

#filein using ($1):(0.02*$5-0.01) t "V(x,t)" w filledcurves y1=-0.01 lc rgb "black" fs transparent solid 0.5 noborder,\
#for[i=1:Morb] "" using 1: (column(nocc(i))*(column(2*Morb+2*(i-1)+10)**2+column(2*Morb+2*(i-1)+11)**2)/(column(4)*column(4)) +(Morb-i+1)*0.75) \
#t orbttl(i) w filledcurves y1=(Morb-i+1)*0.75 lt Morb-i+1  fs trans fs transparent solid 0.8 noborder,\
#"" using ($1):($8) t dnsttl w filledcurves y1=0.0 lc rgb "#00008B" fs transparent solid 0.4 noborder 
