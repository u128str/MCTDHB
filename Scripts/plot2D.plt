#!/usr/bin/gnuplot -persist
print "usage:  gnuplot -e \"fn='$i'\" -e \"fn1='$ii'\" plot2D.plt "
set terminal jpeg medium size 800,600;
set output fn.'.jpg'
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set title 'time='.fn1
set title  offset character 0, 0, 0 font "" norotate
set pm3d implicit at b
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
unset colorbox 
#set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set loadpath 
set fontpath 
set fit noerrorvariables
splot  [-8:8][-6:6][-0.:0.6]  fn  u 1:2:(($6<0.0000001 ? 1/0 :$6)) t "" w l 
#splot  [-8:8][-6:6][-0.:0.6]  "./0.00000000time.dat"  u 1:2:(($6<0.0000001 ? 1/0 :$6)) t "$ii" w l
#    EOF
