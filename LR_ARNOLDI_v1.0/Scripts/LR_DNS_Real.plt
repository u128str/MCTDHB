#!/usr/bin/gnuplot 
print "to use for different Morb you need to execute: gnuplot -e \" N=100; r1=206; r2=207; r3=208; r4=209; r5=210\"  thisScript"
print " Full usage: gnuplot -e \"r1=206; r2=207; r3=208; r4=209; r5=210\"   $mctdhb_dir/Scripts/LR_DNS_Real.plt"
print "where: N= number of particles"
print "     : r1,r2,..r5 LR-roots of interest"
print "     : $mctdhb_dir is the installation directory of the MCTDHB_V3* package"

dirIN="DATA/getLR/"
dirOUT="./media/"

filein=sprintf("%sMC_dens_x.out",dirIN)
#fileout=sprintf("%sLR_DNS_Real.eps",dirOUT) #to see pictures better 
fileout=sprintf("%sLR_DNS_Real.jpeg",dirOUT)
set out fileout
set terminal jpeg enh size 1600,1200  font "arial,24"
set yrange [-0.5:0.5]
set xrange [-4:4]

if (!exists("r1")) r1=218
if (!exists("r2")) r2=219
if (!exists("r3")) r3=220
if (!exists("r4")) r4=221
if (!exists("r5")) r5=222

Npar=10 #Number of Particles 
#scaling factors to see the densities  on the same plot
sc1=1
sc2=5
sc3=5
sc4=5
sc5=1
#phase= Gamma/Abs|Gamma| -- diagonalizers give different phase between real and IM parts of the responce weights (Gamma) so 
# to plot propertly one has to take the respective Gamma {Re,Im} from MC_we.out file!!!!!
#MC_we.out
#       WRITE(198,'(I8,E26.16,1X,F26.16,1X,F26.16,1X,2E16.8,2E16.8,2E16.8,2E16.8,2E16.8)') &
#     & i, abs(wo(i))+Energy,abs(weight_orb(i,1)+weight_CI(i,1)),abs(weight_orb(i,2)+weight_CI(i,2)) &
#     & ,weight_orb(i,1),weight_orb(i,2),weight_CI(i,1),weight_CI(i,2),wo(i)                                    ! 198=MC_we.out
# i Energy, ABS(Gamma_CI+Gamma_Orb))(f1) ABS(Gamma_CI+Gamma_Orb))(f2), ....
#  g=weight_orb(i,1)+weight_CI(i,1) -- complex number !!!!!! because weight_orb(i,1) has Re and Im parts
g1 = {1.0,0.0}
g2 = {1.0,0.0}
g3 = {1.0,0.0}
g4 = {1.0,0.0}
g5 = {1.0,0.0}

#================================== Gnu -stuff
mpl_top    = 0.4 #inch  outer top margin, title goes here
mpl_bot    = 0.7 #inch  outer bottom margin, x label goes here
mpl_left   = 0.9 #inch  outer left margin, y label goes here
mpl_right  = 0.1 #inch  outer right margin, y2 label goes here
mpl_height = 3.0 #inch  height of individual plots
mpl_width  = 3.0 #inch  width of individual plots
mpl_dx     = 0.1 #inch  inter-plot horizontal spacing
mpl_dy     = 0.1 #inch  inter-plot vertical spacing
mpl_ny     = 3   #number of rows
mpl_nx     = 2   #number of columns

# calculate full dimensions
xsize = mpl_left+mpl_right+(mpl_width*mpl_nx)+(mpl_nx-1)*mpl_dx
ysize = mpl_top+mpl_bot+(mpl_ny*mpl_height)+(mpl_ny-1)*mpl_dy

# placement functions
#   rows are numbered from bottom to top
bot(n) = (mpl_bot+(n-1)*mpl_height+(n-1)*mpl_dy)/ysize
top(n)  = 1-((mpl_top+(mpl_ny-n)*(mpl_height+mpl_dy))/ysize)
#   columns are numbered from left to right
left(n) = (mpl_left+(n-1)*mpl_width+(n-1)*mpl_dx)/xsize
right(n)  = 1-((mpl_right+(mpl_nx-n)*(mpl_width+mpl_dx))/xsize)

#set terminal postscript eps enhanced color dl 2.0 size xsize,ysize "Helvetica" 28
#set terminal pdf enhanced color dl 2.0 size xsize,ysize  
set encoding iso_8859_1
set tics scale 1.5

#set output 'LR_DNS_Real.eps'

set offsets
set autoscale fix
set size 1,1
set nokey

# define x-axis settings for all subplots
#set xrange [-4:4]
set xlabel ''
set format x ''
set xtics 4
set mxtics 4

#set yrange [-0.14:0.14]
# start plotting
set multiplot

#-----------------------------------------------
# subplot  1-3
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for third row (top)
set tmargin at screen top(3)
set bmargin at screen bot(3)

set title 'left'
set title ''

#set ylabel "{/Symbol Dr/@\326\140}N" offset 0.5
#set ylabel "Real(DNS)/sqrt(N)" offset 0.5
#set yrange [-1.5:1.5]
set format y "%-2.1f"
set ytics mirror 0.1
set mytics 2

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 1 "state N".r1 at -3,0.1
set label 11 "x".sc1 at 2,0.05
g=g1
plot          \
filein u 1:($2==r1 ? (sc1*real((($3+$5)+{0,1}*($4+$6))*g/abs(g)/$10/$10/sqrt(Npar))):1/0)\
axes x1y1 \
title '' \
with lines lt 1 lc 1 lw 10\
;


#-----------------------------------------------
# subplot  1-2
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for second row (middle)
set tmargin at screen top(2)
set bmargin at screen bot(2)

set title ''

#set ylabel "{/Symbol Dr/@\326\140}N" offset 0.5
set ylabel "Real(DNS)/sqrt(N)" offset 0.5
#set yrange [-1.5:1.5]
set format y "%-1.1f"
set ytics mirror 0.1
set mytics 2

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 2 "state N".r2 at -3,0.1
set label 21 "x".sc2 at 2,0.05
g=g2
plot          \
filein u 1:($2==r2 ? (-sc2*real((($3+$5)+{0,1}*($4+$6))*g/abs(g)/$10/$10/sqrt(Npar))):1/0)\
axes x1y1 \
title '' \
with lines lt 1 lc 2 lw 10\
;

#-----------------------------------------------
# subplot  2-2
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for second row (middle)
set tmargin at screen top(2)
set bmargin at screen bot(2)

set title ''

set ylabel ""             # no label here
#set yrange [-1.5:1.5]
set format y ""           # no tic labels
set ytics mirror 0.1
set mytics 2

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 3 "state N".r3 at -3,0.1
set label 31 "x".sc3 at 2,0.05
g=g3
plot          \
filein u 1:($2==r3 ? (-sc3*real((($3+$5)+{0,1}*($4+$6))*g/abs(g)/$10/$10/sqrt(Npar))):1/0)\
title '' \
with lines lt 1 lc 3 lw 10\
;

#-----------------------------------------------
# subplot  1-1
#  set horizontal margins for first column
set lmargin at screen left(1)
set rmargin at screen right(1)
#  set horizontal margins for first row (bottom)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set title ''

# now set a label and tic marks for the x-axis
set xlabel "x"
#set xtics add ("-{/Symbol p}" -pi, "0" 0, "{/Symbol p}" pi)
set xtics add ("-3" -3, "0" 0, "3" 3)

#set ylabel "{/Symbol Dr/@\326\140}N" offset 0.5
#set yrange [-1.5:1.5]
set format y "%-1.1f"
set ytics mirror 0.1
set mytics 2

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 4 "state N".r4 at -3,0.1
set label 41 "x".sc4 at 2,0.05
g=g4
plot          \
filein u 1:($2==r4 ? (sc4*real((($3+$5)+{0,1}*($4+$6))*g/abs(g)/$10/$10/sqrt(Npar))):1/0)\
title '' \
with lines lt 1 lc 4 lw 10\
;

#-----------------------------------------------
# subplot  2-1
#  set horizontal margins for second column
set lmargin at screen left(2)
set rmargin at screen right(2)
#  set horizontal margins for first row (bottom)
set tmargin at screen top(1)
set bmargin at screen bot(1)

set title ''

set ylabel ""             # no label here
#set yrange [-1.5:1.5]
set format y ""           # no tic labels
set ytics mirror 0.1
set mytics 2

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 5 "state N".r5 at -3,0.1
set label 6 "x".sc5 at 2,0.05
g=g5
plot          \
filein u 1:($2==r5 ? (sc5*real((($3+$5)+{0,1}*($4+$6))*g/abs(g)/$10/$10/sqrt(Npar))):1/0)\
axes x1y1 \
title '' \
with lines lt 1 lc 5 lw 10\
;


unset multiplot
print "Real parts of the Response Densities of the LR excitations are plotted"
print "File name  in: ", filein 
print "File name out: ", fileout
