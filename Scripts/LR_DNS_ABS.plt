#!/usr/bin/gnuplot 
print "to use for different Morb you need to execute: "
print "gnuplot -e \"Xmin=-4; Xmax=4; Morb= 4; Npar=100; r1=206; r2=207; r3=208; r4=209; r5=210; s1=1; s2=1; s3=1; s4=1; s5=1 \"  $mctdhb_dir/Scripts/LR_DNS_Real.plt"
print "where: N= number of particles, Morb= number of orbitals"
print "     : Xmin,Xmax plot region "
print "     : r1,r2,..r5 LR-roots of interest"
print "     : s1,s2,..s5 prefactors to scale the plotted densitis"
print "     : if s_i are not specified all s_i=1"
print "     : if r_i are not specified  r_i=10+2*Morb+1*i"
print "     : first 10+Morb are from -w manifold, next Morb are zeros and the rest ones are from +w manifold"

dirIN="DATA/getLR/"
dirOUT="./media/"

filein=sprintf("%sMC_dens_x.out",dirIN)
fileout=sprintf("%sLR_DNS_ABS.jpeg",dirOUT)
set out fileout
#set terminal jpeg enh size 1600,1200  font "arial,24"
set terminal jpeg enh size 800,600  font "arial,14"
lwv=1
set yrange [-0.0:0.5]
if (!exists("Xmin")) Xmin=-4
if (!exists("Xmax")) Xmax=+4
set xrange [Xmin:Xmax]

if (!exists("Morb"))  {print "Morb must be specified!"; exit gnuplot }
if (!exists("Npar"))  {print "Npar must be specified!"; exit gnuplot }
if (!exists("r1")) r1=10+2*Morb+1
if (!exists("r2")) r2=10+2*Morb+2
if (!exists("r3")) r3=10+2*Morb+3
if (!exists("r4")) r4=10+2*Morb+4
if (!exists("r5")) r5=10+2*Morb+5

print "Roots_i:",r1,r2,r3,r4,r5
if (!exists("s1")) s1=1
if (!exists("s2")) s2=1
if (!exists("s3")) s3=1
if (!exists("s4")) s4=1
if (!exists("s5")) s5=1
#Npar=10 #Number of Particles 
#scaling factors to see the densities  on the same plot
sc1=s1
sc2=s2
sc3=s3
sc4=s4
sc5=s5
print "scalings Scl_i:",sc1,sc2,sc3,sc4,sc5
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

set ytics -0.6,0.2,0.6
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
#set ylabel "Abs(DNS)/sqrt(N)" offset 0.5
#set yrange [-1.5:1.5]
set format y "%-2.1f"
set ytics mirror 0.1
set mytics 2
set ytics -0.6,0.2,0.6

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 2 "Numeration of states as" at 5,0.25
set label 3 "in DATA/getLR/MC\\_anlsplot.out" at 5,0.15
set label 1 "state N".r1 at -3,0.1
#set label 11 "x".sc1 at 2,0.05
s=sprintf("x%3.2f",sc1)
set label 11 s at 2,0.05
plot  [Xmin:Xmax]         \
filein u 1:($2==r1 ? (sc1*abs((($3+$5)+{0,1}*($4+$6))/$10/$10/sqrt(Npar))):1/0)\
axes x1y1 \
title '' \
with lines lt 1 lc 1 lw lwv+3\
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
set ylabel "DNS/sqrt(N)" offset 0.5
#set yrange [-1.5:1.5]
set format y "%-1.1f"
set ytics mirror 0.1
set mytics 2
set ytics -0.6,0.2,0.6

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 2 "state N".r2 at -3,0.1
#set label 21 "x".sc2 at 2,0.05
s=sprintf("x%3.2f",sc2)
set label 21 s at 2,0.05
plot  [Xmin:Xmax]         \
filein u 1:($2==r2 ? (sc2*abs((($3+$5)+{0,1}*($4+$6))/$10/$10/sqrt(Npar))):1/0)\
axes x1y1 \
title '' \
with lines lt 1 lc 2 lw lwv+3\
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
set ytics -0.6,0.2,0.6

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 3 "state N".r3 at -3,0.1
#set label 31 "x".sc3 at 2,0.05
s=sprintf("x%3.2f",sc3)
set label 31 s at 2,0.05
plot  [Xmin:Xmax]         \
filein u 1:($2==r3 ? (sc3*abs((($3+$5)+{0,1}*($4+$6))/$10/$10/sqrt(Npar))):1/0)\
title '' \
with lines lt 1 lc 3 lw lwv+3\
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
set ytics -0.6,0.2,0.6

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 4 "state N".r4 at -3,0.1
#set label 41 "x".sc4 at 2,0.05
s=sprintf("x%3.2f",sc4)
set label 41 s at 2,0.05
plot  [Xmin:Xmax]         \
filein u 1:($2==r4 ? (sc4*abs((($3+$5)+{0,1}*($4+$6))/$10/$10/sqrt(Npar))):1/0)\
title '' \
with lines lt 1 lc 4 lw lwv+3\
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
set ytics -0.6,0.2,0.6

set arrow 1 from graph 0, first 0 rto graph 1,0 nohead lt 1 lw 1 lc 0
set arrow 2 from first 0, graph 0 rto 0, graph 1 nohead lt 1 lw 1 lc 0

unset label
set label 5 "state N".r5 at -3,0.1
#set label 6 "x".sc5 at 2,0.05
s=sprintf("x%3.2f",sc5)
set label 6 s at 2,0.05
plot  [Xmin:Xmax]         \
filein u 1:($2==r5 ? (sc5*abs((($3+$5)+{0,1}*($4+$6))/$10/$10/sqrt(Npar))):1/0)\
axes x1y1 \
title '' \
with lines lt 1 lc 5 lw lwv+3\
;


unset multiplot
print "Abs of the Response Densities of the LR excitations are plotted"
print "File name  in: ", filein 
print "File name out: ", fileout
