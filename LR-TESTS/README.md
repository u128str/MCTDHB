# Arnoldi LR test LR-GP=LR-MCTDHB(1)=BdG results.

Prerequisites: Installed Arnoldi-LR package: https://github.com/u128str/MCTDHB/tree/MCTDHB_V3.3.02

In the dir where all these files are:
1) Run the mctdhb to get GS:
```
$  ../bin/boson_MCTDHB_ifort_MKLFFT 
...
====================================================================================================
 Job->Relax. Forward     Iteration:     200     Time: [      0.000000 ->     19.900000 +       0.050000 ->     20.000000  ]
       Input orbital energy E(t+0):        9.1231859529219328( CI Dim:         1)(ORB Dim:  1*       256=       256)
          OUT  CI  energy E(t+tau):        9.1231859529219328     N =         10   l0*(N-1)= 0.9000000000     kind of W(x-x'):0
                  Delta E     : +/-        0.0000000000000053                  Error due to  dE:        0.0000000000000000
                  Error E     : +/-        0.0000000010000000
      New NO's:   10.00000     | 
 !!!====== Last Point arroaching Tmax=   20.000000000000000     
 !!!====== Time-Step is reduced to   =  -1.4921397450962104E-013
      Itr Time:    0.1553230286     execution time:   31.3989896774
 # This computation has been done in
 # /mnt/worka/alexej/MyGits/MCTDHB-LAB-BUILDER/LR_ARNOLDI_v1.0/Examples_LR_ARNOLDI/1Dmy
 # Date 26/10/2017; time 09:26:43
 #    Current version 3.3.01 l Heidelberg/Kassel Apr (2017)      # 
 # Morb=  1  Npar=      10   Job= Job->Relax. Forward
 # V(x_y_z&t)= x^2                                                                                                 
 # W(R=|r1-r2|&t)= Using Defaults from Get_InterParticle.F                                                             
......
```
 1.1) Compare the outputs with the reference:
```
$ vimdiff  basic_info.out basic_info.out_GP_ref
```


2) Run the OLD lr-mctdhb:
```
$  ../bin/properties_LR_ifort_MKLFFT
```
The results are in DATA/getLR/MC_anlsplot.out:
```
The results are in DATA/getLR/MC_anlsplot.out:
...
#      i,                 Re(w_i),              Re(w_i+Eref),         Re(Norm_Orb),     Re(Norm_CI),   Abs(<PSI_ref|OP1(x)|LR_PSI_i>),  Abs( <PSI_ref|OP2(x)|LR_PSI_i>)
      11       -0.0000000000001121     0.9123185952921821E+01   0.9999999999999999   0.0000000000000000   0.0000000000002047   1.9303668276487767
      13        1.4142135623733161     0.1053739951529525E+02   1.0000000000000000   0.0000000000000000   1.8803015465431998   0.0000000000001247
      14        2.7404920280379570     0.1186367798095989E+02   1.0000000000000004   0.0000000000000000   0.0000000000000975   1.7228203003721720
      15        4.1054294604204893     0.1322861541334242E+02   0.9999999999999994   0.0000000000000000   0.0000000000001021   0.0000000000002507
      16        5.4874832086801710     0.1461066916160210E+02   0.9999999999999999   0.0000000000000000   0.0000000000001161   0.0156027446226764
      17        6.8787705465924178     0.1600195649951435E+02   1.0000000000000000   0.0000000000000000   0.0000000000000352   0.0000000000000472
      18        8.2756260348521788     0.1739881198777411E+02   1.0000000000000009   0.0000000000000000   0.0000000000000257   0.0022610795094706
      19        9.6761281452497272     0.1879931409817166E+02   1.0000000000000004   0.0000000000000000   0.0000000000000122   0.0000000000000192
.....
```
 2.1) Compare the outputs with the reference:
```
vimdiff DATA/getLR/MC_anlsplot.out  MC_anlsplot.out_GP_ref
```

3) Run NEW Arnoldi lr-mctdhb (1.stage MATRIX CONSTRUCTION in __./lr_arnoldi.in__ ``` task=1  ```):
```
$  ../bin/lr_arnoldi_ifort_MKLFFT
```
On terminal you shoulod see:
```
.....
 Upper part constructed in   0.44800000000000001       seconds
 Lower part constructed in   0.44800000000000001       seconds

 LR-MATRIX constructed in  0.69600000000000006      seconds
 It is stored in 'LR_matrix_save.dat'.
...
```

 3.1) Before Runing NEW Arnoldi lr-mctdhb (2.stage MATRIX DIAGONALIZATION ) in __./lr_arnoldi.in__ change the ``` task=1``` to ``` task=2```. Run it:

```
$  ../bin/lr_arnoldi_ifort_MKLFFT
```
On terminal you will see:
```
...
    14             0.8275626035E+001             0.4367763532E-010             0.1346123106E-012
    13            -0.8275626035E+001             0.4370030290E-010             0.4812585667E-013
    12             0.6878770547E+001            -0.1245087281E-012             0.4692088626E-013
    11            -0.6878770547E+001            -0.4673262636E-013             0.2613882852E-013
    10             0.5487483209E+001            -0.5968703257E-010             0.5022654220E-013
     9            -0.5487483209E+001            -0.5971091698E-010             0.8920408817E-013
     8             0.4105429460E+001            -0.7261531083E-012             0.2413928135E-013
     7            -0.4105429460E+001            -0.6507630979E-012             0.2915565890E-012
     6             0.2740492028E+001            -0.6582935154E-011             0.2322927345E-013
     5            -0.2740492028E+001            -0.6505667665E-011             0.1517666073E-012
     4             0.1414213562E+001            -0.9232044955E-012             0.4236521056E-014
     3            -0.1414213562E+001            -0.9060428678E-012             0.1182636316E-012
     2             0.1160684751E-012             0.9144456217E-010             0.2313870903E-010
     1             0.2008573571E-012            -0.8978298577E-010             0.2489769448E-010

 Arnoldi iteration:         135
 ZNAUPD info:           0
 Converged values:           5
 ido:           1
 Time elapsed:   10.464000000000000      seconds
 =======================================
 %-----------------------------------------------%
           5 Ritz values converged
 starting now eigenvector/residual calculation...
 %-----------------------------------------------%
  
 Eigenvectors and residuals calculated.
   
   
# root          Ritz Values (Real)          Ritz Values (Imag)          Residuals
   
     1             0.1160684751E-012             0.9144456217E-010             0.2200436834E-010
     2             0.2008573571E-012            -0.8978298577E-010             0.2358849618E-010
     3             0.1414213562E+001            -0.9232044955E-012             0.6856810816E-011
     4             0.2740492028E+001            -0.6582935154E-011             0.8797614759E-011
     5             0.4105429460E+001            -0.7261531083E-012             0.7361434058E-011
   
   
  Dimension of Krylov subspace           80
  Tolerance    1.0000000133514320E-010
  Mat*vec operations       10800
  Arnoldi restarts         135
  Converged Ritz values           5
  Eigenvectors stored in subfolder "Eigenv/"
  Computation time    10.484000000000000      seconds

```

this reference screen can be find in __./arnoldi_LR.out_ref_GP__


