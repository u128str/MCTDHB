# Arnoldi LR MCTDHB template

## A. How2 get/install/test etc 
[Arnoldi LR MCTDHB](https://github.com/u128str/MCTDHB/blob/master/LR-ARNOLDI.md)

## B.  Arnoldi LR MCTDHB run - __How2 verify/test__
Next stage is to verify/check correctness of the installation. For this copy the LR templates files from $HOME/MCTDHB-master/Templates to $HOME/TEST  to test Arnoldi-LR
Let say you have installed the Arnoldi LR MCTDHB to: ``` HOME/MCTDHB-master ```.
Copy both lr-templates to ~/TEST:
```
cp -r $HOME/MCTDHB-master/Templates/LR-arnoldi-1D-* $HOME/TEST/.
cd $HOME/tmp/LR-arnoldi-1D-contact
```

Now you need 3-steps to get the LR spectrum: 
#### 1- get MCTDHB ground state
#### 2- constraction of the LR-matrix on-top the above GS
#### 3- Diagonalization of the just constructed LR-matrix

* 1.Step: ```$ $HOME/MCTDHB-master/bin/boson_MCTDHB_gnu_FFTW ``` to get the GS with MCTDHB(2)
```
.....
====================================================================================================
 Job->Relax. Forward     Iteration:     400     Time: [      0.000000 ->     19.950000 +       0.050000 ->     20.000000  ]
       Input orbital energy E(t+0):        9.1090785966024903( CI Dim:        11)(ORB Dim:  2*       256=       512)
          OUT  CI  energy E(t+tau):        9.1090785966024868     N =         10   l0*(N-1)= 0.9000000000     kind of W(x-x'):0
                  Delta E     : +/-       -0.0000000000000089                  Error due to  dE:        0.0000000000000000
                  Error E     : +/-        0.0000000010000000
      New NO's:   9.993561     |  0.6438750E-02 |
 !!!====== Last Point arroaching Tmax=   20.000000000000000
 !!!====== Time-Step is reduced to   =  -1.4921397450962104E-013
      Itr Time:    0.1823329926     execution time:   67.5541808605
 # This computation has been done in
 # /root/tmp/LR-arnoldi-1D-contact
 # Date 01/11/2017; time 12:11:09
 #    Current version 3.3.03 l Heidelberg/Kassel Oct (2017)      #
 # Morb=  2  Npar=      10   Job= Job->Relax. Forward
 # V(x_y_z&t)= x^2
 # W(R=|r1-r2|&t)= Using Defaults from Get_InterParticle.F
 # Kind of W== 0 [if ==0 W=delta(R) else W=f(R) see Get_InterParticle.F]  lambda_0=     0.100000 Time-dependent? (T/F):F
```

* 2.Step: ```$ $HOME/MCTDHB-master/bin/lr_arnoldi_ifort_MKLFFT``` to construct the LR-Matrix
```
.....
Lower part constructed in    2.8000000000000247E-002  seconds

 LR-MATRIX constructed in   1.1000000000000001      seconds
                                                                                                                                                                91,1          32%
LR-MATRIX constructed in   1.1000000000000001      seconds
 It is stored in 'LR_bin'.
```
For the third step (diagonalization) you have to open __./lr_arnoldi.in__ and _replace_ in line 5:  __task=1__ to __task=2__ and run lr-exe again:

* 3.Step: ```$ $HOME/MCTDHB-master/bin/lr_arnoldi_ifort_MKLFFT``` to diagonalize the LR-Matrix
```
....
Time elapsed:   102.92000000000000      seconds
 =======================================
 %-----------------------------------------------%
           5 Ritz values converged
 starting now eigenvector/residual calculation...
 %-----------------------------------------------%

 Eigenvectors and residuals calculated.


# root          Ritz Values (Real)          Ritz Values (Imag)          Residuals

     1             0.1414232165E+001            -0.4583126968E-009             0.7052622311E-011
     2             0.2737188981E+001             0.1215855822E-009             0.2224433324E-010
     3             0.2898173091E+001             0.3654130565E-008             0.7254972480E-011
     4             0.3980689542E+001             0.2748625112E-009             0.5737189571E-010
     5             0.4107632252E+001            -0.4882378929E-008             0.2448773957E-009


  Dimension of Krylov subspace           80
  Tolerance    1.0000000133514320E-010
  Mat*vec operations       43840
  Arnoldi restarts         548
  Converged Ritz values           5
  Eigenvectors stored in subfolder "Eigenv/"
  Computation time    102.95200000000000      seconds

```
