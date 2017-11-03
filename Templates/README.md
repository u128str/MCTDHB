# Prerequisites
MCTDHB and/or Arnoldi LR MCTDHB packages installed loccaly or with Docker:
#### [MCTDHB how2 get/install/test etc](https://github.com/u128str/MCTDHB/blob/master/README.md)
#### [Arnoldi LR MCTDHB how2 get/install/test etc](https://github.com/u128str/MCTDHB/blob/master/LR-ARNOLDI.md)

#  Templates
Here are the templates of the input.in and V_W_Psi_string.in files needed to reproduce different data or figures 
from the computations published in different journals. 

For example:
1) create a working directory:  
```
$ mkdir $HOME/TEST 
```
2) copy files from the MCTDHB-master/Templates/PRA_86_063606_Table_1 to HIM_Table_1:
``` 
$ cp -r $HOME/MCTDHB-test/Templates/PRA_86_063606_Table_1 $HOME/TEST/.
```
3) cd to HIM_Table_1
```
$ cd  $HOME/TEST 
```
3) execute the job, typically: 
 ```
 $ mpirun -n 2 $HOME/MCTDHB-master/bin/boson_MCTDHB_gnu_FFTW 
 ``` 
the name of the exe file depends on a compiler used...
4) After successfully done computation compare with reference: 
 ```
vimdiff basic_info.out basic_info.out_Reference
 ```
5) To get other numbers please change M, Npar in input.in  and restart the job

####!!!!  Each directory contains induvidual READ-ME file with more instructions and details
