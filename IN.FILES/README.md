# Here are the examples of the *.in files used by the MCTDHB and/or Arnoldi LR MCTDHB packages 

## Prerequisites
MCTDHB and/or Arnoldi LR MCTDHB packages installed localy or with Docker:
#### [MCTDHB how2 get/install/test etc](https://github.com/u128str/MCTDHB/blob/master/README.md)
#### [Arnoldi LR MCTDHB how2 get/install/test etc](https://github.com/u128str/MCTDHB/blob/master/LR-ARNOLDI.md)

You can use these files to start with...  __BUT it is better__ to download a specific directory with the corresponding template
[for full set of the available templtes see|(https://github.com/u128str/MCTDHB/edit/master/Templates/README.md)
### To download, say, template _LR-arnoldi-1D-contact_:

```
curl  https://codeload.github.com/u128str/MCTDHB/tar.gz/master  | \
tar -xz --strip=2 MCTDHB-master/Templates/LR-arnoldi-1D-contact
```

#  Templates
Every directory contains template files needed to reproduce different data or figures 
from the computations published in different journals.

### There are three basic kinds of computations available within MCTDHB technology: 

| Type of computation       | Input files needed |  Number of steps needed |
| ------------- |:-------------:| :-------------: |
| MCTDHB            | input.in, V_W_Psi_string.in  | 1. standard MCTDHB |
| LR-MCTDHB         | input.in, V_W_Psi_string.in, properties.in | 1. standard MCTDHB  -> 2. LR Analysys|
| Arnoldi LR-MCTDHB | input.in, V_W_Psi_string.in, lr_arnoldi.in | 1. standard MCTDHB  -> 2. Arnoldi LR-Matrix construction -> 3. LR-Matrix diagonalization |

For example, let say you would like to reproduce the results from  __PRA_86_063606_Table_1__.
This is a simple 1.Step MCTDHB computation. 

1) create a working directory TEST and cd to it:
```
mkdir $HOME/TEST \
cd $HOME/TEST 
```
2) Copy the corresponding template directory with all files from the MCTDHB-master/Templates/__PRA_86_063606_Table_1__ to your TEST:
``` 
$ cp -r $HOME/MCTDHB-test/Templates/PRA_86_063606_Table_1 $HOME/TEST/.
```
or download __PRA_86_063606_Table_1__ from web with ___curl___:
```
curl  https://codeload.github.com/u128str/MCTDHB/tar.gz/master  | \
tar -xz --strip=2 MCTDHB-master/Templates/PRA_86_063606_Table_1
```
and cd to the just copied directory: 
```
cd PRA_86_063606_Table_1
```

3) Execute the MCTDHB job: 
 ```
 $ mpirun -n 2 $HOME/MCTDHB-master/bin/boson_MCTDHB_gnu_FFTW 
 ``` 
the name of the exefile depends on a compiler used here it is gnu and FFTW...

4) After successfully done computation compare your result with the reference one: 
 ```
vimdiff basic_info.out basic_info.out_Reference
 ```
5) To get other numbers please change M, Npar in input.in  and restart the job

#### Each directory contains individual README(.md) files with more instructions and output details...
