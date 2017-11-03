#  Templates
Here are the templates of the input.in and V_W_Psi_string.in files needed to reproduce different data or figures 
from the computations published in different journals. 

For example:
1) create a working directory:  
```
$ mkdir TEST 
```
2) copy files from the MCTDHB-master/Templates/PRA_86_063606_Table_1 to HIM_Table_1:
``` 
$ cp MCTDHB-master/Templates/PRA_86_063606_Table_1 TEST/. 
```


3) copy mctdhb-package/bin/* to HIM_Table_1
4) cd to HIM_Table_1 
5) execute the job, typically: mpirun -n 7  ./boson_MCTDHB_gnu_FFTW , the name of the exe file depends on a compiler used...
6) After successfully done computation compare with reference: vimdiff basic_info.out basic_info.out_Reference
7) To get other numbers please change M, Npar in input.in  and restart the job


#            !!!!  Each directory contains induvidual READ-ME file with more instructions and details
