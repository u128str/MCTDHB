# Arnoldi LR MCTDHB 

## A. Arnoldi LR-MCTDHB usage - How2 recompile on Ubintu 16.04
0) In the Makefile file uncomment line: 
 ```mk=./make_systems/ARNOLDI_gcc_mkl.mk ```

ARNOLDI_gcc_mkl.mk file has reference to the DIR_PARPACK with __YOUR__ PARPACK and MKL (it is free) installations:
 ```mk=./make_systems/ARNOLDI_gcc_mkl.mk 
 cmp=ifort
BLAS_LAPACK=intel

FFT=MKLFFT

ARNOLDI=PARPACK
DIR_PARPACK=/home/$USER/ARPACK-NG/lib/
PARPACK=-L $(DIR_PARPACK) -lparpack #-lblas -llapack

mpi_f90=mpif90  -O2 -mkl
mpi_f90=mpif90  -O2  -fopenmp
mpi_f90=mpif90  -g  -fopenmp

#DIR_MKL=$(MKLROOT)/lib/ia32/
DIR_MKL=$(MKLROOT)/lib/intel64/
DIR_MKL=/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64
``` 
  where the respectiv emk-file is:

[mkdnlink]: https://github.com/u128str/MCTDHB/blob/MCTDHB_V3.3.02/make_systems/ARNOLDI_gcc_mkl.mk

{{md make_systems/ARNOLDI_gcc_mkl.mk  }} 

<script src="https://github.com/u128str/MCTDHB/blob/MCTDHB_V3.3.02/make_systems/ARNOLDI_gcc_mkl.mk"></script>

1) ```$ sudo apt-get update && apt-get install -y vim make openmpi-bin libopenmpi-dev fftw3 fftw3-dev libblas-dev liblapack-dev ``` 
2) ```$ cd ```  you are at your $HOME
3) ```$ git clone https://github.com/u128str/MCTDHB.git MCTDHB_V3.3.01```
4) ```$ cd MCTDHB_V3.3.01```
5) ```$ make```
Now copy the templates files to /TEMP reproduce some data from [PRA 86 063606](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063606) [ArXiV:1207.5128](https://arxiv.org/abs/1207.5128) 
6) ```$ cp Templates/PRA_86_063606_Table_1/* /TEMP```
7) ```$ cd /TEMP ```
8)  ```$ $HOME/MCTDHB_V3.3.01/bin/boson_MCTDHB_gnu_FFTW ```
9)  ```$ vimdiff basic_info.out basic_info.out_Reference ```
10) __ENJOY__

__Free__ cross-platform (Mac/Unix/Windows) with GUI: http://www.mctdhb-lab.com

## B. MCTDHB usage with docker (5-steps)
1)   install docker (see e.g. https://www.docker.com/community-edition ) 
2)  ````$ docker pull mctdhb/auto-build````
3)  ````$ docker run --rm -it  mctdhb/auto-build````
and you will be inside the docker container in the TEST directory ```root@97f61e1389e7:/TEST#``` with 
__input.in__ and __V_W_Psi_string.in__ files in it. To reproduce some data from [PRA 86 063606](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063606) [ArXiV:1207.5128](https://arxiv.org/abs/1207.5128)  Table_1 type:
4) ````$ /mctdhb/bin/boson_MCTDHB_gnu_FFTW````
5) Wait ... and compare with reference: 
root@97f61e1389e7:/TEST```$ vimdiff basic_info.out basic_info.out_Reference ```

## C. MCTDHB usage - How2 recompile on Ubintu 16.04
1) ```$ sudo apt-get update && apt-get install -y vim make openmpi-bin libopenmpi-dev fftw3 fftw3-dev libblas-dev liblapack-dev ``` 
2) ```$ cd ```  you are at your $HOME
3) ```$ git clone https://github.com/u128str/MCTDHB.git MCTDHB_V3.3.01```
4) ```$ cd MCTDHB_V3.3.01```
5) ```$ make```
Now copy the templates files to /TEMP reproduce some data from [PRA 86 063606](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063606) [ArXiV:1207.5128](https://arxiv.org/abs/1207.5128) 
6) ```$ cp Templates/PRA_86_063606_Table_1/* /TEMP```
7) ```$ cd /TEMP ```
8)  ```$ $HOME/MCTDHB_V3.3.01/bin/boson_MCTDHB_gnu_FFTW ```
9)  ```$ vimdiff basic_info.out basic_info.out_Reference ```
10) __ENJOY__

```
 #===============================================================#
 #               __  __  ___ _____ ___  _  _ ___                 #
 #  Scientific  |  \/  |/ __|_   _|   \| || | _ ) (2006-present) #
 #              | |\/| | (__  | | | |) | __ | _ \    Germany     #
 #    Package   |_|  |_|\___| |_| |___/|_||_|___/   Heidelberg   #
 #      http://mctdhb.org    http://mctdhb-lab.com/              #
 #===============================================================#
 #   The Multiconfigurational Time-Dependent Hartree For Bosons  #
 #    Current version 3.3.01 l Heidelberg/Kassel Apr (2017)      #
 #===============  BBB: Be superB with the mctdhbB ==============#
 #                        Founders:                              #
 #                                                               #
 #     Alexej I. Streltsov, Ofir E. Alon, Lorenz S. Cederbaum    #
 #===============================================================#
 #            Created, developed and designed by                 #
 #                                                               #
 #           Alexej I. Streltsov (Project Leader)                #
 #          Alexej.Streltsov [At] pci.uni-heidelberg.de          #
 #                   u128str [At] gmail.com                      #
 #                                                               #
 #                     Contributors:                             #
 #            Lorenz S. Cederbaum, Ofir E. Alon                  #
 #      Kaspar Sakmann, Axel U. J. Lode, Julian Grond            #
 #     Oksana I. Streltsova, Shachar Klaiman, Raphael Beinke     #
 #===============================================================#
 #                       Citation:                               #
 #   When citing the MCTDHB Package in the literature,           #
 #   please cite at least one of the papers 1), 2), or 3)        #
 #   as well as the Package 4):                                  #
 #                                                               #
 #    1) A. I. Streltsov, O. E. Alon, and L. S. Cederbaum,       #
 #       Phys. Rev. Lett. 99, 030402 (2007).                     #
 #                                                               #
 #    2) O. E. Alon, A. I. Streltsov, and L. S. Cederbaum,       #
 #       Phys. Rev. A 77, 033613 (2008).                         #
 #                                                               #
 #    3) A. I. Streltsov, O. E. Alon, and L. S. Cederbaum,       #
 #       Phys. Rev. A 81, 022124 (2010).                         #
 #                                                               #
 #    4) The Multiconfigurational Time-Dependent Hartree         #
 #       for Bosons Package, http://mctdhb.org,                  #
 #       A. I. Streltsov,  et al                                 #
 #                                                               #
 #    Current version 3.3.01 l Heidelberg/Kassel Apr (2017)      #
 #===============================================================#
 ____    ____    ______  _________  ______    ____  ____ ______   
|_   \  /   _|.'' ___  ||  _   _  ||_   _ `. |_   ||  _||_   _  \ 
  |   \/   | / .''   \_||_/ | | \_|  | | `. \ | |__| |    | |_)  |
  | |\  /| | | |            | |      | |  | | |  __  |    |  __''.
 _| |_\/_| |_\ `.___.''\   _| |_    _| |_.''/_| |  | |_  _| |__) |
|_____||_____|`.____ .''  |_____|  |______.''|____||____||_______/
        http://mctdhb.org     http://mctdhb-lab.com/              
==================================================================
 #===============================================================#
 #   The Multiconfigurational Time-Dependent Hartree For Bosons  #
 #    Current version 3.3.01 l Heidelberg/Kassel Apr (2017)      #
 #===============  BBB: Be superB with the mctdhbB ==============#
```

