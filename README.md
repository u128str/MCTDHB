# MCTDHB  

#### [Arnoldi LR MCTDHB how2 get/install/test etc](https://github.com/u128str/MCTDHB/blob/master/LR-ARNOLDI.md)

## A. How2 Get the sources of the MCTDHB package?
<details>
<summary> Click here to see how to get ./MCTDHB-master with git, wget or curl:</summary>
a)  Clone the latest version of the MCTDHB package to the directory MCTDHB-master:
<pre><code>
git clone https://github.com/u128str/MCTDHB.git MCTDHB-master
</code></pre>
b)  OR download zip-archive MCTDHB-master.zip:
<pre><code>
wget --no-check-certificate --content-disposition https://github.com/u128str/MCTDHB/archive/master.zip
</code></pre>
<pre><code>
curl -LJO https://github.com/u128str/MCTDHB/archive/master.zip
</code></pre>
Unzip the downloaded archive to the directory MCTDHB-master
<pre><code>
unzip MCTDHB-master.zip
</code></pre>
</details>


## B. MCTDHB usages
There are __three__ ways to do simulations with the MCTDHB package.
The standard and the hardest one is to download on your local computer sources -> recompile -> run.
Second way is to use the popular Docker technology - install docker -> download the mctdhb docker-image -> run it... 
The third and the easiest way is to use the MCTDHB-Laboratory with GUI (installers are available for Win/Unix/Mac). 

### I. MCTDHB with MCTDHB-Laboratory
[__Free__ cross-platform (Mac/Unix/Windows) GUI](http://www.mctdhb-lab.com)

[Mac](http://www.mctdhb-lab.com/images/how2-figs/launch_mac.jpg)
[Windows](http://www.mctdhb-lab.com/images/how2-figs/Win_appearence.jpg)
[Linux](http://www.mctdhb-lab.com/images/how2-figs/launch_linux.jpg)

### II. MCTDHB with Docker technology
#### [How2 install docker] (https://www.docker.com/community-edition)
<details>
<summary> click here to see the steps needed to download and run the mctdhb docker image (you do it once to use afterwards):</summary>
1) To download the latest MCTDHB docker image (532MB) [mctdhb/auto-build](https://hub.docker.com/r/mctdhb/auto-build/) type:
<pre><code>
$ docker pull mctdhb/auto-build
....
$ docker images
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
mctdhb/auto-build   latest              8dad46489fd3        16 minutes ago      532MB
mctdhb/minunix      latest              ff5670deb65e        13 days ago         434MB
</code></pre>
2)  Run mctdhb/auto-build docker:
<pre><code>
$ docker run --hostname mctdhb-user --rm -it -v $(pwd):/tmp mctdhb/auto-build
</code></pre>
At this point you are inside your Docker-Ubuntu system with the MCTDHB package installed in $HOME/MCTDHB-master:
<pre><code>
user@mctdhb-user:~/MCTDHB-master/bin$ ls -ltr
total 2668
-rwxrwxr-x 1 user user 1306232 Nov  2 13:10 boson_MCTDHB_gnu_FFTW
-rwxrwxr-x 1 user user 1423048 Nov  2 13:10 properties_LR_gnu_FFTW
</code></pre>

<details>
<summary> click here to see the steps needed to build the mctdhb docker image (actually you don't need it, but in the case...):</summary>
a)  Get ./MCTDHB-master with the above step A and cd to it:
<pre><code>
$ cd $HOME/MCTDHB-master
</code></pre>
b)  Build (~4 mins) the mctdhb-user Docker-image from available Dockerfile (final image size is about of ~532MB)
<pre><code>
$ docker build --no-cache -f Dockerfile -t mctdhb-user .
....
Successfully tagged mctdhb-user:latest
$ docker images
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
mctdhb-user         latest              266f4b3c721e        33 seconds ago      532MB
mctdhb/auto-build   latest              8dad46489fd3        8 minutes ago       532MB
mctdhb/minunix      latest              ff5670deb65e        13 days ago         434MB
</code></pre>
c)  Run the just-built mctdhb-user Docker-image:
<pre><code>
$ docker run --hostname mctdhb-user --rm -it -v $(pwd):/tmp mctdhb-user
</code></pre>
</details>
</details>

###  III. MCTDHB with Ubuntu - _How2 recompile on Ubuntu 16.04_
<details>
<summary> click here to see how recompile the LR-MCTDHB code on your Ubuntu</summary>

1) ```$ sudo apt-get update && apt-get install -y vim make openmpi-bin libopenmpi-dev fftw3 fftw3-dev libblas-dev liblapack-dev ``` 
2) ```$ cd ```  you are at your $HOME
3) ```$ git clone https://github.com/u128str/MCTDHB.git MCTDHB-master```
4) ```$ cd MCTDHB-master```
5) ```$ make```

At this point the MCTDHB package is installed in your local Ubuntu system at $HOME/MCTDHB-master:

```
MCTDHB-master/bin$ ls -ltr
total 2668
-rwxrwxr-x 1 user user 1306232 Nov  2 13:10 boson_MCTDHB_gnu_FFTW
-rwxrwxr-x 1 user user 1423048 Nov  2 13:10 properties_LR_gnu_FFTW
```
</details>


## C. MCTDHB first run/test
The primary goal now is to verify/check correctness of the installation. 
The secondary goal is to reproduce some data from [PRA 86 063606](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063606) [ArXiV:1207.5128](https://arxiv.org/abs/1207.5128).

##### ([Here you can see how this example works with the MCTDHB-Lab](https://github.com/u128str/MCTDHB/blob/master/Templates/PRA_86_063606_Table_1/README.md))


To run/test your local or Docker installation copy the __PRA_86_063606_Table_1__ templates files from $HOME/MCTDHB-master/Templates to $HOME/TEST. Here we assume that you have installed the MCTDHB package to  ```$HOME/MCTDHB-master```:

1) ```cd && mkdir TEST``` make TEST directory in your $HOME. In the case of usage of the MCTDHB with Docker ./TEST is already exists!
2) ```$ cp -r $HOME/MCTDHB-test/Templates/PRA_86_063606_Table_1 $HOME/TEST/.```
3) ```$ cd $HOME/TEST/PRA_86_063606_Table_1``` All necessary input files, i.e, __input.in__ and __V_W_Psi_string.in__ should be in this directory:
```
user@mctdhb-user:~/TEST/PRA_86_063606_Table_1$ ls -ltrh *.in
-rwxr-xr-x 1 user user 4.9K Nov  2 13:10 input.in
-rwxr-xr-x 1 user user 1.4K Nov  2 13:10 V_W_Psi_string.in
```
4)  ```$ mpirun -n 2 $HOME/MCTDHB-master/bin/boson_MCTDHB_gnu_FFTW ``` run MCTDHB simulation
```
....
====================================================================================================
 Job->Relax. Forward     Iteration:     200     Time: [      0.000000 ->      9.950000 +       0.050000 ->     10.000000  ]
       Input orbital energy E(t+0):        7.0383484153111748( CI Dim:      3003)(ORB Dim:  6*       128=       768)
          OUT  CI  energy E(t+tau):        7.0383484153111668     N =         10   l0*(N-1)= 0.5000000000     kind of W(x-x'):4
                  Delta E     : +/-       -0.0000000000000053                  Error due to  dE:                       NaN
                  Error E     : +/-        0.0000010000000000
      New NO's:   9.968427     |  0.3147304E-01 |  0.9936897E-04 |  0.3137349E-06 |  0.9905159E-09 |  0.3105116E-11 | 
 !!!====== Last Point arroaching Tmax=   10.000000000000000     
 !!!====== Time-Step is reduced to   =  -8.8817841970012523E-015
      Itr Time:    0.2723159790     execution time:   62.0533051491
 # This computation has been done in
 # /home/user/TEST/PRA_86_063606_Table_1
 # Date 02/11/2017; time 14:43:56
 #    Current version 3.3.01 l Heidelberg/Kassel Apr (2017)      # 
 # Morb=  6  Npar=      10   Job= Job->Relax. Forward
 # V(x_y_z&t)= 0.5d0*x^2                                                                                           
 # W(R=|r1-r2|&t)= r^2                                                                                                 
 # Kind of W== 4 [if ==0 W=delta(R) else W=f(R) see Get_InterParticle.F]  lambda_0=     0.055556 Time-dependent? (T/F):F
======================================================================================================================
====================     basic_info.out file (with PATH, time, E, etc) has been written            ===================
======================================================================================================================
  Master               0  is DONE after   64.159628868103027     

```
5)  - to compare with the Reference numbers 
```
$ vimdiff basic_info.out basic_info.out_Reference 
```
6)  __ENJOY__ MCTDHB with more [templates e.g. PRA_86_063606, PRA_88_023606, PRL_106_240401, PRL_99_030402 and others ](https://github.com/u128str/MCTDHB/tree/master/Templates)

```
 #===============================================================#
 #               __  __  ___ _____ ___  _  _ ___                 #
 #  Scientific  |  \/  |/ __|_   _|   \| || | _ ) (2006-present) #
 #              | |\/| | (__  | | | |) | __ | _ \    Germany     #
 #    Package   |_|  |_|\___| |_| |___/|_||_|___/   Heidelberg   #
 #      http://mctdhb.org    http://mctdhb-lab.com/              #
 #===============================================================#
 #   The Multiconfigurational Time-Dependent Hartree For Bosons  #
 #    Current version 3.3.03 l Heidelberg/Kassel Nov (2017)      #
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
 #    Current version 3.3.03 l Heidelberg/Kassel Nov (2017)      #
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
 #    Current version 3.3.03 l Heidelberg/Kassel Nov (2017)      #
 #===============  BBB: Be superB with the mctdhbB ==============#
```

