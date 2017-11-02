# LR Arnoldi MCTDHB:
## A. How2 get it? It is already in your $HOME/MCTDHB-master!
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

## B. Arnoldi LR MCTDHB usage with docker technology (5-steps):
#### [How2 install docker] (https://www.docker.com/community-edition)
<details>
<summary> click here to see the steps needed to build lr-mctdhb-user Docker-image (you do it once to use afterwards):</summary>
1)  Get ./MCTDHB-master with the above A. step and cd to it:
<pre><code>
$ cd $HOME/MCTDHB-master
</code></pre>
4)  Build (~14 mins) the __lr-mctdhb-user__ Docker-image from available Dokerfile.LR.user (Why rebuild locally? Because it  installs/rebuilds MKL+parpack+... final image size is about of ~4.5GB)
<pre><code>
$ docker build --no-cache -f Dockerfile.LR.user -t lr-mctdhb-user .
</code></pre>
You will see:
<pre><code>
MCTDHB-master$ docker build --no-cache -f Dockerfile.LR.user -t lr-mctdhb-user .
Sending build context to Docker daemon  22.44MB
Step 1/11 : FROM mctdhb/minunix:latest
 ---> ff5670deb65e
Step 2/11 : MAINTAINER Alexej I. Streltsov  <u128str@gmail.com>
 ---> Running in 095ef9d2fc04
 ---> 65ebf30fa94a
Removing intermediate container 095ef9d2fc04
Step 3/11 : ENV MKL_PATH /opt/intel
 ---> Running in 6e13ff131f83
 ---> 233ca4ea1ed5
Removing intermediate container 6e13ff131f83
Step 4/11 : RUN apt-get update &&   apt-get install -y man tar wget cpio unzip autoconf sudo
 ---> Running in 53f710f63281
Get:1 http://security.ubuntu.com/ubuntu xenial-security InRelease [102 kB]
Hit:2 http://archive.ubuntu.com/ubuntu xenial InRelease
Get:3 http://archive.ubuntu.com/ubuntu xenial-updates InRelease [102 kB]
.....
Successfully tagged lr-mctdhb-user:latest
</code></pre>
To check the available dockers, type <code> docker images</code>:
<pre><code>
$ docker images
~/MCTDHB-master$ docker images
REPOSITORY          TAG                 IMAGE ID            CREATED              SIZE
lr-mctdhb-user      latest              7284eda37179        About a minute ago   4.3GB
mctdhb/auto-build   latest              8dad46489fd3        44 minutes ago       532MB
mctdhb/minunix      latest              ff5670deb65e        13 days ago          434MB
</code></pre>
</details>

<details>
<summary> click here to see how to run/use the builded above Docker image lr-mctdhb-user</summary>
1)  Run docker:
<pre><code>
$ docker run --hostname lr-mctdhb-user --rm -it lr-mctdhb-user
</code></pre>
and you will be inside the docker container in the TEST directory
<pre><code> 
user@lr-mctdhb-users:~$ ls -ltr
total 11680
-rw-rw-r--  1 user user 11950772 Nov  2 16:53 MCTDHB-MCTDHB_V3.3.03.zip
drwxrwxr-x 17 user user     4096 Nov  2 16:54 MCTDHB-MCTDHB_V3.3.03
drwxrwxr-x  4 user user     4096 Nov  2 16:54 TEST

</code></pre>
</details>


## C. MCTDHB LR Arnoldi - How2 recompile on Ubuntu 16.04

<details>
<summary> click here to see how recompile the LR-MCTDHB code on your Ubuntu</summary>

1) ```$ sudo apt-get update && apt-get install -y man tar wget cpio unzip autoconf vim make openmpi-bin libopenmpi-dev fftw3 fftw3-dev libblas-dev liblapack-dev ``` 
2) ```$ mkdir $HOME/tmp && cd $HOME/tmp ```  you are at your $HOME/tmp
3) ```$ wget -q http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/12070/l_mkl_2018.0.128.tgz  ``` Download MKL install package l_mkl_2018.0.128.tgz
4) ```$ tar -xzf l_mkl_2018.0.128.tgz  && cd l_mkl_2018.0.128 && sed -i 's/ACCEPT_EULA=decline/ACCEPT_EULA=accept/g' silent.cfg ``` 
5) ```$ sudo ./install.sh -s silent.cfg```
6) ```$ wget https://github.com/opencollab/arpack-ng/archive/master.zip && unzip master.zip && cd arpack-ng-master && sh bootstrap && ./configure  --enable-mpi && make ```
7) ```$ sudo make install ```  Installing arpack libs globally to: /usr/local/lib
8) ```$ sudo echo "${MKL_PATH}/mkl/lib/intel64" >> /etc/ld.so.conf.d/intel.conf && ldconfig && echo ". /opt/intel/bin/compilervars.sh intel64" >> /etc/bash.bashrc``` 
9) ```$ . /opt/intel/bin/compilervars.sh intel64``` Defining MKL-related variables
10) ```$ cd ```
11) ```$ wget --no-check-certificate --content-disposition https://github.com/u128str/MCTDHB/archive/MCTDHB_V3.3.03.zip```
12) ```$ unzip MCTDHB-MCTDHB_V3.3.03.zip ```
13) ```$ cd MCTDHB-MCTDHB_V3.3.03```
14) ```$ make mk_file=ARNOLDI_gcc_mkl.mk``` Compilation of the Arnoldi LR MCTDHB
__Congrads!__
At this point the LR-Arnoldi-MCTDHB package is installed in your Ubuntu system at $HOME/MCTDHB-MCTDHB_V3.3.03:

```
~/MCTDHB-MCTDHB_V3.3.03/bin# ls -ltrh
total 7.4M
-rwxr-xr-x 1 root root 2.3M Nov  1 12:04 boson_MCTDHB_ifort_MKLFFT
-rwxr-xr-x 1 root root 2.6M Nov  1 12:04 properties_LR_ifort_MKLFFT
-rwxr-xr-x 1 root root 2.7M Nov  1 12:04 lr_arnoldi_ifort_MKLFFT
```
</details>

## D. MCTDHB LR Arnoldi - How2 verify/test 
Next stage is to verify/check correctness of the installation. For this copy the LR templates files from $HOME/MCTDHB-MCTDHB_V3.3.03/Templates to $HOME/tmp  to test Arnoldi-LR
Let say you have installed the Arnoldi LR MCTDHB to: _``` HOME/MCTDHB-MCTDHB_V3.3.03 ```.
Copy both lr-templates to ~/tmp:

1) ```$ cp -r $HOME/MCTDHB-MCTDHB_V3.3.03/Templates/LR-arnoldi-1D-* $HOME/tmp/.```
2) ```$ cd $HOME/tmp/LR-arnoldi-1D-contact```

Now you need 3-stages to get LR spectrum: 
#### 1- MCTDHB ground state; 
#### 2- constraction of the LR-matrix on-top; 
#### 3- Diagonalization of the LR-matrix

3) ```$ $HOME/MCTDHB-MCTDHB_V3.3.03/bin/boson_MCTDHB_gnu_FFTW ``` #### 1- Stage get Ground MCTDHB(2) State
```
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

4) ```$ $HOME/MCTDHB-MCTDHB_V3.3.03/bin/lr_arnoldi_ifort_MKLFFT``` #### 2 Stage- Construct LR-Matrix

```
Lower part constructed in    2.8000000000000247E-002  seconds

 LR-MATRIX constructed in   1.1000000000000001      seconds
                                                                                                                                                                91,1          32%
LR-MATRIX constructed in   1.1000000000000001      seconds
 It is stored in 'LR_bin'.
```
For the third stage (diagonalization) you have to open __./lr_arnoldi.in__ and _replace in line 5:  __task=1__ to __task=2__ and run again_

5) ```$ $HOME/MCTDHB-MCTDHB_V3.3.03/bin/lr_arnoldi_ifort_MKLFFT``` #### 3 Stage- Diagonalization of the LR-Matrix
```
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

