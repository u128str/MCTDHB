# MCTDHB  

#### [Arnoldi LR MCTDHB how2 get/install/test etc](https://github.com/u128str/MCTDHB/blob/MCTDHB_V3.3.03/LR-ARNOLDI.md)

## A. How2 Get the sources of the MCTDHB package?
<details>
<summary> Click here to see how to get ./MCTDHB-master with git, wget or curl:</summary>
a)  Clone latest iversion of the MCTDHB package to the directory MCTDHB-master:
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


## B. MCTDHB usage
There are __three__ ways to do simulations with the MCTDHB package.
The standard and hardest one is to download the sources, recompile, run.
Second way is to use the popular Docker technology - install docker -> douwnload the mctdhb dockar-image and run it... 
The third and easyest way is to try the MCTDHB-Laboratory (installers avalable for Win/Unix/Mac). 

### I. MCTDHB-Laboratory
[__Free__ cross-platform (Mac/Unix/Windows) GUI](http://www.mctdhb-lab.com)

[Mac](http://www.mctdhb-lab.com/images/how2-figs/launch_mac.jpg)
[Windows](http://www.mctdhb-lab.com/images/how2-figs/Win_appearence.jpg)
[Linux](http://www.mctdhb-lab.com/images/how2-figs/launch_linux.jpg)

### II. MCTDHB with docker
#### [How2 install docker] (https://www.docker.com/community-edition)
<details>
<summary> click here to see the steps needed to download and run the mctdhb docker image (you do it once to use afterwards):</summary>
1) Download the latest MCTDHB docker image (450Mb)
<pre><code>
$ docker pull mctdhb/auto-build
</code></pre>
2)  Run mctdhb/auto-build docker:
<pre><code>
$ docker run --hostname mctdhb-user --rm -it -v $(pwd):/tmp mctdhb/auto-build
</code></pre>

<details>
<summary> click here to see the steps needed to build the mctdhb docker image (you do it once to use afterwards):</summary>
1)  Get ./MCTDHB-master with the above step A and cd to it:
<pre><code>
$ cd $HOME/MCTDHB-master
</code></pre>
2)  Build (~4 mins) the mctdh Docker-image from available Dokerfile (final image size is about of ~450MB)
<pre><code>
$ docker build --no-cache -f Dockerfile -t mctdhb-user .
</code></pre>
2)  Run mctdhb/aouto-build docker:
<pre><code>
$ docker run --hostname mctdhb-user --rm -it -v $(pwd):/tmp mctdhb-user
</code></pre>
</details>
At this point the MCTDHB package is installed in your Docker-Ubuntu system at $HOME/MCTDHB-master:
<pre><code>
user@mctdhb-user:~/MCTDHB-master/bin$ ls -ltr
total 2668
-rwxrwxr-x 1 user user 1306232 Nov  2 13:10 boson_MCTDHB_gnu_FFTW
-rwxrwxr-x 1 user user 1423048 Nov  2 13:10 properties_LR_gnu_FFTW
</code></pre>
</details>



###  III MCTDHB - How2 recompile on Ubuntu 16.04
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


Now copy the templates files to /TEMP reproduce some data from [PRA 86 063606](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063606) [ArXiV:1207.5128](https://arxiv.o
6) ```$ cp Templates/PRA_86_063606_Table_1/* /TEMP```
7) ```$ cd /TEMP ```
8)  ```$ $HOME/MCTDHB-master/bin/boson_MCTDHB_gnu_FFTW ```
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

