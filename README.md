# MCTDHB
## MCTDHB usage with MCTDHB-Lab (Mac/Unix/Windows) http://www.mctdhb-lab.com
## MCTDHB usage with docker (5-steps)
1)   install docker (see e.g. https://www.docker.com/community-edition ) 
2)  ````$ docker pull mctdhb/auto-build````
3)  ````$ docker run --rm -it  mctdhb/auto-build````
and you will be inside the docker container in the TEST directory ```root@97f61e1389e7:/TEST#``` with 
__input.in__ and __V_W_Psi_string.in__ files in it. To reproduce some data from [PRA 86 063606] https://journals.aps.org/pra/abstract/10.1103/PhysRevA.86.063606  Table_1 type:
4) ````$ /mctdhb/bin/boson_MCTDHB_gnu_FFTW````
5) Wait ... and compare with reference: root@97f61e1389e7:/TEST````$ vimdiff basic_info.out basic_info.out_Reference ````

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
 #          Alexej.Streltsov@pci.uni-heidelberg.de               #
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
====================================================================================================
 #===============================================================#
 #   The Multiconfigurational Time-Dependent Hartree For Bosons  #
 #    Current version 3.3.01 l Heidelberg/Kassel Apr (2017)      #
 #===============  BBB: Be superB with the mctdhbB ==============#
```

