#!/bin/bash

# process args
case "$1" in
    clean|cl*)
rm -f ./GetCIJKL*
rm -f ./GetCIJKL1body_Par.f90
rm -f ./GetCIJKL2body_Par.f90
rm -f ./MODULES_CI_SUBR.f90
  exit 0 ;;
    help|-h*)
echo "Makes a CIJKL fortran files for further compilation"
echo "Usage: ./cijkl.bash 2 50 "
echo "where: from M=2 till M=50"
  exit 0 ;;
    "")
echo "Makes a CIJKL fortran files for further compilation"
echo "Usage: ./cijkl.bash 2 50 "
echo "where: from M=2 till M=50"
exit 0 ;
esac

if [ $1 -eq "1" ] ; then
echo "please specify ./cijkl 2 "$2
exit 0 ;
fi


from=$1
till=$2
x=$from
while [ $x -le $till ]
do
#  echo "File GetCIJKL1body_"$x"_OMP.F is created"
#  echo "File GetCIJKL2body_"$x"_OMP.F is created"
#  cat ./part1.txt > File GetCIJKL1body_"$x"_OMP.F 
  x=$(( $x + 1 ))
done

((j=0))
x=$from
while [ $x -le $till ]
do
((j=j+1))
lvarX[j]=${x}
lvarI[j]=i${x}
lvarJ[j]=j${x}
lvar1[j]=GetCIJKL1body_${x}_OMP.F 
#echo ${lvar1[j]}
lvar2[j]=GetCIJKL2body_${x}_OMP.F 
x=$(( $x + 1 ))
done

echo ${lvar1[@]}
echo ${lvar2[@]}
echo ${lvarI[@]}
echo ${lvarJ[@]}
echo ${lvarX[@]}


rm -f ./MODULES_CI_SUBR.f90
#=================== GetCIJKL1body ================================================
echo   "MODULE CI_SUBR " > MODULES_CI_SUBR.f90
echo   "    INTERFACE " >> MODULES_CI_SUBR.f90
echo   "!       SUBROUTINE  GetCIJKL1body_Par(MYID,VIN) " >> MODULES_CI_SUBR.f90
echo   "!        ! CALL  GetCIJKL1body_Par(MYID,VIN) " >> MODULES_CI_SUBR.f90
echo   "!          COMPLEX*16 ::  VIN(:) " >> MODULES_CI_SUBR.f90
echo   "!          INTEGER    :: MYID " >> MODULES_CI_SUBR.f90
echo   "!       END SUBROUTINE  GetCIJKL1body_Par " >> MODULES_CI_SUBR.f90
echo   "!       SUBROUTINE  GetCIJKL2body_Par(MYID,VIN) " >> MODULES_CI_SUBR.f90
echo   "!           ! CALL  GetCIJKL2body_Par(MYID,VIN) " >> MODULES_CI_SUBR.f90
echo   "!          COMPLEX*16 ::  VIN(:) " >> MODULES_CI_SUBR.f90
echo   "!          INTEGER    :: MYID " >> MODULES_CI_SUBR.f90
echo   "!       END SUBROUTINE  GetCIJKL2body_Par " >> MODULES_CI_SUBR.f90
echo   "    SUBROUTINE  PrdCIJKL1body_M_OMP(MYID,VIN)" >> MODULES_CI_SUBR.f90
echo   "       COMPLEX*16 ::  VIN(:)" >> MODULES_CI_SUBR.f90
echo   "       INTEGER    :: MYID" >> MODULES_CI_SUBR.f90
echo   "    END SUBROUTINE  PrdCIJKL1body_M_OMP" >> MODULES_CI_SUBR.f90
echo   "    SUBROUTINE  PrdCIJKL2body_M_OMP(MYID,VIN)" >> MODULES_CI_SUBR.f90
echo   "       COMPLEX*16 ::  VIN(:)" >> MODULES_CI_SUBR.f90
echo   "       INTEGER    :: MYID" >> MODULES_CI_SUBR.f90
echo   "    END SUBROUTINE  PrdCIJKL2body_M_OMP" >> MODULES_CI_SUBR.f90
echo   "    SUBROUTINE  GetCIJKL1body_1(MYID,VIN)" >> MODULES_CI_SUBR.f90
echo   "       COMPLEX*16 ::  VIN(:)" >> MODULES_CI_SUBR.f90
echo   "       INTEGER    :: MYID" >> MODULES_CI_SUBR.f90
echo   "    END SUBROUTINE  GetCIJKL1body_1" >> MODULES_CI_SUBR.f90
echo   "    SUBROUTINE  GetCIJKL2body_1(MYID,VIN)" >> MODULES_CI_SUBR.f90
echo   "       COMPLEX*16 ::  VIN(:)" >> MODULES_CI_SUBR.f90
echo   "       INTEGER    :: MYID" >> MODULES_CI_SUBR.f90
echo   "    END SUBROUTINE  GetCIJKL2body_1" >> MODULES_CI_SUBR.f90
x=2
while [ $x -le $till ]
do
echo   "    SUBROUTINE  GetCIJKL1body_"$x"_OMP(MYID,VIN)" >> MODULES_CI_SUBR.f90
echo   "       COMPLEX*16 ::  VIN(:)" >> MODULES_CI_SUBR.f90
echo   "       INTEGER    :: MYID" >> MODULES_CI_SUBR.f90
echo   "    END SUBROUTINE  GetCIJKL1body_"$x"_OMP" >> MODULES_CI_SUBR.f90
echo   "    SUBROUTINE  GetCIJKL2body_"$x"_OMP(MYID,VIN)" >> MODULES_CI_SUBR.f90
echo   "       COMPLEX*16 ::  VIN(:)" >> MODULES_CI_SUBR.f90
echo   "       INTEGER    :: MYID" >> MODULES_CI_SUBR.f90
echo   "    END SUBROUTINE  GetCIJKL2body_"$x"_OMP" >> MODULES_CI_SUBR.f90
  x=$(( $x + 1 ))
done
echo    "    END INTERFACE" >> MODULES_CI_SUBR.f90
echo    "END MODULE CI_SUBR" >> MODULES_CI_SUBR.f90



rm -f ./GetCIJKL1body_Par.f90
#=================== GetCIJKL1body_Par.f90 ================================================
echo   "      subroutine  GetCIJKL1body_Par(MYID,VIN) " > GetCIJKL1body_Par.f90
echo   "      USE PASS_ARG  " >> GetCIJKL1body_Par.f90
echo   "      USE CI_SUBR " >> GetCIJKL1body_Par.f90
echo   "      USE SHARED_DIMS " >> GetCIJKL1body_Par.f90
echo   "      USE CI_All " >> GetCIJKL1body_Par.f90
echo   "      USE CI_Prod " >> GetCIJKL1body_Par.f90
echo   "      implicit NONE " >> GetCIJKL1body_Par.f90
echo   "!c========================================================= " >> GetCIJKL1body_Par.f90
echo   "      INTEGER ::  MYID " >> GetCIJKL1body_Par.f90
echo   "!c      COMPLEX*16, DIMENSION(Nconf) :: VIN " >> GetCIJKL1body_Par.f90
echo   "      COMPLEX*16 :: VIN(:) " >> GetCIJKL1body_Par.f90
echo   "!c========================================================= " >> GetCIJKL1body_Par.f90
echo   "!     write(6,*)MYID,\"GetCIJKL1body_Par 0\" " >> GetCIJKL1body_Par.f90
echo   "        MorbChoice: SELECT CASE (Morb) " >> GetCIJKL1body_Par.f90
echo   "                    CASE (1) " >> GetCIJKL1body_Par.f90
echo   "                    CALL GetCIJKL1body_1(MYID,VIN) " >> GetCIJKL1body_Par.f90
x=2
while [ $x -le $till ]
do
echo   "                    CASE ("$x") " >> GetCIJKL1body_Par.f90
echo   "      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)" >> GetCIJKL1body_Par.f90
echo   "      IF(CI_Production_1b.eqv..FALSE.) &" >> GetCIJKL1body_Par.f90 
echo   "                       CALL GetCIJKL1body_"$x"_OMP(MYID,VIN) " >> GetCIJKL1body_Par.f90  

x=$(( $x + 1 ))
done
echo   "                          CASE ("$x":100)" >> GetCIJKL1body_Par.f90  
echo   "            write(6,*)\" Still NOT generated!!!!!!!!!\" " >> GetCIJKL1body_Par.f90
echo   "                          CASE DEFAULT " >> GetCIJKL1body_Par.f90
echo   "            write(6,*)\" Something wrong in Number of orbitals!!!!!!!!!\" " >> GetCIJKL1body_Par.f90
echo   "                          END SELECT MorbChoice " >> GetCIJKL1body_Par.f90
echo   "      !     write(6,*) MYID,\"GetCIJKL1body_Par 1\",SUM(VIN) " >> GetCIJKL1body_Par.f90
echo   "            end subroutine GetCIJKL1body_Par " >> GetCIJKL1body_Par.f90


rm -f ./GetCIJKL2body_Par.f90
#=================== GetCIJKL2body_Par.f90 ================================================
echo   "      subroutine  GetCIJKL2body_Par(MYID,VIN) " > GetCIJKL2body_Par.f90
echo   "      USE PASS_ARG  " >> GetCIJKL2body_Par.f90
echo   "      USE CI_SUBR " >> GetCIJKL2body_Par.f90
echo   "      USE SHARED_DIMS " >> GetCIJKL2body_Par.f90
echo   "      USE CI_All " >> GetCIJKL2body_Par.f90
echo   "      USE CI_Prod " >> GetCIJKL2body_Par.f90
echo   "      implicit NONE " >> GetCIJKL2body_Par.f90
echo   "!c========================================================= " >> GetCIJKL2body_Par.f90
echo   "      INTEGER ::  MYID " >> GetCIJKL2body_Par.f90
echo   "!c      COMPLEX*16, DIMENSION(Nconf) :: VIN " >> GetCIJKL2body_Par.f90
echo   "      COMPLEX*16 :: VIN(:) " >> GetCIJKL2body_Par.f90
echo   "!c========================================================= " >> GetCIJKL2body_Par.f90
echo   "!     write(6,*)MYID,\"GetCIJKL2body_Par 0\" " >> GetCIJKL2body_Par.f90
echo   "        MorbChoice: SELECT CASE (Morb) " >> GetCIJKL2body_Par.f90
echo   "                    CASE (1) " >> GetCIJKL2body_Par.f90
echo   "                    CALL GetCIJKL2body_1(MYID,VIN) " >> GetCIJKL2body_Par.f90
x=2
while [ $x -le $till ]
do
echo   "                    CASE ("$x") " >> GetCIJKL2body_Par.f90
echo   "      IF(CI_Production_2b.eqv..TRUE.) CALL PrdCIJKL2body_M_OMP(MYID,VIN)" >> GetCIJKL2body_Par.f90
echo   "      IF(CI_Production_2b.eqv..FALSE.) &" >> GetCIJKL2body_Par.f90 
echo   "                       CALL GetCIJKL2body_"$x"_OMP(MYID,VIN) " >> GetCIJKL2body_Par.f90  

x=$(( $x + 1 ))
done
echo   "                          CASE ("$x":100)" >> GetCIJKL2body_Par.f90  
echo   "            write(6,*)\" Still NOT generated!!!!!!!!!\" " >> GetCIJKL2body_Par.f90
echo   "                          CASE DEFAULT " >> GetCIJKL2body_Par.f90
echo   "            write(6,*)\" Something wrong in Number of orbitals!!!!!!!!!\" " >> GetCIJKL2body_Par.f90
echo   "                          END SELECT MorbChoice " >> GetCIJKL2body_Par.f90
echo   "      !     write(6,*) MYID,\"GetCIJKL2body_Par 1\",SUM(VIN) " >> GetCIJKL2body_Par.f90
echo   "            end subroutine GetCIJKL2body_Par " >> GetCIJKL2body_Par.f90

#=================== GetCIJKL1body ================================================
((j=0))
for i in ${lvar1[@]} ;
do
((j=j+1))
x=${lvarX[$j]}
echo  "        subroutine GetCIJKL1body_"${x}"_OMP(MYID,VIN)" > ${i}
cat ./part1.txt >> ${i}
echo  "       INTEGER, DIMENSION("${lvarX[j]}") :: nvecin,nvecout,Sh_m,Sh_p" >> ${i}
jj=$(($x*($x+1)/2))
echo  "       COMPLEX*16, DIMENSION("${jj}") :: RIJ" >> ${i}

till=${lvarX[$j]}
y=1
while [ $y -le $till ]
do
echo  "        integer :: "i$y,j$y >> ${i}
y=$(( $y + 1 ))
done
cat ./part1_1b.txt >> ${i}
#====================== NESTED DO - LOOPS =======================================
qq="GetCIJKL1body_2_OMP.F"
#echo $qq,$i
echo "Creating file ",$i
if [ $i = $qq ] ; then # Case Morb=2
echo  "        Do i1=Npar+Morb-1-Sh_m(1),Morb-1+Sh_m(2),-1" >> ${i}
echo  "        nvecin(1)=Npar+Morb-1-i1 " >> ${i}
echo  "        nvecin(2)=i1-1 "  >> ${i}
echo  "        K=1+MCNK(i1-1,Morb-1)  "  >> ${i}
echo  "        j1=Npar+Morb-1 - (nvecin(1) - Sh_m(1) + Sh_p(1)) "  >> ${i}
echo  "        J=1+MCNK(j1-1,Morb-1) "  >> ${i}

else                   # Case Morb ne 2
echo  "        Do i1=Npar+Morb-1-Sh_m(1),Morb-1,-1"  >> ${i}
echo  "        nvecin(1)=Npar+Morb-1-i1" >> ${i}
till=$((${lvarX[$j]}-2))
y=2
while [ $y -le $till ]
do
echo  "        Do "i$y"=i"$(($y-1))"-1-Sh_m("$y")   ,Morb-"$y",-1"  >> ${i}
echo  "        nvecin("$y")=i"$(($y-1))"-i"$y"-1" >> ${i}
y=$(( $y + 1 ))
done

last=$((${lvarX[$j]}-1))
echo  "        Do i"$last"=i"$(($last-1))"-1-Sh_m("$last") ,Morb-"$last"+Sh_m("$(($last+1))"),-1"  >> ${i}
echo  "        nvecin("$last")=i"$(($last-1))"-i"$last"-1"  >> ${i}
echo  "        nvecin("$(($last+1))")=i"$last"-1"  >> ${i}

echo  "        K=1+MCNK(i1-1,Morb-1) " >> ${i}
till=$((${lvarX[$j]}-1))
y=2
while [ $y -le $till ]
do
echo  "     $  +MCNK(i"$y"-1,Morb-"$y")   ">> ${i}
y=$(( $y + 1 ))
done

echo  "        j1=Npar+Morb-1 - (nvecin(1) - Sh_m(1) + Sh_p(1))" >> ${i}
till=$((${lvarX[$j]}-1))
y=2
while [ $y -le $till ]
do
echo  "         j"$y"= j"$(($y-1))"-1-(nvecin("$y") - Sh_m("$y") + Sh_p("$y"))">> ${i}
y=$(( $y + 1 ))
done

echo  "        J=1+MCNK(j1-1,Morb-1) " >> ${i}
till=$((${lvarX[$j]}-1))
y=2
while [ $y -le $till ]
do
echo  "     $  +MCNK(j"$y"-1,Morb-"$y")   ">> ${i}
y=$(( $y + 1 ))
done

fi

cat ./part2_1b.txt >> ${i}

till=$((${lvarX[$j]}-1))
y=1
while [ $y -le $till ]
do
echo  "            EndDo">> ${i}
y=$(( $y + 1 ))
done
echo  "                  EndDO">> ${i}

cat ./part3_1b.txt >> ${i}
#====================== NESTED DO - LOOPS DONE =======================================

echo  "        end subroutine GetCIJKL1body_"${x}"_OMP" >> ${i}

done #End DO over every GetCIJKL1body file


#=================== GetCIJKL2body ================================================
((j=0))
for i in ${lvar2[@]} ;
do
((j=j+1))
x=${lvarX[$j]}
echo  "       subroutine GetCIJKL2body_"${x}"_OMP(MYID,VIN)" > ${i}
cat ./part1.txt >> ${i}
echo  "       INTEGER, DIMENSION("${lvarX[j]}") :: nvecin,nvecout,Sh_m,Sh_p" >> ${i}
x=${lvarX[$j]}
#echo "X " ${x}
jj=$(($x*($x+1)/2))
kk=$(($jj*($jj+1)/2))
echo  "       COMPLEX*16, DIMENSION("${kk}") :: RIJKL" >> ${i}
till=${lvarX[$j]}
y=1
while [ $y -le $till ]
do
echo  "       integer :: "i$y,j$y >> ${i}
y=$(( $y + 1 ))
done
cat ./part1_2b.txt >> ${i}
#====================== NESTED DO - LOOPS =======================================
qq="GetCIJKL2body_2_OMP.F"
echo "Creating file ",$i
if [ $i = $qq ] ; then # Case Morb=2
echo  "        Do i1=Npar+Morb-1-Sh_m(1),Morb-1+Sh_m(2),-1" >> ${i}
echo  "        nvecin(1)=Npar+Morb-1-i1 " >> ${i}
echo  "        nvecin(2)=i1-1 "  >> ${i}
echo  "        K=1+MCNK(i1-1,Morb-1)  "  >> ${i}
echo  "        j1=Npar+Morb-1 - (nvecin(1) - Sh_m(1) + Sh_p(1)) "  >> ${i}
echo  "        J=1+MCNK(j1-1,Morb-1) "  >> ${i}
else                   # Case Morb ne 2
echo  "        Do i1=Npar+Morb-1-Sh_m(1),Morb-1,-1"  >> ${i}
echo  "        nvecin(1)=Npar+Morb-1-i1" >> ${i}
till=$((${lvarX[$j]}-2))
y=2
while [ $y -le $till ]
do
echo  "        Do "i$y"=i"$(($y-1))"-1-Sh_m("$y")   ,Morb-"$y",-1"  >> ${i}
echo  "        nvecin("$y")=i"$(($y-1))"-i"$y"-1" >> ${i}
y=$(( $y + 1 ))
done
last=$((${lvarX[$j]}-1))
echo  "        Do i"$last"=i"$(($last-1))"-1-Sh_m("$last") ,Morb-"$last"+Sh_m("$(($last+1))"),-1"  >> ${i}
echo  "        nvecin("$last")=i"$(($last-1))"-i"$last"-1"  >> ${i}
echo  "        nvecin("$(($last+1))")=i"$last"-1"  >> ${i}


echo  "        K=1+MCNK(i1-1,Morb-1) " >> ${i}
till=$((${lvarX[$j]}-1))
y=2
while [ $y -le $till ]
do
echo  "     $  +MCNK(i"$y"-1,Morb-"$y")   ">> ${i}
y=$(( $y + 1 ))
done

echo  "        j1=Npar+Morb-1 - (nvecin(1) - Sh_m(1) + Sh_p(1))" >> ${i}
till=$((${lvarX[$j]}-1))
y=2
while [ $y -le $till ]
do
echo  "         j"$y"= j"$(($y-1))"-1-(nvecin("$y") - Sh_m("$y") + Sh_p("$y"))">> ${i}
y=$(( $y + 1 ))
done

echo  "        J=1+MCNK(j1-1,Morb-1) " >> ${i}
till=$((${lvarX[$j]}-1))
y=2
while [ $y -le $till ]
do
echo  "     $  +MCNK(j"$y"-1,Morb-"$y")   ">> ${i}
y=$(( $y + 1 ))
done

fi

cat ./part2_2b.txt >> ${i}

till=$((${lvarX[$j]}-1))
y=1
while [ $y -le $till ]
do
echo  "            EndDo">> ${i}
y=$(( $y + 1 ))
done
echo  "                  EndDO Iloop">> ${i}

cat ./part3_2b.txt >> ${i}
#====================== NESTED DO - LOOPS DONE =======================================

echo  "        end subroutine GetCIJKL2body_"${x}"_OMP" >> ${i}

done #End DO over every GetCIJKL2body file














