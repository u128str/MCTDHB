#!/bin/bash
if  [ -d ./user_guesslib.ORG ]
then
echo "Dir ./user_guesslib.ORG already exists, please rm or rename it!" 
exit
fi
echo " I make the copy of the current user_guesslib "
cp -r ./user_guesslib ./user_guesslib.ORG

if test "$1" != "all"
then
(k=0)
while test "$1" != "" ; do
        ((k=k+1))
        tst[k]=./test/$1_test.tar
        tstdir[k]=$1_test
        echo "You selected to perform the test:" $1  
        shift
done
elif test "$1" = "all"
then
tst[1]=./test/01_test.tar
tst[2]=./test/02_test.tar
tst[3]=./test/03_test.tar 
tst[4]=./test/04_test.tar
tst[5]=./test/05_test.tar
tst[6]=./test/06_test.tar
tst[7]=./test/07_test.tar
tst[8]=./test/08_test.tar
tst[9]=./test/09_test.tar
tst[10]=./test/10_test.tar
tstdir[1]=01_test
tstdir[2]=02_test
tstdir[3]=03_test 
tstdir[4]=04_test
tstdir[5]=05_test
tstdir[6]=06_test
tstdir[7]=07_test
tstdir[8]=08_test
tstdir[9]=09_test
tstdir[10]=10_test
fi
# collect the tests here
TEST_LIST1=${tst[@]}
# and the dirs there
TEST_LIST2=${tstdir[@]}

echo "TEST_LIST1=" ${TEST_LIST1}

# clean stuff
#make cl
if test -e testing.out; then
rm testing.out
fi
# unpack the tests wanted
for j in ${TEST_LIST1};
do
tar -xvf $j 
done

# perform tests
for j in ${TEST_LIST2}
do
if test "$j" == "10_test"
then
if test -e 09_test/PSI_bin 
then 
cp 09_test/PSI_bin 10_test/
else
echo "Cannot propagate when you do not create the initial guess!!! 10_test needs 09_test!"
exit
fi
if test -e 09_test/CIc_bin 
then
cp 09_test/CIc_bin 10_test/
else
echo "Cannot propagate when you do not create the initial guess!!! 10_test needs 09_test!"
exit
fi
fi
cd test
cp compare.awk ../$j
cd ../$j
cp *.F ../user_guesslib
cd ..
echo "To do the tests I use the makefile 'Makefile',"
echo "if you want to use some other makefile, please"
echo "do 'cp Makefile Makefile.BAK' and"
echo "'cp Makefile.YOUR_SYSTEM Makefile'"

make -f Makefile > crap
cp ./bin/libguess.so $j
cp ./bin/boson* $j
cd $j
echo "I'm doing test $j now"
./boson_MCTDHB_intel > ../$j.`date +%d_%m_%H_%M_%S`.out
#./boson_MCTDHB_intel > ../$j.`date +%d`.out
case $j in
     01_test)
       echo "####################################################################" 
       echo "Performing relaxation in parabolic trap with HO DVR and W=gauss_0.05" 
       echo "####################################################################" 
     ;;
     02_test)
       echo "####################################################################"
       echo "Performing relaxation in parabolic trap with sine DVR and W=gauss_0.05" 
       echo "####################################################################" 
     ;;
     03_test)
       echo "####################################################################"
       echo "Performing relaxation in parabolic trap with FFT and W=gauss_0.05" 
       echo "####################################################################" 
     ;;
     04_test)
       echo "####################################################################"
       echo "Performing relaxation in parabolic trap with EXP DVR and W=gauss_0.05"
       echo "####################################################################"
     ;;
     05_test)
       echo "####################################################################" 
       echo "Performing relaxation in parabolic trap with HO DVR and W=delta" 
       echo "####################################################################" 
     ;;
     06_test)
       echo "####################################################################" 
       echo "Performing relaxation in parabolic trap with sine DVR and W=delta" 
       echo "####################################################################" 
     ;;
     07_test)
       echo "####################################################################"
       echo "Performing relaxation in parabolic trap with FFT DVR and W=delta"
       echo "####################################################################" 
     ;;
     08_test)
       echo "####################################################################" 
       echo "Performing relaxation in parabolic trap with EXP DVR and W=delta" 
       echo "####################################################################" 
     ;;
     09_test)
       echo "####################################################################" 
       echo "Performing relaxation in 2D parabolic trap with sine DVR and W=gauss_0.05" 
       echo "####################################################################"
     ;;
     10_test)
       echo "####################################################################" 
       echo "Performing propagation in 2D DW trap with sine DVR and W=gauss_0.05" 
       echo "####################################################################" 
     ;;
esac
awk -f compare.awk current reference 
awk -f compare.awk current reference >> ../testing.out 
cd ..
chk=`less testing.out | grep problematic`
if test "$chk" != "" 
then 
 echo "TEST $j is problematic!!!"
 echo "`less testing.out | grep problematic`"
else
 echo "TEST $j is OK!!!"
fi
rm crap
done

if test -e ./test-outputs
then
   echo "I copy all the outputs to ./test-outputs"
   mv *_test*.out ./test-outputs/
   ls -ltr ./test-outputs/*
else
   echo "I copy all the outputs to ./test-outputs"
   mkdir test-outputs
   mv *_test*.out ./test-outputs/
   ls -ltr ./test-outputs/*
fi

for j in ${TEST_LIST2}
do
rm -r $j
done
echo " I restore the DIR user_guesslib "
rm -rf user_guesslib
mv user_guesslib.ORG user_guesslib
echo "OK"
exit ;
