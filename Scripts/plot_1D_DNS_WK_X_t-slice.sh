#!/bin/bash
files_to_plot="media/*timeWrkOrb.jpeg \
               media/*timeWrkOrbReIm_3Dview.jpeg "

#for f in $files_to_plot
#do
#	rm -f $f
#	echo "Removing $f"
#done

# process args
case "$1" in
    clean|cl|-cl*)   
for f in $files_to_plot
do
	rm -f $f
	echo "Removing $f"
done
  exit 0 ;;
    help|-h*)   
echo "Plots  at a specified time-slice one-particle density and WORKING orbitals in real space: Normal and 3D views are available"
echo "Usage: $mctdhb_dir/Scripts/plot_1D_DNS_WK_X_t-slice 27.5 4 "
echo "where: Time_slice=27.5 Morb=4"
echo "where: $mctdhb_dir is  the installation directory of the MCTDHB_V3* package"
  exit 0 ;;
    "")   
echo "Sorry, you have to provide more arguments :-)"
echo "Plots  at a specified time-slice one-particle density and WORKING orbitals in real space: Normal and 3D views are available"
echo "Usage: $mctdhb_dir/Scripts/plot_1D_DNS_WK_X_t-slice 27.5 4 "
echo "where: Time_slice=27.5 Morb=4"
echo "where: $mctdhb_dir is  the installation directory of the MCTDHB_V3* package"
exit 0 ;
esac

gnuplot -e "Morb=$2" -e "Time=$1" $mctdhb_dir/Scripts/1D_DNS_work_orb_x.plt
gnuplot -e "Morb=$2" -e "Time=$1" $mctdhb_dir/Scripts/1D_DNS_work_orb_Re_Im_x_3Dview.plt
 

exit 0 
