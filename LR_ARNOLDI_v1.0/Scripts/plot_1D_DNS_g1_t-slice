#!/bin/bash
files_to_plot="media/*x-correlations.jpeg \
               media/*k-correlations.jpeg "

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
echo "Plots  at a specified time-slice g1(x,x'|t) g1(k,k'|t) "
echo "Usage: $mctdhb_dir/Scripts/plot_1D_DNS_g1_t-slice 27.5 4 "
echo "where: Time_slice=27.5 Morb=4"
echo "where: $mctdhb_dir is  the installation directory of the MCTDHB_V3* package"
  exit 0 ;;
    "")   
echo "Sorry, you have to provide more arguments :-)"
echo "Plots  at a specified time-slice g1(x,x'|t) g1(k,k'|t) "
echo "Usage: $mctdhb_dir/Scripts/plot_1D_DNS_NO_X_t-slice 27.5 4 "
echo "where: Time_slice=27.5 Morb=4"
echo "where: $mctdhb_dir is  the installation directory of the MCTDHB_V3* package"
exit 0 ;
esac
# media/27.9000000k-correlations.jpeg
# media/27.9000000x-correlations.jpeg

gnuplot -e "Morb=$2" -e "Time=$1" $mctdhb_dir/Scripts/1D_g1RR.plt
gnuplot -e "Morb=$2" -e "Time=$1" $mctdhb_dir/Scripts/1D_g1KK.plt
 

exit 0 
