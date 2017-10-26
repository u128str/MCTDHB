#!/bin/bash
files_to_plot="media/fig_E_of_t_PLT.jpeg         \
               media/fig_E_of_t_T_V_W_PLT.jpeg   \
               media/fig_OP_of_t_x_X2_PLT.jpeg   \
               media/fig_Err_of_t_PLT.jpeg	 \
               media/fig_Err_of_t_Decomp_PLT.jpeg\
               media/fig_Err_of_t_Log_Decomp_PLT.jpeg \
               media/fig_Nat_occ_of_t_PLT.jpeg   \
               media/fig_Nat_occ_of_t_log_PLT.jpeg "

#for f in $files_to_plot
#do
#	rm -f $f
#	echo "Removing $f"
#done

# process args
case "$1" in
    clean|-cl*)   
for f in $files_to_plot
do
	rm -f $f
	echo "Removing $f"
done
  exit 0 ;;
    help|-h*)   
echo "Plots evolution (Re Im Time) of E, T, V, W, Err, and natural occupations "
echo "Usage: $mctdhb_dir/Scripts/evolution_E_V_T_W_Err_nat_occ  0 15 1 10 "
echo "where: from Time_Bgn=0 till Time_Fnl=15 dt=1 Morb=10"
echo "where: $mctdhb_dir is  the installation directory of the MCTDHB_V3* package"
  exit 0 ;;
    "")   
echo "Usage: $mctdhb_dir/Scripts/evolution_E_V_T_W_Err_nat_occ  0 15 1 10 "
echo "where: from T=0 till T=15 dT=1 Morb=10"
echo "where: $mctdhb_dir is  the installation directory of the MCTDHB_V3* package"
exit 0 ;
esac


gnuplot -e "Time_Bgn=$1; Time_Fnl=$2; dt=$3; Morb=$4"  $mctdhb_dir/Scripts/E_of_t.plt
gnuplot -e "Time_Bgn=$1; Time_Fnl=$2; dt=$3; Morb=$4"  $mctdhb_dir/Scripts/Err_of_t.plt
gnuplot -e "Time_Bgn=$1; Time_Fnl=$2; dt=$3; Morb=$4"  $mctdhb_dir/Scripts/Nat_occ_of_t.plt
 exit 0 
