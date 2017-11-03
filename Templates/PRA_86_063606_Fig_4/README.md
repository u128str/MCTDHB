# MCTDHB template

Here are the input.in and V_W_Psi_string.in files needed to reproduce M=7 N=10 energy from Fig.4 of PRA_86_063606 (2012)
with calibration of the MCTDHB method on time-dependent HIM model.
For comparisions one has to plot MB results with the time-dependent trap and time-dependent inter-particle interaction
Hence there are 3 computational steps:
.step1 - to get initial state for many-body -- it is a ground stat  for N=10 with lambda_0=0.5
          make a working directory step1
         copy the input.in and V_W_Psi_string.in and exe -files
         execute  the job:       mpirun -n 7  ./boson_MCTDHB_gnu_FFTW
         compare with reference: vimdiff basic_info.out basic_info.out_Reference
.step2 - to propagate the MCTDHB equation with initial wave-function obtained at step1
          make a working directory step2
 BEFORE running second job one has to copy this wave-function stored in  CIc_bin and PSI_bin files
         if you are already in step2 type:
         copy the respective input.in and V_W_Psi_string.in exe -files and basic_info.out_Reference
         cp ../step1/*bin .
         run the propagation:  mpirun -n 7  ./boson_MCTDHB_gnu_FFTW  
         compare with reference: vimdiff basic_info.out basic_info.out_Reference
.step3 - to solve one-particle problem in time-dependent trap potential 
          make a working directory step3
          copy the respective input.in and V_W_Psi_string.in exe -files and basic_info.out_Reference from the Template
         run the propagation:  mpirun -n 1  ./boson_MCTDHB_gnu_FFTW  #ATTN: Only one processor is needed
         compare with reference: vimdiff basic_info.out basic_info.out_Reference
.   Plot with gnuplot:
         plot "step2/NO_PR.out" u 1:9, "step3/NO_PR.out" u 1:3 w l 
         plot "step2/NO_PR.out" u 1:($9-14.924811557), "step3/NO_PR.out" u 1:($3-14.924811557) w l
         Here 14.924811557=D/2*(N-1)\delta_N shift of the total energy see PRA_86_063606 for details
