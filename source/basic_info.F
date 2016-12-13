c      write a few numbers (energy and occupations of the last timestep 
c      writes in output file 'basic_info.out' some basic information about computation 
c      it can be used to compare computations done on diferent platforms...
!      subroutine basic_info(time,exec_time,nprocs,nthreads)

      subroutine basic_info(time,exec_tm,nprocs,nthreads)
      USE SHARED_DIMS
      USE rR_hW
      USE W_INTERPARTICLE
      implicit NONE

      integer i,nprocs,nthreads,ind
      integer jj,DVRX,DVRY,DVRZ
      Real*8 time,exec_tm
      character rhoname
         DVRX= Time_DVRMETHODX
         DVRY= Time_DVRMETHODY
         DVRZ= Time_DVRMETHODZ
       
      write(6,'(125a1)')("=",jj=1,118)
      write(6,'(A24,a63,a31)') '====================    ',
     ." basic_info.out file (with PATH, time, E, etc) has been written",
     .' ==================='
      write(6,'(125a1)')("=",jj=1,118)
      open(UNIT=102,FILE='basic_info.out',STATUS='UNKNOWN')
      call wr_info(102) 
      write(102,'(A8,I4,A6,I4,A6,I4,A12,I4,A8,I4,A8,I4,A8,I4)') 
     .' # NDX= ',NDX,' NDY= ',NDY,' NDZ= ',NDZ,' DIM-MCTDH= ',dim_mctdhb
     .,' DVR_X= ',DVRX,' DVR_Y= ',DVRY,' DVR_Z= ',DVRZ
      write(102,'(A43,F16.12,A12,F16.12,A23,F12.6)') 
     ." # MCTDHB-EOM have been solved from T_bgn= ", Time_Bgn,
     ." till T_fnl= ",time," with Time_print_step= ",Time_print_step
      write(102,'(A10,F20.4,A18,I4,A16,I4 )') 
     .' # T_cpu= ', exec_tm,' Num_OMP_threads= ',nthreads,
     .' Num_MPI_proc = ',nprocs

      write(102,'(A10,F25.16)') ' Energy = ', energy
      do i=1,morb
          ind=Morb-i+1
            write(102,'(a4,I2,a4,F25.16)')'  n(',i,') = ',Nocc(ind)
      end do


      close(102)
            
      end subroutine basic_info 

c=================================================
c============ Writes some info to output file specfied by fln 
       subroutine  wr_info(fln)
 
       USE    PASS_ARG
       USE   SHARED_DIMS
       USE   DVR_ALL
       USE   CI_ALL
       USE   W_INTERPARTICLE
 
       CHARACTER(len=255) :: pwd 
       character*20 jtp 
       integer fln 
       integer, dimension(3) :: today,now
 

       PreFac=JOB_PreFac
c========================================================================
      IF(PreFac.eq.(-1.0d0,0d0))write(jtp,*)"Job->Relax. Forward"
      IF(PreFac.eq.(+1.0d0,0d0))write(jtp,*)"Job->Relax.Backward"
      IF(PreFac.eq.(0d0,-1.0d0))write(jtp,*)"Job->Prop.  Forward"
      IF(PreFac.eq.(0d0,+1.0d0))write(jtp,*)"Job->Prop. Backward"
           CALL getcwd(pwd)
           call idate(today)   ! today(1)=day, (2)=month, (3)=year
           call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
!           print *, "Workinf dir:",pwd
       WRITE(fln,'(a36,a255)') " # This computation has been done in: "
       WRITE(fln,*)"# ",trim(pwd)
       WRITE(fln, 1 )  today(1), today(2), today(3), now 
       WRITE(fln,'(a67)') adjustl( banner(9))
       write(fln,'(a8,i3,a7,i8,a7,a20)')
     . " # Morb=",Morb," Npar=",Npar," Job=",jtp
       write(fln,'(a15,a100)')" # V(x_y_z&t)= ",adjustl(stringV)
       write(fln,'(a19,a100)')" # W(R=|r1-r2|&t)= ",adjustl(stringW)
       write(fln,'(a14,i2,a56,a12,F12.6,a23,l)')
     ." # Kind of W==", Wxx_TYPE,
     .  " [if ==0 W=delta(R) else W=f(R) see Get_InterParticle.F]"
     . ," lambda_0= ", xlambda0," Time-dependent? (T/F):",WTD
        return 
 1    format ( ' # Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
       end subroutine wr_info
c=================================================

