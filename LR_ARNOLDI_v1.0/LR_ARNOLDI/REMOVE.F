       subroutine remove(path)
c=================================================
!        CHARACTER delimiter
!        character, allocatable :: path
        CHARACTER*50 path
        CHARACTER*56 command
        logical :: L_EXISTS
           
!       CALL getenv('DELIMITER',delimiter)
!        path='DATA/orb_1'
!        write(command,'(A56)') 'rm -f ',path
         command='rm -f '//path
!          write(6,*) command
!            pause
!gfortran does not support it        inquire(directory=//trim(path), exist=L_EXISTS)
!        if(.not.L_EXISTS) then
!gfortran does not support it        print *,"The directory ",C_DIRSPEC," does not exist"
        CALL system(command)
!gfortran does not support it        else
!gfortran does not support it        write(6,*),"The directory ",path," ALREADY exist!"
!gfortran does not support it        endif
        return 
       end subroutine remove
c=================================================
