       subroutine remove(pathNew)
c=================================================
!        CHARACTER delimiter
!        character, allocatable :: path
        CHARACTER*50 pathNew
        CHARACTER*56 command
        logical :: L_EXISTS
        integer:: stat
           
!       CALL getenv('DELIMITER',delimiter)
!        path='DATA/orb_1'
!         command='rm -f '//pathNew
!        write(command,'(A56)') 'rm -f ',path
!          write(6,*) command
!            pause
!gfortran does not support it        inquire(directory=//trim(path), exist=L_EXISTS)
!        if(.not.L_EXISTS) then
!gfortran does not support it        print *,"The directory ",C_DIRSPEC," does not exist"
!        CALL system(command)
!gfortran does not support it        else
!gfortran does not support it        write(6,*),"The directory ",path," ALREADY exist!"
!gfortran does not support it        endif
      open(unit=1234, iostat=stat, file=trim(pathNew), status='old')
      if (stat == 0) then
       write(6,*) trim(pathNew), " file was deleted"
       close(1234, status='delete') 
      endif
        return 
       end subroutine remove
c=================================================
