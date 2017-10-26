        subroutine READ_init_LR_Arnoldi
        use LR_ARNOLDI_MOD
        implicit none
        integer :: nargs,ist,j,k
        CHARACTER(LEN=255) :: string,string0,srting1  ! for command line arguments
        CHARACTER(LEN=255) :: message           ! for I/O error messages


!========================================================================        
         nargs=command_argument_count()
         write(6,*)" Number of command line argumenst", nargs



!========================================================================        
11       continue
       if(nargs==0)then
        write(6,*) "Trying to Read properties.in file"
        open(8,file="lr_arnoldi.in",
     &         status='OLD',recl=80, delim='APOSTROPHE')
c        read(8,NML=ZERO_body) 
        read(8,NML= LR)
        write(6,*) "Read of properties.in file is OK"
        close(8)
        endif

        end subroutine  READ_init_LR_Arnoldi


