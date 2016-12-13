!       subroutine read_parse_V_W_Psi(stringV,stringW,stringPSI)
       subroutine read_parse_V_W_Psi
       USE   DVR_ALL
       use interpreter
!       USE   SHARED_DIMS
!       USE   rR_hW
       USE   CI_All
!       USE   W_INTERPARTICLE
       IMPLICIT NONE
!       INTEGER, DIMENSION(MaxUserCnf) :: UserCNfN !Added to read from V_W_Psi_string.in file
!       real*8, DIMENSION(MaxUserCnf) :: UserCNfW  !Added to read from V_W_Psi_string.in file
       integer ::  i,j,k,l,NFock,ind
       character*18 lname
       character*10 path
       character*11 var2
       character*1255  VAr1,Value1,READ_TMP
       integer, dimension(100) :: noccwork,noccwork1
       character (len=*) , parameter   ::    end_of_line = &     
                 char(13)//char(11)//char(0) 
!       character arrayValue1*(*)
        real*8  :: ww
          character*1255 stringCNF
          real*8  DNRM2
          external DNRM2
!DVR_ALL       character*255, DIMENSION(100)  :: stringPSI
!=================================================
!DVR_ALL       character(len = 10),  dimension(10) :: variables
!DVR_ALL       real*8,               dimension(10) :: variablesvalues
!DVR_ALL       character (len = 5)  :: statusflag
!=================================================
          path='.'
            lname="V_W_Psi_string.in"
!DVR_ALL          stringV=''
!DVR_ALL          stringW=''
!DVR_ALL          stringPSI=''
      open(unit=112,file=trim(path)//'/'//lname,form='formatted',status='old',err=102)
                    ind=1
            Do i=1,100
       read(112,"(A1255)",err=102,end=101) READ_TMP
!               j=VAr1,Value1
              j=index(READ_TMP,"#")
              if(j.eq.1) cycle
!          write(6,*)"I read", trim(VAr1), " and value ",trim(Value1)
              j=index(READ_TMP,"=",BACK = .TRUE.)
                Var1=READ_TMP(1:j)
              Value1=READ_TMP(j+1:len(READ_TMP))
          write(var2,'(a11)') trim(VAr1)
!          write(6,*) trim(VAr1)
!          write(6,*) var2
          select case(Var2)
              case('V(x_y_z&t)=')
                 stringV=trim(Value1)
      call init (stringV, variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed V(x,y,z)=',trim(stringV)
        write(6,*)'V(x,y,z) is from VTRAP_EXT_TD.F'
        stringV='Using defaults from VTRAP_EXT_TD.F'
      else
        write(6,*)statusflag,'in parsed V(x,y,z)=',trim(stringV)
      endif
      call destroyfunc()
              case('W(R=|r1-r2|')
                 stringW=Value1
      call init (stringW, variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed W(r1-r2)=',trim(stringW)
        write(6,*)'W(R) is from Get_InterParticle.F'
        stringW='Using Defaults from Get_InterParticle.F'
      else
        write(6,*)statusflag,'in parsed W(R)=',trim(stringW)
      endif
      call destroyfunc()
              case('Psi_1(x_y_z')
                 stringPSI(1)=Value1
      call init (stringPSI(1), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_1',trim(stringPSI(1))
        write(6,*)'Psi_1 is from  Guess_PSI.F'
        stringPSI(1)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_1=',trim(stringPSI(1))
      endif
      call destroyfunc()
              case('Psi_2(x_y_z')
                 stringPSI(2)=Value1
      call init (stringPSI(2), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_2',trim(stringPSI(2))
        write(6,*)'Psi_2 is from  Guess_PSI.F'
        stringPSI(2)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_2=',trim(stringPSI(2))
      endif
      call destroyfunc()
              case('Psi_3(x_y_z')
                 stringPSI(3)=Value1
      call init (stringPSI(3), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_3',trim(stringPSI(3))
        write(6,*)'Psi_3 is from  Guess_PSI.F'
        stringPSI(3)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_3=',trim(stringPSI(3))
      endif
      call destroyfunc()
              case('Psi_4(x_y_z')
                 stringPSI(4)=Value1
      call init (stringPSI(4), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_4',trim(stringPSI(4))
        write(6,*)'Psi_4 is from  Guess_PSI.F'
        stringPSI(4)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_4=',trim(stringPSI(4))
      endif
      call destroyfunc()
              case('Psi_5(x_y_z')
                 stringPSI(5)=Value1
      call init (stringPSI(5), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_5',trim(stringPSI(5))
        write(6,*)'Psi_5 is from  Guess_PSI.F'
        stringPSI(5)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_5=',trim(stringPSI(5))
      endif
      call destroyfunc()
              case('Psi_6(x_y_z')
                 stringPSI(6)=Value1
      call init (stringPSI(6), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_6',trim(stringPSI(6))
        write(6,*)'Psi_6 is from  Guess_PSI.F'
        stringPSI(6)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_6=',trim(stringPSI(6))
      endif
      call destroyfunc()
              case('Psi_7(x_y_z')
                 stringPSI(7)=Value1
      call init (stringPSI(7), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_7',trim(stringPSI(7))
        write(6,*)'Psi_7 is from  Guess_PSI.F'
        stringPSI(7)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_7=',trim(stringPSI(7))
      endif
      call destroyfunc()
              case('Psi_8(x_y_z')
                 stringPSI(8)=Value1
      call init (stringPSI(8), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_8',trim(stringPSI(8))
        write(6,*)'Psi_8 is from  Guess_PSI.F'
        stringPSI(8)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_8=',trim(stringPSI(8))
      endif
      call destroyfunc()
              case('Psi_9(x_y_z')
                 stringPSI(9)=Value1
      call init (stringPSI(9), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_9',trim(stringPSI(9))
        write(6,*)'Psi_9 is from  Guess_PSI.F'
        stringPSI(9)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_9=',trim(stringPSI(9))
      endif
      call destroyfunc()
              case('Psi_10(x_y_')
                 stringPSI(10)=Value1
      call init (stringPSI(10), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_10',trim(stringPSI(10))
        write(6,*)'Psi_10 is from  Guess_PSI.F'
        stringPSI(10)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_10=',trim(stringPSI(10))
      endif
      call destroyfunc()
              case('Psi_11(x_y_')
                 stringPSI(11)=Value1
      call init (stringPSI(11), variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =Psi_11',trim(stringPSI(11))
        write(6,*)'Psi_11 is from  Guess_PSI.F'
        stringPSI(10)='Using Defaults from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed Psi_10=',trim(stringPSI(11))
      endif
      call destroyfunc()
              case('Imprint_MOM')
                 stringMOM=Value1
      call init (stringMOM, variables, statusflag)
      if(statusflag.ne.'ok') then 
        write(6,*)statusflag,' detected in parsed =ImprintMOM',trim(stringMOM)
        write(6,*)'Imprint default momentum from Guess_PSI.F'
        stringMOM='Imprint defaults K from Guess_PSI.F'
      else
        write(6,*)statusflag,'in parsed ImprintMOM=',trim(stringMOM)
      endif
      call destroyfunc()
              case('Df_cnf_Fock')
              j=index(Value1,"|",BACK = .TRUE.)
!           print * ,"in Cnf | at pos:",j
          if(j.ne.0) Value1=Value1(j:len(Value1))
                stringCNF=Value1
              j=index(Value1,"|")
!           print * ,"in Cnf next | at pos:",j
          if(j.ne.1) cycle 
        write(6,*)'CI CNF detected in parsed:',trim(Value1)
              call strFock(Value1,NFock,ww)
                 UserCNfN(ind)=Nfock
                 UserCNfW(ind)=ww
         write(6,'(a25,i8)')"In Fock space its number:",NFock
!       call         GET_Ind_vs_ii(NFock,Npar,Morb,noccwork,noccwork1)
!            write(6,*) "I parsed as ",noccwork1(1:Morb),ww 
!            write(6,*) ind,"I parsed as ",UserCNfN(ind),UserCNfW(ind)
!            stop "fhjfsdhjkh CI"
                  ind=ind+1
              case default
     write(6,*) Var1,'- i cannot recognize it in V_W_Psi_string.in'
!                   STOP "I stop in V_W_Psi_string.in"
!                   return
           end select
                 enddo
101    write(6,*)"I have reached the end-of-file V_W_Psi_string.in"
             close(112)
            lname="V_W_Psi_string.out"
       open(unit=113,file=trim(path)//'/'//lname,form='formatted')
       write(113,'(a55,a1250)')"Trap V(x_y_z&t)== ",adjustL(stringV)
       write(113,'(a55,a1250)') &
      "Interaction W(R=|r1-r2|&t)== ",adjustL(stringW)
       write(113,'(a55,a1250)') &
      "Imprinted Momentum Psi(j)*exp i*(k1x+k2y+k3z)= ",adjustL(stringMOM)
           Do i=1,Morb
       write(113,'(a51,i2,a2,a1250)')"Psi_",i,"= ",adjustL(stringPSI(i))
           endDo
        UserCnfW=UserCnfW/DNRM2(MaxUserCnf,UserCnfW,1) !Normalization of the input weights
           Do i=1,ind-1
                 Nfock=UserCnfN(i)
                 ww=UserCnfW(i)
       call      GET_Ind_vs_ii(NFock,Npar,Morb,noccwork,noccwork1)
       write(113,'(a5,i10,F10.5,a3,25i4)')" cnf ",NFock,ww,"*|",noccwork1(1:Morb)
           endDo
       close(113)
!==================================================================================
       write(6,*)"OK parsing is done"
       write(6,*)"parse_V_W_Psi.in file is in  used i.e, V,W,Psi are from it!!!"
!              stop "OK parsing"
             close(112)
               return
102    write(6,*)"parse_V_W_Psi.in file is not used : either empty or corrupted"
             close(112)
!              stop ("not OK parsing")
                return
          END SUBROUTINE read_parse_V_W_Psi
         subroutine strFock(st,cind,W)
         USE CI_All
         IMPLICIT NONE
        integer      ::     i,strlen,j,k,l,ind,cind
        integer      ::     nn(100)
        character    ::     st*(*)
        character*1255    ::    stORG,strW
        real*8  :: W 
            nn=0
            stORG=st
        strlen = len(st)
              j=index(st,"|")
         if(j.eq.0) stop " CNF Wrong format! should istart from |"
         if(st(1:1).ne."|") stop "CNF Wrong format! should start from |"
              st=st(j+1:strlen)
              l=index(st,">*")
             ind=1
        do while (l.ne.0)
             k=index(st,":")
         if(k.ne.0) then 
                    read(st(1:k-1),*) nn(ind)
!         write(6,'(a5,i4)')"n =" ,nn(1)
                    st=st(k+1:strlen)
                ind=ind+1
                    else 
                    read(st(k+1:l-1),*) nn(ind)
!         write(6,*)"w =" ,st(l+2:strlen)
                    read(st(l+2:strlen),*)  w
                    st=st(1:1)
!         write(6,'(a5,i4,a10,f16.8)')"n =" ,n(1), "weight=",w
         endif 
              l=index(st,">*")
!              j=0
         !write(6,*)"k=",k,"L=",l
         enddo
              nn=abs(nn)
!             write(6,*) "Npar=", sum(nn),Npar
!             write(6,*) "Morb=", ind,Morb
!             write(6,*) "PARSE:",nn(1:ind), w
         if(Npar.ne.sum(nn)) then
        stop "Wrong Npar!=n_1+n_2+...+n_M |n_1:n_2:...:n_M>"
          endif
         if(Morb.ne.ind)  then 
        stop "Wrong Morb!=M |n_1:n_2:n_3...:n_M>"
          endif
!          Write(6,*)"Before GetInd= ", Npar,Morb,nn(1:Morb)
          call GetInd(Npar,Morb,nn,cind)
!          Write(6,*)"ConfN= ",cind
!          stop "jhljhjklH"
        end subroutine strFock
