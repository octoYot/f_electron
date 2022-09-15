PROGRAM f_electrons !   LAST MODIFIED:  1.53.5  15/09/22
USE f_e_data
USE f_e_fit
USE f_e_input
USE f_e_wigner
USE f_e_calculate
USE f_e_group
USE f_e_edipole
USE f_e_parameters

      IMPLICIT none
      REAL*4 Time(11)
      CHARACTER(LEN=4), PARAMETER :: command_keyw(43)=                          &
     (/"AOM ","ALTP","BASE","BLOC","CF  ","CONF",  "CONS","ECHO","EDIP","END ", &
       "EXP1","EXPB","EXPE","EXPG","FEXT","FIT ",  "FITO","GEXS","GRID","INTR", &
       "JUDO","LF1E","LINK","MDIP","MFLD","OFFS",  "OPTN","OUTP","PLT1","PLT2", &
       "PRED","PRNG","REDF","ROTL","SPIN","SYML",  "TIME","TITL","VARI","VECT", &
       "WARN","WORD","XREF"/)
      CHARACTER(LEN=4), PARAMETER :: command_keyw3(10)=                          &
     (/"AOMX","CCF ","EXPL","EXPM","EXPT", "EXSO","FAST","LANC","PRT4","RORB"/)
      CHARACTER(LEN=4), PARAMETER :: debug_keyw(7)=                             &
     (/"CHK1","CHK2","FNCR","PRT1","PRT2",  "PRT3","TEST"/)
      CHARACTER(LEN=4), PARAMETER :: parameter_keyw1(20)=                       &
     (/"ALPH","BETA","GAMM","EAVE","F2  ",  "F4  ","F6  ","M0  ","M2  ","M4  ", &
       "P2  ","P4  ","P6  ","T2  ","T3  ",  "T4  ","T6  ","T7  ","T8  ","ZETA"/)
      CHARACTER(LEN=4), PARAMETER :: parameter_keywCFR(15)=                     &
     (/"B20 ","B21 ","B22 ","B40 ","B41 ",  "B42 ","B43 ","B44 ","B60 ","B61 ", &
       "B62 ","B63 ","B64 ","B65 ","B66 "/)
      CHARACTER(LEN=4), PARAMETER :: parameter_keywCFI(15)=                     &
     (/"    ","B21'","B22'","    ","B41'",  "B42'","B43'","B44'","    ","B61'", &
       "B62'","B63'","B64'","B65'","B66'"/)
      CHARACTER(LEN=6) parameter_keywAOM(25)
      data parameter_keywAOM    &
         /"Esig1 ","Esig2 ","Esig3 ","  :   ","Esig20",  &
          "Epi1  ","Epi2  ","Epi3  ","  :   ","Epi20 ",  &
          "theta3","theta2","theta3","  :   ","thet20",  &
          "phi1  ","phi2  ","phi3  ","  :   ","phi20 ",  &
          "chi1  ","chi2  ","chi3  ","  :   ","chi20 "/
!
!   f_electrons 1.5x  
!   ----------------
!
!  This program calculates the energy of the d-d or f-f transitions.
!  Reduced matrix elements are read from the files "d_electron.dat" and "f_electron.dat".
!  Point group and symmetry operation information are read from: "group.dat"
!  
!  The ligand field can be specified in terms of crystal field parameters
!  or Angular Overlap Model (AOM) parameters.
!
!  NOTE: The program is dimensioned for any d- or f-electron configuration.
!
!  Uses: CPU_TIME                  F95 Standard 
!        GET_ENVIRONMENT_VARIABLE  Fortran 2003 standard
!        GET_COMMAND_ARGUMENT        "
!        GETCWD                    GNU extension 
!        SIGNAL                    GNU extension 
!
      INTEGER*4 i,k,icalc,ivar,minNv,n,n1,n2,values(8),UserLevel/3/  ! 1,2,3 = Basic,Advanced,Expert
      CHARACTER*32 ARG,ARG1,ARG2,ARGS(4)
      character xval*60,yval*60,yval1*15,CWD*100
      LOGICAL  noGrpDat, plt2/.false./,test

!      external endProg
!
      do i=1,11; Time(i)=0.0; enddo
      call CPU_TIME(Time(1))
      call DATE_AND_TIME(VALUES=values)
      write( *,'(/," F-ELECTRON  version:",A6,"; compile date:",A,"; this run:(",         &
                                                   I2,"/",I2.2,"/",I4,") ",I2,":",I2.2)') &
                 version,trim(CompileDate),values(3),values(2),values(1),values(5),values(6)
      INFILE="INPUT.DAT"
      OUTFIL="OUTPUT.DAT"      
      call GET_ENVIRONMENT_VARIABLE("F_E_PATH",F_E_PATH)
      call GETCWD(CWD)
      write(*,'(" F_E_PATH=",A)') trim(F_E_PATH)
      write(*,'(" Current Working Directory=",A)') trim(CWD)
!
      call SIGNAL(2,endProg)     !  minGW GNU extension  2=SIGINT 
!
      DO i=1,4
        args(i)=" "
!        write(*,'(" Command line",I2)') i
        CALL get_command_argument(i, arg)
!        write(*,'(" Argument:",A32)') arg
        IF (LEN_TRIM(arg) == 0) EXIT
        args(i)=TRIM(arg)
!        write(*,'(" Command line",I2,":",A32)') i,args(i)
      END DO !  i
!      write(*,'(" TRIM(args(1)):",A32)') adjustl(args(1))
!      write(*,'(" LEN_TRIM(args(1)):",I2)') LEN_TRIM(args(1))
!      write(*,'(" TRIM(args(2)):",A32)') adjustl(args(2))
!      write(*,'(" LEN_TRIM(args(2)):",I2)') LEN_TRIM(args(2))
      IF (LEN_TRIM(args(1))/=0) then
        arg1=args(1); CALL TO_UPPER(arg1)
        if (adjustl(arg1)=="VERSION" .or. adjustl(arg1)=="VER" .or. adjustl(arg1)=="V" .or.  &
            adjustl(arg1)=="-VERSION" .or. adjustl(arg1)=="-VER" .or. adjustl(arg1)=="-V" ) then
          STOP
!*        elseif (adjustl(arg1)=="LICENSE" .or. adjustl(arg1)=="LIC") then
!*          call checkLic(1)
!*          STOP
        elseif (adjustl(arg1)=="HELP") then
          arg2=args(2); CALL TO_UPPER(arg2)
          if (LEN_TRIM(arg2)==0) then
!             write(*,'(" about to call HELP with:",A32)') adjustl(arg1)
             call HELP(adjustl(arg1))
          else 
!             write(*,'(" about to call HELP with:",A32)') adjustl(arg2)
            call HELP(adjustl(arg2))
          endif
          STOP
        elseif (adjustl(arg1)=="HIDE" .and. UserLevel.gt.2) then
          arg2=args(2); CALL TO_UPPER(arg2)
          if (LEN_TRIM(arg2)==0) then
             call HIDE("HIDE")
          else 
            call HIDE(adjustl(arg2))
          endif
          STOP
        else
          INFILE=adjustl(args(1))
          if (LEN_TRIM(args(2))/=0) OUTFIL=adjustl(args(2))
        endif
      endif        
      OPEN(IO,FILE=OUTFIL,STATUS='UNKNOWN')
      OPEN(IN,FILE=INFILE,STATUS='OLD',ERR=900)    
!      
      write(io,'(" F-ELECTRON  version:",A6,"; compile date:",A,"; this run:(",   &
                                                   I2,"/",I2.2,"/",I4,") ",I2,":",I2.2)') &
                 version,trim(CompileDate),values(3),values(2),values(1),values(5),values(6)
!      call checkLic(0)
      call input(INFILE)
      call CPU_TIME(Time(2))
!
      icalc=0
      Time(3)=Time(2);
 !**** Judd Ofelt calculations     
      if (JuddOfelt) then
        CALL calcEDipJO()
        Call output(0)
        call CPU_TIME(Time(4))
        if (nfit.gt.0) then
          call prepareFitJO()
          CALL doFit()
          Call output(0)
          call CPU_TIME(Time(4))
          write( *,'(13x,"J Multiplet Fit Calc:",F8.2," secs")') Time(4)-Time(2)
        else 
          write( *,'(13x,"J Multiplet calculation:",F8.2," secs")') Time(4)-Time(2)
        endif
        goto 990   ! skip main calculation
      endif  
!**** |SLJ> basis calculations
      if (.not.FullMJBasis .and. IMD.ne.0) then  ! ie BASE "J" and MDIP commands
        CALL calcMDipJBasis() ! Calculation in |SLJ> basis.
        Call output(0) 
        call CPU_TIME(Time(4))
        write( *,'(13x,"J Multiplet Mag.Dipole calculation:",F8.2," secs")') Time(4)-Time(2)
        goto 990   ! skip main calculation
      endif    
!   
!      if (ESO2J.ne.0) then  
!        CALL calc_ESO() ! Do Extended Stevens operators in |J> basis.
!        if (ESOskip.eq.1) return
!      endif
!
!**** PRED prep calculation
      if (nPRED.gt.0) then
!        write(*,'("About to: doPRED1")')
        call doPRED1()  ! Atomic parameters only using the full |LSJMJ> basis; save nPRED eigenvectors
      endif
!
!**** Fit calculation
      if (nFit.gt.0) then
!        if (fitType.eq.4) then
!          CALL calc_ESO() ! One calc so ESOcmat(i,j) is loaded. 
!        endif  
        call prepareFit()
        CALL doFit()
        call PrintParameters(2)
! FitType
!    1   EXPE    Energies/Intensities  NEXP
!    2   EXPE    Judd-Ofelt            NEXP  <--ToDo, this will never happen goto 990 above
!    3   EXPB    Bkq data              NEXP
!    4   EXPM    ESO matrix elements    -
!    5   EXPG    g-values              NexpG
!    6   EXP1    one electron LF MEs   1/2*N*(N+1)   N = 2L+1
        if (fitType.eq.1 .or. fitType.eq.2 .or. fitType.eq.5) then 
          call prepCalc(0)
          call doCalc(1) !  do one last calc with fitted parameters
        elseif (fitType.eq.6) then
          call OUTPUT_EXP1()        
        endif
        goto 990 
      else ! nFit=0
        if (fitType.eq.6) then
          call OUTPUT_EXP1()        
          goto 990 
        endif
      endif

!**** Multiple nvar calculations
      if (nvar.ne.0) then
!        write( io,'("nvar.ne.0")')
        if (grid) then ! Make grid of all possible combination of variables. 
                       ! Number of calculations will be n1*n2*n3..  
!          do ivar=1,nvar;  P(400+ivar)=allVarVals(ivar,1);  enddo
          do ivar=1,2
!            write(io,'("ivar=",i4,";  varNum(ivar)=",I4,";  nVarVals(ivar)=",i4)') ivar,varNum(ivar),nVarVals(ivar)
!            write(io,'("allVarVals(varNum(ivar),n)=",12F8.2)') (allVarVals(varNum(ivar),n),n=1,nVarVals(varNum(ivar)))
            if (varNum(ivar).lt.1 .or. varNum(ivar).gt.nvar) then
              write( *,'(" Error in GRID command")')
              write(io,'(" Error in GRID command")')
              stop
            endif
          enddo  
          do n1=1,nVarVals(varNum(1))
            do n2=1,nVarVals(varNum(2))
              icalc=icalc+1
              P(401)=allVarVals(varNum(1),n1)
              P(402)=allVarVals(varNum(2),n2)
              WRITE(IO,'(/,"Calculation number:",I3)') icalc       
              write(io,'("n1,n2=",2i4,";  P(401),P(402)=",2F8.2)') n1,n2,P(401),P(402)
              call prepCalc(0)      
              if (OUTP(8).and.icalc.ne.1) call PrintParameters(1)
              call doCalc(icalc)
            enddo
          enddo
        else   !  Step variables together; number of calculation will be min(n1,n2,n3..)
          minNv=Nvarvals(1); do i=1,nVar; minNv=min(Nvarvals(i),minNv); enddo
          do n=1,minNv
            icalc=icalc+1
            do ivar=1,nvar;  P(400+ivar)=allVarVals(ivar,n);  enddo
            WRITE(IO,'(/,"Calculation number:",I3,"; Variables:",10(A,":",F10.2,", "))') icalc,  &
                       (trim(variables(ivar)),allVarVals(ivar,n),ivar=1,nvar) 
!      write(*,'("About to call prepCalc")') 
            call prepCalc(0)
            if (OUTP(8).and.icalc.ne.1) call PrintParameters(1)
            call doCalc(icalc)
            do k=1,10; if (PLOTout2(k)) plt2=.true.; enddo
            if (plt2) then
              yval=" "
              WRITE(xval,'(10E12.6)') (allVarVals(ivar,n),ivar=1,nvar)
              do i=1,10
                if (PlotOut2(i)) then
                  WRITE(yval1,'(E12.6)') Yplot2(i) 
!                  write(io,'("yval1=",A)') yval1
                  yval=trim(yval)//" "//trim(yval1)
                endif
!                write(io,'("yval=",A)') yval
              enddo 
              write(ip2,'(A,2x,A)') trim(xval),trim(yval)
            endif  
!      write(*,'("Back from call doCalc")') 
          enddo                           
        endif
      else  
!        write( io,'("ELSE")')
!        write( io,'("About to: doCalc")')
!**** One-off calculation      
        call doCalc(0)   ! nvar = 0; no Variables, do single calculation
      endif
      goto 990
!      
  900 write(*, '("***FATAL: Input file not found:",A40)') INFILE
      write(IO,'("***FATAL: Input file not found:",A40)') INFILE     
      STOP
!
  990 call CPU_TIME(Time(11))
!
      test=Option(2).or.Option(3).or.Option(4).or.Option(6).or.Option(7)                          
      write( *,'(65("-"),/,          "      Timing:            input:",F8.2," secs")') Time(2)-Time(1)
      if (icalc.eq.0) then
                          write( *,'("               diagonalisation:",F8.2," secs")') Time(3)-Time(2)
        if (test)         write( *,'("   WF spin & Free ion analysis:",F8.2," secs")') Time(4)-Time(3)
        if (idg.ne.0)     write( *,'("               Symmetry labels:",F8.2," secs")') Time(5)-Time(4)
        if (GEXS(1).ne.0) write( *,'("                      g-values:",F8.2," secs")') Time(6)-Time(5)
        if (IMD.ne.0)     write( *,'("      Mag. Dipole calculations:",F8.2," secs")') Time(7)-Time(6)
        IF (IED.ne.0)     write( *,'("      Elec Dipole calculations:",F8.2," secs")') Time(8)-Time(7)
        IF (n_expect.ne.0)write( *,'("Expectation value calculations:",F8.2," secs")') Time(9)-Time(8)
      elseif (icalc.gt.0) then
        write( *,'("calculations:",i3,"   diagonalise:",F8.2" secs",/,  &
                   "                          post:",F8.2," secs")') icalc,Time(3)-Time(2),Time(11)-Time(3) 
      endif 
      write( *,'(65("-"),/,          "                total CPU time:",F8.2," secs")') Time(11)-Time(1) 
!
      if (iTime) then
        write(io,'(/,115("-"),/,       "  Timing:                input:",F8.2," secs")') Time(2)-Time(1)
        if (icalc.eq.0) then
                            write(io,'("                   diagonalise:",F8.2," secs")') Time(3)-Time(2)
          if (test)         write(io,'("   WF spin & Free ion analysis:",F8.2," secs")') Time(4)-Time(3)
          if (idg.ne.0)     write(io,'("               Symmetry labels:",F8.2," secs")') Time(5)-Time(4)
          if (GEXS(1).ne.0) write(io,'("                      g-values:",F8.2," secs")') Time(6)-Time(5)
          if (IMD.ne.0)     write(io,'("      Mag. Dipole calculations:",F8.2," secs")') Time(7)-Time(6)
          IF (IED.ne.0)     write(io,'("      Elec Dipole calculations:",F8.2," secs")') Time(8)-Time(7)
          IF (n_expect.ne.0)write(io,'("Expectation value calculations:",F8.2," secs")') Time(9)-Time(8)
        endif 
      endif
      write(io,'(16x,"total CPU time:",F8.2," secs",/,115("-"))') Time(11)-Time(1)
!             
      STOP
      CONTAINS
!
!-----------------------------------------------------------------------
!
      subroutine checkLic(ptime) 
      IMPLICIT NONE
      integer i,j,VALUES(8),ptime,OUT1(64),s1(64),KEY(64),dd,mm,yy
      integer*8 T,date_now,date_lic,times
      CHARACTER DATE1*8, TIME1*10, ZONE1*5, PASS*8, PASS16*16 
      CHARACTER C1(8)*1, C2(8)*2
      LOGICAL ex

      inquire(file=trim(F_E_PATH)//'\license.dat', exist=ex)
      if (.not.ex) then; write(io,'("License file ''license.dat'' missing")')
                         write(*,'("License file ''license.dat'' missing")'); stop; endif
                         
      OPEN(Ilic,FILE=trim(F_E_PATH)//'\license.dat',STATUS='OLD')
      READ(Ilic,'(8A2)') (C2(i),i=1,8)
!      write(*,'("   read from file:  C2=",8A2)')(C2(i),i=1,8)
      call CONVERT(T,S1,C2,C1,4)    ! C2->S1
      call CONVERT(T,S1,C2,C1,2)    ! S1->T
!      write(*,'("   read from file:   T=",i24)') T
      do i=1,64; KEY(i)=0; enddo; KEY(33)=1; KEY(14)=1; KEY(4)=1;
      CALL DES(S1,KEY,1,OUT1)    !  decryption  S1 -> OUT1
      call CONVERT(T,OUT1,C2,C1,2)  ! OUT1->T
!      write(*,'("decoded from file:   T=",i24)') T
      times=T/100000000; date_lic=T-times*100000000; 
!      write(*,'("times=",I24,";  license date:=",i24)') times,date_lic
!
      call DATE_AND_TIME(DATE1, TIME1, ZONE1, VALUES)
              !  year           month          day    
      date_now=VALUES(1)*10000+VALUES(2)*100+VALUES(3)
!
      if (ptime.eq.0) then
        i=date_lic-date_now
        if (i.lt.0) then
          write(*,'("****FATAL: Invalid license.")')
          write(*,'("License date:",I8,"; Current date:",I8)')date_lic,date_now
          stop
        elseif (i.gt.10600) then  ! this is 18 months
          write(*,'("****FATAL: License has expired.")')
          write(*,'("License date:",I8,"; Current date:",I8)')date_lic,date_now
          stop
        elseif (times.lt.0 .or. times.ge.2000) then
          write(*,'("****FATAL: Invalid license number.")')
          write(*,'("times used:",I8,"; Current date:",I8)') times
          stop
        else
          yy=date_lic/10000; mm=(date_lic-yy*10000)/100; dd=mod(date_lic,100)
          write(*,'("Licensed until:",i2,"/",I2,"/",I4,"; program has run:",i3," times" )') dd,mm,yy, times
          times=times+1
!          times=1
          T=date_lic+times*100000000 
          call CONVERT(T,S1,C2,C1,1)   ! T->S1
          CALL DES(S1,KEY,0,OUT1)   !  encryption  S1 -> OUT1
          call CONVERT(T,OUT1,C2,C1,3) ! OUT1->C2
          rewind(Ilic)
          write(Ilic,'(8A2)') (C2(i),i=1,8)
        endif
      ELSE
 10     write(*,'("enter password")')
        read(*,'(A8)') pass
        pass=adjustr(pass)
        do i=1,8; C1(i)=pass(i:i); enddo
        call CONVERT(T,S1,C2,C1,5)   ! C1->C2
        call CONVERT(T,S1,C2,C1,4)   ! C2->S1
        CALL DES(S1,KEY,0,OUT1)    !  encryption  S1 -> OUT1
        call CONVERT(T,OUT1,C2,C1,3)   ! OUT1->C2
        PASS16=C2(1)//C2(2)//C2(3)//C2(4)//C2(5)//C2(6)//C2(7)//C2(8)
        if (PASS16.ne."5D202DEA12DE24F5") then
          write(*,'("invalid password")')
          goto 10
        endif  
        write(*,'("License date:",I8,"; Current date:",I8)')date_lic,date_now
 20     write(*,'("license to dd,mm,yyyy?")')
        read(*,*) dd,mm,yy 
        date_lic=yy*10000+mm*100+dd  
        if (date_lic.lt.date_now) then
          write(*,'("License date:",I12," is less than current date:",I8)')date_lic,date_now
          goto 20 
        endif
        date_lic=date_lic+100000000  ! reset to use 1
        call CONVERT(date_lic,S1,C2,C1,1)   ! date_lic->S1
        CALL DES(S1,KEY,0,OUT1)    !  encryption  S1 -> OUT1
        call CONVERT(T,OUT1,C2,C1,2)   ! OUT1->T
        call CONVERT(T,OUT1,C2,C1,3)   ! OUT1->C2
!        write(*,'("****NEW LICENSE(before encrypt):",I24)') date_lic
!        write(*,'("****NEW LICENSE (after encrypt):",I24)') T
        write(*,'("****NEW LICENSE (hex):",8A2)') (C2(i),i=1,8)
        rewind(Ilic)
        write(Ilic,'(8A2)') (C2(i),i=1,8)
      endif
      close(Ilic)
!  
      return
      end subroutine checkLic
!     
!-----------------------------------------------------------------------
!
      subroutine doCalc(icalc) 
      IMPLICIT NONE
      integer i,j,neng,icalc
!     
      if (option(8)) then ! writing matrix skipped and diagonalisation skipped
        return 
      endif   
      if (nPRED.eq.0) then
        call checkLFcomplex()
        do i=1,nBlocks
          CALL buildAndDiagonalise(i,0)
        enddo
!      call calcEDip()  ! TODO delete temp
        if (option(9)) return  ! diagonalisation skipped
        if (nBlocks.ne.1) call orderE_V(calcVecs)
!  Offset
        if (offsetType.eq.2) then
          DO i=1,n_matrix
            do j=1,offsetN
             if (i.ge.ioffset(j,1) .and. i.le.ioffset(j,2)) Engs(i)=Engs(i)+Roffset(j)
            enddo
          enddo  
          if (nBlocks.ne.1) call orderE_V(calcVecs)
        endif 

      else  ! nPRED.ne.0
!        write( io,'("About to: doPRED2")')
        call doPRED2()  ! Ligand Field parameters only using the full |LSJMJ> basis; save LF 
        call buildAndDiagonalisePRED()
      endif

      call CPU_TIME(Time(3))
      
      neng=n_matrix
      if (nPRED.ne.0) neng = nPRED
      if (fastMat(2) .and. .not.fastMat(3)) neng=neng/2 ! One componant of Kramers doublet
      if (fastMat(2) .and.      fastMat(3)) call DoubleKramers()  !
      if (calcVecs) then
        if (Option(2).or.(IDG.ne.0 .and. ISO.eq.0)) call calcSpinPR(neng)  ! spin projections 
        if (SpAllow(2)-SpAllow(1).gt.0) call CalcSpinAllowed(neng)
        If (Option(3) .or. Option(4) .or. Option(6) .or. Option(7)) call FreeIonPR1(neng) ! free ion %  or MJ %
        call CPU_TIME(Time(4))
        if (idg.gt.0)     then; call GROUPL(neng); call CPU_TIME(Time(5)); else; Time(5)=Time(4); endif
        if (GEXS(1).ne.0) then; call calcGval();   call CPU_TIME(Time(6)); else; Time(6)=Time(5); endif
        if (IMD.ne.0)     then; call calcMDip();   call CPU_TIME(Time(7)); else; Time(7)=Time(6); endif
        IF (IED.ne.0)     then; call calcEDip();   call CPU_TIME(Time(8)); else; Time(8)=Time(7); endif
        IF (n_expect.ne.0)then; call calcExpect(); call CPU_TIME(Time(9)); else; Time(9)=Time(8); endif
      endif
      ! ToDo the following 2 CALLs will never happen because of the GOTO 990 for JuddOfelt 
      if (JuddOfelt) CALL calcEDipJO() ! Judd-Ofelt calculation doesn't need eigenvectors of full calculation
      if (.not.FullMJBasis) CALL calcMDipJBasis() ! Calculation in |SLJ> basis.
      Call output(neng)
!  
      return
      end subroutine doCalc
!     
!-----------------------------------------------------------------------
!
      SUBROUTINE help(key)
      IMPLICIT NONE
      integer i,j,n1,lamda,t,p
      CHARACTER*4 KEY 
!      write(*,'("Subroutine help, key:",A4)') key
      write(*,'(/)')
      select case (KEY)
      case ("HELP")
        write(*,'(/,"Use:",/,">f_electron ",/,">f_electron infile outfile",/,           &
                    ">f_electron help",/,">f_electron help keyword")')
        write(*,'(/,"If infile or outfile are not specified, the program assumes the",/, &
                  " default names: ''input.dat'' and ''output.dat''.",//,                &
                  "The structure of the input.dat file should be as follows:",/,         &
                  "keyword   value1  value2 ..",/,                                       &
                  "The ''keyword'' must appear in the first 4 characters and ''values'' after this.",/, &
                  "For this reason, it is best to avoid the use of tabs.",/,             &
                  "The remainder of the line contains one or more values, depending on the keyword.",/, &
                  "The keywords can be Commands or Parameters",//,                       &
                  "Commands:",/,(10(A4,4X)))') (command_keyw(i),i=1,43)
        write(*,'("Parameters: (all in units of cm-1)",/,(10(A4,4X)))') (parameter_keyw1(i),i=1,20)  ! ,(parameter_keywCF(i),i=1,27)
        write(*,'(6X,"To get the ligand field parameters type:",/,                       &
                        "    ''f_electron help CF'' when using Bkq parameters",/,        &
                        "or: ''f_electron help AOM'' when using AOM parameters")')
        write(*,'(/,"Type: >''f_electron help kword'' for more information about a particular keyword ''kword''.")')
!
! DEBUGGING help
!    (/"CHK1","CHK2","FNCR","PRT1","PRT2",  "PRT3","TEST"/)
!  plus: "DEBU"

      CASE ("DEBU")
        write(*,'("Debugging keywords:",/,(10(A4,4X)))') (debug_keyw(i),i=1,7)
        write(*,'("Type: >''f_electron help kword'' for more information about a particular  keyword ''kword''.")')      
      CASE ("CHK1") ! 1-11: Basis,Un,V11,EE,abg, MnPn,Tn,AOM,Eigen,Pt Groups 
        write(*,'("CHK1  L1-L10",/,6X,"T/F for checks to be made on:",/,6X,"1) Basis",/,                         &
                  6X,"2) Un matrices",/,   6X,"3) V11 matrices",/, 6X,"4) EE matrices",/,6X,"5) abg matrices",/, &
                  6X,"6) Mn/Pn matrices",/,6X,"7) Tn matrices",/,  6X,"8) AOM matrices",/,6X,"9) Eigenvalues",/, &
                  6X,"10) Point Groups")') 
        write(*,'(6X,"The extra information is put into the file: ''debug.dat''")') 
      CASE ("CHK2") ! 1-2: skip Diag, check Kramers
        write(*,'("CHK2  L1-L10",/,6X,"T/F for checks to be made on:",/,                        &
                  6X,"1) Only write matrices, skip the diagonalisation.",/,                     &
                  6X,"2) Check the wavefunctions form correct Kramers doublets (odd e- only)",/,&
                  6X,"3) Check the L,S matrices give L.S SOC matrix",/,                         &
                  6X,"4) Check the conversion between Altp and Blki parameters",/,              &
                  6X,"5-10) Unused")') 
        write(*,'(6X,"The extra information is put into the file: ''debug.dat''")') 
      CASE ("FNCR") ! 1-8: use errors from: Un  V11  Fn  abg  MnPn  Tn T2  Fnmp
        write(*,'("FNCR L1-L8",/,6X,"T/F for using the errors in the fncross files for the matrices:")')
        write(*,'(6X,"1) Un",/,6X,"2) V11",/,6X,"3) Fn",/,6X,"4) alpha,beta,gamma",/,6x,"5) Mn/Pn",/,  &
                  6X,"6) Tn",/,6X,"7) bad T2",/,6X,"8) Fnmp file used")') 
        write(*,'(6X,"This is useful to compare with published data")') 
        write(*,'(6X,"''badT2'' is the wrongly using the -T2 for the complementary configuration")') 
        write(*,'(6X,"''Fnmp file used'' does not have all the fncross errors")') 
      CASE ("PRT1") ! 1-10: Uk,V11,Fn,Mn,Pn,Sn,Tn in the |SL> basis
        write(*,'("PRT1 L1-L10",/,  &
                  6X,"T/F for printing the matrices in the |SL> basis to the file ""matrices.txt""")')
        write(*,'(6X,"1) Uk matrix",/,  6X,"2) V11 matrix",/, 6X,"3) Fn matrices",/, 6X,"4) Mn matrices",/,                &
                  6X,"5) Pn matrices",/,6X,"6) Sn matrices",/,6X,"7) Tn matrices",/, 6X,"8) Gik matrices",/,               &
                  6X,"9) Unused",/,6X,"10) Selection of matrix printed. (in which case a line i1,i2, j1,j2 must follow)",/, &
                  9X,"rows i1-i2 and columns j1-j2 are printed")') 
      CASE ("PRT2") ! 1-10: matrix,(*/.)matrix,Uk(JO)    in the |SLJ> basis
        write(*,'("PRT2 L1-L10",/,  &
                  6X,"T/F for printing the matrices in the |SLJ> basis to the file ""matrices.txt""",/, &
                  6X,"Unless otherwise stated the maximum size printed is (20x20)")')
        write(*,'(6X,"1) Full Matrix",/, 6X,"2) Full matrix as array of (*/.) (max 100x100)",/,                     &
                  6X,"3) Uk matrix",/,   6X,"4) <L>/<S> matrices",/,  6X,"5) <L+2S> matrix",/,6X,"6)-9) Unused",/,  & 
                  6X,"10) Selection of matrix printed. (in which case a line i1,i2, j1,j2 must follow)",/,           &
                  9X,"rows i1-i2 and columns j1-j2 are printed")') 
      CASE ("PRT3") ! 1-10: matrix,(*/.)matrix,Ukq,L/S,kL+gS  in the |SLJMJ> basis
        write(*,'("PRT3 L1-L10",/,  &
                  6X,"T/F for printing the matrices in the |SLJMJ> basis to the file ""matrices.txt""",/, &
                  6X,"Unless otherwise stated the maximum size printed is (20x20)")')
        write(*,'(6X,"1) Full Matrix",/,                           6X,"2) Full matrix as array of (*/.) (max 100x100)",/, &
                  6X,"3) <|Ukq|> matrices (non-zero, max 100)",/,  6X,"4) Ukq matrices (full, max 14x14)",/,              & 
                  6X,"5) <|ED|> matrices",/,                       6X,"6) <|L|> and <|S|> matrices",/,                    &
                  6X,"7) <|kL+gS|> matrices",/,                    6X,"8) LF matrix in |SLJMJ> and PRED basis",/,         &
                  6X,"9) Unused",/,6X,"10) Selection of matrix printed. (in which case a line i1,i2, j1,j2 must follow)",/,&
                  9X,"rows i1-i2 and columns j1-j2 are printed")') 
      CASE ("TEST") ! 
        write(*,'("TEST",/,  & 
                  6X,"If given, the fit to experiment GOF=sqrt(delta^2/N) is appended to the file test.out.",/, &
                  6X,"This command is intended for use with a series of calculations run from a batch job,",/,  &
                  6X,"in which the results are kept as a record. Also records the version, compile date, ",/,   &
                  6X,"run date and input filename. Can be used to test whether a new program version still ",/, &
                  6X,"gives the same results as previous versions.",/,                                          &
                  6X,"If the GOF is > TestCrit, then this is flagged in the output.",/,                         &
                  6X,"(Default value of TestCrit = 1.0, but can be changed using the INTR command)")')
! PARAMETERS help
      CASE ("ALTP")
        write(*,'("ALTP  N1 N2",/,                                                     &
             6X,"When this command is given, the intensity parameters follow.",/,      &
             6X,"N1=1 Altp parameters on the next 9 lines immediately following.",/,   &
             6X,"   2 Blkr parameters on the next 15 lines",/,                         &
             6X,"   3 Tki  parameters on the next NLIGAND lines",/,                    &
             6X,"   4 Cki  parameters on the next NLIGAND lines",/,                    &
             6X,"N2=1,2,3  parameters: pure real,pure imag,complex")')
        write(*,'("Enter N1:")')
        read(*,*) N1
        if (N1.eq.1) then     
          write(*,'(6X,"These Altp parameters are in units of 10^-11 cm and are as defined by",/,      &
                    6X,"Reid&Richardson,J.Phys.Chem.,88,3579,(1984).")')
          write(*,'(10X,2(2X,"A21",i1))') (i-1,i=1,2)
          write(*,'(10X,3(2X,"A22",i1))') (i-1,i=1,3)
          write(*,'(10X,4(2X,"A23",i1))') (i-1,i=1,4)
          write(*,'(10X,4(2X,"A43",i1))') (i-1,i=1,4)
          write(*,'(10X,5(2X,"A44",i1))') (i-1,i=1,5)
          write(*,'(10X,6(2X,"A45",i1))') (i-1,i=1,6)
          write(*,'(10X,6(2X,"A65",i1))') (i-1,i=1,6)
          write(*,'(10X,7(2X,"A66",i1))') (i-1,i=1,7)
          write(*,'(10X,8(2X,"A67",i1))') (i-1,i=1,8)
        else if (N1.eq.2) then
          write(*,'(6X,"See Burdick et al, Phys.Rev.B, 59, R7789,(1999).",/,   &
                    6X,"for the relationship between these and the Altp parameters.")')
          do i=2,6,2 
            do j=0,i
              write(*,'(10X,"B",2I1,"x,  B",2i1,"y,  B",2i1,"z")') i,j,i,j,i,j
            enddo
          enddo
        else if (N1.eq.3) then
          write(*,'(6X,"See Brown et al, Mol.Phys., 64, 771,(1988).",/,   &
                    6X,"for a discussion on these parameters with respect to d-d transitions.")')
          write(*,'(10X,"There will be NLIGAND lines i=1,NLIGANDS given by:")') 
          write(*,'(10X,"PTiSig  FTiSig  RTiSig   PTiPiX  FTiPiX  RTiPiX   PTiPiY  FTiPiY  RTiPiY (for d-d transitions)")') 
          write(*,'(10X,"DTiSig  GTiSig  RTiSig   DTiPiX  GTiPiX  RTiPiX   DTiPiY  GTiPiY  RTiPiY (for f-f transitions)")') 
        else if (N1.eq.4) then
          write(*,'(6X,"See Brown et al, Mol.Phys., 64, 771,(1988).",/,   &
                    6X,"for a discussion on these parameters with respect to d-d transitions.")')
          write(*,'(10X,"There will be NLIGAND lines i=1,NLIGANDS given by:")') 
          write(*,'(10X,"Ci1x  Ci2x  Ci1y  Ci2y  Ci1z  Ci2z  Ci3z (for d-d transitions)")') 
          write(*,'(10X,"Ci1x  Ci2x  Ci1y  Ci2y  Ci1z  Ci2z  Ci3z (TODO)(for f-f transitions)")') 
        endif
      CASE ("EAVE"); write(*,'("EAVE   If OPTN(1)=T, calculated energies are relevant to the lowest level",/, &
                               "       (The lowest level will then be 0.0)",/,                                &
      -                        "       If OPTN(1)=F, calculated energies are - EAVE",/,                       &
                               "       If you wish the absolute energies, use OPTN(1)=F and EAVE=0")')
      CASE ("F2  ","F4  ","F6  ")
        write(*,'("F2, F4, F6   are the electron repulsion parameters.")')
        write(*,'("For d-electron calculations the Racah parameters B,C can be used instead of F2,F4.")')        
      CASE ("ZETA"); write(*,'("ZETA is the spin orbit coupling parameter")')
      CASE ("ALPH","BETA","GAMM"); write(*,'("ALPH,BETA,GAMM are the alpha,beta,gamma parameters.")')
      CASE ("M0  ","M2  ","M4  "); write(*,'("M0, M2, M4 parameters.")')
      CASE ("P2  ","P4  ","P6  "); write(*,'("P2, P4, P6 parameters.")')
      CASE ("T2  ","T3  ","T4  ","T6  ","T7  ","T8  "); write(*,'("Tn parameters.")')
!      
      CASE ("AOMP")
        write(*,'("The parameter ''standard order'' when AOM ligand field parameterization is used:")')
        write(*,'(5(4X,I2," : ",A4,2X))') (i,parameter_keyw1(i),i=1,20)
        do i=1,3
        write(*,'(4(4X,I2," : ",A6))') 20+i,parameter_keywAOM(i),   40+i,parameter_keywAOM(5+i),  &
          60+i,parameter_keywAOM(10+i),80+i,parameter_keywAOM(15+i),80+i,parameter_keywAOM(20+i)
        enddo
        write(*,'(4(4X,2X," : ",A6))')      parameter_keywAOM(4),        parameter_keywAOM(5+4),  &
               parameter_keywAOM(10+4),     parameter_keywAOM(15+4),     parameter_keywAOM(20+4)
        write(*,'(4(4X,I2," : ",A6))') 35+i,parameter_keywAOM(5),   55+i,parameter_keywAOM(5+5),  &
          75+i,parameter_keywAOM(10+5),95+i,parameter_keywAOM(15+5),95+i,parameter_keywAOM(20+5)
      CASE ("CFPA")
        write(*,'("The parameter ''standard order'' when Bkq ligand field parameterization is used:")')
        write(*,'(5(4X,I2," : ",A4,2X))') (i,parameter_keyw1(i),i=1,20)
        do i=1,15
          write(*,'(2(4X,I2," : ",A4,2X))') 20+i,parameter_keywCFR(i),40+i,parameter_keywCFI(i)
        enddo
      CASE ("INTP")
        write(*,'("The parameter ''standard order'' when Altp Intensity parameters are used:")')
        write(*,'(8(3X,I4,2X,"A21",i1))') (200+i,i-1,i=1,2)
        write(*,'(8(3X,I4,2X,"A22",i1))') (202+i,i-1,i=1,3)
        write(*,'(8(3X,I4,2X,"A23",i1))') (205+i,i-1,i=1,4)
        write(*,'(8(3X,I4,2X,"A43",i1))') (209+i,i-1,i=1,4)
        write(*,'(8(3X,I4,2X,"A44",i1))') (213+i,i-1,i=1,5)
        write(*,'(8(3X,I4,2X,"A45",i1))') (218+i,i-1,i=1,6)
        write(*,'(8(3X,I4,2X,"A65",i1))') (224+i,i-1,i=1,6)
        write(*,'(8(3X,I4,2X,"A66",i1))') (230+i,i-1,i=1,7)
        write(*,'(8(3X,I4,2X,"A67",i1))') (237+i,i-1,i=1,8)
!
!  COMMANDS help
!     (/"AOM ","ALTP","BLOC","CF  ","CONF",  "CONS","ECHO","EDIP","END ","EXP1", &
!       "EXPB","EXPE","EXPG","FEXT","FIT ",  "FITO","GEXS","GRID","INTR","JUDO", &
!       "LF1E","LINK","MDIP","MFLD","OFFS",  "OPTN","OUTP","PLOT","PRED","PRNG", &
!       "REDF","ROTL","SPIN","SYML","TIME",  "TITL","VARI","VECT","WARN","WORD", &
!       "XREF"/)
! plus "FORM","NLIG"
      CASE ("AOM")
        write(*,'("AOM  N1",/,   &
             6X,"If this command is given the ligand field is specified by the Angular Overlap Model (AOM).",/)')
        write(*,'(6X,"|N1|: is the number of ligands",/, & 
             6X,"      there will be |N1| lines following this command specifying each ligand position.",/, &
             6X,"      If N1>0, then each ligand is specified by:   Name, Esig,EpiX,EpiY, phi,theta,chi",/, &
             6X,"      If N1<0, then each ligand is specified by:   Name, Esig,EpiX,EpiY,   x,    y,  z")')
        write(*,'(6X,"Name: only first 4 characters are used.",/,                                     &
             6X,"Esig,EpiX,EpiY: are the AOM bonding parameters (units cm-1).",/,                     &
             6X,"x,y,z: are cartesian coordinates of the ligands (units can be anything,",/,          &
             6X,"       used to calculate angles).",/,                                                &
             6X,"theta,phi,chi: are the position specified by the AOM angles (degrees);",/,           &
             6X,"               (Note: chi need not be given, a default value of 0 will be used)",/,  &
             6X,"               see Schaffer: Struct&Bond,1,1968.")')
      CASE ("BASE") 
        write(*,'("BASE T1",/,                                                                             &
             6X,"If T1 = ""J"" then the reduced |SLJ> basis is used instead of the full |SLJMJ> basis.",/, &
             6X,"(If T1!= ""J"" then the full basis is used. )")')
      CASE ("BLOC") 
        write(*,'("BLOC [N1]",/,                                                                       &
             6X,"Symmetry blocks the matrix to be diagonalised.",/,                                    &
             6X,"N1 is optional and if not equal 0, the crystal quantum numbers are printed for each level",/, &
             6X,"The blocking is according to the smallest q>0 in the non-zero Bkq values.",/,         &
             6X,"Beware if the magnetic field is non-zero.",/,                                         &
             6X,"Eigenvectors are reconstructed in the original basis before any subsequent calculations.")')
      CASE ("CF  ")
        write(*,'("CF",/,6X,"The ligand field is to be specified by the crystal field parameters Bkq (units cm-1).")')
        write(*,'(6X,"These Bkq parameters must immediately follow the CF command.")')
        write(*,'(6X,"The Bkq follow the definitions of Goeller-Walrand,",/,                       &
                  6X,"''Handbook of Physics & Chemistry of Rare Earths'',Vol 23, Ch155, 1996")')
        write(*,'(6X,"There will be 27 Bkq values for Bkq k=2,4,6; q=0,k; (note: No B00) given in 15 lines:")')
        write(*,'(16X,"B20",/,                                                                     &
             16X,"B21   B21''",/,16X,"B22   B22''",/,                                              &
             16X,"B40",/,                                                                          &
             16X,"B41   B41''",/,16X,"B43   B42''",/,16X,"B43   B43''",/,16X,"B44   B44''",/,      &
             16X,"B60",/,                                                                          &
             16X,"B61   B61''",/,16X,"B62   B62''",/,16X,"B63   B63''",/,16X,"B64   B64''",/,      &
             16X,"B65   B65''",/,16X,"B66   B66''")')
      CASE ("CONF") 
        write(*,'("CONF  N1  N2",/,                             &
             6X,"N1=2/3 for d/f electron calculation.",/,       &
             6X,"N2=number of electrons, 1-9 for d and 1-13 for f electrons.")')
      CASE ("CONS") 
        write(*,'("CONS cName T1_R1",/,                                                                  &
             6X,"A constant called cName can be either text (defined by an expression)",/,               &
                " or a number (defined by a value).",/,                                                  &
             6X,"It can also contained other constants, as long as they have been defined previously.",/,&
             6X,"Multiple CONS commands can be given.")')
        write(*,'(6X,"Examples:",/,"CONS  sr2    sqrt(2.)",/,"CONS  x1    3.22)",/,"CONS  pi   4.0*atan(1.0)",/,  &
                     "CONS   trigonal  acos(1/sqrt(3))*180/pi")')        
      CASE ("ECHO") 
        write(*,'("ECHO L1",/,6X,"T/F for the input file to be echoed to the output file as each line is read.")')
      CASE ("EDIP")
        write(*,'("EDIP N1 N2 N3 N4 N5 ",/,                                                          &
             6X,"Calculation of the intensities of electric dipole transitions ",/,                  &
             6X,"in units of Debye^2 between a set of initial and final levels.",//,                 &
             6X,"n1,n2 is the range of the initial levels (inclusive)",/,                            &
             6X,"n3,n4 is the range of the final levels (inclusive)",/,                              &
             6X,"There is an automatic summation over any initial or final state degeneracies.",/,   &
             6X,"Note: n1-n4 refer to the actual (not summed) energy levels",//,                     &
             6X,"  If n5=1, the calculation is made in cartesian coordinates.",/,                    &
             6X,"  If n5=2, the calculation is also made for left (AL) and right (AR) circularly",/, & 
             6X,"           polarised light, MCD (dA= AL – AR) and MLD  (dA= Ax – Ay).",/,           &
             6X,"           Note it is assumed that the direction of the incident light is ",/,      &
             6X,"           k||z; the magnetic field is H||z and H||x for MCD and MLD respectively.",/,&
             6X,"  If n5=3, the magnetic dipole amplitudes (not intensities) are given, ",/,         &
             6X,"           including relative phases. (Note that in this case there is no ",/,      &
             6X,"           summing over degeneracies.)",/,                                          &
             6X,"The EDIP command requires that the Altp intensity parameters to be set. ",/,        &
             6X,"Type: ''f_electron help ALPT'' for more details")')
      CASE ("END "); write(*,'("END",/,6X,"Starts the calculation. Any commands after this will be ignored.")')
      CASE ("EXP1") 
        write(*,'("EXP1  ",/,                                                                              &
                  6X,"Flags that calculated matrix elements of the one-electron LF matrix are fitted ",/,  &
                  6X,"against a provided set. The complete lower triangle must given."                     &
                  6X,"This will be 7(5) lines with a total of 28(15) MEs for f(d) systems resp.",/,        &
                  6X,"The weight should be given after each value.",/,                                     &
                  6X,"f11 w11 ",/,"f21 w21  f22 w22",/,":",/,                                              &
                  6X,"f71 w71  f72 w72  f73 w73  f74 w74  f75 w75  f76 w76  f77 w77 ")')                  
      CASE ("EXPB") 
        write(*,'("EXPB  N1",/,6X,"There will be N1 lines of crystal field Bqk values following this command")')
        write(*,'(6X,"Each line will contain Value, Bkq",/,                                      &
             6X,"Value: the value to be fitted (cm-1)",/,                                        &
             6X,"Bkq:   the Bkq or Bkq'' CF coefficient to be fitted (as a 3 or 4 char text).")')
      CASE ("EXPE") 
        write(*,'("EXPE  N1 N2 [N3]",/,                                                          &
             6X,"There will be N1 lines of experimental energies following this command",/,      &
             6X,"They will contain Energy (N2=1) or Energy and Intensity (N2=2) data.",/,        &
             6X,"The optional N3 is for fitting the energy levels using additional criteria.")')
        write(*,'(/,"N2=1; Each line will contain Eng, Wgt, Nass1,Nass2...",/,                   &
             6X,"Eng:  the experimental energy (cm-1)",/,                                        &
             6X,"Wgt:  a weighting factor for Eng (cannot be negative)",/,                       &
             6X,"Nass(1-20): the number of the calculated levels that Eng is to be assigned to.")')
        write(*,'(/,"N2=2; Each line will contain Eng, WgtE, Int, WgtI, Nass1,Nass2...",/,       &
             6X,"Eng:  the experimental energy (cm-1)",/,                                        &
             6X,"WgtE: a weighting factor for Eng (cannot be negative)",/,                       &
             6X,"Int:  the experimental intensity",/,                                            &
             6X,"WgtI: a weighting factor for Int (cannot be negative)",/,                       &
             6X,"Nass(1-20): the numbers of the calculated levels that the experimental",/,      &
             6X,"            transition (Eng,Int) is to be assigned to.",/,                      &
             6X,"The intensities will be in the units specified in either the JUDO or EDIP commands",/, &
             6X,"Up to 20 calculated levels can be assigned to a single experimental level Eng.")')
        write(*,'(/,"N3=0; (Default) Use energy ordering to make assignments.",/,                &
             6X,"The given level will be assigned to the Nass1,Nass2,... levels.",/,             &
             "If N3 not equal to 0, then Nass1,Nass2 must be given for each level",/,            &
             "N3=1; Use crystal quantum numbers to make assignments.",/,                         &
             6X,"The given level will be assigned to the Nass1 level of the Nass2 block.",/,     &
             "N3=2; Use Spin (S) to make assignments (spin-orbit coupling must be zero)",/,      &
             6X,"The level will be assigned to the Nass1 level of spin 2*S+1=Nass2.",/,          &
             6X,"Note: This only makes sense for calculations on Transition Metal ions.",/,      &
             "For N3=1/2, only a single calculated level can be assigned to the given level.")')
      CASE ("EXPG") 
        write(*,'("EXPG  N1",/,6X,"There will be N1 lines of g-values following this command")')
        write(*,'(6X,"Each line will contain g1,w1, g2,w2, g3,w3, Nassign",/,                    &
             6X,"g1,g2,g3: the g-values in ascending order",/,                                   &
             6X,"w1,w2,w3: the relative weighting of these values in the fit.",/,                &
             6X,"Nassign: the lower calculated level of the Kramers doublet that these g-values are assigned to.")')
      CASE ("FEXT")
        write(*,'("FEXT fext")')
        write(*,'(6X,"fext: text (<22 chars) that is appended to the name of debug, matrices, plot files.",/,   &
             6X,"      This is to ensure these files are not overwritten in a batch job.")')
      CASE ("FIT ")
        write(*,'("FIT  Nfit maxIts FitTol")')
        write(*,'(6X,"Nfit:   The number of parameters to be varied in the fit.",/,               &
             6X,"maxIts: The maximum number of iterations in the fit.",/,                         &
             6X,"FitTol: The fit will stop when maxIts reached or the fitting function (Chi^2) falls below FitTol.")')
        write(*,'(6X,"There must be Nfit lines following this command which should be of the form:",/,  &
             6x,"   Npar   Pmin  Pmax ",/,                                                        &
             6x,"Npar: the number of the parameter to be fitted (in the standard order).",/,      &
             6x,"[Pmin,Pmax]: The range of values that this parameter is restricted to.")')
        write(*,'(6X,"To get the parameter standard order type: ''f_electron help CFPA'' when using Bkq parameters")')
        write(*,'(6X,"                                      or: ''f_electron help AOMP'' when using AOM parameters")')
        write(*,'(6X,"                                      or: ''f_electron help INTP'' when using Altp intensity parameters")')
      CASE ("FITO")
        write(*,'("FITO L1,..L8")')
        write(*,'(6X,"L1: The RESHAPing fit matrices printed.",/,                                 &
                  6X,"L2: The Covariance matrix will be printed.",/,                              &
                  6X,"L3: The Minimized function and parameters printed on each iteration.",/,    &
                  6X,"L4: The fitted parameters are written to the file ''fit.out''",/,           &
                  6X,"L5: The Fit matrix written to the file ''matrices.dat'' each iteration.",/, &
                  6X,"L6: The parameter limits will be treated as hard; otherwise ",              &
                         "100*(P-Plimit)^2 is added to the penalty function.",/,                  &
                  6X,"L7: The fit will sum over all the degeneracies.",/,                         &
                  6X,"L8: The fit will also give the mean deviation: Sum(|Eexp-Ecal|)/Nexp.")')
      CASE ("FORM")
        write(*,'("Formulae can be used to specify the values of any of the parameters.",/,       & 
                  6X,"These formula can take the form of an algebraic expression built from",/,   &
                  6X,"  Operators: +,-,*,/,^,unary(-)",/,                                         &
                  6X,"  Functions: sin,cos,tan,asin,acos,atan,sinh,cosh,tanh,",/,                 &
                  6X,"             sind,cosd,tand,log,log10,nint,anint,aint,exp,sqrt,abs",/,      &  
                  6X,"  Numbers:",/,                                                              &
                  6X,"  Variables:  defined with the command ''VARI''",/,                         &
                  "Note: at the time that the calculation occurs (after the ''END'' command is given),",/, &
                  "all the Variables used in any formula must be defined, otherwise an error occurs.")')
      CASE ("GEXS"); write(*,'("GEXS N1 N2 N3",/,                                                 &
                               6X,"N1: if #0, the principle axes of the g-tensor will be calculated.",/,   &
                               6X,"Calculates the g-values between pairs of states in the range N2:N3",/,  &
                               6X,"Must be an even number of values (N3-N2+1 even) given for this command.")')
      CASE ("GRID")
        write(*,'("GRID N1 N2    ")')
        write(*,'(6X,"The variables number N1 and N2 are varied in multiple calculations",/,                                      &
                  6X,"Consider:",/,"VARI variable1  R11 R12 R13 R14",/,"VARI variable2  R21 R22 R23 R24",/,                       &
                  6X,"If GRID 1 2 command is given, then each variable can vary independently, and there are 16 calculations.",/, &
                  6X,"If GRID 1 2 has nor been given, then each variable steps together, and there are 4 calculations.",/,        &
                  6X,"Note if GRID is not given, there must be the same number of each of the variables.",/,                      &
                  6X,"If the numbers of variables are different the lower number is taken.")')
      CASE ("INTR")
        write(*,'("INTR T1 R1    ")')
        write(*,'(6X,"Used to change the value of the internal parameter given by T1 from the default value to R1",/,       &
                  6X,"Valid values for T1(default value):")')
        write(*,'(6X,"FImax   (",F7.4,"): Threshold fraction of a Free Ion term to appear in the WF output (OPTN(4,5))",/,  &
                  6X,"MJmax   (",F7.4,"): Threshold fraction of an MJ component to appear in the WF output (OPTN(6,7))",/,  &
                  6X,"CFbit   (",F7.4,"): Threshold belwo which a CF element is considered zero.",/,                        &
                  6X,"Edegen  (",F7.4,"): Threshold for two energy levels to be considered degenerate.",/,                  &
                  6X,"PFactor (",F7.4,"): The factor by which the output to matrix.dat is multiplied.",/,                   &
                  6X,19X,                "(Only needed if one of PRTn is true)",/,                                          &
                  6X,"TestCrit(",F7.4,"): The criteria that will be used to flag when a test is not met.")')                &
                                                                   FImax,MJmax,CFbit,Edegen,PFactor,TestCrit
      CASE ("JUDO")
        write(*,'("JUDO N1 N2 N3 N4 N5 Dielectric Omega2 Omega4 Omega6",/,   &
             6X,"Calculation of the electric dipole induced intensities of the Multiplets",/, &
             6X,"(not the individual levels) using the the Judd-Ofelt model.",/,              &
             6X,"n1,n2:  the range of the initial levels (inclusive)",/,                      &
             6X,"n3,n4:  the range of the final levels (inclusive).",/,                       &
             6X,"n5: (1-4) Calculate intensities in units of f, Debye, Area, all 3.",/,       &
             6X,"Note: that this range refers to the levels number of the multiplets",/,      &
             6X,"      (level number after summing over the degeneracies",//,                 & 
             6X,"Dielectric is the Dielectric factor and ",/,                                 &
             6X,"Omega2,Omega4,Omega6 are the 3 J-O parameters:",/,                           &
             6X,"The calculation is between multiplets averaged over all polarisations.",/,   &
             6X,"There is an automatic summation over any initial or final state degeneracies.")')
      CASE ("LF1E") 
        write(*,'("LF1e",/,6X,"There will be 5 (7) lines following this command for d(f) electrons",/,     &
                  6X,"Each line will contain a line of the lower triangle of the 1-e LF matrix. ",/,       &
                  6X,"in terms of the real d(f) orbitals.",/, 6X,"d11",/,6X,"d21  d22",/,                  &
                  6X,"d31  d32  d33",/,6X,"d41  d42  d43  d44",/,6X,"d51  d52  d53  d54  d55",/,           &
                  6X,"The standard order is that given in: Schaffer: Struct&Bond,1,1968.")')
      CASE ("LINK")
        write(*,'("LINK  N1 N2 A B")')
        write(*,'(6X,"The parameter number N1 is related to the parameter number N2 according to: ")')
        write(*,'(6X,"       P(N1) = A*P(N2) + B")')
        write(*,'(6X,"Note: 1) Multiple LINK commands can be given.")')
        write(*,'(6X,"      2) The ''LINKing'' will be carried out in the order given.")')
        write(*,'(6X,"To get the parameter standard order type: ''f_electron help CFPA'' when using Bkq parameters")')
        write(*,'(6X,"                                      or: ''f_electron help AOMP'' when using AOM parameters")')
      CASE ("MDIP")
        write(*,'("MDIP N1 N2 N3 N4 N5 N6 R1",/, &
             6X,"Calculation of the intensities of magnetic dipole transitions ",/,                  &
             6X,"between a set of initial and final levels.",//,                                     &
             6X,"N1,N2 is the range of the initial levels (inclusive)",/,                            &
             6X,"N3,N4 is the range of the final levels (inclusive)",/,                              &
             6X,"There is an automatic summation over any initial or final state degeneracies.",/,   &
             6X,"Note: n1-n4 refer to the actual (not summed) energy levels",/,                      &
             6X,"N5: Calculation type",/,                                                            &
             6X,"N5=1, the calculation is made in cartesian coordinates.",/,                         &
             6X,"N5=2, the calculation is also made for left (AL) and right (AR) circularly",/,      & 
             6X,"      polarised light, MCD (dA= AL – AR) and MLD  (dA= Ax – Ay).",/,                &
             6X,"      Note it is assumed that the direction of the incident light is ",/,           &
             6X,"      k||z; the magnetic field is H||z and H||x for MCD and MLD respectively.",/,   &
             6X,"N5=3, the magnetic dipole amplitudes (not intensities) are given, ",/,              &
             6X,"      including relative phases. (Note that in this case there is no ",/,           &
             6X,"      summing over degeneracies.)",//,                                              &
             6X,"N6: Calculation units",/,                                                           & 
             6X,"N6=1-3: Bohr magneton^2, ×10-7 Debye^2, ×10-24 cm^2 ",/,                            &
             6X,"        (if N5=3): Bohr magneton, ×10-3 Debye, ×10-12 cm ",/,                       &
             6X,"R1:   A multiplication factor to account for solvent effects.",/,                   &
             6X,"      (if N5=3) Set to 1.0.",/,                                                     &
             6X,"The orbital reduction parameters, k, will influence these calculations.")')
      CASE ("MFLD"); write(*,'("MFLD Hx Hy Hz",/,                              &
                               6X,"An applied magnetic field (in Tesla).",/,   &
                               6X,"Note that a non-zero Hy implies a longer calculation (complex matrix)." )')
      CASE ("NLIG")
        write(*,'("The comamnd NLIG has been changed to AOM in version 1.48")')
      CASE ("OFFS")
        write(*,'("OFFS  OffsetType OffsetN",/,                                                               &
                  6X,"OffsetType = 1 for diagonal atomic basis |LSJ>. PRED must also be used)",/,             &
                  6X,"OffsetN = the number of multiplets (N1=1) or individual levels (N1=2) to be offset.",/, &
                  6X,"There will be OffsetN lines of offset values following this command")')
        write(*,'(/,6X,"For OffsetType=1, Each line will contain: N1 Eng ",/,                                 &
                  6X,"N1 in the multiplet in the diagonalized atom only |LSJ> basis in energy order",/,       &
                  6X,"For OffsetType=2, Each line will contain: N1 N2 Eng ",/,                                &
                  6X,"N1,N2 is the range of levels in the full diagonalization (in energy order) to be",/,    &
                  6X,"shifted by the same energy Eng.")')
      CASE ("OPTN") ! 1-10: EAVE?,exp data?,exp int?,spin proj?,free-ion%?,  Hss?,MJ%,skip main?,g-prin.axes?,NU
        write(*,'("OPTN  L1-L10",/6X,"Each T/F to the following:",/,        &
             6X,"1) Calculated levels relative to lowest energy?",/,        &
             6X,"2) Calculate spin projections?",/,                         &
             6X,"3) Calculate free ion %? (short)",/,                       &
             6X,"4) Calculate free ion %? (long)",/,                        &
             6X,"5) Do NOT include the Hss spin-spin interactions?(use M0,M2,M4)",/,   &
             6X,"6) Calculate MJ%? (short)",/,                              &
             6X,"7) Calculate MJ%? (long)",/,                               &
             6X,"8) Skip the main calculate?",/,                            &
             6X,"9) g-value principal axes?",/,                             &
             6X,"10) not used")')
        write(*,'(6X,"Note that L2,3,4,6,7,9 requires more time (eigenvectors must be calculated).")')
      CASE ("OUTP") ! L1-L10: Ligand Dist(AOM),Equiv.CF(AOM),AOM mat,Print eigenvectors,Before/after ROTL,Ignore Deg,Size info,All Par,NU,NU
        write(*,'("OUTP  L1-L10",/6X,"Each T/F to the following:",/,               &
             6X,"1) Print ligand distance/angle matrix (AOM parameterization)",/,  &
             6X,"2) Print the equivalent CF parameters (AOM parameterization)",/,  &
             6X,"3) Print one electrom AOM matrix (AOM parameterization)",/,       &
             6X,"4) Print equivalent intensity parameters, Blki if Altp given, etc.?",/,  &
             6X,"5) If ROTL used, print CF before & after rotation",/,             &
             6X,"6) Print all energies, ignoring the degeneracies?",/,             &
             6X,"7) Print info on matrix size, non-zero elements?",/,              &
             6X,"8) Print all parameters for each calculation in multiple calcs",/,&               
             6X,"9) Print AOM orbital energies & eigenvectors",/,                  &               
             6X,"10)  currently not used.")')
      CASE ("PLOT") 
        write(*,'("PLOT  This command removed in ver 1.52, use PLT1 or PLT2 instead.")')
      CASE ("PLT1") ! L1-L10: PLT1: Eng,Spin fract,Spin Int,MJ fract,MDIP,  EDip, Fit,LF par,g-val,OrbE.
        write(*,'("PLT1  L1-L10",/6X,"The following is written for each T/F:",/,  &
             6X,"1) Energies ",/,                                   &
             6X,"2) Spin fractions ",/,                             &
             6X,"2) Spin-allowed intensities ",/,                   &
             6X,"4) MJ fractions ",/,                               &
             6X,"5) Magnetic Dipole intensities ",/,                &
             6X,"6) Elect.Dipole (JO) intensities ",/,              &
             6X,"7) Parameters (during a Fit) ",/,                  &
             6X,"8) Ligand Field Bkq parameters ",/,                &
             6X,"9) g-values ",/,                                   &
             6X,"10) Orbital Energies.")')
      CASE ("PLT2") ! L1-L10: PLT2: MinFunc,(L2-10 not used).
        write(*,'("PLT2  L1-L10",/6X,"The following is written for each T/F:",/,  &
             6X,"1) Minimized function ",/,                                       &
             6X,"2)-10) Not used")')
      CASE ("PRED")
        write(*,'("PRED  N1",/,6X,"The ligand field will be evaluated in the basis of the N1 lowest states",/,    &
                 6X," found in an atomic only calculation. The Eigenvectors found from the atomic calculation",/, &
                 6X," are used as basis functions to create an N1 x N1 block of the ligand field matrix elements.")')
      CASE ("PRNG")
        write(*,'("PRNG  N1 N2 N3",/,6X,"Control how many levels are printed.",/,   &
             6X,"N1: 1/2 Energy(cm-1)/number range.",/,                             &
             6X,"If N1=1, levels with energies in range [N2,N3] are printed.",/,    &
             6X,"If N1=2, levels in range [N2,N3] are printed (N2 must be < N3)")')
      CASE ("REDF")
        write(*,'("REDF  Kx Ky Kz",/,6X,"The orbital reduction factors. Must be 0<=k<=1.",/,  &
             6X,"Used for d calculations. Determines the extent that the L is reduced in ",/, &
             6X,"calculation of the g-values or the magnetic dipole transitions.",/,          &
             6X,"Note: It does not reduce the spin-orbit coupling constant in the calculation.")')
      CASE ("ROTL")
        write(*,'("ROTL  R1 R2 R3",/,6X,"The rotation of the ligand field. R1,R2,R3 re the Euler angles.",/,  &
             6X,"Since all the other parameters are spherically symmetric, rotating the ligand field",/,      &
             6X,"is the same as rotating the molecule. Note that the global axes do not change.",/,           &
             6X,"(The global axes are used to assign the symmetry labels (SYML), and polarization directions.")')
      CASE ("SPIN")
        write(*,'("SPIN  N1,N2",/,6X,"The initial states ([N1-N2] inclusive) from which a calculation for",/, &
             6X,"spin-allowed character is made. The [N1, N2] states themselves will be set to 1.0",/,        &
             6X,"The excited state degeneracies will be summed.")')
      CASE ("SYML")
        write(*,'("SYML  T1 [N2] [N3]")')
        write(*,'(6X,"T1 is the label of the point group. Valid upper or lower case labels given below")')
        write(*,'(6X,"N2 = 0/1 for Mulliken/Gamma notation (appropriate for spin-orbit coupling zero/non-zero)")')
        write(*,'(6X,"N3 = 0/1 if projections onto each character is made")')
        write(*,'(6X,"N2,N3 are both optional; if left off N2 = 0/1 for spin-orbit coupling zero or non-zero,")')
        write(*,'(6X,"                                     N3 = 0")')
        write(*,'(6(3X,I2,2X,A3,3X))') (I,Groups(I),I=1,37)
        write(*,'(6X,"The character tables of Koster, Dimmock, Wheeler& Schatz are used:",/,  &
             6X,"Properties of the 32 Point groups, MIT Press, (1963).",/,                    &
             6X,"Note that using this can be very time consuming for some f configurations.")')
      CASE ("TIME")
        write(*,'("TIME",/,6X," If given more detailed CPU timing of calculation is given.")')
      CASE ("TITL")
        write(*,'("TITL  Title of calculation (<120 characters)",/,             &
             6X,"This line can be repeated multiple times.",/,                  &
             6X,"It is intended to be used for comments or documentation",/,    &
             6X,"It is echo into the output file." )')
      CASE ("VARI")
        write(*,'("VARI  variablename R1, R2, R3,...R20")')
        write(*,'(6X,"A variable called variablename is defined and takes the values: R1, R2, R3....",/,  &
             6X,"A calculation can be made for each different variable value.",/,                    &
             6X,"If a parameter (or geometry) is defined in terms of this variable,",/,              &
             6X,"this can be used to generate a Tanabe-Sugano type diagram.",/,                      &
             6X,"Multiple commands of VARI can be given to define multiple variables, or multiple",/,&
             6X,"commands of VARI with the same variablename can be given to define more values.",/, &
             6X,"The total number of values a particular variable can have is 100.",/,               &
             6X,"An alternative to giving each variable value explicitly, you can use:",/,           &
             "VARI  variablename FROM R1 TO R2 STEP R3",/,                                           &
             "VARI  variablename FROM R1 STEP R2 FOR N1")')
      CASE ("VECT")
        write(*,'("VECT  Nlow,Nhigh,N",/,                                                       &
             6X,"Print the [1:N] coefficients of the eigenvectors in the range[Nlow:Nhigh] ",/, &
             6X,"The coefficients are printed in the basis ""standard order"". ",/,             &
             6X,"If N is not specified, all coefficients will be printed.")')
      CASE ("WARN") ! 1-10: check Bkq, check Alpt, 3-10 NU
        write(*,'("WARN L1-L10",/,6X,"T/F for warning messages to be given.")')
        write(*,'(6X,"1) Check Bkq and that they are consistent with the point group.",/,    & 
             6X,"2) check Alltp and that they are consistent with the point group.",/,       &
             6X,"3)-10) are currently not used.",/,                                          &
             6X,"For 1) & 2) a point group must be specified (SYML)")')
      CASE ("WORD") ! 1-10: check Bkq, check Alpt, 3-10 NU
        write(*,'("WORD L1 ",/,6X,"T/F for warning messages to be given.")')
        write(*,'(6X,"If L1 true, prints extra (wordy) output.")')
      CASE ("XREF")
        write(*,'("XREF  N1,N2")')
        write(*,'(6X,"Defines a molcule axis system determined by the 2 ligand atoms given by N1,N2.",/, &
              6X,"The metal is assumed at the origin. The Z axis will be in the 0->N1 direction.",/,     &
              6X,"The Y axis will be perpendicular to the plane defined by atom 0,N1,N2.",/,             &
              6X,"The X axis will be perpendicular to Y & Z axes, and form a right-handed system.",/,    &
              6X,"Note: An error will be generated if the atoms given by 0,N1,N2 are co-linear.")')
      CASE default
        Write(*,'("***WARNING: Unrecognized KEY=",A4)') KEY
      END SELECT 
!
      END subroutine help
!
!-----------------------------------------------------------------------
!
      SUBROUTINE hide(key)
! help for hidden commands. (untested, incomplete commands)
!
      IMPLICIT NONE
      integer i
      CHARACTER*4 KEY 
!      write(*,'("Subroutine hide, key:",A4)') key
      write(*,'(/)')
      select case (KEY)
      case ("HIDE")
        write(*,'(/,"Use:",/,">f_electron ",/,"The hidden keywords (untested, incomplete) are as follows:",/,  &
                  "Commands:",/,(10(A4,4X)))') (command_keyw3(i),i=1,10)
        write(*,'(/,"Type: >''f_electron hide kword'' for more information about a particular keyword ''kword''.")')
!
!  COMMANDS
!     (/"AOMX","CCF ","EXPL","EXPM","EXPT", "EXSO","FAST","LANC","PRT4","RORB"/)

! CRNG deleted

      CASE ("AOMX") 
        write(*,'("AOMX NL",/,                                                                             &
                  6X,"Uses the extended AOM parameters. ",/,                                               &
                  6X," NL is the number of ligands, and must equal the NL given in the AOM command.",/,    &
                  6X," An AOM command must also be given.",/,                                              &
                  6X," NL lines must immediately follow this command, each containing 4 parameters:",/,    &
                  6X," delS, delC, phiS, phiC (see : Urland,Chem.Phys,14,393,(1976) for more details)")')
      CASE ("CCF ")
        write(*,'("CCF    N1  N2",/,                                                                       &
             6X,"Correlation Crystal field parameters.",                                                   &
             6X,"N1 gives the CCF model (1: delta-function, 2: SCCF, 3: general)"                          &
             6X,"N2 is the number of CCF parameters (limit of 60).")')
        write(*,'(6X,"N2 lines must immediately follow the CCF command, with the following format:",/,     &  
             6X,"k,q   del(k,q) (N1=1,2)",/,6X,"i,k,q, Bikq (N1=3)",/,                                     &
             6X,"For more info see: Reid, J.Chem.Phys.,87,2875,(1987).")')
      CASE ("CRNG") 
        write(*,'("CRNG This command has been removed")')
      CASE ("EXPL") 
        write(*,'("EXPL N1-N2   (or N1 N2 N3 ..)",/,                                                       &
                  6X,"Either the range of |SLJ> states N1-N2 are used in the basis",/,                     &
                  6X,"or the explicit list N1 N2 N3 N4 ... |SLJ> states are used.",/,                      &
                  6X,"All |L S J MJ> states for the given |SLJ> multiplets are included.",/,               &
                  6X,"The hyphen is important: for a f^6 calculation ",/,                                  &
		  6X,"EXPL 1-7 uses 49 basis funtions corresponding to 7FJ, J=0,1,2,3,4,5,6",/,            &
		  6X,"EXPL 1,7 uses 14 basis funtions corresponding to 7F0 and 7F6",/,                     &
                  6X,"The |SLJ> basis functions are in the standard order of Neilson & Koster.",/,         &
                  6X,"Use CHK1(1) to see the |SLJ> basis functions printed to ""debug.dat""." )')
      CASE ("EXPM") 
        write(*,'("EXPM  ",/,                                                                              &
                  6X,"Flags that matrix elements of the lowest multiplet are fitted against the",/,        &
                  6X,"     matrix elements of the ESO (Extended Stevens Operators).",/,                    &
                  6X,"An EXSO command must have been given to provide the ESO matrix.",/,                  &
                  6X,"The J value of the lowest multiplet must equal that of the ESO matrix.",/,           &
                  6X,"For the calculation to be chemically reasonable, PRED command should also have ",/,  &
                  6X,"    been given so the lowest multiplet is diagonal in the atomic parameters.")') 
      CASE ("EXPT") 
        write(*,'("EXPT  N1, N2, N3 ..",/,                                                                 &
                  6X,"Expectation values of the parameters given by the numbers N1,N2,.. in their standard order.")') 
      CASE ("EXSO") 
        write(*,'("EXSO  N1 N2 N3 N4",/,                                                                   &
                  6X,"The (2J+1) energy levels of a multiplet will calulated defined in terms of ",/,      &
                  6X,"ESO (Extended Stevens Operators).",/,                                                &
                  6X,"N1: The value of 2J",/,                                                              &
                  6X,"N2: The value of kmax (must be even)",/,                                             &
                  6X,"N3: (0/1) No normalisation / Kmn factor used.",/,                                    &
                  6X,"N4: (0/1) The main calculation is done/skipped.",/,                                  &
                  6X," There will be N2/2 lines following this command.",/,                                &
                  6X,"If you cant fit all values on one line, multiple lines can be used,",/,              &
                  6X,"but a new line must be used for each new value of k.")')                  
      CASE ("FAST") 
        write(*,'("FAST L1-L4  ",/,                                                                        &
                  6X,"Options that may speed up a calculation in some cases",/,                            &
                  6X,"L1: Leave out the higher MJ values? (TODO)",/,                                       &
                  6X,"L2: Only calculate half of a Kramers doublet.",/,                                    &
                  6X,"    (must have symmetry such that Bkq for odd q is zero)",/,                         &
                  6X,"L3: If L2, Regenerate the other Kramers component.")')                  
      CASE ("LANC") 
        write(*,'("LANC N1,RL,RU",/,                                                                       &
                  6X,"Lanczos diagonalisation method used.",/,                                             &
                  6X,"N1: The number of Lanczos iterations to use (try ~number of levels you expect). ",/, &
                  6X,"RL,RU: The upper and lower bound of energies to be searched for energy levels")')
      CASE ("PRT4") ! 1-10: ESO info, ESO Mat,8*NU    in the |MJ> basis
        write(*,'("PRT4 L1-L10",/,                                                                      &
                  6X,"T/F for printing the matrices in the |MJ> basis to the file ""matrices.txt""",/,  &
                  6X,"1) ESO info in the (2J+1) basis of the ground state multiplet",/,                 &
                  6X,"2) Full ESO matrix and calculated multiplet matrices",/,                          &
                  6X,"3) Print matrix in decoupled |SLMlMs> basis.",/,                                  &
                  6X,"4) Print matrix in decoupled |SLMlMs=S> basis.",/,                                &
                  6X,"5) Print full matrix in |Sfi(i=1,7)> where fi is the occupancy (0,1,2) of the",/, &
                  9X," real f-orbitals. The definition and order of these are given by: ",/,            &
                  9X," Harung,Struct&Bond,12,201,(1972).",/,                                            &
                  9X," The number of basis functions will be the same as in the |SLMlMs=S> basis.",/,   &
                  6X,"6) Print Lx,Ly,LZ matrices in decoupled |SLMlMs=S> basis.",/,                     &
                  6X,"7)-9) Unused.",/,                                                                 &
                  6X,"10) Selection of matrix printed.",/,                                              &
                  9X,"In this case a line i1,i2, j1,j2 must follow; rows i1-i2 and columns j1-j2 are printed")')
      CASE ("RORB")
        write(*,'("RORB  OrbitalsMl,OrbitalsReal ",/,                                           &
                  6X,"OrbitalsMl = True (TODO) calculate the population of the ml orbitals",/,  &
                  6X,"OrbitalsReal = True (TODO) calculate the population of the real d/f orbitals")')
      CASE default
        Write(*,'("***WARNING: Unrecognized KEY=",A4," in HIDE")') KEY
      END SELECT 
!
      END subroutine hide
!
!-----------------------------------------------------------------------
!      
      integer function endProg()
!  Exit the program gracefully...
!  Gives the possibility of closing files to save data in a premature termination.     
      character*1 ans
      logical op
      write(*,'("Exit?.....E<return> to stop, C<return> to continue")')

      endProg=1
      read(*,*) ans
      write(*,*) ans
      if (ans.eq.'c' .or. ans.eq.'C') then
        write(*,'("......continuing")')
        return
      endif
      inquire (unit=IN,     opened=op); if (op) CLOSE(IN)
      inquire (unit=IO,     opened=op); if (op) CLOSE(IO)
      inquire (unit=IG,     opened=op); if (op) CLOSE(IG)
      inquire (unit=Idata,  opened=op); if (op) CLOSE(Idata)
      inquire (unit=IP1,    opened=op); if (op) CLOSE(IP1)
      inquire (unit=IP2,    opened=op); if (op) CLOSE(IP2)
      inquire (unit=Idebug, opened=op); if (op) CLOSE(Idebug)
      inquire (unit=Imatrix,opened=op); if (op) CLOSE(Imatrix)
! If one array allocated in a fit, they all are:
      if (ALLOCATED(ARfit))  DEallocate(ARfit,FitMatR,iIndexR,jIndexR)
      if (ALLOCATED(AIfit))  DEallocate(AIfit,FitMatI,iIndexI,jIndexI)
      if (ALLOCATED(Ukq)) DEallocate(Ukq,Ikq,Jkq,EDTrans)
! If one array allocated in a PRED, they all are:
      if (ALLOCATED(PRvec)) DEallocate(PRvec,LFMAT,CLFMAT)
!      
      stop
      end function endProg
!
!------------------------------------------------------------------------------
!  
      SUBROUTINE CONVERT(LI8,B64,C2,C1,convert_type)
!  Does various conversions.
!  convert_type=1  INTEGER*8 LI8 into INTEGER*4 B64(64) array of bits
!               2  reverse of 1
!               3  INTEGER*4 B64(64) array of bits into CHARACTER*2 C2(8) of hex numbers
!               4  reverse of 3
!               5  Char*1 C1(8) array of characters into a CHARACTER*2 C2(8) array of hex numbers each writtten as pair of characters
!               6  reverse of 5
!  no input is changed
      IMPLICIT NONE
      integer*8 LI8,I8,j8 
      integer*4 I,J,B64,convert_type,io,II,ios
      CHARACTER C1(8)*1, C2(8)*2
      DIMENSION B64(64)  
!
      SELECT CASE(convert_type)
        CASE (1) ! LI8 -> B64(64)
          i8=LI8 
          do i=1,64; b64(i)=0; enddo
          do i=63,1,-1
            j8=int8(2)**(i-1)
            if (i8.ge.j8) then; b64(I)=1; i8=i8-j8; endif
          enddo
        case(2)  !  B64 -> LI8  Note for this bit(64) (sign bit) is ignored.
          LI8=0; do i=1,63; LI8=LI8+int8(2)**(i-1)*b64(i); enddo 
        case(3)  !  B64(64) -> C2(8)
          do i=1,8
            II=0; do j=1,8; II=II+2**(j-1)*b64((i-1)*8+j); enddo 
            write(C2(i),'(Z2.2)') II
          enddo
        case(4) ! C2(8) -> B64(64)  
          do i=1,64; b64(i)=0; enddo
          do j=1,8
            READ(C2(j),'(Z2.2)') II
            do i=8,1,-1
              j8=int8(2)**(i-1)
              if (ii.ge.j8) then; b64((j-1)*8+I)=1; ii=ii-j8; endif
            enddo
          enddo
        case(5) !  C1(8) -> C2(8)  
          do i=1,8
            WRITE(C2(i),'(Z2.2)', IOSTAT=IOS) IACHAR(C1(i))  !  reads a character from array C1 to an hexadecimal ascii value stored as a character in array C2.
            if (IOS .ne. 0) then; write(*,'("Error in subroutine CONVERT converting C1 -> C2")'); STOP; endif
          enddo
        case(6) !  C2(8) -> C1(8)  
          do i=1,8
            READ(C2(i),'(Z2.2)', IOSTAT=IOS) II  !  reads a hexadecimal value from an element in C2 into II
            if (IOS .ne. 0) then; write(*,'("Error in subroutine CONVERT converting C2 -> C1")'); STOP; endif
            C1(i)=achar(ii)
          enddo
        case default
          WRITE(io,'("Unrecognised convert_type")') convert_type 
      end select
      RETURN
      END SUBROUTINE CONVERT
!
!------------------------------------------------------------------------------
!  
      SUBROUTINE DES(INPUT,KEY,ISW,JOTPUT)
!  DES encryption standard. Encrypts 64 bits, stored 1 bit per word, in array INPUT
!  into JOTPUT using KEY. Set ISW=0 for encryption, =1 for decryption
      IMPLICIT NONE
      integer*4 I,J,II,IC, ISW
      integer*4 INPUT,KEY,JOTPUT,ITMP,IP,IPM,ICF,KN,KNS
      DIMENSION INPUT(64),KEY(64),JOTPUT(64),ITMP(64),IP(64),IPM(64),ICF(32),KN(48),KNS(48,16)
      DATA IP /58,50,42,34,26,18,10, 2,60,52, 44,36,28,20,12, 4,62,54,46,38,      &
               30,22,14, 6,64,56,48,40,32,24, 16, 8,57,49,41,33,25,17, 9, 1,      &
               59,51,43,35,27,19,11, 3,61,53, 45,37,29,21,13, 5,63,55,47,39, 31,23, 15, 7/
      DATA IPM/40, 8,48,16,56,24,64,32,39, 7, 47,15,55,23,63,31,38, 6,46,14,      &
               54,22,62,30,37, 5,45,13,53,21, 61,29,36, 4,44,12,52,20,60,28,      &
               35, 3,43,11,51,19,59,27,34, 2, 42,10,50,18,58,26,33, 1,41, 9, 49,17,57,25/
! Recalculated KNS each time subroutine called
        DO I=1,16
          CALL KS(KEY,I,KN)  !  Get the 16 sub-master keys from the master key
          do j=1,48
            KNS(J,I)=KN(J)
          enddo 
        enddo 
      
      DO J=1,64                    ! Initial permutation.
        ITMP(J)=INPUT(IP(J))
      enddo 
!
      DO I=1,16                    ! 16 stages of encryption,
        II=I
        IF (ISW.EQ.1) II=17-I      ! Use the sub-master keys in reverse order for decryption.
        CALL CYFUN(ITMP(33),KNS(1,II),ICF) ! get cipher function
        DO J=1,32                
          IC=ICF(J)+ITMP(J)        ! Pass one half-word through unchanged, while encrpting the other 
          ITMP(J)=ITMP(J+32)       !  half-word and exchanging the two half-word output. 
          ITMP(J+32)=IAND(IC,1)
!         ITMP(J+32)=MOD(MOD(IC,2)+2,2)  ! Use if you don't have IAND
        enddo
      enddo  ! 16 stages
      DO J=1,32     ! Final exchange of 2 half-words.
        IC=ITMP(J)
        ITMP(J)=ITMP(J+32)
        ITMP(J+32)=IC
      enddo 
      DO J=1,64     !  Final permutation
        JOTPUT(J)=ITMP(IPM(J))
      enddo 
      RETURN
      END SUBROUTINE DES
!
!------------------------------------------------------------------------------
!  
      SUBROUTINE KS(KEY,N,KN)
!  Key schedule calculation, returns KN given Key and N=1,2..16; Must be called with N in that order.
      IMPLICIT NONE
      integer*4 N,J,I,IT,IC,ID
      integer*4 KEY(64),KN(48),ICD(56),IPC1(56),IPC2(48)
      DATA IPC1/57,49,41,33,25,17, 9, 1,58,50,  42,34,26,18,10, 2,59,51,43,35,  &
                27,19,11, 3,60,52,44,36,63,55,  47,39,31,23,15, 7,62,54,46,38,  &
                30,22,14, 6,61,53,45,37,29,21,  13, 5,28,20,12, 4/
      DATA IPC2/14,17,11,24, 1, 5, 3,28,15, 6,  21,10,23,19,12, 4,26, 8,16, 7,  &
                27,20,13, 2,41,52,31,37,47,55,  30,40,51,45,33,48,44,49,39,56,  &
                34,53,46,42,50,36,29,32/
      IF (N.EQ.1) THEN   !  Initial selection and permutation.
        DO J=1,56
          ICD(J)=KEY(IPC1(J))
        enddo 
      ENDIF
      IT=2   ! For most N perform 2 shifts
      IF (N.EQ.1.OR.N.EQ.2.OR.N.EQ.9.OR.N.EQ.16) IT=1  ! Excepts for these
      DO I=1,IT    ! Circular left-shifts of the 2 halves of the array ICD.
        IC=ICD(1)
        ID=ICD(29)
        DO J=1,27
          ICD(J)=ICD(J+1)
          ICD(J+28)=ICD(J+29)
        enddo 
        ICD(28)=IC
        ICD(56)=ID
      enddo   ! Finish shifts
      DO J=1,48
        KN(J)=ICD(IPC2(J))  ! Sub-master key is a selection of bits from shifted ICD
      enddo 
      RETURN
      END SUBROUTINE KS
!
!------------------------------------------------------------------------------
!  
      SUBROUTINE CYFUN(IR,K,IOUT)
!  Returns the cipher function of IR and K in IOUT
      IMPLICIT NONE
      integer*4 J,JJ,KK,KI,ISS,IROW,ICOL
      INTEGER*4 IR(32),K(48),IOUT(32),IE(48),IET(48),IP(32),ITMP(32),IS(16,4,8),IBIN(4,16)
      DATA IET/32, 1, 2, 3, 4, 5, 4, 5, 6, 7,   8, 9, 8, 9,10,11,12,13,12,13,  &
               14,15,16,17,16,17,18,19,20,21,  20,21,22,23,24,25,24,25,26,27,  &
               28,29,28,29,30,31,32,1/
      DATA IP /16, 7,20,21,29,12,28,17, 1,15,  23,26, 5,18,31,10, 2, 8,24,14,  &
               32,27, 3, 9,19,13,30, 6,22,11,   4,25/
!  The S-box
      DATA IS /14, 4,13, 1, 2,15,11, 8, 3,10, 6,12, 5, 9, 0, 7,  &
                0,15, 7, 4,14, 2,13, 1,10, 6,12,11, 9, 5, 3, 8,  &
                4, 1,14, 8,13, 6, 2,11,15,12, 9, 7, 3,10, 5, 0,  &
               15,12, 8, 2, 4, 9, 1, 7, 5,11, 3,14,10, 0, 6,13,  &
               15, 1, 8,14, 6,11, 3, 4, 9, 7, 2,13,12, 0, 5,10,  &
                3,13, 4, 7,15, 2, 8,14,12, 0, 1,10, 6, 9,11, 5,  &
                0,14, 7,11,10, 4,13, 1, 5, 8,12, 6, 9, 3, 2,15,  &
               13, 8,10, 1, 3,15, 4, 2,11, 6, 7,12, 0, 5,14, 9,  &
               10, 0, 9,14, 6, 3,15, 5, 1,13,12, 7,11, 4, 2, 8,  &
               13, 7, 0, 9, 3, 4, 6,10, 2, 8, 5,14,12,11,15, 1,  &
               13, 6, 4, 9, 8,15, 3, 0,11, 1, 2,12, 5,10,14, 7,  &
                1,10,13, 0, 6, 9, 8, 7, 4,15,14, 3,11, 5, 2,12,  &
                7,13,14, 3, 0, 6, 9,10, 1, 2, 8, 5,11,12, 4,15,  &
               13, 8,11, 5, 6,15, 0, 3, 4, 7, 2,12, 1,10,14, 9,  &
               10, 6, 9, 0,12,11, 7,13,15, 1, 3,14, 5, 2, 8, 4,  &
                3,15, 0, 6,10, 1,13, 8, 9, 4, 5,11,12, 7, 2,14,  &
                2,12, 4, 1, 7,10,11, 6, 8, 5, 3,15,13, 0,14, 9,  &
               14,11, 2,12, 4, 7,13, 1, 5, 0,15,10, 3, 9, 8, 6,  &
                4, 2, 1,11,10,13, 7, 8,15, 9,12, 5, 6, 3, 0,14,  &
               11, 8,12, 7, 1,14, 2,13, 6,15, 0, 9,10, 4, 5, 3,  &
               12, 1,10,15, 9, 2, 6, 8, 0,13, 3, 4,14, 7, 5,11,  &
               10,15, 4, 2, 7,12, 9, 5, 6, 1,13,14, 0,11, 3, 8,  &
                9,14,15, 5, 2, 8,12, 3, 7, 0, 4,10, 1,13,11, 6,  &
                4, 3, 2,12, 9, 5,15,10,11,14, 1, 7, 6, 0, 8,13,  &
                4,11, 2,14,15, 0, 8,13, 3,12, 9, 7, 5,10, 6, 1,  &
               13, 0,11, 7, 4, 9, 1,10,14, 3, 5,12, 2,15, 8, 6,  &
                1, 4,11,13,12, 3, 7,14,10,15, 6, 8, 0, 5, 9, 2,  &
                6,11,13, 8, 1, 4,10, 7, 9, 5, 0,15,14, 2, 3,12,  &
               13, 2, 8, 4, 6,15,11, 1,10, 9, 3,14, 5, 0,12, 7,  &
                1,15,13, 8,10, 3, 7, 4,12, 5, 6,11, 0,14, 9, 2,  &
                7,11, 4, 1, 9,12,14, 2, 0, 6,10,13,15, 3, 5, 8,  &
                2, 1,14, 7, 4,10, 8,13,15,12, 9, 0, 3, 5, 6,11/
      DATA IBIN/ 0,0,0,0,0,0,0,1,0,0,  1,0,0,0,1,1,0,1,0,0,  &
                 0,1,0,1,0,1,1,0,0,1,  1,1,1,0,0,0,1,0,0,1,  &
                 1,0,1,0,1,0,1,1,1,1,  0,0,1,1,0,1,1,1,1,0, 1,1,1,1/
      DO J=1,48   ! Expand IR to 48 bits and combine it with K
        IE(J)=IAND(IR(IET(J))+K(J),1)
!       IE(J)=MOD(MOD(IR(IET(J))+K(J),2)+2,2)   !  If you don't have IAND function
      ENDDO 
      DO JJ=1,8   !  Loop over 8 groups of 6 bits.
        J=6*JJ-5
        IROW=IOR(IE(J+5),ISHFT(IE(J),1))  ! Find place in S-box table.
        ICOL=IOR(IE(J+4),ISHFT(IOR(IE(J+3),ISHFT(IOR(IE(J+2),ISHFT(IE(J+1),1)),1)),1))
!        IROW=2*IE(J)+IE(J+5)   !  If you don't have the above bit functions
!       ICOL=8*IE(J+1)+4*IE(J+2)+2*IE(J+3)+IE(J+4)
        ISS=IS(ICOL+1,IROW+1,JJ)  ! Look up the number in the S-box table
        KK=4*(JJ-1)
        DO KI=1,4
          ITMP(KK+KI)=IBIN(KI,ISS+1)  ! ..and plug it's bits into the output.
        enddo 
      enddo 
      DO J=1,32  ! Final permutation
        IOUT(J)=ITMP(IP(J))
      enddo 
      RETURN
      END SUBROUTINE CYFUN
!      
!-----------------------------------------------------------------------
!
      END PROGRAM f_electrons
!   1355