MODULE f_e_input 

USE f_e_block  ! calls BlockBasis
USE f_e_check  ! calls checkAltp,checkAssignments,checkDiagOpt,checkESO,checkLinkFitP,checkPrint,checkValidERanges
USE f_e_data   ! GLOBALS
USE f_e_ESO    ! calls calc_ESO
USE f_e_interpreter   ! calls init(), destroyfunc
USE f_e_LF            ! calls convertToAOM,AOMmatrixD,AOMmatrixF,RotateLF,writeIntPara,A_B_Mat,zeroE
USE f_e_parameters    ! calls parLinearity, unloadP,checkLFcomplex
USE f_e_group         ! calls GETGRP()
USE f_e_fncrosserrors ! calls FncrossErrors
USE f_e_magnetics     ! calls calcMagMom
USE f_e_readFile  ! calls getLSBasis,getFullBasis,readUMat,readVMat,readMnPnSnMat,readTnMat,readEEMat,readABGMat,
                  !       t2MatComp, testTnMat
!
!  This module reads the input file.
!
! TODO
! Writes output to plotfile.
!    TODO: spin & MJ for IRREPS
! TODO: spin project not working for degeneracy of 4 states.

IMPLICIT NONE

! PRIVATE

PUBLIC 

CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE input(INFILE)
      IMPLICIT NONE
      integer nline,ninput,i,j,j1,k,q,ii,iL,NTitle,N,n1,n2
      INTEGER nr, ni, i1,i2,i3
      logical echoOn,UseFncrossErrors(8),test,test2,magFld,CONFcom
      logical Fn, UseDefaultIRREPS
      CHARACTER KEY*4, COMMAND*200, TITLE*200/""/, INFILE*80, txt*20
      
      logical newVar
      integer nV,minNv,gridNv
      character START*20,STOP*20,STEP*20, nfor*20, vname*10, noNames(1)*10, func*100
      CHARACTER lab*3, lab1*8,lab2*8,lab3*8  !  Temporary variables
      real*8 vStart,vStop,vStep,noValues(1)
      character (len = 5)  :: statusflag
      
      integer nExplicitbasis, nExBasisNums(max_Jstates)  ! total number of f7 JMJ states
      
      real*8 A(100),D0,D1,D2,bit
      real*8 DistMat(max_ligands,max_ligands),AngMat(max_ligands,max_ligands),DistBonds(max_ligands),r1,r2
      character*1 DistLab(max_ligands,max_ligands)
      character*30 fnamed,fnamem,fnamep,Fext
      parameter(D0=0.0d+00,D1=1.0d+00,D2=2.0d+00, bit=1.0d-12) 
! local variables:
      ninput=0            ! line number read in the input file.
      echoOn=.false.      ! Don't echo input commands.
      CONFcom=.false.     ! make sure that a CONF command has been given.
      NTitle=0            ! Number of title lines read.
      do k=1,8; UseFncrossErrors(k)=.false.; enddo   ! don't use (incorrect) fncross values.
      Fn=.false.; UseDefaultIRREPS=.false.; Fext=""
      nExplicitbasis=0    ! No explicit basis
      
      noNames(1)="noName"; noValues(1)=d0
      do i=1,maxCon        ! Constants
        ConstantV(i) = d0  !   value
        Constants(i) = ""  !   name
        cEq(i) = "0"       !   formula (eg sqrt(3)) 
      enddo
!      nCon=1; Constants(1) = "Pi"; cEq(1) = "atan(1.0,1.0)*4"
! Every parameter can be given as:
!  i) a number,  
!  ii) a formula that can be evaluated immediately, (eg sqrt(3)) 
!  iii) a formula that depends on a variable that may be varied. 
! Initially all parameters zero.
      do i=1,maxTotPar
        P(i) = d0         
        pEq(i) = "0"      
        pEqF(i) = .false.  ! True if this parameter is a formula containing at least one variable
      enddo
      p(310)=D1;p(311)=D1;p(312)=D1; pEq(310)="1";pEq(311)="1";pEq(312)="1"  ! Reduction factors default 1
      do i=1,maxTotVar
        variables(i) = ""
!        variablesvalues(i) = d0
        do j=1,maxTotVarVals
          allVarVals(i,j)=d0  ! Each variable can have maxTotVarVals values.
        enddo
      enddo
      
 10   ninput=ninput+1
      read(IN,"(A4,A120)",END=902) KEY,COMMAND
      CALL TO_UPPER(KEY)
      if (echoOn) WRITE(IO,'(">",A4,"|",A)') KEY,trim(COMMAND)
      
      if (key.eq."    " .or. index(key,"!").ne.0) goto 10  ! skip comments or blank lines      
      if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
        i=index(COMMAND,"!"); COMMAND(i:)=""
      endif
      select case (KEY)
!
!----------------------------------------------------------------------------------------------------------        
!   Commands:  
      CASE("BASE")  ! whether the full MJ or J basis to be used.
        COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); lab=COMMAND(:i); COMMAND=COMMAND(i+1:)
!        write(io,'("Rest of command:",A)') COMMAND
        FullMJBasis=.true.
        if (trim(lab).eq."J") FullMJBasis=.false.
      CASE("BLOC")  ! 
        block=.true.; 
        call DECODE(A,N,COMMAND,KEY)
!        read(COMMAND,*,err=930,END=930) i
        if (n.eq.1 .and. A(1).ne.0) PrintCQN=.true.
      CASE("CONF")
        READ(COMMAND,*,err=930,END=930) Lvalue,Nelectrons 
        CONFcom=.true.
        if (Lvalue.eq.2) then
          e_lab="d"; n_states=16
          OPEN(Idata,FILE=trim(F_E_PATH)//'\d_electron.dat',STATUS='OLD',ERR=917)
        elseif (Lvalue.eq.3) then
          e_lab="f"; n_states=119;
          OPEN(Idata,FILE=trim(F_E_PATH)//'\f_electron.dat',STATUS='OLD',ERR=918)
        else
          write(IO,'("CONF Invalid input, Lvalue=",I4)') Lvalue 
          goto 900
        endif
        If (Nelectrons.lt.1 .or. Nelectrons.gt.(4*Lvalue+1)) then
          write(IO,'("CONF Invalid input, Nelectrons=",I4," Must be <=",I2)') Nelectrons,2*Lvalue+1 
          goto 900
        endif
        N_odd=.true.; Complementary=.false. 
        if (MOD(Nelectrons,2).eq.0) N_odd=.false. 
        IF (Nelectrons.gt.(2*Lvalue+1)) then  ! Nelectrons is always 0<Nelectrons<=7
          Nelectrons = 2*(2*Lvalue+1)-Nelectrons
          Complementary=.true.
          if (idebug.gt.0) write(idebug,'("Number of electrons redefined to complementary number:",I2)') Nelectrons 
        ENDIF 
        CALL getLSBasis()   ! loads d or f LS basis
      CASE("EDIP")
        read(COMMAND,*,err=930,END=930) ED(1,1),ED(2,1),ED(1,2),ED(2,2),IED,EDunits,EDconst ! electric dipole calcs
        ED_IntMech="Full"
      CASE("EXP1")   ! Flag that the one electron LF matrix elements are to be used in a fit.
        if (FitType.ne.0) then; write(IO,'("***WARNING: more than one of EXPB,EXPE,EXPG,EXPM,EXP1.")'); endif 
        FitType=6; j=0
        Do k=1,2*Lvalue+1   !  
          ninput=ninput+1
          READ(IN,'(A120)',ERR=906) COMMAND
          if (echoOn) write(io,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").ne.0) then;  ii=index(COMMAND,"!"); COMMAND(ii:)="";  endif ! skip comments at end of line      
          call DECODE(A,N,COMMAND,KEY)
 !         write(io,'("N,COMMAND=",i1,2x,A)') N,COMMAND
 !         write(io,'("A(N)=",7(2X,2F8.1))') (A(I),I=1,N)
          if (N.ne.2*k) then 
            write(IO,'("***FATAL: line after ENG1 Expecting",I3,": (values + weight) for each ME",I4)') 2*k
            STOP
          endif
          do i=1,k
            j=j+1
            EXP_E(j)=A(2*i-1); WGT_E(j)=A(2*i); WGTchar(j)=" "
            if (WGT_E(j).lt.1.0d+00) then; WGTchar(j)="x"; wgtNote=.true.; endif
            if (WGT_E(j).gt.1.0d+00) then; WGTchar(j)="+"; wgtNote=.true.; endif
          enddo
        enddo
        nexp=j ! should be 15 (d) or 28 (f)
      CASE("EXPB")   ! Bkq data to be used in a fit.
        if (FitType.ne.0) then; write(IO,'("***WARNING: more than one of EXPB,EXPE,EXPG,EXPM,EXP1")'); endif 
        FitType=3
        read(COMMAND,*) NEXP
        if (Nexp.lt.0 .or. e_lab.eq."d".and.Nexp.gt.14 .or. e_lab.eq."f".and.Nexp.gt.27) then
          write(io,'("***FATAL: Invalid value for Nexp=",I4)') NEXP; stop
        endif
        do i=1,nexp
 20       ninput=ninput+1          
          read(IN,'(A120)',ERR=925) COMMAND
          COMMAND=adjustl(COMMAND)
          if (echoOn) WRITE(IO,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").eq.1) goto 20 ! skip comment lines
          if (index(COMMAND,"!").ne.0) then; ii=index(COMMAND,"!"); COMMAND(ii:)=""; endif   ! skip comments at end of line
          read(COMMAND,*,ERR=925) EXP_E(I),WGT_E(I),txt
          WGTchar(I)=" "; 
          if (WGT_E(I).lt.1.0d+00) then; WGTchar(I)="x"; wgtNote=.true.; endif
          if (WGT_E(I).gt.1.0d+00) then; WGTchar(I)="+"; wgtNote=.true.; endif
          txt=ADJUSTL(txt)
          NAssign(i,1)=0
          do j=1,16; if (BkqR_label(j).eq.txt(1:4)) NAssign(i,1)= j; enddo
          do j=1,16; if (BkqI_label(j).eq.txt(1:4)) NAssign(i,1)=-j; enddo
          if (echoOn) write(io,'(I3," EXP_E(I),WGT_E(I),txt,NAssign(i)=",2F10.2,2X,A20,I4)') i,EXP_E(I),WGT_E(I),txt,NAssign(i,1)
          if (NAssign(i,1).eq.0) then
            write(io,'("****FATAL: Crystal field coefficient:",i2,2X,A4," not found")') i,txt(1:4)
            stop
          endif  
        enddo
      CASE("EXPE")  ! Energies/Intensities to be used in a fit.
        if (FitType.ne.0) then; write(IO,'("***WARNING: more than one of EXPB,EXPE,EXPG,EXPM,EXP1.")'); endif 
        FitType=1
        if (JuddOfelt) FitType=2
        call DECODE(A,N,COMMAND,KEY)
        if (n.lt.2 .or. n.gt.3) then  ! First two numbers compulsory, 3rd optional        
          write(io,'("***FATAL: Invalid number of values in EXPE command:")');
          write(io,'("EXPE",A)') COMMAND; stop
        endif
        NEXP=INT(A(1)); EngInt=INT(A(2)); if (n.eq.3) assignMethod=INT(A(3))
!        read(COMMAND,*,ERR=921,end=921) NEXP,EngInt,assignMethod
!        Write(io,'("NEXP,EngInt,assignMethod=",3I4)') NEXP,EngInt,assCQN
        if (Nexp.lt.0 .or. Nexp.gt.max_exp) then; write(io,'("***FATAL: Invalid value for Nexp=",I4)') NEXP; 
                                                  write(io,'(10X,"Must be between 0 and",I4)') max_exp; stop; endif
        if (EngInt.lt.1 .or. EngInt.gt.2) then; write(io,'("***FATAL: Invalid value for EngInt=",I4," in EXPE command")') EngInt; 
                                                write(io,'(10X,"Must be either 1(Engs) or 2(Eng&Ints)")'); stop; endif
        if (assignMethod.lt.0 .or. assignMethod.gt.2) then
          write(io,'("***FATAL: Invalid value for assignMethod=",I4," in EXPE command")') assignMethod 
          write(io,'(10X,"Must be either 0, 1 or 2 (assign with energy order,CQNs,spin)")'); stop
        endif
        do i=1,nexp
 30       ninput=ninput+1          
          read(IN,'(A120)',ERR=920) COMMAND
          COMMAND=adjustl(COMMAND)
          if (echoOn) WRITE(IO,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").eq.1) goto 30 ! skip comment lines
          if (index(COMMAND,"!").ne.0) then; ii=index(COMMAND,"!"); COMMAND(ii:)=""; endif   ! skip comments at end of line      
          call DECODE(A,N,COMMAND,KEY)
          if (N.eq.0) then 
            write(IO,'("***FATAL: line after EXPE Expecting nexp=",I4," levels, found only:",I4)') nexp,i
            STOP
          endif
!          write(io,'("COMMAND=",A120)') COMMAND
!          write(io,'("N,A(N)=",I3,20i3)') N,(INT(A(j)),j=1,N)
          EXP_E(I)=A(1); WGT_E(I)=A(2); WGTchar(I)=" ";
          if (WGT_E(I).lt.1.0d+00) then; WGTchar(I)="x"; wgtNote=.true.; endif
          if (WGT_E(I).gt.1.0d+00) then; WGTchar(I)="+"; wgtNote=.true.; endif
          if (WGT_E(i).lt.0.0d+00) write(io,'("****Warning: WGT_E(",i3,")=",F10.2,"<0.0")') i,WGT_E(i)
          if (i.gt.1 .and. assignMethod.eq.0) then
            if (EXP_E(i-1).gt.EXP_E(i)) write(io,'("****Warning: Exp(",i3,")=",F10.2,">exp(",i3,")=",F10.2)') i-1,exp_E(i-1), &
                                                                                                                i,exp_E(i)
          endif
          if (EngInt.eq.2) then !  energies & intensities
            if (N.lt.5 .or. N.gt.25) then; write(IO,'("***FATAL: line after EXPE  N=",I4," must be >4 & <26.")') N; STOP; endif !
            Nass(i)=N-4
            EXP_I(I)=A(3);WGT_I(I)=A(4);
            if (EXP_I(i).lt.0.0d+00) write(io,'("****Warning: Int(",i3,")=",F10.2,"<0.0")') i,exp_I(i)
            if (WGT_I(i).lt.0.0d+00) write(io,'("****Warning: WGT_I(",i3,")=",F10.2,"<0.0")') i,WGT_I(i)
            do j=1,Nass(i); NAssign(i,j)=NInt(A(4+j)); enddo
          else ! just energies
            if (N.lt.3 .or. N.gt.23) then; write(IO,'("***FATAL: line after EXPE N=",I4," must be >2 & <24.")') N; STOP; endif !
            Nass(i)=N-2
            do j=1,Nass(i); NAssign(i,j)=NInt(A(2+j)); enddo         
          endif
          if (assignMethod.eq.1 .and. Nass(i).gt.2) write(io,'("****Warning: Only 2 assignment numbers required for fitting ", &
                            "the crystal quantum numbers",i3," given, extra ones ignored.")') i
          if (assignMethod.eq.2 .and. Nass(i).gt.2) write(io,'("****Warning: Only 2 assignment numbers required for fitting ", &
                            "the energy levels by spin",i3," given, extra ones ignored.")') i
          if (echoOn) then
            if (.not.JuddOfelt) write(io,'(I4,"; EXP_E(I),WGT_E(I),NAssign(i,j)=",2F10.2,20I4)') i,EXP_E(I),WGT_E(I),  &
                                                                                                 (NAssign(i,j),j=1,Nass(i))
            if (JuddOfelt) write(io,'(I4,"; EXP_E,WGT_E,EXP_I,WGT_I,NAssign=",4F10.2,20I4)') i,EXP_E(I),WGT_E(I),  &
                                                                              EXP_I(I),WGT_I(I), (NAssign(i,j),j=1,Nass(i))
          endif
        enddo
      CASE("EXPG")   !  g-values to be used in a fit.
        read(COMMAND,*,ERR=922,end=922) NexpG 
        if (FitType.ne.0) then; write(IO,'("***WARNING: more than one of EXPB,EXPE,EXPG,EXPM,EXP1.")'); endif 
        if (NexpG.lt.0 .or. NexpG.gt.max_exp) then; write(io,'("***FATAL: Invalid value for NexpG=",I4)') NexpG; stop; endif
        FitType=5
        do i=1,NexpG
 40       ninput=ninput+1          
          read(IN,'(A120)',ERR=920) COMMAND
          COMMAND=adjustl(COMMAND)
          if (echoOn) WRITE(IO,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").eq.1) goto 40 ! skip comment lines
          if (index(COMMAND,"!").ne.0) then; ii=index(COMMAND,"!"); COMMAND(ii:)=""; endif   ! skip comments at end of line      
          call DECODE(A,N,COMMAND,KEY)
          if (N.ne.7) then; write(IO,'("***FATAL: Expecting ",I4," lines of g1,w1, g2,w2, g3,w3, assign.")')NexpG; STOP; endif !
!          write(io,'("COMMAND=",A120)') COMMAND
!          write(io,'("N,A(N)=",I3,20i3)') N,(INT(A(j)),j=1,N)
          EXP_G(1,I)=A(1); WGT_G(1,I)=A(2); EXP_G(2,I)=A(3); WGT_G(2,I)=A(4); EXP_G(3,I)=A(5); WGT_G(3,I)=A(6); NAssignG(i)=A(7)
          do j=1,3; if (WGT_G(j,i).lt.0.0d+00) write(io,'("****Warning: WGT_G(",2i3,")=",F10.2,"<0.0")') j,i,WGT_G(j,i); enddo
          if (echoOn) then
            write(io,'(I4,"; G1,W1,G2,W2,G3,W3,NAssign(i)=",3(F10.4,F6.2),I4)') i,EXP_G(1,I),WGT_G(1,I),  &
                                                          EXP_G(2,I),WGT_G(2,I),EXP_G(3,I),WGT_G(3,I),NAssignG(i)
          endif
        enddo
      CASE("EXPM")   ! Flag that the ESO matrix elements are to be used in a fit.
        if (FitType.ne.0) then; write(IO,'("***WARNING: more than one of EXPB,EXPE,EXPG,EXPM,EXP1.")'); endif 
        FitType=4
      CASE("EXPT"); call DECODE(A,N,COMMAND,KEY)   ! Expectation values of parameters
        if (N.le.0) then; write(io,'("***FATAL: no values found for EXPT command")'); stop; endif
        if (N.gt.20) then; write(io,'("***FATAL: Too many values found for EXPT command (max=20)")'); stop; endif
        N_expect=n
        do i=1,N_expect
          P_expect(i)=int(A(i)); 
          if (P_expect(i).le.0 .or. P_expect(i).gt.140) then
            write(io,'("***FATAL: invalid value for",I3,"th parameter:",I4)') i,P_expect(i); stop
          endif
        enddo
      CASE("EXPL")  ! Explicit basis
         do while (trim(COMMAND).ne."")   ! delimiters can be a mixture of spaces and commas
          COMMAND=ADJUSTL(COMMAND);  i=scan(COMMAND,", "); lab1=COMMAND(:i-1); 
          COMMAND=ADJUSTL(COMMAND(i+1:)); if (COMMAND(1:1).eq.',') COMMAND=COMMAND(2:)
 !         write(io,'("i=",I2,", lab1:",A8,", COMMAND:",A)') i,lab1,COMMAND
          if (scan(lab1,"-").ne.0) then
            i=scan(lab1,"-")
 !           write(io,'("A - was found at position",I2," of chunk:",A8)') i,lab1
            lab2=lab1(:i-1); j=scan(lab1,", "); lab3=lab1(i+1:j-1)
 !           write(io,'("i,j=",2I2," Numbers to range from:",A8," to ",A8)') i,j,lab2,lab3
            read(lab2,*) i2; READ(lab3,*) i3
 !           write(io,'("Numbers in range from:",i3," to ",i3)') i2,i3
 ! no check to see if basis numbers are monatomic, or valid range for particular f/d configuration, or contain duplicates.
            if (i2.gt.0 .and. i2.le.max_Jstates .and. i3.gt.0 .and. i3.le.max_Jstates .and. i3-i2.gt.0) then
              do i=i2,i3
                nExplicitbasis=nExplicitbasis+1
                nExBasisNums(nExplicitbasis)=i
              enddo
            else
              write(io,'("***FATAL: Invalid values for range:",i3,"-",i3," in ",A)') i2,i3,COMMAND
              STOP
            endif            
          else
 !           write(io,'("Numbers to add:",A8)') lab1
            read(lab1,*) i1; ! write(io,'("Numbers to add:",i3)') i1
            if (i1.gt.0 .and. i1.le.max_Jstates) then  
              nExplicitbasis=nExplicitbasis+1
              nExBasisNums(nExplicitbasis)=i1
            else  
              write(io,'("***FATAL: Invalid values for basis function:",i3," in ",A)') i1,COMMAND
              stop
            endif
          endif
!          write(io,'("Rest of command:",A)') COMMAND
        enddo  
 !       write(io,'(i3," explicit basis functions to be used.")') nExplicitbasis
 !       write(io,'((20i3))') (nExBasisNums(i),i=1,nExplicitbasis)
      CASE("FAST") !  1-10: MJ,Kramers,Regen 
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for FAST")'); STOP; endif            
          fastMat(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("GEXS"); call DECODE(A,N,COMMAND,KEY)
        if (A(1).ne.0) g_axes=.true.
        if (a(2).gt.0) then
          if (N.ne.3) then
            write(io,'("***FATAL: GEXS N1 N2 N3 expected; Number of values=",I3)') N
            stop
          endif
          if (mod(int(A(3)-A(2))+1,2).ne.0) then
            write(io,'("***FATAL: Must have an even number of levels specified. (N3-N2+1)=",I4)') INT(A(3)-A(2)+1)
            stop
          endif
          ngexs=int(A(3)-A(2)+1)/2
          do i=1,ngexs 
            GEXS(2*i-1)=int(ABS(A(2)))+2*(i-1)
            GEXS(2*i  )=int(ABS(A(2)))+2*(i-1)+1
          enddo
        else
          ngexs=n-1
          if (mod(ngexs,2).ne.0) then
            write(io,'("***FATAL: Must have an even number of levels specified. (N)=",I4)') ngexs
            stop
          endif
          ngexs=ngexs/2
          do i=1,ngexs 
            GEXS(2*i-1)=ABS(int(A(2*i  )))
            GEXS(2*i  )=ABS(int(A(2*i+1))) 
          enddo
        endif
        if (ngexs.gt.200) then
          write(io,'("***FATAL: GEXS you can only have a maximum of 200 g-values. N=",I4)') ngexs
          stop
        endif
!        write(*,'("ngexs=",i3)') ngexs
!        write(*,'("gvals=",5(2I3,2X))') (GEXS(2*i-1),GEXS(2*i),i=1,ngexs)
      CASE("JUDO") ! JO multiplet calcs
        read(COMMAND,*,err=940,END=940) ED(1,1),ED(2,1),ED(1,2),ED(2,2),EDunits,Dielectric,JO_omega(1),JO_omega(2),JO_omega(3)
!        write(IO,'("ED(1,1),ED(2,1),ED(1,2),ED(2,2),EDunits=",5I3)') ED(1,1),ED(2,1),ED(1,2),ED(2,2),EDunits
!        write(IO,'("Dielec,Omega(1),Omega(2),Omega(3)  =",4F8.4)') Dielectric,JO_omega(1),JO_omega(2),JO_omega(3)
        if (EDunits.lt.1 .or. EDunits.gt.4) then; write(IO,'("***FATAL: EDunits=",I3," must be >=1 & <=4")') EDunits; STOP; endif 
        write(COMMAND,'(F10.6)') JO_omega(1); read(COMMAND,*) pEq(301)
        write(COMMAND,'(F10.6)') JO_omega(2); read(COMMAND,*) pEq(302)
        write(COMMAND,'(F10.6)') JO_omega(3); read(COMMAND,*) pEq(303)
        JuddOfelt=.true.
        FullMJBasis=.false.
        ED_IntMech="J-O"
      CASE("MDIP"); read(COMMAND,*,err=930,END=930) MD(1,1),MD(2,1),MD(1,2),MD(2,2),IMD,MDunits,MDconst ! magnetic dipole calcs
      CASE("MFLD"); 
        do i=1,3   
          COMMAND=ADJUSTL(COMMAND);  ii = index(COMMAND," "); if (ii.eq.0) goto 911; pEq(303+i)=COMMAND(:ii); 
          COMMAND=COMMAND(ii+1:)
        enddo  
! MAGF(1),MAGF(2),MAGF(3) are equations in pEq(304-306) and are Hx, Hy, Hz of magnetic field in Tesla.
!
      CASE("OPTN")  ! (T/F)  E0?, spin proj?, free ion % short?, free ion % long?,No Hss?, MJ % short?, MJ % long?, Skip write&diag?, Skip diag, NU
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for OPTN")'); STOP; endif            
          Option(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("OUTP")  ! (L1-L10) Print:  Lig dist(aom),Equiv CF(aom), AOM mat, Equiv Int.Par,CF before/after rotation,
                    !                  No degeneracy sum,matrix info,all parameters in multiple calcs, NU*,NU*
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for OUTP")'); STOP; endif            
          OUTP(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("PLT1")  ! (L1-L10) PLOT: Energies,+Spin fract,+SpinAllowedInt,+MJ fract,MDIP,  EDIP,Params,Bkq,g-val,d-orbEngs  
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for PLOT")'); STOP; endif            
          PLOTout1(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("PLT2")  ! (L1-L10) PLOT: MinVal,..  
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for PLOT")'); STOP; endif            
          PLOTout2(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("PRNG");  call DECODE(A,N,COMMAND,KEY)
        if (N.ne.3) then; write(io,'("***FATAL: There must be 3 parameters for PRNG command")'); stop; endif
        i=INT(A(1));  if (i.lt.1.or.i.gt.2) then; write(io,'("***FATAL: N1 must be 1/2 in PRNG command")'); stop; endif
        if (i.eq.2) then; NrangeL=int(A(2)); NrangeH=int(A(3)); else; PrangeL=A(2); PrangeH=A(3); endif
        if (PrangeL.gt.PrangeH .or. NrangeL.gt.NrangeH) then
          if (i.eq.1) write(io,'("***FATAL: Invalid value for range of energies to be printed=",2F12.2)') PrangeL,PrangeH
          if (i.eq.2) write(io,'("***FATAL: Invalid value for range of energies to be printed=",2i4)') NrangeL,NrangeH
          stop
        endif
        if (PrangeL.eq.PrangeH) then; PrangeL=D0; PrangeH=D0; endif
      CASE("REDF"); read(COMMAND,*,err=970,END=970) RK(1),RK(2),RK(3)  !  The orbital reduction factors.
        do i=1,3
          if (RK(i).lt.0.0 .or. RK(i).GT.1.0d+00) then
            write(io,'("***FATAL: Invalid value for Orbital Redction Factor RK(",i1,")=",F8.4)') RK(i)
          endif
        enddo
      CASE("SPIN"); read(COMMAND,*,err=970,END=970) SpAllow(1),SpAllow(2)  !  The initial states to calculate 
                                                                           !  spin-allowed transitions from.
        if (SpAllow(1).lt.1 .or. SpAllow(2).lt.1 .or. SpAllow(1).gt.SpAllow(2)) then
          write(io,'("***FATAL: Invalid values for initial values for a SPIN command:N1,N2=",2I6)') SpAllow(1),SpAllow(2)
          stop
        endif
      CASE("SYML")   ! can either be: "SYML GROUP" or "SYML GROUP ISO IGP"
        OPEN(IG,FILE=trim(F_E_PATH)//'\group.dat',STATUS='OLD',ERR=914)
!        COMMAND=adjustl(COMMAND); GROUP=COMMAND(1:3); COMMAND=COMMAND(4:120)
        COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); GROUP=COMMAND(:i); COMMAND=COMMAND(i+1:)
        call DECODE(A,N,COMMAND,KEY)
        if (N.eq.0) UseDefaultIRREPS=.true.
        if (N.eq.1) ISO=int(A(1))
        if (N.eq.2) then; ISO=int(A(1)); IGP=int(A(2)); endif
!        read(COMMAND,*,err=930,END=930) GROUP,ISO,IGP
        If (ISO.NE.0 .and. ISO.ne.1) then
          write(IO,'("SYML: Invalid input, ISO=",I4," Must be 0 or 1")') ISO 
          goto 900
        endif 
        If (IGP.NE.0 .and. IGP.NE.1) then
          write(IO,'("SYML: Invalid input, IGP=",I4," Must be 0 or 1")') IGP 
          goto 900
        endif
        CALL TO_UPPER(GROUP)
        IDG=-1
      CASE("TITL")
        NTitle=NTitle+1
        read(COMMAND,'(A200)') TITLE
        if (NTitle.eq.1) WRITE(IO,'(" Title:",A)') trim(TITLE)
        if (NTitle.gt.1) WRITE(IO,'(7X,A)') (trim(TITLE))
      CASE("VECT")
        call DECODE(A,N,COMMAND,KEY)
        if (n.lt.2 .or. n.gt.3) goto 930
        vectLow=int(A(1))
        vectHigh=int(A(2))
        if (n.eq.3) vectN=int(A(3))  ! last value optional
      CASE("WORD"); read(COMMAND,*,err=930,END=930) wordy  !  Extra output if true. 
!
!----------------------------------------------------------------------------------------------------------        
!  Fitting Commands:


!----------------------------------------------------------------------------------------------------------        
!  Extra Commands:      
      CASE("VARI")  ! variables can be defined on multiple lines(?) 
        COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); vname=COMMAND(:i); COMMAND=COMMAND(i+1:)
        call TO_LOWER(vname)
        newVar=.true.
        do i=1,nvar
          if (vname==variables(i)) then
            newVar=.false.
            nv=i
          endif
        enddo
        if (newVar) then
          nVar=nVar+1
          nV=nVar
          variables(nV)=vname
        endif
        call TO_UPPER(COMMAND)
        if (index(COMMAND, "TO") /= 0) then
          ii=index(COMMAND,"FROM"); n1=index(COMMAND,"TO"); n2=index(COMMAND,"STEP")
          if (n1.le.ii .or. n2.le.n1) then
            write(io,'("***FATAL: Invalid command VARI ",A10," ",A)') variables(nV),COMMAND
            write(io,'(10x,"Correct format: VARI varname FROM r1 TO R2 STEP r3")')
            stop
          endif  
          start=COMMAND(ii+4:n1-1); stop=COMMAND(n1+2:n2-1); step=COMMAND(n2+4:)
          read(start,*,err=930,END=930) vstart
          read(stop,*,err=930,END=930) vstop
          read(step,*,err=930,END=930) vstep
          if (abs(vStep).lt.1.0d-20) then
            write(io,'("***FATAL: Invalid value for STEP of variable:",A10,", STEP=",E16.8)') variables(nV),vStep
            stop
          endif  
          n=Int((vStop-vStart)/vStep) + 1
          if (n.lt.1 .or. nVarVals(nV)+n.gt.maxTotVarVals) then
            write(io,73) variables(nV),n, maxTotVarVals
 73   FORMAT("***FATAL: Invalid number of values for variables:",A10,", nVarVals=",I4,"; Minimum allowed=1, Maximum allowed:",I4)
            stop
          else
            do i=1,n; allVarVals(nV,nVarVals(nV)+i)=vStart+dble(i-1)*vStep; enddo
            nVarVals(nV)=nVarVals(nV)+N
          endif
          
        else if (index(COMMAND, "FOR") /= 0) then   ! allows zero step size, Stop=start
          ii=index(COMMAND,"FROM"); n1=index(COMMAND,"STEP"); n2=index(COMMAND,"FOR")
          if (n1.le.ii .or. n2.le.n1) then
            write(io,'("***FATAL: Invalid command VARI ",A10," ",A)') variables(nV),COMMAND
            write(io,'(10x,"Correct format: VARI varname FROM r1 STEP R2 FOR n1")')
            stop
          endif  
          start=COMMAND(ii+4:n1-1); step=COMMAND(n1+4:n2-1); nfor=COMMAND(n2+3:)
          read(start,*,err=930,END=930) vstart
          read(step,*,err=930,END=930) vstep
          read(nfor,*,err=930,END=930) n
!          write(io,'("COMMAND:",A120)') COMMAND
!
          if (n.lt.1 .or. nVarVals(nV)+n.gt.maxTotVarVals) then
            write(io,73) variables(nV),n, maxTotVarVals
            stop
          else
            do i=1,n; allVarVals(nV,nVarVals(nV)+i)=vStart+dble(i-1)*vStep; enddo
            nVarVals(nV)=nVarVals(nV)+N
          endif
  
        ELSE  ! .not. TO or FROM
!        
          IF (CheckVariableIsNumberOnly(COMMAND)) then
            write(io,'("***FATAL: Variable: ",A," contains non-numerical values:",A)') trim(vname),trim(COMMAND)
            stop
          endif
          call DECODE(A,N,COMMAND,KEY)
          if (n.lt.1 .or. nVarVals(nV)+n.gt.maxTotVarVals) then
            write(io,73) variables(nV),n, maxTotVarVals
            stop
          else
            do i=1,n; allVarVals(nV,nVarVals(nV)+i)=A(i); enddo
            nVarVals(nV)=nVarVals(nV)+N
          endif
        ENDIF
      CASE("CONS")  ! constants  COMMAND  *120 
        COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); vname=COMMAND(:i); COMMAND=COMMAND(i+1:)
        call TO_LOWER(vname)
        if (ncon.gt.0) then  
          do i=1,nCon
            if (vname.eq.Constants(i)) then
              write(io,'("***FATAL: Can only define a Constant once: CONS ",A10," ",A)') vname,COMMAND
              stop
           endif
          enddo
        endif
        if (nCon.gt.maxCon) then; write(io,'("***FATAL: Too many constants defined: nCon=",I4)') nCon; STOP; endif
        func = COMMAND
        if (len_trim(func) == 0) then
          write(io,'(/,"****FATAL: There is no function defined for the constant:",A10,  &
                     /,"Make sure there is a space as a delimiter after the constant name")') vname
          STOP
        end if
        if (index(func,'=') /= 0) then
          write(io,'(/,"****FATAL: Do not define a constant using an ""="" sign:",A)') func
          STOP
        end if
        call init(func, Constants, statusflag)
!        call init(func, noNames, statusflag)
        if(statusflag == "ok") Then
          nCon=nCon+1
          ceq(nCon)=func
          Constants(nCon) = vname
          ConstantV(nCon) = evaluate(ConstantV)
!          write(io,'("The Constant ",I2," ",A9," is evaluated to:",F12.4," by the formula:",A)')   &
!                     nCon,constants(nCon),ConstantV(nCon),trim(ceq(nCon))
        else
          write(io,'("***WARNING: Constant ",A10," has been removed, as it could not be evaluated by expression:",A)') vname,func
        endif
        call destroyfunc()
      CASE("FEXT"); read(COMMAND,*,err=930,END=930) Fext
        if (LEN_TRIM(Fext).gt.22) then; write(io,'("***FATAL: FEXT is too long: length=",I4)') LEN_TRIM(Fext); stop; ENDIF       
!        
      CASE("FIT "); read(COMMAND,*,err=935,END=935) nFit, maxIts, FitTol  
! Nfit: The number of parameters to fit. 
! maxIts: The maximum number of iterations
! FitTol: Fit will stop when number of iterations reaches maxIts, or the function being minimized becomes less than FitTol
        if (nFit.lt.0 .or. nFit.gt.maxFit) then; write(io,'("***FATAL: Invalid value for nFit=",I4)') nFit; stop; endif
        if (maxIts.le.0 .or. FitTol.le.D0) then
          write(io,'("***FATAL: Parameters for FIT must be >0")'); stop; 
        endif
        do i=1,nFit
 110      ninput=ninput+1          
          read(IN,'(A120)',ERR=925) COMMAND
          COMMAND=adjustl(COMMAND)
          if (echoOn) WRITE(IO,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").eq.1) goto 110 ! skip comment lines
          if (index(COMMAND,"!").ne.0) then; ii=index(COMMAND,"!"); COMMAND(ii:)=""; endif   ! skip comments at end of line      
          read(COMMAND,*,ERR=960,END=960) FitN(i),FitL(i),FitH(i)
          if (echoOn) write(io,'("i=",I2,", FitN(i),FitL(i),FitH(i)=",I3,2F12.4)') i,FitN(i),FitL(i),FitH(i)
          if (FitN(i).lt.1 .or. FitN(i).gt.maxTotPar) then
            write(io,'("***FATAL: Invalid value for FitN(i) parameter number to fit =",I4, &
                       "; maxTotPar=",i3)') FitN(i),maxTotPar; stop; 
          endif
          V1(FitN(i))="*"
        enddo
      CASE("FITO")
! FitOpt(1) print fitmatrices.
! FitOpt(2) print covariance matrix
! FitOpt(3) print minimized function and parameters on each iteration.
! FitOpt(4) output fitted parameters to file "fit.out".
! FitOpt(5) The fit matrix to be diagonalised is printed to file "matrices.dat" each iteration.
! FitOpt(6) If true the parameters limits will be hard, otherwise 100*(P-Plimit)^2 is added to the penalty function.
! FitOpt(7) If true the fit will sum over the degeneracies.
! FitOpt(8) The fit will also give the mean deviation: Sum(|Eexp-Ecal|)/Nexp. 
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for FITO")'); STOP; endif            
          FitOpt(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
!
      CASE("GRID"); read(COMMAND,*,err=930,END=930) n1,n2 ! 
        varNum(1)=n1; varNum(2)=n2
        grid=.true.
      CASE("INTR"); read(COMMAND,*,err=930,END=930) txt,r1 ! 
        CALL TO_UPPER(txt)
        if (txt.eq."FIMAX") then
          if (R1.lt.0.0d+00 .or. R1.gt.1.0d+00) then 
            write(io,'("***FATAL: Invalid value for R1 in INTR FImax R1=",F8.4)') r1; stop
          endif
          if (r1.ne.FImax) write(io,'("***Warning: Internal parameter FImax reset from:",F8.4," to ",F8.4)') FImax,r1
          FImax=R1
        elseif (txt.eq."MJMAX") then
          if (R1.lt.0.0d+00 .or. R1.gt.1.0d+00) then 
            write(io,'("***FATAL: Invalid value for R1 in INTR MJmax R1=",F8.4)') r1; stop
          endif
          if (r1.ne.MJmax) write(io,'("***Warning: Internal parameter MJmax reset from:",F8.4," to ",F8.4)') MJmax,r1
          MJmax=R1
        elseif (txt.eq."CFBIT") then
          if (R1.lt.z .or. R1.gt.10.0d+00) then
            write(io,'("***FATAL: Invalid value for INTR CFbit=",F8.5)') r1; stop
          endif
          if (r1.ne.CFbit) write(io,'("***Warning: Internal parameter CFbit reset from:",F8.5," to ",F8.5)') CFbit,r1
          CFbit = r1
        elseif (txt.eq."EDEGEN") then
          if (R1.lt.z .or.  R1.gt.d1) then
            write(io,'("***FATAL: Invalid value for INTR Edegen=",F8.5)') r1; stop
          endif
          if (r1.ne.Edegen) write(io,'("***Warning: Internal parameter Edegen reset from:",F8.5," to ",F8.5)') Edegen,r1
          Edegen=r1
        elseif (txt.eq."NDECPTS") then
          if (int(R1).lt.0 .or.  int(R1).gt.4) then
            write(io,'("***FATAL: Invalid value for INTR Ndecpts=",i3)') int(r1); stop
          endif
          if (int(r1).ne.Ndecpts) write(io,'("***Warning: Internal parameter Ndecpts reset from:",I2," to ",I2)') Ndecpts,int(r1)
          Ndecpts=int(r1)
        elseif (txt.eq."PFACTOR") then
          if (R1.le.z) then; write(io,'("***FATAL: Invalid value for INTR Pfactor=",F9.4)') r1; stop; endif
          if (r1.ne.Pfactor) write(io,'("***Warning: Internal parameter Pfactor reset from:",F9.4," to ",F9.4)') Pfactor,r1
          Pfactor=r1
        elseif (txt.eq."TESTCRIT") then
          if (R1.lt.z) then; write(io,'("***FATAL: Invalid value for INTR TestCrit=",F8.5)') r1; stop; endif
          if (r1.ne.TestCrit) write(io,'("***Warning: Internal parameter TestCrit reset from:",F8.5," to ",F8.5)') TestCrit,r1
          TestCrit=r1
        ELSE
          write(io,'("***FATAL: Invalid value for T1 in INTR T1 R2=",A20)') txt; stop
        endif
      CASE("LINK"); read(COMMAND,*,err=950,END=950) n1,n2,r1,r2
        if (n1.lt.1 .or. n1.gt.maxTotPar) then; write(io,'("***FATAL: Invalid value for n1=",I2)') n1; stop; endif
        if (n2.lt.1 .or. n2.gt.maxTotPar) then; write(io,'("***FATAL: Invalid value for n2=",I2)') n2; stop; endif
        nlinks=nlinks+1  !  The number of parameters to be linked to another.
        if (nlinks.lt.0 .or. nlinks.gt.50) then; write(io,'("***FATAL: Invalid value for nlinks=",I4)') nlinks; stop; endif
        nlink(1,nlinks)=n1; nlink(2,nlinks)=n2
        rlink(1,nlinks)=r1; rlink(2,nlinks)=r2
      case("EXSO"); read(COMMAND,*,err=930,END=930) ESO2J,ESOkmax,ESONorm,ESOskip
        if (ESO2J.le.0 .or. ESO2J.gt.21) then
          write(io,'("***FATAL: Invalid value for EXSO N1 N2 N3 N4; N1=",I4," must be positive & < 22")') ESO2J; stop
        endif
        call checkESO()
        if (ESOkmax.le.0 .or. ESOkmax.ne.2*INT(ESOkmax/2)) then
          write(io,'("***FATAL: Invalid value for EXSO N1 N2 N3 N4; N2=",I4," must be positive & even")') ESOkmax; stop
        endif
        if (ESONorm.lt.0 .or. ESONorm.gt.1) then
          write(io,'("***FATAL: Invalid value for EXSO N1 N2 N3 N4; N3=",I4," must be 0/1")') ESONorm; stop
        endif
        if (ESOskip.lt.0 .or. ESOskip.gt.1) then
          write(io,'("***FATAL: Invalid value for EXSO N1 N2 N3 N4; N4=",I4," must be 0/1")') ESOskip; stop
        endif
        do k=1,ESOkmax/2
          ii=0
 54       ninput=ninput+1          
          read(IN,'(A200)',ERR=920) COMMAND
          call DECODE(A,N,COMMAND,KEY)
          do q=1,n; ESOBkq(k,ii+q)=A(q); enddo;
          if (echoOn) write(io,'("> k=",I2,",q=",i3," to ",i2," Bkq(ESO)=",(13E12.4))') 2*k, -2*k+ii, -2*k+ii+n-1, &
                                                                                     (ESOBkq(k,ii+q),q=1,n)
          if (n+ii.ne.4*k+1) then
            ii=ii+n
            goto 54
          endif
        enddo  
      case("PRED"); read(COMMAND,*,err=930,END=930) nPRED
        if (nPRED.lt.0 ) then
          write(io,'("***FATAL: Invalid value for PRED N1=",I4)') nPRED; stop
        endif
!      CASE("RORB"); read(COMMAND,*,err=930,END=930) OrbitalsMl,OrbitalsReal  ! TODO calculate the population of ml or real d/f orbitals      
      CASE("ROTL")  !  
! RotLF(1),RotLF(2),RotLF(3) are equations in pEq(307-309) and are the Euler angles alpha, beta, gamma for the rotation of the ligand field.
        do i=1,3   
          COMMAND=ADJUSTL(COMMAND);  ii = index(COMMAND," "); if (ii.eq.0) goto 912; pEq(306+i)=COMMAND(:ii); 
          COMMAND=COMMAND(ii+1:)
        enddo      
      CASE("XREF"); read(COMMAND,*,err=930,END=930) Nref(1),Nref(2)      !  The definition of the coordinate system.
        xref=.true.
!
!----------------------------------------------------------------------------------------------------------        
! Debugging Commands:      
      CASE("CHCK")  ! Obsolete
        write(io,'("****FATAL: CHCK command obsolete; Use new Commands: CHK1, CHK2")'); STOP 
      CASE("CHK1")  !  1-10: basis,Umat,V11,EEmat,abgMat,MnPnmat,Tn,AOM_LF,Eigenvalues,group 
                    !  Goes into "debug.dat" file. Warning, must place before CONF command to do some checks.
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for CHK1")'); STOP; endif            
          check1(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("CHK2")  !  1-10: skip diag, Check WF are Kramers doublets (odd e- systems only), testLS
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for CHK2")'); STOP; endif            
          check2(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("ECHO")
        read(COMMAND,*,err=125,END=125) echoOn
        goto 126
 125    write(IO,'("Invalid input, ECHO must be T or F")')  
 126    if (echoOn) WRITE(IO,'(">Input Echo:")')
      CASE("FNCR") ! Deliberate errors introduced for: Un,V11,Fn,abg,MnPn,Tn,whole t2
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.8) then; write(io,'("****FATAL: Too much input for FNCR")'); STOP; endif            
          UseFncrossErrors(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("PRNT")  
        write(io,'("****FATAL: PRNT command obsolete; Use new Commands: PRT1, PRT2, PRT3")'); STOP 
      CASE("PRT1")  !  1-10: Uk,V11,Fn,Mn,Pn, Sn,Tn,gik,NU,section  in |SL> basis
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for PRT1")'); STOP; endif            
          printMat1(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
        if (printMat1(10)) then; ninput=ninput+1; read(IN,*,err=982) (PrintRng1(i),i=1,4); endif
      CASE("PRT2")  !  1-10: matrix,(*/.)matrix,<Uk>,<L/S>,<kL+gS>,4*NU,Section in |SLJ> basis
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for PRT2")'); STOP; endif            
          printMat2(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
        if (printMat2(10)) then; ninput=ninput+1; read(IN,*,err=982) (PrintRng2(i),i=1,4); endif
      CASE("PRT3")  !  1-10: matrix,(*/.)matrix,<Uk>list,<Uk>,<ED>,<L/S>,<kL+gS>,LForig+PRED,NU,Section in |SLJMJ> basis
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for PRT3")'); STOP; endif            
          printMat3(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
        if (printMat3(10)) then; ninput=ninput+1; read(IN,*,err=982) (PrintRng3(i),i=1,4); endif
      CASE("PRT4")  !  1-10: ESO info, ESO mat,8*NU in |MJ> basis
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for PRT4")'); STOP; endif            
          printMat4(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
         if (printMat4(10)) then; ninput=ninput+1; read(IN,*,err=982) (PrintRng4(i),i=1,4); endif
      CASE("WARN")  !  1-10: testBkq,testAlpt,3-10 NU
        i=0; ii=scan(COMMAND,"TFtf")
        DO while (ii.ne.0)   !
          i=i+1
          if (i.gt.10) then; write(io,'("****FATAL: Too much input for WARN")'); STOP; endif            
          warnings(i)=(COMMAND(ii:ii).eq."T".or.COMMAND(ii:ii).eq."t"); COMMAND(ii:ii)=" "
          ii=scan(COMMAND,"TFtf") 
        enddo
      CASE("TEST")
        TestOutput=.true.
        OPEN(iTest,FILE='test.out',STATUS='OLD',ACTION='WRITE',POSITION='APPEND',Err=80)
        goto 10
  80    write(io,'("****FATAL: The file ""test.out"" must exist for the command TEST")')
        write(io,'("****FATAL: Create an empty ""test.out"" file if running for the first time")'); stop
      CASE("TIME")
        iTime=.True.
!
!----------------------------------------------------------------------------------------------------------        
! Hidden Commands:      
      CASE("ATOM")    ! Prediagonalisation with atomic basiswith different parameters
        doAtomic=.true.
        do i=1,20; read(COMMAND,*,err=952,END=952) Atomic(i); enddo
      CASE("LANC")    ! Lanczos method of diagonalisation
        read(COMMAND,*,err=951,END=951) lancIt,RLB,RUB
        if (lancIt.le.0) then 
          write(io,'("****FATAL: Invalid number of Lanczos iteration in LANC command: lancIt=",i4)') lancIt; stop
        endif  
        if (RUB.le.RLB) then 
          write(io,'("****FATAL: Invalid energy range in LANC command: RLB,RUB=",2F10.1)') RLB,RUB; stop
        endif  
!        
!---------------------------------------------------------------------------------------------------------- 
! Parameters:
!     CASE("EAVE");  read(COMMAND,*) EAVE
!     CASE("F2  ");  read(COMMAND,*) F2; Fn=.true.
!     CASE("F4  ");  read(COMMAND,*) F4; Fn=.true.
!     CASE("F6  ");  read(COMMAND,*) F6; Fn=.true.
!     CASE("B   ");  read(COMMAND,*) Bracah; Racah=.true.
!     CASE("C   ");  read(COMMAND,*) Cracah; Racah=.true.      
!     CASE("ALPH");  read(COMMAND,*) ALPHA
!     CASE("BETA");  read(COMMAND,*) BETA
!     CASE("GAMM");  read(COMMAND,*) GAMMA
!     CASE("ZETA");  read(COMMAND,*) ZETA
!     CASE("M0  ");  read(COMMAND,*) M0
!     CASE("M2  ");  read(COMMAND,*) M2
!     CASE("M4  ");  read(COMMAND,*) M4
!     CASE("P2  ");  read(COMMAND,*) P2
!     CASE("P4  ");  read(COMMAND,*) P4
!     CASE("P6  ");  read(COMMAND,*) P6
!     CASE("T2  ");  read(COMMAND,*) T2
!     CASE("T3  ");  read(COMMAND,*) T3
!     CASE("T4  ");  read(COMMAND,*) T4
!     CASE("T6  ");  read(COMMAND,*) T6
!     CASE("T7  ");  read(COMMAND,*) T7
!     CASE("T8  ");  read(COMMAND,*) T8

      CASE("EAVE");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(1)
      CASE("F2  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(2); Fn=.true.
      CASE("F4  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(3); Fn=.true.
      CASE("F6  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(4); Fn=.true.
      CASE("B   ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(2); Racah=.true.; ! write(io,'("pEq(2)=",A)') pEq(2)
      CASE("C   ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(3); Racah=.true.; ! write(io,'("pEq(3)=",A)') pEq(3)      
      CASE("ZETA");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(5)
      CASE("ALPH");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(6)
      CASE("BETA");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(7)
      CASE("GAMM");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(8)
      CASE("M0  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(9)
      CASE("M2  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(10)
      CASE("M4  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(11)
      CASE("P2  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(12)
      CASE("P4  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(13)
      CASE("P6  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(14)
      CASE("T2  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(15)
      CASE("T3  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(16)
      CASE("T4  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(17)
      CASE("T6  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(18)
      CASE("T7  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(19)
      CASE("T8  ");  read(COMMAND,'(A100)',ERR=990,END=990) pEq(20)
      CASE("CF  ")  ! Following lines will be crystal field Bkq coefficents given explicitly
        if (.not.CONFcom) then; write(io,'("The CONF command must be specified before the CF command")'); stop; endif 
        if (LFtype.ne."    ") then; write(io,'("Cannot use CF, ligand field is already type:",A4)') LFtype; stop; endif  
        LFtype="CF  "
        Do k=2,2*Lvalue,2   ! don't read B00
          do q=0,k
            ninput=ninput+1
            i=Bkq_index(k,q)
!              if (q.eq.0) READ(IN,*,ERR=910) BKQR(i)
!              if (q.gt.0) READ(IN,*,ERR=910) BKQR(i),BKQI(i)
!              if (echoOn .and. q.eq.0) write(io,'("B",I1,"0=",A)') k,COMMAND
!              if (echoOn .and. q.gt.0) write(io,'("B",2I1,"=",A)') k,q,COMMAND
            READ(IN,'(A120)',ERR=910) COMMAND
            if (echoOn) write(io,'(">",A)') trim(COMMAND)
            if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
              ii=index(COMMAND,"!"); COMMAND(ii:)=""
            endif
            COMMAND=ADJUSTL(COMMAND);  ii = index(COMMAND," "); if (ii.eq.1) goto 910; pEq(19+i)=COMMAND(:ii); 
            COMMAND=COMMAND(ii+1:)
            if (q.gt.0) then
              COMMAND=ADJUSTL(COMMAND);  ii = index(COMMAND," "); ! write(io,'("k,q,i,ii=",4I3)') k,q,i,ii; 
              if (ii.eq.1) goto 910; pEq(39+i)=COMMAND(:ii); ! write(io,'("pEq("I2,")=",A)')39+i,pEq(39+i)
            endif
          enddo
        enddo
      CASE("AOM")
        if (LFtype.ne."    ") then; write(io,'("Cannot use AOM, ligand field is already type:",A4)') LFtype; stop; endif  
        read(COMMAND,*) NLIGANDS     ! N: -ve (cart ooordinates), +ve (AOM angles)
        if (NLIGANDS.ne.0) then 
          LFtype="AOM "
!          WRITE(*,'("About to read ",I4," ligands")') NLIGANDS
          do iL=1,abs(NLIGANDS)
            ninput=ninput+1
            if (NLIGANDS.gt.0) then
              AOM_angles=.true.
 140          READ(IN,'(A120)',ERR=905) COMMAND
              COMMAND=adjustl(COMMAND)
              if (echoOn) write(io,'(">",A)') trim(COMMAND)
              if (index(COMMAND,"!").eq.1) goto 140 ! skip comment lines
              if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
                ii=index(COMMAND,"!"); COMMAND(ii:)=""
              endif
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; Lname(iL)=COMMAND(:i); COMMAND=COMMAND(i+1:)
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(20+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! Esig
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(40+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! EpiX
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(60+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! EpiY
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(80+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! theta
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.le.1) goto 905; pEq(100+iL)=COMMAND(:i);COMMAND=COMMAND(i+1:) ! phi
              pEq(120+iL)="0"   ! default value for Chi angles.         i.eq.1 means that no space found or space found in first position.
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.le.1) goto 150; pEq(120+iL)=COMMAND(:i);COMMAND=COMMAND(i+1:) ! chi
!
!              READ(IN,*,ERR=905) Lname(i),eSigma(i),ePiX(i),ePiY(i),theta(i),phi(i)
!              if (echoOn) WRITE(IO,'(A4,4X,3(2X,F8.2),4X,3(2X,F6.2))') Lname(i),eSigma(i),ePiX(i),ePiY(i),theta(i),phi(i),chi(i)
 150          if (echoOn) WRITE(IO,'(A4,4X,6(2X,A))') Lname(iL),trim(pEq(20+iL)),trim(pEq(40+iL)),trim(pEq(60+iL)),  &
                                                                trim(pEq(80+iL)),trim(pEq(100+iL)),trim(pEq(120+iL))
            ELSE  ! cartesian coordinates cannot be written as formula
              AOM_angles=.false.
 155          READ(IN,'(A120)',ERR=905) COMMAND
              COMMAND=adjustl(COMMAND)
              if (echoOn) write(io,'(">",A)') trim(COMMAND)
              if (index(COMMAND,"!").eq.1) goto 155 ! skip comment lines
              if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
                ii=index(COMMAND,"!"); COMMAND(ii:)=""
              endif
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; Lname(iL)=COMMAND(:i); COMMAND=COMMAND(i+1:)
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(20+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! Esig
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(40+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! EpiX
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(60+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! EpiY
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(80+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! x
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.le.1) goto 905; pEq(100+iL)=COMMAND(:i);COMMAND=COMMAND(i+1:) ! y
              COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.le.1) goto 905; pEq(120+iL)=COMMAND(:i);COMMAND=COMMAND(i+1:) ! z
              if (echoOn) WRITE(IO,'(A4,4X,6(2X,A),3F12.6)') Lname(iL),trim(pEq(20+iL)),trim(pEq(40+iL)),trim(pEq(60+iL)), &
                                                                trim(pEq(80+iL)),trim(pEq(100+iL)),trim(pEq(120+iL))
!              READ(IN,*,err=905) Lname(i),eSigma(iL),ePiX(i),ePiY(i),xLig(iL),yLig(iL),zLig(iL)
!              if (echoOn) WRITE(IO,'(A4,4X,3(2X,F8.2),4X,3(2X,F6.2))') Lname(iL),eSigma(iL),ePiX(i),ePiY(i),xLig(iL),yLig(iL),zLig(iL)
            endif
          enddo !  iL 
          NLIGANDS=abs(NLIGANDS)
        ENDIF
      CASE("AOMX")  ! extension to AOM: delS,delC, PsiS,PsiC parameters
        if (LFtype.ne."AOM ") then; write(io,'("Must have AOM command before AOMX")'); stop; endif  
        read(COMMAND,*) N     !
        if (N.ne.NLIGANDS)  then; write(io,'("Must have same number of ligands in AOM and AOMX commands")'); stop; endif
        if (NLIGANDS.gt.10) then; write(io,'("Maximum of 10 ligands with AOMX commands")'); stop; endif
        AOMX=.true.
        
!       WRITE(*,'("About to read ",I4," ligands")') NLIGANDS
        do iL=1,NLIGANDS
          ninput=ninput+1
   156    READ(IN,'(A120)',ERR=905) COMMAND
          COMMAND=adjustl(COMMAND)
          if (echoOn) write(io,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").eq.1) goto 156 ! skip comment lines
          if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
            ii=index(COMMAND,"!"); COMMAND(ii:)=""
          endif
          COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(140+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! EDelS
          COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(150+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! EDelC
          COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(160+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! EPsiS
          COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 905; pEq(170+iL)=COMMAND(:i); COMMAND=COMMAND(i+1:) ! EPsiC
          if (echoOn) WRITE(IO,'(A4,4X,6(2X,A))') Lname(iL),trim(pEq(140+iL)),trim(pEq(150+iL)),trim(pEq(160+iL)),trim(pEq(170+iL))
        enddo !  iL 
      CASE("LF1E")  ! Following lines will be the one electron LF matrix elements in terms of real orbitals
        if (LFtype.ne."    ") then; write(io,'("Cannot use LF1E, ligand field is already type:",A4)') LFtype; stop; endif  
        if (.not.CONFcom) then; write(io,'("The CONF command must be specified before the LF1E command")'); stop; endif 
        LFtype="LF1E"; ii=20
        Do k=1,2*Lvalue+1   !  
          ninput=ninput+1
          READ(IN,'(A120)',ERR=906) COMMAND
          if (echoOn) write(io,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
            ii=index(COMMAND,"!"); COMMAND(ii:)=""
          endif
          do j=1,k
            ii=ii+1
           COMMAND=ADJUSTL(COMMAND);  i = index(COMMAND," "); if (i.eq.0) goto 906; pEq(ii)=COMMAND(:i); COMMAND=COMMAND(i+1:)
          enddo
        enddo
      CASE("CCF ")  ! Correlated Crystal field specified in the following line.
        if (.not.CONFcom) then; write(io,'("The CONF command must be specified before the CF command")'); stop; endif ! need the number of electrons
 !       if (LFtype.ne."    ") then; write(io,'("Cannot use CF, ligand field is already type:",A4)') LFtype; stop; endif  
        read(COMMAND,*) CCFtype,N_CCF     ! CCFtype=1 (delta-function model),2 (spin-CCF),3 (general CCF) N_CCF number of CCF parameters specified. 
        if (CCFtype.lt.1 .or. CCFtype.gt.3) then
          write(io,'("***FATAL: Invalid value for CCFtype=",I3," in CCF command")') CCFtype; stop; endif  
        if (N_CCF.lt.1 .or. N_CCF.gt.60) then
          write(io,'("***FATAL: Invalid value for N_CCF=",I3," in CCF command (Max=60)")') N_CCF; stop; endif  
        do ii=1,N_CCF
          ninput=ninput+1
          READ(IN,'(A120)',ERR=915) COMMAND
          if (echoOn) write(io,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
            i=index(COMMAND,"!"); COMMAND(i:)=""
          endif                             !    i             k            q
          if (CCFtype.eq.1) then  ! delta-function model
                                            CCFindex(ii,1)=1
                            read(COMMAND,*)                CCFindex(ii,2),CCFindex(ii,3),pEq(140+ii) 
          endif                  
          if (CCFtype.eq.3) read(COMMAND,*) CCFindex(ii,1),CCFindex(ii,2),CCFindex(ii,3),pEq(140+ii)
        enddo
 
        OPEN(iccf,FILE=trim(F_E_PATH)//'\ccf.dat',STATUS='OLD',ERR=916)
      CASE("ALTP") 
        if (JuddOfelt) then; write(IO,'("***FATAL: cannot have both JUDO and ALTP commands.")'); STOP; endif 
        read(COMMAND,*) IntN1,IntN2 ! IntN1=1/2/3/4 for Altp/Blki/Tli/Cki; IntN2=1(pure real)/2(pure imag)/3(complex)
        if ((IntN1.eq.3 .or. IntN1.eq.4) .and. NLIGANDS.eq.0) then
          write(IO,'("***FATAL: Must specifiy the number of ligands in the AOM command ",/, &
                     " before the ALTP command for IntN1=3 (AOM intensities)")'); STOP
        endif 
        If (IntN1.lt.1 .or. IntN1.gt.4) then; write(io,'("Invalid value IntN1 (=",I2,") in ALPT command")') IntN1; stop; endif 
        If (IntN2.lt.1 .or. IntN2.gt.3) then; write(io,'("Invalid value IntN2 (=",I2,") in ALPT command")') IntN2; stop; endif 
        do i=1,3; do j=1,3; do k=1,8; AltpR(i,j,k)=z; AltpI(i,j,k)=z; enddo; enddo; enddo
        do i=1,15; do j=1,3; BlkiR(i,j)=z; BlkiI(i,j)=z; enddo; enddo
        if (IntN1.eq.1) then
          do k=1,3
            do i=1,3
              ninput=ninput+1
              read(IN,"(A120)",END=900) COMMAND
              if (echoOn) WRITE(IO,'(">",A)') trim(COMMAND)
              if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
                ii=index(COMMAND,"!"); COMMAND(ii:)=""
              endif
!              read(COMMAND,*,ERR=945) (AltpR(k,i,j),AltpI(k,i,j),j=1,i+1+(k-1)*2)
              if (k.eq.1) then
                if (i.eq.1) then; i1=200; j1=2; endif
                if (i.eq.2) then; i1=204; j1=3; endif
                if (i.eq.3) then; i1=210; j1=4; endif
              elseif (k.eq.2) then   
                if (i.eq.1) then; i1=218; j1=4; endif
                if (i.eq.2) then; i1=226; j1=5; endif
                if (i.eq.3) then; i1=236; j1=6; endif
              elseif (k.eq.3) then   
                if (i.eq.1) then; i1=248; j1=6; endif
                if (i.eq.2) then; i1=260; j1=7; endif
                if (i.eq.3) then; i1=274; j1=8; endif
              endif
              if (IntN2.eq.1) then
                do j=1,2*j1
                  COMMAND=ADJUSTL(COMMAND); ii = index(COMMAND," "); if (ii.eq.1) goto 945; 
                  pEq(i1+2*j-1)=COMMAND(:ii); COMMAND=COMMAND(ii+1:)
                enddo
!                if (k.eq.1 .and. i.eq.1) read(COMMAND,'(A100)',ERR=945) (pEq(200+(j-1)*2+1),j=1,2)
!                if (k.eq.1 .and. i.eq.2) read(COMMAND,'(A100)',ERR=945) (pEq(204+(j-1)*2+1),j=1,3)
!                if (k.eq.1 .and. i.eq.3) read(COMMAND,'(A100)',ERR=945) (pEq(210+(j-1)*2+1),j=1,4)
!                if (k.eq.2 .and. i.eq.1) read(COMMAND,'(A100)',ERR=945) (pEq(218+(j-1)*2+1),j=1,4)
!                if (k.eq.2 .and. i.eq.2) read(COMMAND,'(A100)',ERR=945) (pEq(226+(j-1)*2+1),j=1,5)
!                if (k.eq.2 .and. i.eq.3) read(COMMAND,'(A100)',ERR=945) (pEq(236+(j-1)*2+1),j=1,6)
!                if (k.eq.3 .and. i.eq.1) read(COMMAND,'(A100)',ERR=945) (pEq(248+(j-1)*2+1),j=1,6)
!                if (k.eq.3 .and. i.eq.2) read(COMMAND,'(A100)',ERR=945) (pEq(260+(j-1)*2+1),j=1,7)
!                if (k.eq.3 .and. i.eq.3) read(COMMAND,'(A100)',ERR=945) (pEq(274+(j-1)*2+1),j=1,8)
              elseif (IntN2.eq.2) then
                do j=1,2*j1
                  COMMAND=ADJUSTL(COMMAND); ii = index(COMMAND," "); if (ii.eq.1) goto 945; 
                  pEq(i1+2*j)=COMMAND(:ii); COMMAND=COMMAND(ii+1:)
                enddo
!                if (k.eq.1 .and. i.eq.1) read(COMMAND,'(A100)',ERR=945) (pEq(200+j*2),j=1,2)
!                if (k.eq.1 .and. i.eq.2) read(COMMAND,'(A100)',ERR=945) (pEq(204+j*2),j=1,3)
!                if (k.eq.1 .and. i.eq.3) read(COMMAND,'(A100)',ERR=945) (pEq(210+j*2),j=1,4)
!                if (k.eq.2 .and. i.eq.1) read(COMMAND,'(A100)',ERR=945) (pEq(218+j*2),j=1,4)
!                if (k.eq.2 .and. i.eq.2) read(COMMAND,'(A100)',ERR=945) (pEq(226+j*2),j=1,5)
!                if (k.eq.2 .and. i.eq.3) read(COMMAND,'(A100)',ERR=945) (pEq(236+j*2),j=1,6)
!                if (k.eq.3 .and. i.eq.1) read(COMMAND,'(A100)',ERR=945) (pEq(248+j*2),j=1,6)
!                if (k.eq.3 .and. i.eq.2) read(COMMAND,'(A100)',ERR=945) (pEq(260+j*2),j=1,7)
!                if (k.eq.3 .and. i.eq.3) read(COMMAND,'(A100)',ERR=945) (pEq(274+j*2),j=1,8)
              elseif (IntN2.eq.3) then
                do j=1,2*j1
                  COMMAND=ADJUSTL(COMMAND); ii = index(COMMAND," "); if (ii.eq.1) goto 945; 
                  pEq(i1+j)=COMMAND(:ii); COMMAND=COMMAND(ii+1:)
                enddo
              endif ! IntN2              
            enddo ! i
          enddo ! k 
        else if (IntN1.eq.2) then
          do i=1,15
            ninput=ninput+1
            read(IN,"(A120)",END=900) COMMAND
            if (echoOn) WRITE(IO,'(">",A)') trim(COMMAND)
            if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
              ii=index(COMMAND,"!"); COMMAND(ii:)=""
            endif
            do j=1,3
              if (IntN2.eq.1) read(COMMAND,*,ERR=945) (pEq(200+(i-1)*6+(k-1)*2+1),k=1,3)
              if (IntN2.eq.2) read(COMMAND,*,ERR=945) (pEq(200+(i-1)*6+(k-1)*2+2),k=1,3)
              if (IntN2.eq.3) read(COMMAND,*,ERR=945) (pEq(200+(i-1)*6+k),k=1,6)
            enddo ! j
          enddo ! i 
        else if (IntN1.eq.3) then
          do i=1,NLIGANDS
            ninput=ninput+1
            read(IN,"(A120)",END=900) COMMAND
            if (echoOn) WRITE(IO,'(">",A)') trim(COMMAND)
            if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
              ii=index(COMMAND,"!"); COMMAND(ii:)=""
            endif
            read(COMMAND,*,ERR=945) (pEq(200+(i-1)*9+k),k=1,9)
          enddo ! i 
        else if (IntN1.eq.4) then
          do i=1,NLIGANDS
            ninput=ninput+1
            read(IN,"(A120)",END=900) COMMAND
            if (echoOn) WRITE(IO,'(">",A)') trim(COMMAND)
            if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
              ii=index(COMMAND,"!"); COMMAND(ii:)=""
            endif
            read(COMMAND,*,ERR=945) (pEq(200+(i-1)*7+k),k=1,7)
          enddo ! i 
        ENDIF ! IntN1
      CASE("OFFS")  ! Following lines will be offset the basis diagonal in |LSJ> or full |LSJMJ>.
        if (.not.CONFcom) then; write(io,'("The CONF command must be specified before the OFFS command")'); stop; endif 
        read(COMMAND,*) OffsetType,OffsetN
        if (OffsetType.lt.1 .or. OffsetType.gt.2) then
          write(io,'("The invalid number for OffsetType =",I2)') OffsetType; stop
        endif
        if (OffsetN.lt.1 .or. OffsetN.gt.50) then
          write(io,'("The invalid number for OffsetN =",I4)') OffsetN; stop
        endif   
        Do i=1,OffsetN   
          ninput=ninput+1
          READ(IN,'(A120)',ERR=904) COMMAND
          if (echoOn) write(io,'(">",A)') trim(COMMAND)
          if (index(COMMAND,"!").ne.0) then  ! skip comments at end of line      
            ii=index(COMMAND,"!"); COMMAND(ii:)=""
          endif
          if (OffsetType.eq.1) then
            read(COMMAND,*,ERR=904,END=904) iOffset(i,1), pEq(350+i)
!            write(io,'("iOffset:",I4,5X,A)') iOffset(i,1), pEq(350+i)
         else
            read(COMMAND,*,err=904,END=904) iOffset(i,1), iOffset(i,2), pEq(350+i)
!            write(io,'("iOffset:",2I4,5X,A)') iOffset(i,1),iOffset(i,2), pEq(350+i)
          endif            
!           write(io,'("Parameter label:",A9)') PLab(350+i)         
        enddo
!---------------------------------------------------------------------------------------------------------- 
! End of input:
      CASE("END ")
        if (.not.CONFcom) goto 902
!
        WRITE(IO,'(" Input from file:",A,/)') trim(INFILE)
        if (.not.Complementary) WRITE(IO,'(" Calculation: ",A1,"^",I2," electrons")') e_lab,Nelectrons
        if (     Complementary) WRITE(IO,'(" Calculation: ",A1,"^",I2," electrons; complementary ",   &
                                A1,"^",I2," configuration used.")') e_lab,2*(2*Lvalue+1)-Nelectrons, e_lab,Nelectrons
        if (nPRED.ne.0) then
          WRITE(IO,'(" The calculation will do a pre-diagonalisation using only the atomic parameters:",/,  &
                     " The eigenvectors from this calculation will then be used when the ligand field is included.")')
          if (nPRED.ne.TS_full(Nelectrons)) WRITE(IO,'(" A reduced (from full: ",i4,") Pre-Diagonalised basis of",i4, &
                                                       " functions will be used.",/)') TS_full(Nelectrons),nPRED
          if (nPRED.eq.TS_full(Nelectrons)) WRITE(IO,'(" The full: ",i4," Pre-Diagonalised basis will be used.",/)') &
                                                       nPRED
        endif
                               
        if (JuddOfelt) option(8)=.true.  ! J-O; skip main calculation.
!        if (ESO2J.ne.0) option(8)=.true. ! Calculation using Extended Stevens Operators 
        if (.not.FullMJBasis) then
          option(8)=.true.                ! skip main calculation.
          WRITE(IO,'(" A reduced |SLJ> basis of size:",I3," will be used instead of the full ", &
                                         "|LSJMJ> basis of size:",I4)') TSLJ_f(Nelectrons),TS_full(Nelectrons)
          if (LFtype.ne."    ") WRITE(IO,'(" Warning: the ligand/crystal field will not be used.")')
          write(io,'(/)')
        endif
!        
        test=.false.; do k=1,10; if (Option(k)) test=.true.; enddo
        if (test) Write(IO,'(" Calculation Options: OPTN(1-",I2,") = ",4(5L2,2X))') 10,(Option(I),I=1,10)
        if (option(3) .and. option(4)) then
          WRITE(IO,'(5X,"****WARNING: Cannot have both long and short % free ion labels. option(4) set to false.")')
          option(4)=.false. 
        endif
        if (option(6) .and. option(7)) then
          WRITE(IO,'(5X,"****WARNING: Cannot have both long and short % MJ labels. option(7) set to false.")')
          option(7)=.false. 
        endif
        if (option(1)) WRITE(IO,'(5X,"Energies relative to lowest energy=0.")') 
        if (option(1) .and. EAVE.ne.D0) WRITE(IO,'("****WARNING: EAVE will have no effect when OPTN(1)=T.")') 
        if (option(2)) WRITE(IO,'(5X,"Calculate spin projections.")') 
        if (option(3)) WRITE(IO,'(5X,"Calculate % of free ion terms (short 40 char).")') 
        if (option(4)) WRITE(IO,'(5X,"Calculate % of free ion terms (long 80 char).")')
        if (Lvalue.eq.3.and.Nelectrons.ne.1.and.Nelectrons.ne.13) then        
          if (option(5)) WRITE(IO,'(5X,"The Hss spin-spin interactions (also Mn parameters) are NOT included.")') 
          if (.not.option(5)) WRITE(IO,'(5X,"The Hss spin-spin interactions (using the Mn parameters) were included.")') 
        endif
        if (option(6))  WRITE(IO,'(5X,"Calculate % of MJ (short 40 char).")') 
        if (option(7))  WRITE(IO,'(5X,"Calculate % of MJ (long 80 char).")') 
        if (option(8))  WRITE(IO,'(5X,"Skip Writing matrix & diagonalisation **********")') 
        if (option(9))  WRITE(IO,'(5X,"Write matrix; Skip diagonalisation **********")')
        if (option(10)) write(IO,'("****WARNING: option(10) not used")')                 
!
        call checkPrint() 
!
        test=.false.; do k=1,10; if (OUTP(k)) test=.true.; enddo
        if (test) Write(IO,'(" Output Options:      OUTP(1-",I2,") = ",4(5L2,2X))') 10,(OUTP(I),I=1,10)
        if (OUTP(1)) write(IO,'(5X,"The ligand distance/angle matrix is printed (for AOM).")') 
        if (OUTP(2)) write(IO,'(5X,"Print equivalent crystal field (Bkq parameters).")') 
        if (OUTP(3)) write(IO,'(5X,"The one electron AOM matrix and energies will be printed.")') 
        if (OUTP(4)) write(IO,'(5X,"Print equivalent intensity parameters, Blki if Altp given, etc.")') 
        if (OUTP(5)) write(IO,'(5X,"If ROTL used, print CF parameters before & after rotation.")') 
        if (OUTP(6)) write(IO,'(5X,"Degenerate energies are printed.")') 
        if (OUTP(7)) write(IO,'(5X,"Print info on matrix size.")') 
        if (OUTP(8)) write(IO,'(5X,"All parameters are printed each time for multiple calculations.")') 
        if (OUTP(9)) write(IO,'(5X,"The one electron AOM energies and eigenvectores will be printed.")') 
!
        test=.false.; do k=1,10; if (fastMat(k)) test=.true.; enddo
        if (test) Write(IO,'(" FAST Options:        FAST(1-",I2,") = ",4(5L2,2X))') 10,(fastMat(I),I=1,10)
        if (fastMat(1)) write(IO,'(5X,"Higher MJ values are removed from basis.")') 
        if (fastMat(2)) write(IO,'(5X,"Time-reversal symmetry exploited.")') 
        if (fastMat(3)) write(IO,'(5X,"Other Kramers component regenerated.")') 
!
        test=.false.; do k=1,10; if (PLOTout1(k)) test=.true.; enddo
        if (test) then
          Write(IO,'(" Plot Output Options: PLT1(1-",I2,") = ",4(5L2,2X))') 10,(PLOTout1(I),I=1,10)
          fnamep="plot1"//TRIM(Fext)//".dat"  ! will be plot1.dat if FExt not set.
          write(IO,'(5X,"Data will be written to a plotfile:",A30)') fnamep 
          OPEN(IP1,FILE=fnamep,STATUS='UNKNOWN')
          write(IP1,'("Plotfile made from calculation made using the input file:",A)') trim(INFILE)
          write(IP1,'(4(5L2,2X)," !  Eng,S,Sint,MJ,MD,  ED,Fit,Bkq,g-val,OrbE")') (PLOTout1(I),I=1,10)
        endif
        if (PLOTout1(1)) write(IO,'(5X,"Energies written.")') 
        if (PLOTout1(2)) write(IO,'(5X,"+Spin fractions written.")') 
        if (PLOTout1(3)) write(IO,'(5X,"+spin-allowed intensities written.")') 
        if (PLOTout1(4)) write(IO,'(5X,"+MJ fractions written.")') 
        if (PLOTout1(5)) write(IO,'(5X,"Mag.Dipole transitions written.")') 
        if (PLOTout1(6)) write(IO,'(5X,"Elect.Dipole transitions written.")') 
        if (PLOTout1(7)) write(IO,'(5X,"Fit Output written.")') 
        if (PLOTout1(8)) write(IO,'(5X,"Bkq parameters written.")') 
        if (PLOTout1(9)) write(IO,'(5X,"g-values written.")') 
        if (PLOTout1(10))write(IO,'(5X,"d/f-orbital energies.")') 
        if (PLOTout1(2)) option(2)=.true.  ! Calculate spin fraction
        if (PLOTout1(4)) option(7)=.true.  ! Calculate MJ fraction
!
        test=.false.; do k=1,10; if (PLOTout2(k)) test=.true.; enddo
        if (test) then
          Write(IO,'(" Plot Output Options: PLT2(1-",I2,") = ",4(5L2,2X))') 10,(PLOTout2(I),I=1,10)
          fnamep="plot2"//TRIM(Fext)//".dat"  ! will be plot.dat if Fext not set.
          write(IO,'(5X,"Data will be written to a plotfile:",A30)') fnamep 
          OPEN(IP2,FILE=fnamep,STATUS='UNKNOWN')
          write(IP2,'("Plotfile made from calculation made using the input file:",A)') trim(INFILE)
          write(IP2,'(4(5L2,2X)," !  MinFunc,9*NU")') (PLOTout2(I),I=1,10)
        endif
        if (PLOTout2(1)) write(IO,'(5X,"Minimized function written.")') 
!
        test=.false.; do k=1,10; if (printMat1(k).or.printMat2(k).or.printMat3(k).or.printMat4(k).or.FitOpt(5)) test=.true.; enddo
        if (test) then
          fnamem="matrices"//TRIM(Fext)//".dat"
          OPEN(Imatrix,FILE=fnamem,STATUS='UNKNOWN')
          write(Imatrix,'("Matrix file made from calculation using the input file:",A)') trim(INFILE)
        endif  
        test=.false.; do k=1,10; if (printMat1(k)) test=.true.; enddo
        if (test) then
          Write(IO,'(" Matrix Options:      PRT1(1-",I2,") = ",4(5L2,2X))') 10,(printMat1(I),I=1,10)
          write(IO,'(5X,"The following matrices in |SL> basis will be printed to the file ",A30)') fnamem
          write(Imatrix,'("This file contains the following matrices from PRT1:")') 
          if (printMat1(1)) then; write(Imatrix,'("Uk  matrix.")');  write(IO,'(5X,"Uk  matrix.")'); endif
          if (printMat1(2)) then; write(Imatrix,'("V11 matrix.")');  write(IO,'(5X,"V11 matrix.")'); endif
          if (printMat1(3)) then; write(Imatrix,'("Fn matrix.")');   write(IO,'(5X,"Fn matrix.")'); endif
          if (printMat1(4)) then; write(Imatrix,'("Mn matrix.")');   write(IO,'(5X,"Mn matrix.")'); endif
          if (printMat1(5)) then; write(Imatrix,'("Pn matrix.")');   write(IO,'(5X,"Pn matrix.")'); endif
          if (printMat1(6)) then; write(Imatrix,'("Sn matrix.")');   write(IO,'(5X,"Sn matrix.")'); endif
          if (printMat1(7)) then; write(Imatrix,'("Tn matrix.")');   write(IO,'(5X,"Tn matrix.")'); endif
          if (printMat1(8)) then; write(Imatrix,'("gik matrix.")');   write(IO,'(5X,"gik matrix.")'); endif
          if (printMat1(9)) write(IO,'("****WARNING: printMat1(",i2,") not used")') k                 
          if (printMat1(10)) then; write(Imatrix,'("matrix subset printed; i1,i2, j1,j2:":4i3)') (printRng1(i),i=1,4);
                                   write(IO,'("matrix subset printed; i1,i2, j1,j2:":4i3)') (printRng1(i),i=1,4); endif
        endif  
!
        test=.false.; do k=1,10; if (printMat2(k)) test=.true.; enddo
        if (test) then
          Write(IO,'(" Matrix Options:      PRT2(1-",I2,") = ",4(5L2,2X))') 10,(printMat2(I),I=1,10)
          write(IO,'(5X,"The following matrices in |SLJ> basis will be printed to the file ",A30)') fnamem
          write(Imatrix,'("This file contains the following matrices from PRT2:")') 
          if (printMat2(1)) then; write(Imatrix,'("Full matrix in |SLJ> basis.")')
                                  write(IO,  '(5X,"Full matrix in |SLJ> basis.")');   endif
          if (printMat2(2)) then; write(Imatrix,'("*/. matrix in |SLJ> basis.")') 
                                  write(IO,  '(5X,"*/. matrix in |SLJ> basis.")');    endif
          if (printMat2(3)) then; write(Imatrix,'("Uk matrices in |SLJ> basis.")')
                                  write(IO,  '(5X,"Uk matrices in |SLJ> basis.")'); endif
          if (printMat2(4)) then; write(Imatrix,'("L,S matrices in |SLJ> basis.")') 
                                  write(IO,  '(5X,"L,S matrices in |SLJ> basis.")'); endif
          if (printMat2(5)) then; write(Imatrix,'("mu matrices in |SLJ> basis.")') 
                                  write(IO,  '(5X,"mu matrices in |SLJ> basis.")'); endif
          do k=6,9; if (printMat2(k)) write(IO,'("****WARNING: printMat2(",i2,") not used")') k; enddo 
          if (printMat2(10)) then; write(Imatrix,'("matrix subset printed; i1,i2, j1,j2:":4i3)') (printRng2(i),i=1,4);
                                   write(IO,'("matrix subset printed; i1,i2, j1,j2:":4i3)') (printRng2(i),i=1,4); endif
        endif  
!        
        test=.false.; do k=1,10; if (printMat3(k)) test=.true.; enddo
        if (test) then
          Write(IO,'(" Matrix Options:      PRT3(1-",I2,") = ",4(5L2,2X))') 10,(printMat3(I),I=1,10)
          write(IO,'(5X,"The following matrices in |SLJMJ> basis will be printed to the file ",A30)') fnamem
          write(Imatrix,'("This file contains the following matrices from PRT3:")') 
          if (printMat3(1)) then; write(Imatrix,'("Full matrix.")'); write(IO,'(5X,"Full matrix.")'); endif
          if (printMat3(2)) then; write(Imatrix,'("*/. matrix.")');  write(IO,'(5X,"*/. matrix.")'); endif
          if (printMat3(3)) then; write(Imatrix,'("Ukq matrices in |SLJMJ> basis (non-zero MEs only, max 100).")'); 
                                  write(IO,  '(5X,"Ukq matrices in |SLJMJ> basis (non-zero MEs only, max 100).")'); endif
          if (printMat3(4)) then; write(Imatrix,'("Ukq matrices in |SLJMJ> basis (full matrix, max 14x14).")'); 
                                  write(IO,  '(5X,"Ukq matrices in |SLJMJ> basis (full matrix, max 14x14).")'); endif
          if (printMat3(5)) then; write(Imatrix,'("<|ED|> MEs in |SLJMJ> basis.")'); 
                                  write(IO,  '(5X,"<|ED|> MEs in |SLJMJ> basis.")'); endif
          if (printMat3(6)) then; write(Imatrix,'("L,S matrices in |SLJMJ> basis.")'); 
                                  write(IO,  '(5X,"L,S matrices in |SLJMJ> basis.")'); endif
          if (printMat3(7)) then; write(Imatrix,'("<|mu|> matrices in |SLJMJ> basis.")'); 
                                  write(IO,  '(5X,"<|mu|> matrices in |SLJMJ> basis.")'); endif
          if (printMat3(8)) then; write(Imatrix,'("LF matrix in original |SLJMJ> and PRED basis.")'); 
                                  write(IO,  '(5X,"LF matrix in original |SLJMJ> and PRED basis.")'); endif
          if (printMat3(9)) write(IO,'("****WARNING: printMat3(",i2,") not used")') k
          if (printMat3(10)) then; write(Imatrix,'("matrix subset printed; i1,i2, j1,j2:":4i3)') (printRng3(i),i=1,4);
                                   write(IO,'("matrix subset printed; i1,i2, j1,j2:":4i3)') (printRng3(i),i=1,4); endif
        endif          

        test=.false.; do k=1,10; if (printMat4(k)) test=.true.; enddo
        if (test) then
          Write(IO,'(" Matrix Options:      PRT4(1-",I2,") = ",4(5L2,2X))') 10,(printMat4(I),I=1,10)
          write(IO,'(5X,"The following matrices in |SLJ> basis will be printed to the file ",A30)') fnamem
          write(Imatrix,'("This file contains the following matrices from PRT4:")') 
          if (printMat4(1)) then; write(Imatrix,'("Info on ESOs in a |MJ> basis.")') 
                                  write(IO,  '(5X,"Info on ESOs in a |MJ> basis.")'); endif
          if (printMat4(2)) then; write(Imatrix,'("Half ESO matrices in |MJ> basis.")') 
                                  write(IO,  '(5X,"Half ESO matrices in |MJ> basis.")'); endif
          if (printMat4(3)) then; write(Imatrix,'("Matrix in uncoupled |SLMlMs> basis.")') 
                                  write(IO,  '(5X,"Matrix in uncoupled |SLMlMs> basis.")'); endif
          if (printMat4(4)) then; write(Imatrix,'("Matrix in uncoupled |SLMlMs=S> basis.")') 
                                  write(IO,  '(5X,"Matrix in uncoupled |SLMlMs=S> basis.")'); endif
          if (printMat4(5)) then; write(Imatrix,'("Matrix in uncoupled |Sfi(i=1,7)> basis.")') 
                                  write(IO,  '(5X,"Matrix in uncoupled |Sfi(i=1,7)> basis.")'); endif
          if (printMat4(6)) then; write(Imatrix,'("Lx,Ly,Lz matrices in uncoupled |SLMlMs=S> basis.")') 
                                  write(IO,  '(5X,"Lx,Ly,Lz matrices in uncoupled |SLMlMs=S> basis.")'); endif
          do k=7,9; if (printMat4(k)) write(IO,'("****WARNING: printMat4(",i2,") is not used")') k; enddo 
          if (printMat4(10)) then; write(Imatrix,'("matrix subset printed; i1,i2, j1,j2:":4i3)') (printRng4(i),i=1,4);
                                   write(IO,'("matrix subset printed; i1,i2, j1,j2:":4i3)') (printRng4(i),i=1,4); endif
        endif  
!
        test=.false.; idebug=0; do k=1,10; if (check1(k).or.check2(k)) test=.true.; enddo
        IF (nExplicitbasis.gt.0) test=.true.  ! writes explicit basis into debug file.
        if (test) then
          Write(IO,'(" Check Options:       CHK1(1-",I2,") = ",4(5L2,2X))') 10,(check1(I),I=1,10)
          Write(IO,'(" Check Options:       CHK2(1-",I2,") = ",4(5L2,2X))') 10,(check2(I),I=1,10)
          idebug=4; fnamed="debug"//TRIM(Fext)//".dat"
          OPEN(Idebug,FILE=fnamed,STATUS='UNKNOWN')
          write(Idebug,'("Debug file made from calculation made using the input file:",A)') trim(INFILE)
          write(Idebug,'("CHK1 = ",4(5L2,2X))') (check1(I),I=1,10)
          write(Idebug,'("CHK2 = ",4(5L2,2X))') (check2(I),I=1,10)
          if (check1(1)) then; write(Idebug,'("Basis matrices.")'); write(IO,'(5X,"Basis matrices.")'); call writeLSBasis(); endif 
          if (check1(2)) then; write(Idebug,'("UMat  matrices.")'); write(IO,'(5X,"UMat  matrices.")'); endif 
          if (check1(3)) then; write(Idebug,'("V11   matrices.")'); write(IO,'(5X,"V11   matrices.")'); endif 
          if (check1(4)) then; write(Idebug,'("EE    matrices.")'); write(IO,'(5X,"EE    matrices.")'); endif 
          if (check1(5)) then; write(Idebug,'("abg   matrices.")'); write(IO,'(5X,"abg   matrices.")'); endif 
          if (check1(6)) then; write(Idebug,'("MnPnSn matrices.")');write(IO,'(5X,"MnPnSn matrices.")');endif 
          if (check1(7)) then; write(Idebug,'("Tn    matrices.")'); write(IO,'(5X,"Tn    matrices.")'); endif 
          if (check1(8)) then; write(Idebug,'("AOM   matrices.")'); write(IO,'(5X,"AOM   matrices.")'); endif 
          if (check1(9)) then; write(Idebug,'("Debugging the diagonalisation.")');
                            write(IO,'(5X,"Debugging the diagonalisation.")'); endif 
          if (check1(10))then; write(Idebug,'("Debugging the group theory module.")');
                            write(IO,'(5X,"Debugging the group theory module.")'); endif 
                            
          if (check2(1)) write(IO,'("****WARNING: Check2(1) not used")')               
          if (check2(2))then; write(Idebug,'("Check the WF are Kramers doublets.")'); 
                            write(IO,'(5X,"Check the WF are Kramers doublets.")'); endif
          if (check2(3))then; write(Idebug,'("Check the L,S matrices give L.S SOC matrix")'); 
                            write(IO,'(5X,"Check the L,S matrices give L.S SOC matrix")'); endif
          if (check2(4))then; write(Idebug,'("Check the conversion between Altp and Blki parameters")'); 
                            write(IO,'(5X,"Check the conversion between Altp and Blki parameters")'); endif
          if (check2(5))then; write(Idebug,'("Check the CCF matrix elements")'); 
                            write(IO,'(5X,"Check the CCF matrix elements")'); endif
          do k=6,10; if (check2(k)) write(IO,'("****WARNING: Check2(",i2,") not used")') k; enddo                  
          write(Idebug,'(105("-"))')
        endif
!
        test=.false.; do k=1,10; if (warnings(k)) test=.true.; enddo
        if (test) then
          Write(IO,'(" Warning Options:     WARN(1-",I2,") = ",4(5L2,2X))') 10,(warnings(I),I=1,10)
          if (warnings(1)) write(IO,'(5X,"Check of Bkq values will be made.")') 
          if (warnings(2)) write(IO,'(5X,"Check of Altp values will be made.")') 
        endif

        test=.false.; do k=1,8; if (fitOpt(k)) test=.true.; enddo
        if (test) Write(IO,'(" Fit Options:         FITO(1-",I2,") = ",4(5L2,2X))') 8,(FitOpt(I),I=1,8)
        if (FitOpt(1)) write(IO,'(5X,"RESHAPing fit matrices printed")') 
        if (FitOpt(2)) write(IO,'(5X,"Covariance matrix will be printed")') 
        if (FitOpt(3)) write(IO,'(5X,"Minimized function and parameters printed on each iteration.")') 
        if (FitOpt(4)) write(IO,'(5X,"Fitted parameters are written to the file ""fit.out""")') 
        if (FitOpt(5)) write(IO,'(5X,"Fit matrix written to the file ""matrices.dat"" each iteration.")') 
        if (FitOpt(6)) write(IO,'(5X,"The parameters limits will be hard.")') 
        if (FitOpt(7)) write(IO,'(5X,"The fit will sum over all the degeneracies.")') 
        if (FitOpt(8)) write(IO,'(5X,"The fit will also give the mean deviation: Sum(|Eexp-Ecal|)/Nexp.")') 
        if (test .and. nFit.eq.0) write(IO,'(5X,"****Warning: A fit must be made for option in ""FITO"" to come into effect.")') 
        
        if (fitType.ne.0) then
          write(io,'(/," Experimental data provided:")')
          if (fitType.eq.1) write(IO,'(5X,"A comparison of calculated to given energies will be made.")')
          if (assignMethod.eq.1) write(IO,'(5X,"Crystal Quantum numbers will be used for energy assignments.")') 
          if (assignMethod.eq.2) write(IO,'(5X,"Spin will be used for energy assignments.")') 
          if (fitType.eq.2) write(IO,'(5X,"A comparison of multiplet calculation to experimental multiplets will be made.")')
          if (fitType.eq.3) write(IO,'(5X,"A comparison of calculated to given Bkq values will be made.")')
          if (fitType.eq.4) write(IO,'(5X,"A comparison of calculated and given Extend Stevens Operators will be made.")')
          if (fitType.eq.5) write(IO,'(5X,"A comparison of calculated and given g-values will be made.")')
          if (fitType.eq.6) write(IO,'(5X,"A comparison of calculated and given 1-electron LF matrix will be made.")')
        endif
!

        write(IO,'(115("-"))')
!
!        CALL getLSBasis()   ! loads d or f LS basis
        if (IDG.eq.-1) call GETGRP()   ! Sets IDG non-zero
        CALL getFullBasis(0,nExBasisNums) ! gets full SLJMj basis (nExplicitbasis,nExBasisNums ignored in this call)
        nline=1
        Call readUMat(nline)
        Call readVMat(nline)   
        CALL readMnPnSnMat(nline)
        CALL readTnMat(nline)
        call readEEMat(nline)
        CALL readABGMat()
        if (CCFtype.ne.0) Call readCCFMat()   
!
        test=.false.; do k=1,8; if (UseFncrossErrors(k)) test=.true.; enddo
        if (test) then
          WRITE(IO,'(80("*"),/,"***Warning; deliberate errors have been introduced to reproduce erroneous", &
            " published results.",/,80("*"),/,  &
            "                      FNCR(1-",I2,") = ",4(5L2,2X))') 8,(UseFncrossErrors(k),k=1,8)
          if (UseFncrossErrors(8)) then; write(io,'("Errors are taken from the appropriate ''fnmp file'' (Crosswhite programs)")') 
            ELSE; write(io,'("Errors are taken from the appropriate ''fncross'' file (Reid programs)")'); ENDIF    
          call FncrossErrors(UseFncrossErrors)
        endif
        if (Complementary) call t2MatComp(UseFncrossErrors(7)) ! Must call EEmat before t2MatComp, as we need the e3 values.  
        if (check1(7)) call testTnMat()      ! Must call EEmat before testTnMat, as we need the e3 values.  
!        
!*********************         
!  All parameters and constants should be set before here. Can now load into array P.
   !     call loadP()
  ! No; only constants set
!*********************  
!        
!  Handle the variables:
        if (nVar.gt.0) then
          minNv=Nvarvals(1); gridNv=1
          do i=1,nVar
            minNv=min(Nvarvals(i),minNv); gridNv=gridNv*Nvarvals(i)
            variablesvalues(i)=allVarVals(i,1); P(400+i)=allVarVals(i,1)  ! initial values of variables
          enddo
          if (nVar.gt.1) then
            if (grid.and.gridNv.gt.1) then
              write(io,'(" A series of calculations will be made on a grid of variable values. Total calcs=",I3)') gridNv
              if (gridNv.gt.1000) then
                write(io,'("***FATAL: Too many calculations, gridNv=",I5)') gridNv
                stop
              endif
            else  
              if (minNv.gt.1) write(io,'(" A series of calculations will be made, stepping the variable values together.", &
                                         " Total calcs=",I3)') minNv
            endif
          else
            if (grid) write(io,'("***WARNING: GRID command has no effect when only one variable exists.")') 
          endif          
        endif
        test=.false.; do k=1,10; if (PLOTout1(k)) test=.true.; enddo
        if (test) then
          Write(IO,'(" Plot Options:        PLOT(1-",I2,") = ",4(5L2,2X))') 10,(PlotOut1(I),I=1,10)
          write(IP1,'(4I4,"   ! Nplots, Nvar, LValue, ID")') minNv,Nvar,LValue,0  ! ID=0
          if (nvar.gt.0) then
            do i=1,nvar; WRITE(IP1,'(A10,I4,/,(5E16.8))') variables(i),minNv,(allVarVals(i,n),n=1,minNv); enddo
          endif
        endif
!
!  Handle the constants:
        if (Ncon.gt.0) then
          do i=1,nCon
            variables(nVar+i)=Constants(i)
            variablesvalues(nVar+i)=ConstantV(i)
          enddo
!          write(io,'("      variables(",i2,")=",10(2X,A10))') nVar+nCon,(variables(i),i=1,nVar+nCon)
!          write(io,'("variablesvalues(",i2,")=",10F12.4)') nVar+nCon,(variablesvalues(i),i=1,nVar+nCon)
        endif
!  Run through once without variables set, determine which parameters are numbers (or can be evaluated immediately, ie: sqrt(2)),
!  and which are dependent on the variable values.         
        atLeastOneVarForm=.false.
        do i=1,maxTotPar  ! Don't go to maxTotPar, or will reset variables P(401)-P(450) to zero
          if (i.gt.400 .and. i.le.450) goto 160
          func=pEq(i)
!          do j=1,nCon
!            if (index(func,trim(Constants(j))).ne.0) write(io,'("Found ",A10," in ",A100)') Constants(j),func
!          enddo 
          call init(func, Constants, statusflag)
!          call init(func, noNames, statusflag)
          if(statusflag == "ok") Then
            P(i) = evaluate(ConstantV)
!            P(i) = evaluate(noValues)
            pEqF(i) = .false.
          else
            if (i.gt.400 .and. i.le.450) then 
              write(io,'("The variable ",A9," cannot be specified by a formula:",A)') pLab(i),trim(pEq(i))
              stop
            endif
            atLeastOneVarForm=.true.
            pEqF(i) = .true.  ! requires variables
            if (V1(i).eq."*") then
              write(io,'("***FATAL: Cannot have a parameter:",A9," defined by an equation and varied in a fit.")') Plab(i)
              stop
            endif 
            V1(i)="F"   ! parameter determined by formula involving variables, (otherwise can be evaluated immediately).
!            write(IO,'("func: ",A)') trim(func)
!            write(IO,'("Parameter:",i3," ",A9," will be determined by variables: ",A)') i,pLab(i),trim(pEq(i))
          endif  
!          if (P(I).ne.0.0d+00) write(IO,'("Parameter:",i3," ",L1," ",A9," P(i)=",F12.4)') i,pEqF(i),pLab(i),P(i)
          call destroyfunc()
 160      continue
        enddo  ! i
        call parLinearity()   !  Checks parameters/variables being varied in a fit for linearity (sets V2). 
! Set the variables and make sure that they are all there        
        if (atLeastOneVarForm) then
          do i=1,maxTotPar
            if (pEqF(i)) then
              func=pEq(i)
              call init(func, variables, statusflag)
              if(statusflag == "ok") Then
                P(i) = evaluate(variablesvalues)
              else
                write(IO,'("***FATAL: Parameter:",i3," ",A9," contains unknown variables: ",A)') i,pLab(i),trim(func)
                stop              
              endif  
              call destroyfunc()
            endif            
          enddo
        endif  
!        write(io,'("Contents of P:")'); write(io,'((10F12.3))') (P(i),i=1,330)
!        write(io,'(" # 2:",3F12.3)') P(310),P(311),P(312)
        call unloadP()  ! Empties array P back into the parameters: EAVE,F2,F4,..
!        write(io,'(" # 3:",3F12.3)') P(310),P(311),P(312)
!-------------
! Parameter values now set to initial values.

        magFld=.false.; if (abs(MAGF(1))+abs(MAGF(2))+abs(MAGF(3)).gt.bit) magFld=.true.
        if (FullMJBasis .and. (GEXS(1).ne.0 .or. magFld .or. IMD.ne.0)) call calcMagMom()
        if (magFld) WRITE(IO,'(" Applied magnetic field: (Hx,Hy,Hz)=",3F8.4," Tesla.")') MAGF(1),MAGF(2),MAGF(3)
!
        if (RK(1).ne.D1 .or. RK(2).ne.D1 .or. RK(3).ne.D1) then
          WRITE(IO,'(" Reduction factors (Kx,Ky,Kz)=",3F8.4," used for g-values and mag.dip. transitions.")') RK(1),RK(2),RK(3)
        endif
        
        if (fitType.eq.4) then
          if (mod(ESO2J,2).eq.0) then    ! 2J even; 2J+1 odd
            nr=(ESO2J/2+1)**2;           ni=(ESO2J/2+1)*ESO2J/2   ! real and imag MEs
          else  ! 2J odd; 2J+1 even
            nr=(ESO2J+1)/2*(ESO2J+3)/2;  ni=((ESO2J+1)/2)**2  ! real and imag MEs
          endif 
          nexp = nr + ni
          write(io,'("Fitting ESOs:",i3," MEs will be fitted: ",i2," real and ",i2," imaginary")') nexp,nr,ni 
        endif
        
        call checkValidERanges(Fn,nExplicitbasis,nExBasisNums)  !  Check limit & parameter consistency; can make changes to parameter values.
                                                            !  Must be called after getLSBasis(), to have the TS_full values set.  
        if (ED_IntMech.eq."Full") call checkAltp() 

        if (UseDefaultIRREPS) then
          if (zeta.eq.0.0d+00) then; ISO=0; else; ISO=1; endif
        endif
! getFullBasis() gets the SLJMj basis that have the |SLJ> basis numbers explicitly given in nExBasisNums
        if (nExplicitbasis.ne.0) CALL getFullBasis(nExplicitbasis,nExBasisNums) 
!        
!  Handle the links:     
        if (nlinks.gt.0) then
          WRITE(IO,'(/,"Linked parameters:",/,115("-"))'); 
          do i=1,nlinks
            if (V1(nlink(1,i)).eq."*") then
              write(io,'("***FATAL: Cannot have a parameter:",A9," defined by an LINK and varied in a fit.")') Plab(nlink(1,i))
              stop 
            endif
            V1(nlink(1,i))="L"
            write(IO,'(1X,A9,"=",F12.6," * ",A9,"+",F12.6)')  PLab(nlink(1,i)),Rlink(1,i),PLab(nlink(2,i)),Rlink(2,i)
          enddo
        endif        
        if (nlinks.gt.0 .or. nfit.gt.0) call checkLinkFitP()
        if (nlinks.gt.0) then
          call makeLinks()  ! apply links before writing parameters.
          call unloadP()  ! Empties array P back into the parameters: EAVE,F2,F4,..
        endif
!
!  Handle ligand field:
        if (LFtype.eq."AOM ") then   ! AOM parameters used, converted to Bkq values, but these are not now parameters.
          call convertToAOM(DistMat,AngMat,DistLab,Distbonds)
          if (OUTP(1)) then
            WRITE(IO,'(/,"Ligand distance matrix")'); WRITE(IO,'(115("-"))')
            WRITE(IO,'(7X,20(5X,A4))') (Lname(i),i=1,NLIGANDS)
            WRITE(IO,'(" Bonds  |",20(F8.4,1X))') (Distbonds(i),i=1,NLIGANDS)
            do i=1,NLIGANDS
              write(IO,'(4X,A4,"|",20(F8.4,A1))') Lname(i),(DistMat(i,j),DistLab(i,j),j=1,i)
            enddo
            WRITE(IO,'(/,"Ligand angle matrix")'); WRITE(IO,'(115("-"))')
            WRITE(IO,'(7X,20(5X,A4))') (Lname(i),i=1,NLIGANDS)
            do i=1,NLIGANDS
              write(IO,'(4X,A4,"|",20F9.4)') Lname(i),(AngMat(i,j),j=1,i)
            enddo
          endif  
        endif  
         if (LFtype.eq."AOM ".or. LFtype.eq."LF1E") then 
           if (Lvalue.eq.2) call AOMmatrixD()
           if (Lvalue.eq.3) call AOMmatrixF()
        endif  
        if (abs(RotLF(1))+abs(RotLF(2))+abs(RotLF(3)).gt.bit) then
          WRITE(IO,'(" Crystal Field will be rotated by the Euler angles:",3F8.2)') (RotLF(i),i=1,3)
          call RotateLF()
        endif
        call checkLFcomplex() ! tests BKQI(i), & set to 0 if < CFbit   ! Sets CFcomplex
        if (fastMat(2)) call checkTimeRev()   ! tests to see if there is Time reversal blocking.
        call BlockBasis(0)    ! needs the Bkq values first; Call even if no symmetry blocking formally requested.
        call checkAssignments()
!
!********* OUTPUT *********
        call PrintParameters(0)  ! Output Constants, Variables, Variable expressions
        call PrintParameters(1) 
!********* OUTPUT *********
        MFcomplex=.false.; if (abs(MAGF(2)).gt.bit) MFcomplex=.true.
        MatComplex=CFcomplex.or.MFcomplex
        call checkDiagOpts()
!          spin proj?, free ion % short?, free ion % long?, MJ % short?, MJ % long?
        test=Option(2).or.Option(3).or.Option(4).or.Option(6).or.Option(7)
        if (test .or. IDG.ne.0 .or. gexs(1).ne.0 .or. IMD.ne.0 .or. IED.ne.0) calcVecs=.true. 
        if (N_expect.ne.0 .or. vectLow.ne.0) calcVecs=.true. 
        close(idata)  ! d_electron.dat or f_electron.dat
        call flush(IO)  ! Allows output to be seen immediately.

        return
!---------------------------------------------------------------------------------------------------------- 

      CASE default
        Write(io,'("***FATAL: Unrecognized KEY=",A4," in line:",I4,5X,A80)') KEY,ninput,COMMAND; STOP
      END SELECT 
      goto 10
! 
 990  write(io,'("***FATAL: Error reading a parameter:",A4,", remember you need at least one value")') KEY; stop
 982  write(io,'("***FATAL: For PRT*(10) true, a line with i1, i2, j1, j2 must immediately follow",I3)'); stop
 980  write(io,'("***FATAL: Error reading an OFFS parameter:",I3)') i; stop
 970  write(io,'("***FATAL: Error reading the REDF command (3 values required)")'); stop
 960  write(io,'("***FATAL: Error reading a FIT parameter")'); stop
 952  write(io,'("***FATAL: Error reading the ATOM command")'); stop
 951  write(io,'("***FATAL: Error reading the LANC command")'); stop
 950  write(io,'("***FATAL: Error reading the LINK command after NLinks=",I2)') NLinks; stop
 945  write(io,'("***FATAL: Error reading the Alpt intensity parameters, they must be on 9 separate lines")'); stop
 940  write(io,'("***FATAL: Error reading the JUDO command, (expect N1,N2,N3,N4,ID,Om2,Om4,Om6:",A80)') COMMAND; stop
 935  write(io,'("***FATAL: Error in command: FIT nFit,maxIts,FitTol; more input expected:",A80)') COMMAND; stop
 930  write(io,'("***FATAL: Error in command online:",I4,/,A4,A80)') ninput,KEY,COMMAND; stop
 925  WRITE(IO,'("***FATAL: Error reading the ",I4," experimental CF values")') NEXP;  stop
 922  WRITE(IO,'("***FATAL: Remember 1 parameters: EXPG  Nexp")'); stop
 921  WRITE(IO,'("***FATAL: Remember 2 parameters: EXPE  Nexp  EngInt")'); stop
 920  WRITE(IO,'("***FATAL: Error reading the ",I4," experimental energies")') NEXP;  stop
 918  write(io,'("***FATAL: f-electron calculation and could not find file:""f_electron.dat"" in directory:",A20)') F_E_PATH; stop
 917  write(io,'("***FATAL: d-electron calculation and could not find file:""d_electron.dat"" in directory:",A20)') F_E_PATH; stop
 916  write(io,'("***FATAL: Commmand CCF and could not find file:""ccf.dat"" in directory:",A20)') F_E_PATH; stop
 915  WRITE(IO,'("***FATAL: Error reading the correlated crystal field Bikq parameters on line:",I2,/,   &
                 " CCF  command has CCFtype="I3,"; N_CCF=",i3)') CCFtype,N_CCF;  stop
 914  write(io,'("***FATAL: Commmand SYML and could not find file:""group.dat"" in directory:",A20)') F_E_PATH; stop
 912  WRITE(IO,'("***FATAL: Error reading the euler angles for ROTL; remember you need 3 values")');  stop
 911  WRITE(IO,'("***FATAL: Error reading the magnetic field for MFLD; remember you need 3 values")');  stop
 910  WRITE(IO,'("***FATAL: Error reading the crystal field Bkq parameters on line:",I2,/,   &
                 " CF  command must be followed by 16 lines of Bkq")') ninput;  stop
 906  WRITE(IO,'("***FATAL: Error reading the LF1E data on line:",I2,/)') ninput;  STOP
 905  WRITE(IO,'("***FATAL: Error reading the AOM data for a ligand on line:",I2,/,   &
                 " The lines after AOM must be in the format: ""Name,Esig,EpiX,EpiY,theta,phi,chi"" if NLIG>0",/, &
                 " (Optionally ""chi"" can be left out, and it will be set to zero.)",/, &
                 " or ""Name, Esig, EpiX, EpiY, x, y, z"" if NLIG<0")') ninput;  STOP
 904  WRITE(IO,'("***FATAL: Error reading the OFFS command on line:",I2,/,   &
                 " Type ""> f_electron help OFFS"" to see expected format")') ninput; STOP
 902  write(io,'("***FATAL: Error reading input; remember you need a CONF command and an END command.")');  stop
 900  write(io,'("***FATAL: Error reading input.")');  stop
!
      end subroutine input
!
!-----------------------------------------------------------------------

      SUBROUTINE AssignCalcJ(assNum,J,Jhigh,KK)
!
!  Assigns the calculated Energy level (J) to an experimental one.
!  assignMethod=0 Normal energy ordering 
!  assignMethod=1 CQN 
!  assignMethod=2 Spin 
!
      IMPLICIT NONE
      integer j,jhigh,jj,kk,assNum(max_exp), k,k1
!      
      kk=0
      do k=1,nexp
        if (assignMethod.eq.0) then 
          do k1=1,Nass(k)
            if (NAssign(k,k1).ge.j .and. NAssign(k,k1).le.jhigh) kk = k
          enddo ! k1
        elseif (assignMethod.eq.1 .or. assignMethod.eq.2) then
          do jj=j,jhigh
            if (assNum(k).eq.jj) kk = k
          enddo ! jj
        endif ! assignMethod  
      enddo ! k
       
!      
      end subroutine AssignCalcJ
!
!-----------------------------------------------------------------------
!
      SUBROUTINE AssignNum(assNum,neng)
!
!  For assignMethod=1 (CQN) or 2 (Spin), returns in assNum(k) the k assignment values 
!
      IMPLICIT NONE
      integer i,j,j1,jj,k,assNum(max_exp),neng
      REAL*8 D1/1.0D00/, D2/2.0D00/, bit/1.0d-08/ 
      if (assignMethod.eq.1) then
        do i=1,Nexp
          k=0
          do j=1,neng
            if (iCQN(j).eq.NAssign(i,2)) then
              k=k+1;jj=j
            endif  
            if (k.eq.NAssign(i,1)) goto 105
          enddo
          write(io,'("***FATAL: AssignNum: level",i3," of CQN block",i2," not found.")') NAssign(i,1),NAssign(i,2)
          stop
 105      assNum(i)=jj
!          write(io,'("i,NAssign(i,1),NAssign(i,2),assNum(i)=",4i4)') i,NAssign(i,1),NAssign(i,2),assNum(i)
        enddo ! i
      elseif (assignMethod.eq.2) then
        do i=1,Nexp
          k=0
          do j=1,neng
!            write(io,'("i,j,2*SPIN(j,1)+1,2*SPIN(j,2)+1=",4i4)') i,j,INT(SPIN(j,1)*S_mult(1)+bit),INT(SPIN(j,2)*S_mult(2)+bit)
            do j1=1,NSPIN  ! SPIN(J,K),K=1,NSPIN
              if (INT(SPIN(j,j1)*S_mult(j1)+bit).eq.NAssign(i,2)) then
                k=k+1; jj=j
              endif  
              if (k.eq.NAssign(i,1)) goto 106
            enddo  ! j1
          enddo  ! j
          write(io,'("***FATAL: AssignNum: level",i3," of Spin 2*S+1=",i2," not found.")') NAssign(i,1),NAssign(i,2)
          stop
 106      assNum(i)=jj          
!          write(io,'("i,NAssign(i,1),NAssign(i,2),assNum(i)=",4i4)') i,NAssign(i,1),NAssign(i,2),assNum(i)
        enddo ! i
      endif  ! assignMethod
!      
      end subroutine AssignNum
!
!-----------------------------------------------------------------------
!
      SUBROUTINE PrintParameters(icall)
! icall=0 Ouput that need not be repeated every calcaulation
!       1 Parameters as read from input.
!       2 Parameters at the end of a fit.
! 
      IMPLICIT NONE
      integer i,j,ii,k,l,q,icall,j1,j2
      real*8 LF_str,D0,D1,D2,bit,xx,sum1
      logical someSmall
      parameter(D0=0.0d+00,D1=1.0d+00,D2=2.0d+00, bit=1.0d-12) 
      character lab*8,cL(maxTotPar)*8,lab10*10
             
      if (icall.eq.0) Then
       if (nCon.gt.0) then
          WRITE(IO,'(/,"Constants:  Value        Formula",/,115("-"))'); 
          do i=1,nCon    !  variablesvalues(nVar+i) instead of ConstantV(i), as ConstantV(i) wiped
            write(io,'(A10,F12.6,3X,A75)') Constants(i),variablesvalues(nVar+i),cEq(i)
          enddo
        endif 
        if (nvar.gt.0) then
          WRITE(IO,'(/,"Variables:",/,115("-"))');
          do i=1,nvar 
            if (Nvarvals(i).eq.1) write(IO,'(I4,1X,A10,2A1,F12.4)') 400+i,trim(variables(i)),v1(400+i),v2(400+i),P(400+i)
            if (Nvarvals(i).gt.1) write(IO,'(I4,1X,A10,8F12.4,/,(15X,8F12.4))') 400+i,trim(variables(i)),   &
                                                                                (allVarVals(i,j),j=1,Nvarvals(i))
          enddo 
          WRITE(IO,'(115("-"))')
        endif
        if (atLeastOneVarForm) then
          WRITE(IO,'("Parameters defined as a function of variables:")'); 
          write(IO,'("Parameter    Init.Val.    Function    ")')
          do i=1,maxTotPar
            if (pEqF(i)) write(IO,'(2(i3,2X,A9,F10.4,3X,A40))') i,pLab(i),P(i),ADJUSTL((pEq(i)))
          enddo
          WRITE(IO,'(115("-"))'); 
        endif  
      endif
!
      if (icall.eq.1) Then
        iF (ESO2J.ne.0) then
          write(io,'(" A calculation using Extended Stevens Operators")') 
          if (ESO2J.eq.2*INT(ESO2J/2)) then  ! even
            if (ESONorm.eq.0) then 
              write(io,'(" Basis: J=",I2,";  k(max)=",i2,";   normalisation using Kkq constants")') ESO2J/2,ESOkmax
            else   
              write(io,'(" Basis: J=",I2,";  k(max)=",i2,";   No normalisation used")') ESO2J/2,ESOkmax
            endif  
          else
            if (ESONorm.eq.0) then 
              write(io,'(" Basis: J=",I2,"/2;  k(max)=",i2,";   normalisation using Kkq constants")') ESO2J,ESOkmax
            else 
              write(io,'(" Basis: J=",I2,"/2;  k(max)=",i2,";   no normalisation used")') ESO2J,ESOkmax 
            endif  
          endif
          write(io,'("  k    q             Bkq")') 
          do k=1,ESOkmax/2
            write(io,'(i3," [",i3,":",i2,"]",6E16.8,/,(12X,6E16.8))') 2*k,-2*k,2*k,(ESOBkq(k,q),q=1,4*k+1)
          enddo ! k  
          CALL calc_ESO() ! Do Extended Stevens operators in |J> basis.
          if (ESOskip.eq.1) return
        endif
!      
        if (Nfit.gt.0) then
          WRITE(IO,'("Fitting Options: nFit:",I2,",   maxIts:",I4,",   FitTol:",E12.6,/,115("-"))') Nfit,maxIts,FitTol 
          if (.not.linearFit) then 
            write(io,'("The fitting is a nonlinear function of parameters")') 
            write(io,'("This will require a new matrix to be calculated every iteration")') 
          endif
          WRITE(IO,'(/,"Initial values for parameters: (*l/*n/L/F indicates Fit linear/Fit nonlinear/Linked/Formula)",/,115("-"))') 
        else
          WRITE(IO,'(/,"Atomic Parameters:  (L:LINKed; F:formula)",/,115("-"))'); 
        endif
        WRITE(IO,22) "           Offset:", 1," EAVE",EAVE, V1(1),V2(1)
        if (racah) then
        WRITE(IO,21) "E-E    parameters:", 2,"    B",Bracah,V1(2),V2(2), 3,"    C",Cracah,V1(3),V2(3) 
        else
        WRITE(IO,20) "E-E    parameters:", 2,"   F2",  F2, V1(2),V2(2),  3,"   F4",  F4, V1(3),V2(3),  4,"   F6",  F6,V1(4),V2(4)
        endif
        WRITE(IO,22) "SOC    parameter: ", 5," ZETA",ZETA, V1(5),V2(5)
        WRITE(IO,20) "E-E CI parameters:", 6,"ALPHA",ALPHA,V1(6),V2(6),  7," BETA",BETA, V1(7),V2(7),  8,"GAMMA",GAMMA,V1(8),V2(8)
        WRITE(IO,20) "soo    parameters:", 9,"   M0",  M0, V1(9),V2(9), 10,"   M2",  M2,V1(10),V2(10),11,"   M4", M4,V1(11),V2(11)
        WRITE(IO,20) "ec-so  parameters:",12,"   P2",  P2,V1(12),V2(12),13,"   P4",  P4,V1(13),V2(13),14,"   P6", P6,V1(14),V2(14)
        WRITE(IO,20) "3 body parameters:",15,"   T2",  T2,V1(15),V2(15),16,"   T3",  T3,V1(16),V2(16),17,"   T4", T4,V1(17),V2(17)
        WRITE(IO,20) "                  ",18,"   T6",  T6,V1(18),V2(18),19,"   T7",  T7,V1(19),V2(19),20,"   T8", T8,V1(20),V2(20)
 20     format(A18,i3,1x,A5,"=",f12.4,2A1,2(i3,1x,A5,"=",f12.4,2A1),"|")     
 21     format(A18,i3,1x,A5,"=",f12.4,2A1,i3,1x,A5,"=",f12.4,2A1,24X,"|")     
 22     format(A18,i3,1x,A5,"=",f12.4,2A1,48x,"|")     
        WRITE(IO,'(" Ligand Field Parameters:",/,115("-"))') 
        if (LFtype.eq."AOM ") then 
          if (AOM_angles) WRITE(IO,'(" AOM Ligand field of ",I2," ligands. Angles given, X,Y,Z generated.")') NLIGANDS
          if (.not.AOM_angles) WRITE(IO,'(" AOM Ligand field of ",I2," ligands. X,Y,Z given, angles generated.")') NLIGANDS
          WRITE(IO,'(115("-"),/,"   Name   |   E(sigma)        E(piX)        E(piY)   |",            &
                                "     Theta        Phi        Chi       |       X           Y           Z")')
          do i=1,NLIGANDS
            write(IO,'(4X,A4,2X,"|",3(1X,I2,F9.2,2A1),"|",3(1X,I3,F7.2,2A1),"|",2X,3(1X,F11.8))')    &
                  Lname(i),20+i,eSigma(i),V1(20+i),V2(20+i),40+i,ePiX(i),V1(40+i),V2(40+i), 60+i,ePiY(i),V1(60+i),V2(60+i),    &
                           80+i,theta(i),V1(80+i),V2(80+i),100+i,phi(i),V1(100+i),V2(100+i),120+i,chi(i),V1(120+i),V2(120+i),  &
                          xLig(i),yLig(i),zLig(i)
          enddo
          WRITE(IO,'(115("-"))')
          if (AOMX) then 
            WRITE(IO,'(" Extended AOM Ligand field")')
            WRITE(IO,'(115("-"),/,"   Name   |   E(DeltaS)     E(DeltaC)        E(PsiS)      E(PsiC)  |")')
            do i=1,NLIGANDS
              write(IO,'(4X,A4,2X,"|",4(1X,I2,F9.2,2A1),2X,"|")') Lname(i),  &
                  140+i,eDelS(i),V1(140+i),V2(140+i),150+i,eDelC(i),V1(150+i),V2(150+i), 160+i,ePsiS(i),V1(160+i),V2(160+i),    &
                  170+i,ePsiC(i),V1(170+i),V2(170+i)
            enddo
            WRITE(IO,'(115("-"))')
          endif
          if (OUTP(2)) WRITE(IO,'(" Equivalent crystal field:")')
        endif  ! LFtype.eq."AOM "
        
        if (LFtype.eq."LF1E") then 
          if (Lvalue.eq.2) then
            WRITE(IO,'(115("-"))')
            WRITE(IO,'(14x,"| z2  >",10x,"| yz  >",10x,"| xz  >",10x,"| xy  >",8x,"|x2-y2> ")')
            WRITE(IO,'(" < z2 | ",5(I2,":" F12.4,2A1))') (20+i,LF1emat(1,i),V1(20+i), V2(20+i),i=1,1) 
            WRITE(IO,'(" < yz | ",5(I2,":" F12.4,2A1))') (21+i,LF1emat(2,i),V1(21+i), V2(21+i),i=1,2) 
            WRITE(IO,'(" < xz | ",5(I2,":" F12.4,2A1))') (23+i,LF1emat(3,i),V1(23+i), V2(23+i),i=1,3) 
            WRITE(IO,'(" < xy | ",5(I2,":" F12.4,2A1))') (26+i,LF1emat(4,i),V1(26+i), V2(26+i),i=1,4) 
            WRITE(IO,'(" <x2y2| ",5(I2,":" F12.4,2A1))') (30+i,LF1emat(5,i),V1(30+i), V2(30+i),i=1,5)
            WRITE(IO,'(115("-"))')
          else if (Lvalue.eq.3) then
            WRITE(IO,'(127("-"))')
            WRITE(IO,'(16x,"|  z3  >",9X,"| yz2  >",9X,"| xz2  >",9x,"| xyz  >",9x,"|z(x2-y2)>",5x,"|y(x2-3y2)>", &
                                                                                                5X,"|x(3x2-y2)>")')
            WRITE(IO,'("  <  z3   | ",7(I2,":" F12.4,2A1))') (20+i,LF1emat(1,i),V1(20+i), V2(20+i),i=1,1) 
            WRITE(IO,'("  <  yz2  | ",7(I2,":" F12.4,2A1))') (21+i,LF1emat(2,i),V1(21+i), V2(21+i),i=1,2) 
            WRITE(IO,'("  <  xz2  | ",7(I2,":" F12.4,2A1))') (23+i,LF1emat(3,i),V1(23+i), V2(23+i),i=1,3) 
            WRITE(IO,'("  <  xyz  | ",7(I2,":" F12.4,2A1))') (26+i,LF1emat(4,i),V1(26+i), V2(26+i),i=1,4) 
            WRITE(IO,'(" <z(x2-y2)| ",7(I2,":" F12.4,2A1))') (30+i,LF1emat(5,i),V1(30+i), V2(30+i),i=1,5) 
            WRITE(IO,'("<y(x2-3y2)| ",7(I2,":" F12.4,2A1))') (35+i,LF1emat(6,i),V1(35+i), V2(35+i),i=1,6) 
            WRITE(IO,'("<x(3x2-y2)| ",7(I2,":" F12.4,2A1))') (41+i,LF1emat(7,i),V1(41+i), V2(41+i),i=1,7) 
            WRITE(IO,'(127("-"))')
          endif
          if (OUTP(2)) WRITE(IO,'(" Equivalent crystal field:")')
        endif  ! LFtype.eq."LF1E"

        ii=20; LF_str=0; someSmall=.false.
        do k=2,2*Lvalue,2
          do q=0,k
            i=Bkq_index(k,q)
            if (q.eq.0) then 
              ii=ii+1
              if (LFtype.eq."CF  ") WRITE(IO,'(6X,I2,1X,A4,F14.6,2A1)') ii,BkqR_label(i),BkqR(i),v1(ii),v2(ii)
              if (OUTP(2).and. LFtype.ne."CF  ") WRITE(IO,'(6X,A4,F14.6)') BkqR_label(i),BkqR(i)
              LF_str=LF_str+BkqR(i)**2/(D2*k+D1)
              xx=Abs(BkqR(i)); if (xx.gt.bit.and.xx.lt.CFbit) then
!                  write(io,'("Small:",A4,E12.4)') BkqR_label(i),BkqR(i)
                someSmall=.true.
              endif
            else if (q.ne.0) then
              ii=ii+1 
              if (LFtype.eq."CF  ") WRITE(IO,'(2(6X,I2,1X,A4,F14.6,2A1))') ii,BkqR_label(i),BkqR(i),v1(II),v2(II),  &
                                                                      ii+20,BkqI_label(i),BkqI(i),v1(ii+20),v2(ii+20) 
              if (OUTP(2).and. LFtype.ne."CF  ") WRITE(IO,'(2(6X,A4,F14.6))') BkqR_label(i),BkqR(i), &
                                                                               BkqI_label(i),BkqI(i) 
              LF_str=LF_str+d2*(BkqR(i)**2+BkqI(i)**2)/(D2*k+D1)
              xx=Abs(BkqR(i))+abs(BkqI(i)); if (xx.gt.bit.and.xx.lt.CFbit) then
!                  write(io,'("Small:",A4,E12.4,", ",A4,E12.4)') BkqR_label(i),BkqR(i),BkqI_label(i),BkqI(i) 
                  someSmall=.true.
                endif
            endif
          enddo
        enddo  
        LF_str=sqrt(LF_str)
        if (CFcomplex) write(IO,'("Nv/sqrt(4pi)=",F10.2,"   The ligand field is complex sum|Bkq''| =",E12.6)') LF_str,sumCLF
        if (.not.(CFcomplex)) write(IO,'("Nv/sqrt(4pi)=",F10.2,"   The ligand field is real")') LF_str
        if (someSmall) write(IO,'("Parameters where bit=",E8.2," < Bkq  < CFbit=",F6.4," will be set to zero")') bit,CFbit
! CCF        
        IF (CCFtype.ne.0) then
          lab10=" "
          if (CCFtype.eq.1) lab10="delta-func"; if (CCFtype.eq.2) lab10=" spin-CCF "; if (CCFtype.eq.3) lab10="  general "; 
          WRITE(IO,'(/," Correlated Crystal Field Parameters: (",A10," model)",/,115("-"))') lab10
          if (CCFtype.lt.3) then
            WRITE(IO,'(10x,"k  q",12X,"G(k,q)",/,115("-"))') 
            do i=1,N_CCF
              WRITE(IO,'(5X,i3,(2I3),7X,F12.4)') 140+i,(CCFindex(i,j+1),j=1,2),CCF(i)
            enddo
          else
            WRITE(IO,'(10x,"i  k  q",8X,"B(i,k,q)",/,115("-"))') 
            do i=1,N_CCF
              WRITE(IO,'(5X,i3,(3I3),7X,F12.4)') 140+i,(CCFindex(i,j),j=1,3),CCF(i)
            enddo
          endif
        endif
        
        WRITE(IO,'(115("-"))')
        if (ED_IntMech.eq."Full") then
          call writeIntPara(IntN1)
          if (OUTP(4)) then
            WRITE(IO,'(" Equivalent intensity parameters:")')
            call A_B_Mat()
            if (IntN1.eq.1) call writeIntPara(2)
            if (IntN1.eq.2) call writeIntPara(1)
          endif  
       endif
        if (offsetN.gt.0) then
          if (offsetType.eq.1) then
            WRITE(IO,'("Offset Parameters for the multiplets of the diagonalised |LSJ> atomic basis:")')
            write(IO,'((3(3X,2I3,1X,A9,F10.2,5X)))') (350+j,ioffset(j,1),Plab(350+j),Roffset(j),j=1,offsetN)
          elseif (offsetType.eq.2) then
            WRITE(IO,'("Offset Parameters for the diagonalised  energies:")')
            write(IO,'((4(3X,3I3,1X,A9,F10.2,5X)))') (350+j,ioffset(j,1),ioffset(j,2),Plab(350+j),Roffset(j),j=1,offsetN)
          endif
          WRITE(IO,'(115("-"))')
        endif ! offsetN.gt.0
      endif ! icall.eq.1
!      
      if (icall.eq.2) then
        do i=1,maxTotPar
          cL(i)="        "
        enddo
        if (Nfit.gt.0) then
          WRITE(IO,'(/,"Final values for parameters: (*l/*n/L/F indicates Fit linear/Fit nonlinear/Linked/Formula)",/,115("-"))') 
          do i=1,nFit
            write(lab,'("(",F6.1,")")') cov(i,i)
            cL(fitN(i))=lab
          enddo
        endif
        WRITE(IO,32) "           Offset:", 1," EAVE",EAVE, V1(1),V2(1),cL(1)
        if (racah) then
        WRITE(IO,31) "E-E    parameters:", 2,"    B",Bracah,V1(2),V2(2),cL(2),  3,"    C",Cracah,V1(3),V2(3),cL(3)
        else
        WRITE(IO,30) "E-E    parameters:", 2,"   F2",  F2, V1(2),V2(2),cL(2),   3,"   F4",   F4, V1(3),V2(3),cL(3), &
                                           4,"   F6",  F6, V1(4),V2(4),cL(4)
        endif        
        WRITE(IO,32) "SOC    parameter: ", 5," ZETA",ZETA, V1(5),V2(5),cL(5)
        WRITE(IO,30) "E-E CI parameters:", 6,"ALPHA",ALPHA,V1(6),V2(6),cL(6),   7," BETA", BETA, V1(7),V2(7),cL(7),  &
                                           8,"GAMMA",GAMMA,V1(8),V2(8),cL(8)
        WRITE(IO,30) "soo    parameters:", 9,"   M0",  M0, V1(9),V2(9),cL(9),  10,"   M2",  M2,V1(10),V2(10),cL(10), &
                                          11,"   M4",  M4,V1(11),V2(11),cL(11)
        WRITE(IO,30) "ec-so  parameters:",12,"   P2",  P2,V1(12),V2(12),cL(12),13,"   P4",  P4,V1(13),V2(13),cL(13), &
                                          14,"   P6",  P6,V1(14),V2(14),cL(14)
        WRITE(IO,30) "3 body parameters:",15,"   T2",  T2,V1(15),V2(15),cL(15),16,"   T3",  T3,V1(16),V2(16),cL(16), &
                                          17,"   T4",  T4,V1(17),V2(17),cL(17)
        WRITE(IO,30) "                  ",18,"   T6",  T6,V1(18),V2(18),cL(18),19,"   T7",  T7,V1(19),V2(19),cL(19), &
                                          20,"   T8",  T8,V1(20),V2(20),cL(20)
 30     format(A18,i3,1x,A5,"=",f12.4,2A1,A8,2(i3,1x,A5,"=",f12.4,2A1,A8),"|")     
 31     format(A18,i3,1x,A5,"=",f12.4,2A1,A8,i3,1x,A5,"=",f12.4,2A1,A8,32X,"|")     
 32     format(A18,i3,1x,A5,"=",f12.4,2A1,A8,64X,"|")     
!23456789*123456789*123456789*123456789*123456789*123456789*123456789*123456789*123456789*123456789*123456789*12345
        WRITE(IO,'(" Ligand Field Parameters:",89x,"|")') 
        if (LFtype.eq."AOM ") then 
          if (AOM_angles)      WRITE(IO,'(" AOM Ligand field of ",I2," ligands. Angles given, X,Y,Z generated.",51x,"|")') NLIGANDS
          if (.not.AOM_angles) WRITE(IO,'(" AOM Ligand field of ",I2," ligands. X,Y,Z given, angles generated.",51x,"|")') NLIGANDS
          WRITE(IO,'(115("-"),/,  &
  "   Name   |   E(sigma)",16X,"E(piX)",10X,"E(piY)",14X,"|",2x,"Theta",15x,"Phi"16x,"Chi",16X,"|",7x,"X       Y       Z   |")')
          do i=1,NLIGANDS
            write(IO,'(4X,A4,2X,"|",3(I3,F8.2,2A1,A8),"|",3(I3,F7.2,2A1,A8),"|",2X,3(1X,F7.4))') &
                  Lname(i),20+i,eSigma(i),V1(20+i), V2(20+i), cL(20+i),      &
                           40+i,  ePiX(i),V1(40+i), V2(40+i), cL(40+i),      &
                           60+i,  ePiY(i),V1(60+i), V2(60+i), cL(60+i),      &
                           80+i, theta(i),V1(80+i), V2(80+i), cL(80+i),      &
                          100+i,   phi(i),V1(100+i),V2(100+i),cL(100+i),     &
                          120+i,   chi(i),V1(120+i),V2(120+i),cL(120+i),xLig(i),yLig(i),zLig(i)
          enddo
          WRITE(IO,'(115("-"))')
          if (AOMX) then 
            WRITE(IO,'(" Extended AOM Ligand field")')
            WRITE(IO,'(115("-"),/,"   Name   |   E(DeltaS)     E(DeltaC)        E(PsiS)      E(PsiC)  |")')
            do i=1,NLIGANDS
              write(IO,'(4X,A4,2X,"|",4(1X,I2,F9.2,2A1,A8),2X,"|")') Lname(i),  &
                  140+i,eDelS(i),V1(140+i),V2(140+i),cL(140+i), 150+i,eDelC(i),V1(150+i),V2(150+i),cL(150+i),   &
                  160+i,ePsiS(i),V1(160+i),V2(160+i),cL(160+i), 170+i,ePsiC(i),V1(170+i),V2(170+i),cL(170+i)
            enddo
            WRITE(IO,'(115("-"))')
          endif
          if (OUTP(2)) WRITE(IO,'(" Equivalent crystal field:")')
        endif  ! LFtype.eq."AOM "

        ii=20; LF_str=0; someSmall=.false.
        if (FitType.ne.3) then ! The CF coefficients are parameters
          do k=2,2*Lvalue,2
            do q=0,k
              i=Bkq_index(k,q)
              if (q.eq.0) then 
                ii=ii+1
                if (LFtype.eq."CF  ") WRITE(IO,'(6X,I2,1X,A4,F14.6,2A1,A8,77x,"|")') ii,BkqR_label(i),BkqR(i),v1(ii),v2(ii),cL(ii)
                if (OUTP(2).and. LFtype.ne."CF  ") WRITE(IO,'(6X,A4,F14.6,90x,"|")') BkqR_label(i),BkqR(i)
                LF_str=LF_str+BkqR(i)**2/(D2*k+D1)
                xx=Abs(BkqR(i)); if (xx.gt.bit.and.xx.lt.CFbit) someSmall=.true.
              else if (q.ne.0) then
                ii=ii+1 
                if (LFtype.eq."CF  ") WRITE(IO,'(2(6X,I2,1X,A4,F14.6,2A1,A8),40x,"|")') ii,BkqR_label(i),BkqR(i),v1(ii),v2(ii), &
                                                              cL(ii),ii+20,BkqI_label(i),BkqI(i),v1(ii+20),v2(ii+20),cL(ii+20) 
                if (OUTP(2).and. LFtype.ne."CF  ") WRITE(IO,'(2(6X,A4,F14.6),66x,"|")') BkqR_label(i),BkqR(i), &
                                                                                 BkqI_label(i),BkqI(i) 
                LF_str=LF_str+d2*(BkqR(i)**2+BkqI(i)**2)/(D2*k+D1)
                xx=Abs(BkqR(i))+abs(BkqI(i)); if (xx.gt.bit.and.xx.lt.CFbit) someSmall=.true.
              endif
            enddo
          enddo
        else if (FitType.eq.3) then  ! The CF coefficients are fitted:
          WRITE(IO,'(" Crystal Field Coefficients fitted:",79X,"|")')
          sum1=d0
          do k=2,2*Lvalue,2
            do q=0,k
              i=Bkq_index(k,q)
              if (NAssign(i,1).lt.0) sum1=sum1+Wgt_E(i)*(EXP_E(i)-BKQI(abs(NAssign(i,1))))**2
              if (NAssign(i,1).gt.0) sum1=sum1+Wgt_E(i)*(EXP_E(i)-BKQR(NAssign(i,1)))**2
              if (q.eq.0) then 
                ii=ii+1; j1=0
                do j=1,nexp; if (Nassign(j,1).eq.i) j1=j; enddo
                if (j1.eq.0) then
                  WRITE(IO,'(6X,I2,1X,A4,F14.6,87x,"|")') ii,BkqR_label(i),BkqR(i)
                else
                  WRITE(IO,'(6X,I2,1X,A4,F14.6," F ",F11.4,73x,"|")') ii,BkqR_label(i),BkqR(i),EXP_E(j1)
                endif 
                LF_str=LF_str+BkqR(i)**2/(D2*k+D1)
                xx=Abs(BkqR(i)); if (xx.gt.bit.and.xx.lt.CFbit) someSmall=.true.
              else if (q.ne.0) then
                ii=ii+1; j1=0; j2=0 
                do j=1,nexp; if (Nassign(j,1).eq.i) j1=j; if (Nassign(j,1).eq.-i) j2=j; enddo
                if (j1.eq.0 .and. j2.eq.0) then
                  WRITE(IO,'(6X,I2,1X,A4,F14.6,17X,           I2,1X,A4,F14.6,49x,"|")') & 
                                                  ii,BkqR_label(i),BkqR(i),20+ii,BkqI_label(i),BkqI(i) 
                else if (j1.eq.0 .and. j2.ne.0) then
                  WRITE(IO,'(6X,I2,1X,A4,F14.6,17X,           I2,1X,A4,F14.6," F ",F11.4,35x,"|")') &
                                                  ii,BkqR_label(i),BkqR(i),20+ii,BkqI_label(i),BkqI(i),EXP_E(j2) 
                else if (j1.ne.0 .and. j2.eq.0) then
                  WRITE(IO,'(6X,I2,1X,A4,F14.6," F ",F11.4,3X,I2,1X,A4,F14.6,49x,"|")')  &
                                                  ii,BkqR_label(i),BkqR(i),EXP_E(j1),20+ii,BkqI_label(i),BkqI(i) 
                else if (j1.ne.0 .and. j2.ne.0) then
                  WRITE(IO,'(6X,I2,1X,A4,F14.6," F ",F11.4,3X,I2,1X,A4,F14.6," F ",F11.4,35x,"|")') & 
                                                  ii,BkqR_label(i),BkqR(i),EXP_E(j1),20+ii,BkqI_label(i),BkqI(i),EXP_E(j2) 
                endif                
                LF_str=LF_str+d2*(BkqR(i)**2+BkqI(i)**2)/(D2*k+D1)
                xx=Abs(BkqR(i))+abs(BkqI(i)); if (xx.gt.bit.and.xx.lt.CFbit) someSmall=.true.
              endif
            enddo  ! q
          enddo  ! k
          if (Nexp.ne.0) sum1=sqrt(sum1)/Nexp
          if (nFit.eq.0) WRITE(IO,'("sqrt(sum(Bkq-Bex)^2/Nex)=",E12.5,72X," |")') SUM1 
          if (nFit.ne.0) WRITE(IO,'("sqrt(sum(Bkq-Bex)^2/Nex)=",E12.5,5X,   & 
                                  "(with initial parameters=",E12.5,")",34x,"|")') SUM1,fnorm_init 

        endif        
        LF_str=sqrt(LF_str)
        if (CFcomplex) write(IO,'("Nv/sqrt(4pi)=",F10.2,"   The ligand field is complex sum|Bkq''| =",E12.6,37x,"|")') LF_str,sumCLF
        if (.not.(CFcomplex)) write(IO,'("Nv/sqrt(4pi)=",F10.2,"   The ligand field is real",64x,"|")') LF_str
        if (someSmall) write(IO,'("Parameters where bit=",E8.2," < Bkq  < CFbit=",F6.4," will be set to zero",43x,"|")') bit,CFbit
        WRITE(IO,'(115("-"))')
        
! CCF        
        IF (CCFtype.ne.0) then
          lab10=" "
          if (CCFtype.eq.1) lab10="delta-func"; if (CCFtype.eq.2) lab10=" spin-CCF "; if (CCFtype.eq.3) lab10="  general "; 
          WRITE(IO,'(" Correlated Crystal Field Parameters: (",A10," model)",/,115("-"))') lab10
          if (CCFtype.lt.3) then
            WRITE(IO,'(10x,"k  q",12X,"G(k,q)",/,115("-"))') 
            do i=1,N_CCF
              WRITE(IO,'(5X,i3,(2I3),7X,F12.4,2A1,A8)') 140+i,(CCFindex(i,j+1),j=1,2),CCF(i),v1(140+i),v2(140+i),cL(140+i)
            enddo
          else
            WRITE(IO,'(10x,"i  k  q",8X,"B(i,k,q)",/,115("-"))') 
            do i=1,N_CCF
              WRITE(IO,'(5X,i3,(3I3),7X,F12.4,2A1,A8)') 140+i,(CCFindex(i,j),j=1,3),CCF(i),v1(140+i),v2(140+i),cL(140+i)
            enddo
          endif
        endif

        if (offsetN.gt.0) then
          WRITE(IO,'("Offset Parameters:")')
          write(IO,'((3(3X,I3,1X,A9,F10.2,2A1,A8,2X)))') (350+i,Plab(350+i),Roffset(i),v1(350+i),v2(350+i),cL(350+i),i=1,offsetN)
          WRITE(IO,'(115("-"))')
        endif
       if (nvar.gt.0) then
          WRITE(IO,'("Variables:")'); 
          write(IO,'(3(I4,1X,A10,F10.4,2A1,A8,2X))') (400+i,variables(i),P(400+i),v1(400+i),v2(400+i),cL(400+i),i=1,nVar)
          WRITE(IO,'(115("-"))')
        endif
      
        if (FitOpt(2)) then
          write(io,'(/,"Covariance Matrix:")')
          write(io,'(14X,10(A9,3X))') (Plab(fitN(i)),i=1,nFit)
          do i=1,nFit
            write(io,'(A9,3X,10(F10.4,2X))') Plab(fitN(i)),(cov(i,j),j=1,nFit)        
          enddo  
          WRITE(IO,'(115("-"))')
        endif
        if (FitOpt(4)) then
          OPEN(ifit,FILE="fit.out",STATUS='UNKNOWN')
          write(ifit,'(/,"Fitted parameters")')
          do i=1,nFit
            write(ifit,'(A9,3X,F18.8)') Plab(fitN(i)),P(fitN(i))        
          enddo 
          close(ifit)          
        endif
        
      endif  ! icall.eq.2
!
      return
      end subroutine PrintParameters
!
!-----------------------------------------------------------------------
!
      SUBROUTINE DECODE(A,N,COMMAND,COMD)
!
!  This subroutine decodes the character string in COMMAND and returns 
!  the N (real) numerical elements in the array A. 
!  A maximum of 100 elements may be returned in A.
!
      IMPLICIT NONE
      REAL*8 A(100)
      INTEGER*4 i,j,n,n1,n2
      logical found
      CHARACTER COMMAND*200, cmd*20, COMD*4
      CHARACTER*1 M(12)/'0','1','2','3','4','5','6','7','8','9','.','-'/
!
!      write(io,'("decoding:",A120)') COMMAND
      N=0; N1=1
  10  continue
      found=.false.  
      do i=N1,200
        DO J=1,12
          IF (Command(N1:N1).eq.M(J)) found=.true.
        enddo
        if (found) exit        
        N1=N1+1
      enddo !i  
      N2=N1+1
      if (N2.gt.200) return      
      do I=N1+1,200
        if (Command(N2:N2).eq." ".or.Command(N2:N2).eq.",") exit
        N2=N2+1
        if (N2.gt.200) then; write(io,'("****FATAL: Input line too long:",A200)') Command; stop; endif
      enddo ! I
      if (N1.lt.N2) then
        cmd=Command(N1:N2)
        N=N+1
        Read(cmd,*,ERR=20) A(N)
!        write(io,'("N1=",I4,", N2=",I4,", cmd=",A20,", N=",I4,", A(N)=",F8.2)') N1,N2,cmd,N,A(N)
        N1=N2+1
      endif  
      goto 10
      
      RETURN
 20   Write(IO,'("***FATAL: Error in DECODE, reading Command: ",A4)') COMD
      STOP
      END SUBROUTINE DECODE
!
!-----------------------------------------------------------------------
!
      logical function CheckVariableIsNumberOnly(COMMAND)
!
!  This function checks to see if the character string in COMMAND is composed of only numbers 
!
      IMPLICIT NONE
      INTEGER*4 i
      logical notNumber
      CHARACTER COMMAND*200
      CHARACTER*13 number/'0123456789.- '/
!
      notNumber=.false.  
      do i=1,200
        IF (scan(Number,COMMAND(i:i)).eq.0) notNumber=.true.
!        if (notNumber) write(io,'("****WARNING: Non-numerical character found:",A1)') COMMAND(i:i)
        if (notNumber) exit        
      enddo !i  
      CheckVariableIsNumberOnly = notNumber
!      
      RETURN
      END function CheckVariableIsNumberOnly
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OUTPUT(neng)
!
!  Called with neng=0 for JO calculation (option(8)).
!
!  TODO does not give S & MJ output in plotfile for IRREPS
! 
      IMPLICIT none
      REAL*8 AMV(max_N),E0,EE,DE,sum1,sum2,SumInt,DelSqE,DelSqI,F1,TotI
      INTEGER*4 NCLASS1,I,J,K,kk,k1,II,IDEG,IWRITE,JHIGH,neng,sumn,IG,IE
      integer*4 TwoJ,LS,indexAMV(max_N),iblock,values(8),n1,n2,n3,i3,assNum(max_exp) 
      logical nonKramers,last,notLast,test,test1
      character lab*3, label*13, symlab*6, Format1a*40,Format1b*40,Format1c*40 
      CHARACTER FMT5*61,FMT6*61,bit*81 
      CHARACTER LINE1*150,LINE2*120,thisRun*18
      CHARACTER ExpEng*10,DelEng*10,ExpInt*12,CalcInt*12,DelInt*12, WgtEng*8,WgtInt*8
!
      if (option(8)) goto 100  ! skip writing matrix & diagonalisation
      sumn=0
      write(FMT5,'("(",I1,"F6.2,"" |"")")') NSPIN
      write(FMT6,'("(",I1,"I6,"" |"")")') NSPIN
!
      if (option(1)) E0=Engs(1)   !  Order relative to lowest energy
      call zeroE(Engs,neng)
!  Plotfile    
      test=.false.; do k=1,10; if (PLOTout1(k)) test=.true.; if (PLOTout2(k)) test=.true.; enddo
      if (test) call PlotfileOUT(neng)
!
      if (wordy) then
        if (PrangeL.ne.PrangeH) write(io,'(4X,"*wordy: Only levels with energies in range: [",F12.2,",",F12.2, &
                                                                              "] will be printed.")') PrangeL,PrangeH
        if (NrangeL.ne.NrangeH) write(io,'(4X,"*wordy: Only levels: [",I4,",",I4,"] will be printed.")') NrangeL,NrangeH
      endif 
!      
      IF (Option(1)) write(io,'(4X,"Energy relative to lowest calculated E0=",F11.1)') E0
      line1=" "
      line2="  Level    Rel.Eng. (deg)|"
      if (EngInt.eq.1) line2=trim(line2)//"   - Exp     del    |" 
      IF (IDG.ne.0) then
        if (EngInt.eq.0) write(line1,'(25X,"| ",A4," group |")') GROUP
        if (EngInt.ne.0) write(line1,'(46X,"| ",A4," group |")') GROUP
        line2=trim(line2)//"   Symmetry |"
      endif
      if (PrintCQN) line2=trim(line2)//" CrQN|"
      IF (Option(2))  then
        if (EngInt.eq.0 .and. IDG.eq.0 .and. .not.PrintCQN) write(line1,'(25X,"|  2S+1 ")')
        if (EngInt.eq.0 .and. IDG.eq.0 .and.  PrintCQN)     write(line1,'(31X,"|  2S+1 ")')
        if (EngInt.ne.0 .and. IDG.eq.0 .and. .not.PrintCQN) write(line1,'(46X,"|  2S+1 ")')
        if (EngInt.ne.0 .and. IDG.eq.0 .and. PrintCQN)      write(line1,'(50X,"|  2S+1 ")')
        if (EngInt.ne.0 .and. IDG.ne.0) line1=trim(line1)//"  2S+1 "
        write(bit,FMT6) (S_mult(K),K=1,NSPIN)
        line2=trim(line2)//bit
      endif
      if (SpAllow(2)-SpAllow(1).gt.0) line2=trim(line2)//"SpAllow|"   ! spin-allowedness
      !                                  123456789*123456789*123456789*123456789*
      if (Option(6)) line2=trim(line2)//"      MJ projections                   |"
      !                                  123456789*123456789*123456789*123456789*123456789*123456789*123456789*123456789* 
      if (Option(7)) line2=trim(line2)//"      MJ projections                                                           |"
      IF (Option(3).or.Option(4)) line2=trim(line2)//"     Free Ion %"
      write(io,'(A)') trim(line1)
      write(io,'(A)') trim(line2)
      WRITE(IO,'(115("-"))')
!      
      IDEG=1; sum1=0.0d+00; sum2=0.0d+00; 
      write(Format1a,'("(I4,5X,F11.",i1,",5X,""|"")")') Ndecpts                    ! single energy levels
      write(format1b,'("(I4,""-"",I4,F11.",i1,",1X,""("",I2,"")|"")")') Ndecpts    ! summed over deg energy levels
      write(format1c,'("("" -"",F10.",i1,",F7.",i1,",a1,""|"")")') Ndecpts,Ndecpts
!      write(io,*) Format1a
!      write(io,*) Format1b
!      write(io,*) Format1c
      if (assignMethod.eq.1 .or. assignMethod.eq.2) call AssignNum(assNum,neng)
!
      DO 200 I=1,neng
        IWRITE=1
        J=n_matrix-I+1
        IF (neng.NE.n_matrix) J=neng-I+1
        IF (IDEG.EQ.1) JHIGH=J
        EE=engs(J)
!        test1 = PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                    &
!         .or. PrangeL.ne.PrangeH .and. (EE.gt.PrangeL .and. EE.lt.PrangeH)   & 
!         .or. NrangeL.ne.NrangeH .and. ( J.ge.NrangeL .and.  J.le.NrangeH) 
!        write(io,'("J=",I4,"; NrangeL,NrangeH:",2I4,"; test=",L1)') J,NrangeL,NrangeH, test1
        if  ( PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                    &
         .or. PrangeL.ne.PrangeH .and. (EE.gt.PrangeL .and. EE.lt.PrangeH)   & 
         .or. NrangeL.ne.NrangeH .and. ( J.ge.NrangeL .and.  J.le.NrangeH)  ) then
          IF (J.GT.1) THEN
            IF (.not.OUTP(6) .AND. ABS(EE-engs(J-1)).LT.Edegen) THEN
              IDEG=IDEG+1
              IWRITE=0
            ENDIF
          ENDIF
          IF (IWRITE.EQ.1) THEN
            IF (IDEG.EQ.1) WRITE(line1,Format1a) J,      engs(J)
            IF (IDEG.GT.1) WRITE(line1,Format1b) J,JHIGH,engs(J),IDEG
            if (EngInt.eq.1 ) then
              write(bit,'(19X," |")') ! experimental
              call AssignCalcJ(assNum,J,Jhigh,KK)  ! assigns calculated E(J) to Exp_E(KK), KK=0 otherwise
              if (kk.ne.0) then
                write(bit,Format1c) exp_E(kk),engs(j)-exp_E(kk),wgtchar(kk)
                if (FitOpt(7)) then  ! sum over degeneracies
                  sum1=sum1+ideg*wgt_e(kk)*(engs(j)-exp_E(kk))**2 
                  sum2=sum2+ideg*wgt_e(kk)*abs(engs(j)-exp_E(kk)) 
                  sumn=sumn+ideg
                else  
                  sum1=sum1+wgt_e(kk)*(engs(j)-exp_E(kk))**2 
                  sum2=sum2+wgt_e(kk)*abs(engs(j)-exp_E(kk)) 
                  sumn=nexp
                endif  
              endif
              line1=trim(line1)//bit
            endif
            if (IDG.ne.0) then  ! symmetry labels
              write(bit,'(3A4,"|")') IRREP(J-1+IDEG,3),IRREP(J-1+IDEG,1),IRREP(J-1+IDEG,2)
              line1=trim(line1)//bit
            endif
            if (PrintCQN) then
              write(bit,'(A5,"|")') allCQN(J)
              line1=trim(line1)//bit
            endif
            if (Option(2)) then  ! spin projections
              write(bit,FMT5) (SPIN(J,K),K=1,NSPIN)
              line1=trim(line1)//bit
            endif
            if (SpAllow(1)+SpAllow(2).gt.0) then  ! spin-allowedness
              write(bit,'(F6.2," |")') SpinAllowed(J)
              line1=trim(line1)//bit
            endif
            if (Option(6)) then  ! MJ %
              write(bit,'(A40,"|")') cAllMJ(J)
              line1=trim(line1)//bit
            endif
            if (Option(7)) then  ! MJ %
              write(bit,'(A80,"|")') cAllMJ(J)
              line1=trim(line1)//bit
            endif
            if (Option(3).or.Option(4)) then  ! free ion %
              line1=trim(line1)//cFreeIon(J)
            endif
            write(io,'(A)') trim(line1)
!
            IDEG=1
          endif  ! IF (IWRITE.EQ.1) 
        ENDIF  ! if (EE.gt.PrangeL .and. EE.lt.PrangeH) or (i.gt.PrangeL .and. i.lt.PrangeH) 
 200  CONTINUE ! i
      WRITE(IO,'(115("-"))')
      if (EngInt.ne.0) then
        if (sumn.ne.0) then; sum1=sqrt(SUM1/sumn); sum2=SUM2/sumn; endif
        if (nFit.eq.0) WRITE(IO,'("sqrt(sum(Ec-Eex)^2/Nex)= |     ",E12.5,"   |")') SUM1 
        if (nFit.ne.0) WRITE(IO,'("sqrt(sum(Ec-Eex)^2/Nex)= |     ",E12.5,"   | With initial parameters=",E12.5)') SUM1,fnorm_init
        if (FitOpt(8)) WRITE(IO,'("        sum|Ec-Eex|/Nex= |     ",E12.5,"   | ")') sum2 
        if (wgtNote) WRITE(IO,'("Some levels not weighted(x), or weighted more(+) in comparison with experiment")')
        WRITE(IO,'(115("-"))') 
        if (PlotOut2(1)) Yplot2(1)=SUM1
        if (TestOutput)  then
          lab="   "; if (sum1.gt.TestCrit) lab="***" 
          call DATE_AND_TIME(VALUES=values)
          write(thisRun,'("(",I2,"/",I2.2,"/",I4,") ",I2,":",I2.2)') &
                             values(3),values(2),values(1),values(5),values(6)
          write(iTest,'(A3,"ver:",A6,"; compile:",A,"; this run:",A18,"; sqrt(delE^2/N)=",F10.4,"; input file:",A)')  &
                             lab,version,trim(CompileDate),thisRun,sum1,trim(infile)
        endif                     
      endif  
!
!  g-values
      if (GEXS(1).ne.0) then
        if (NexpG.eq.0) then
          if (g_axes) WRITE(IO,'(/,"g-Values",60x,"principal axes",/,115("-"),/,  &
                       " States   Energy        g1      g2      g3          1",22x,"2",22x,"3"/,115("-"))')
          if (.not.g_axes) WRITE(IO,'(/,"g-Values",/,115("-"),/,  &
                       " States   Energy        g1      g2      g3"/,115("-"))')
          do i=1,ngexs  ! (GEXS(2)-GEXS(1)+1)/2
            nonKramers=.false.;          
            if (abs(engs(GEXS(2*i))-engs(GEXS(2*i-1))).gt.Edegen) nonKramers=.true.; 
            if (g_axes) write(io,'(i3,"-",i3,F9.2,4X,3F8.4,3(2X,3F6.3))') gexs(2*i-1),gexs(2*i),engs(gexs(2*i)),   &
                                    (gval(j,i),j=1,3),(g_prin(1,j,i),g_prin(2,j,i),g_prin(3,j,i),j=1,3)
            if (.not.g_axes) write(io,'(i3,"-",i3,F9.2,4X,3F8.4)') gexs(2*i-1),gexs(2*i),engs(gexs(2*i)),(gval(j,i),j=1,3) 
            if (nonKramers) write(io,'("****Warning: not a Kramers doublet: deltaE >",F8.4)') Edegen                                 
          enddo ! ! i
          WRITE(IO,'(115("-"))')
        else if (NexpG.gt.0) then
          sum1=z; sumn=0
          WRITE(IO,'(/,"g-Values",/,115("-"),/,  &
          " States   Energy |   g1(calc)  g1(exp)  diff  |   g2(calc)  g2(exp)  diff  |   g3(calc)  g3(exp)  diff  |", &
                       /,115("-"))')
          do i=1,ngexs  ! (GEXS(2)-GEXS(1)+1)/2
            kk=0
            line1=" "; line2=" "
            do k=1,nexpG
              if (NAssignG(k).eq.GEXS(2*i).or.NAssignG(k).eq.GEXS(2*i-1)) then
                kk=kk+1
                write(line1,'(i3,"-",i3,F8.1,"  |")') gexs(2*i-1),gexs(2*i),engs(gexs(2*i))
                if (wgt_g(1,i).gt.z) then; write(line2,'(3F9.4," |")') gval(1,i),exp_G(1,k),gval(1,i)-exp_G(1,k)
                else; write(line2,'(F9.4,18X," |")') gval(1,i); endif
                line1=trim(line1)//line2
                if (wgt_g(2,i).gt.z) then; write(line2,'(3F9.4," |")') gval(2,i),exp_G(2,k),gval(2,i)-exp_G(2,k)
                else; write(line2,'(F9.4,18X," |")') gval(2,i); endif
                line1=trim(line1)//line2
                if (wgt_g(3,i).gt.z) then; write(line2,'(3F9.4," |")') gval(3,i),exp_G(3,k),gval(3,i)-exp_G(3,k)
                else; write(line2,'(F9.4,18X," |")') gval(3,i); endif 
                line1=trim(line1)//line2
                write(io,'(A)') trim(line1)
                sum1 = sum1 + wgt_g(1,i)*(gval(1,i)-exp_G(1,k))**2 + wgt_g(2,i)*(gval(2,i)-exp_G(2,k))**2  &
                                                                   + wgt_g(3,i)*(gval(3,i)-exp_G(3,k))**2
                if (wgt_g(1,i).ne.z) sumn=sumn+1; if (wgt_g(2,i).ne.z) sumn=sumn+1; if (wgt_g(3,i).ne.z) sumn=sumn+1
              endif
            enddo ! k
            if (kk.eq.0) write(io,'(i3,"-",i3,F8.1,"  | ",F8.4,18X," | ",F8.4,18X," | ",F8.4,18X," |")')  &
                                           gexs(2*i-1),gexs(2*i),engs(gexs(2*i)), gval(1,i),gval(2,i),gval(3,i)
          enddo ! i
          WRITE(IO,'(115("-"))')
          sum1=sqrt(sum1/sumn)
          if (nFit.eq.0) WRITE(IO,'("sqrt(sum(Gc-Gex)^2/Nex)= |     ",E12.5," |")') SUM1 
          if (nFit.ne.0) WRITE(IO,'("sqrt(sum(Gc-Gex)^2/Nex)= |     ",E12.5,   & 
                                  " |  (with initial parameters value=",E12.5,")")') SUM1,fnorm_init 
          WRITE(IO,'(115("-"))')
          if (TestOutput)  then
            lab="   "; if (sum1.gt.TestCrit) lab="***" 
            call DATE_AND_TIME(VALUES=values)
            write(thisRun,'("(",I2,"/",I2.2,"/",I4,") ",I2,":",I2.2)') &
                               values(3),values(2),values(1),values(5),values(6)
            write(iTest,'(A3,"ver:",A6,"; compile:",A,"; this run:",A18,"; sqrt(delG^2/N)=",F10.4,"; input file:",A)')  &
                               lab,version,trim(CompileDate),thisRun,sum1,trim(infile)
          endif                     
        endif ! NexpG.gt.0
      endif
!
!  Magnetic Dipole Transitions
      if (IMD.NE.0) then
        F1=1.0d+0; if (MDunits.eq.2) F1=9.274009d-03*1.0d+03     ! 1 BM = 9.274009d-3 (10^-3 Debye)
                   if (MDunits.eq.3) F1=9.274009d-03/4.80319d-04  ! 1 BM = 9.274009d-3/4.80319d04 (10^-12 cm)
!                   if (MDunits.eq.3) F1=4.80319d04/9.274009d-03  ! 1 BM = 9.274009d-3/4.80319d04 (10^-12 cm)
        F2=1.0d+0; if (MDunits.eq.2) F2=(9.274009d-03)**2*1.0d+07     ! 1 BM^2 = (9.274009d-3)^2 (10^-7 Debye^2)
                   if (MDunits.eq.3) F2=(9.274009d-03/4.80319d-04)**2  ! 1 BM^2 =  (10^-24 cm^2)
!                   if (MDunits.eq.3) F2=(4.80319d04/9.274009d-03)**2  ! 1 BM^2 =  (10^-24 cm^2)
        if (MDconst.eq.1.0d+00) WRITE(IO,'(/," Magnetic Dipole Transitions")')
        if (MDconst.ne.1.0d+00 .and. IMD.ne.3) WRITE(IO,'(/," Magnetic Dipole Transitions   (multiplied by:",F10.4,")")')MDconst
          WRITE(IO,'(115("-"),/," kx=",F6.4,"  ky=",F6.4,"  kz=",F6.4)')  (RK(i),i=1,3)
          if (IMD.eq.1) then
            if (MDunits.eq.1) WRITE(IO,'("    Energy(deg.)    Dipole strength(Bohr Mag**2)")')
            if (MDunits.eq.2) WRITE(IO,'("    Energy(deg.)    Dipole strength(x10-7 Debye**2)")')
            if (MDunits.eq.3) WRITE(IO,'("    Energy(deg.)    Dipole strength(x10-24 cm**2)")')
            WRITE(IO,'(20X,"deltaE       <Mx>**2     <My>**2     <Mz>**2       Sum/3",/,115("-"))')
          elseif (IMD.eq.2) then
            if (MDunits.eq.1) WRITE(IO,'("    Energy(deg.)    Dipole strength(Bohr Mag**2)")')  
            if (MDunits.eq.2) WRITE(IO,'("    Energy(deg.)    Dipole strength(x10-7 Debye**2)")')  
            if (MDunits.eq.3) WRITE(IO,'("    Energy(deg.)    Dipole strength(x10-24 cm**2)")')
            WRITE(IO,'(20X,"deltaE     <Mx>**2   <My>**2   <Mz>**2   <M+>**2   <M->**2",3X,"DA(MCD)   DA(MLD)",/,115("-"))')
          elseif (IMD.eq.3) then
            if (MDunits.eq.1) WRITE(IO,'("    Energy(deg.)",10X,"Mag.Dipole Amplitude(Bohr Mag.)")') 
            if (MDunits.eq.2) WRITE(IO,'("    Energy(deg.)",10X,"Mag.Dipole Amplitude(x10-3 Debye)")') 
            if (MDunits.eq.3) WRITE(IO,'("    Energy(deg.)",10X,"Mag.Dipole Amplitude(x10-12 cm)")')
            WRITE(IO,'(21X,"deltaE",9X,"<Mx>",16X,"<My>/i",14X,"<Mz>",/,115("-"))')
          endif
!
          if (MDG.le.0 .or. MDE.le.0) then; write(io,'("****FATAL No magnetic dipole transitions calculated.")'); stop; endif
          do I=1,MDG
            WRITE(IO,'(F9.1," (",I2,") -->")') RENGG(I),IDEGG(I)
            DO J=1,MDE
              EE=ABS(RENGE(J)-RENGG(I)) 
              TotI=F2*MDconst*(RMD(I,J,1)+RMD(I,J,2)+RMD(I,J,3))/3.0d+00 
              IF(IMD.EQ.1) WRITE(IO,'(F12.1," (",I2,")",F10.1,3E12.4,2X,E12.4)') RENGE(J),IDEGE(J),EE, &
                                                                       (F2*MDconst*RMD(I,J,K),K=1,3),TotI
!              IF(IMD.EQ.1) WRITE(IO,'(F12.1," (",I2,")",F10.1,3F11.5,F13.5)') RENGE(J),IDEGE(J),EE, &
!                                                                       (F2*MDconst*RMD(I,J,K),K=1,3),TotI
              IF(IMD.EQ.2) WRITE(IO,'(F12.1," (",I2,")",F10.1,7F10.5)') RENGE(J),IDEGE(J),EE,  &
                                       (MDconst*F2*RMD(I,J,K),K=1,5),MDconst*F2*(RMD(I,J,5)-RMD(I,J,4)), &
                                        MDconst*F2*(RMD(I,J,1)-RMD(I,J,2))
              IF(IMD.EQ.3) WRITE(IO,'(F12.1," (",I2,")",F10.1,3(2X,2F8.4))') RENGE(J),IDEGE(J),EE,(F1*CMD(I,J,K),K=1,3)
            enddo ! J
          enddo ! I
          WRITE(IO,'(115("-"))')
      endif ! IMD
!
!  Electric Dipole Transitions
      if (IED.NE.0) then
        if (EDconst.eq.1.0d+00) WRITE(IO,'(/," Electric Dipole Transitions")')
        if (EDconst.ne.1.0d+00 .and. IED.ne.3) WRITE(IO,'(/," Electric Dipole Transitions   (multiplied by:",F10.4,")")')EDconst
        WRITE(IO,'(115("-"))')
        if (IED.eq.1) then
          if (EDunits.eq.2) WRITE(IO,'("    Energy(deg.)    Dipole strength(x10^-7 D^2)")')
          if (EDunits.eq.3) WRITE(IO,'("    Energy(deg.)    Dipole strength(x10^-24 cm^2)")')
          WRITE(IO,'(19X,"deltaE     <Ex>**2   <Ey>**2   <Ez>**2     Sum/3",/,115("-"))')
        elseif (IED.eq.2) then
          if (EDunits.eq.2) WRITE(IO,'("    Energy(deg.)    Dipole strength(x10^-7 D^2)")')
          if (EDunits.eq.3) WRITE(IO,'("    Energy(deg.)    Dipole strength(x10^-24 cm^2)")')
          WRITE(IO,'(19X,"deltaE     <Ex>**2   <Ey>**2   <Ez>**2   <E+>**2   <E->**2",3X,"DA(MCD)   DA(MLD)",/,115("-"))')
        elseif (IED.eq.3) then
          if (EDunits.eq.2) WRITE(IO,'("    Energy(deg.)",10X,"Dipole Amplitude(x10^-7 D)")')
          if (EDunits.eq.3) WRITE(IO,'("    Energy(deg.)",10X,"Dipole Amplitude(x10^-12 cm)")')
          WRITE(IO,'(21X,"deltaE",9X,"<Ex>",16X,"<Ey>/i",14X,"<Ez>",/,115("-"))')
        endif
!
! Debye defined as 10^-19/c  C cm   where c is speed of light in m/s.
! The multiplication factor below is from multiplying by electron charge (in C), speed of light (in m/s)
!  1.60217602D-19 x (2.99792458D+08/10^-19) x 1.0D-12, last factor due to Altp parameters in unit 10^-12 cm 
!  = 4.80319d-04

        F1=1.0d00;                                                     ! ToDo 
        if (EDunits.eq.2) F1=-1.60217602d00*2.99792458d-04*1.0d+07     ! F1 converts to (x10-7 D)
        if (EDunits.eq.3) F1=1.0d00                                    ! F1 converts to (x10-12 cm)
        F2=1.0d00;                                                     ! ToDo
        if (EDunits.eq.2) F2=(1.60217602d00*2.99792458d-04)**2*1.0d+07 ! F2 converts to (x10^-7 D^2)
        if (EDunits.eq.3) F2=1.0d00                                    ! F2 converts to (x10-24 cm^2)
                                        
        do I=1,EDG
          WRITE(IO,'(F9.1," (",I2,") -->")') RENGG(I),IDEGG(I)
          DO J=1,EDE
            EE=ABS(RENGE(J)-RENGG(I)) 
            TotI=F2*EDconst*(RED(I,J,1)+RED(I,J,2)+RED(I,J,3)) 
            IF(IED.EQ.1) WRITE(IO,'(F12.1," (",I2,")",F10.1,3F10.3,F12.3)') RENGE(J),IDEGE(J),EE,  &
                                                                            (F2*EDconst*RED(I,J,K),K=1,3),TotI
            IF(IED.EQ.2) WRITE(IO,'(F12.1," (",I2,")",F10.1,7F10.3)') RENGE(J),IDEGE(J),EE,  &
               (F1*EDconst*RED(I,J,K),K=1,5),F1*EDconst*(RED(I,J,5)-RED(I,J,4)),F2*EDconst*(RED(I,J,1)-RED(I,J,2))
            IF(IED.EQ.3) WRITE(IO,'(F12.1," (",I2,")",F10.1,3(2X,2F9.1))') RENGE(J),IDEGE(J),EE,(F1*CED(I,J,K),K=1,3)
          enddo ! J
        enddo ! I
        WRITE(IO,'(115("-"))')
      endif ! IED
!
! Skip main calculation
 100  continue
!  Magnetic Dipole Transitions in |LSJ> basis
      if (.not.FullMJBasis .and. IMD.ne.0) then
        write(io,'("****MDG=",I4,",  MDE=",i4)') MDG,MDE
        if (MDG.le.0 .or. MDE.le.0) then; write(io,'("****FATAL No magnetic dipole transitions calculated.")'); stop; endif
        WRITE(IO,'(/," Magnetic Dipole Transitions")')
        if (IMD.EQ.3) WRITE(IO,'(10X,"SL(J)",7X,"Energy    <L+2S> (B.M.)       <L+2S> (x10-3 Debye)")')
        if (IMD.NE.3) WRITE(IO,'(10X,"SL(J)",7X,"Energy   |<L+2S>|^2 (B.M.^2) |<L+2S>|^2 (x10-7 Debye^2)")')
        DO IG=1,MDG 
          if (wfg(ig).le.0) then; write(io,'("****FATAL wfg(ig)=",i2)') wfg(ig); stop; endif
          TwoJ=jBasis(3,wfg(ig)); LS=jBasis(4,wfg(ig))  !  2J, Place in the tables of Nielson&Koster
          if (mod(TwoJ,2).eq.0) then
            write(label,'(I4,1X,A3,"(",I2,")")') IG,TS_labels(LS,nelectrons),TwoJ/2  
          else
            write(label,'(I2,1X,A3,"(",I2,"/2)")') IG,TS_labels(LS,nelectrons),TwoJ 
          endif
          WRITE(IO,'(115("-"),/,A12,2X,F9.2," -->")') label,engs(IG)-engs(1)
          DO IE=1,MDE
!            WRITE(IO,'("IG,IE,EngJO,IntJO=",2I3,F10.2,f10.4)') IG,IE,EngJO(IG,IE),IntJO(IG,IE)
!            WRITE(IO,'("PrangeL,PrangeH=",2F10.2)') PrangeL,PrangeH
            if (PrangeL.ne.PrangeH .and. (EngJO(IG,IE).lt.PrangeL .or. EngJO(IG,IE).gt.PrangeH)) goto 105
            TwoJ=jBasis(3,wfe(ie)); LS=jBasis(4,wfe(ie))  !  2J, Place in the tables of Nielson&Koster
            if (mod(TwoJ,2).eq.0) then
              write(label,'(I5,1X,A3,"(",I2,")")') IE,TS_labels(LS,nelectrons),TwoJ/2 
            else
              write(label,'(I3,1X,A3,"(",I2,"/2)")') IE,TS_labels(LS,nelectrons),TwoJ 
            endif
            if (engJO(IG,IE).gt.1.0d-10) then ! stop counting 0->0 transition
              if (IMD.EQ.3) WRITE(IO,'(4X,A13,2X,F12.2,2f13.4)') label,EngJO(IG,IE),IntJO(IG,IE),IntJO(IG,IE)*9.274009d+0
              if (IMD.ne.3) WRITE(IO,'(4X,A13,2X,F12.2,2f13.4)') label,EngJO(IG,IE),IntJO(IG,IE)**2,   &
                                                                                   (IntJO(IG,IE)*9.274009D-03)**2*1.0d07
            endif ! engJO(IG,IE).gt.1.0d-10
 105        continue
          enddo ! IE
        enddo ! IG
        WRITE(IO,'(115("-"))')
      endif  ! .not.FullMJBasis   
!  Electric Dipole Transitions: Judd-Ofelt
      if (JuddOfelt) then
        WRITE(IO,'(/," Electric Dipole Transitions",/,115("-"))')
        if (nfit.gt.0) then 
          WRITE(IO,'(" Judd-Ofelt Intensity Fit:  Dielectric constant=",F6.4,/,  &
                     " Omega(2)=",E10.4," Omega(4)=",E10.4," Omega(6)=",E10.4," x10^-20 cm^2")') Dielectric,(JO_omega(i),i=1,3)
          write(io,'(/,"Fitted values:")')
          write(io,'("  Number   Parameter   Value  ")') 
          do i=1,nFit; write(io,'(I2,2X,I3,2X,A9,F12.4)') i,FitN(i),Plab(FitN(i)),P(FitN(i)); enddo
          IF (EDunits.EQ.1) WRITE(IO,'(23X," E(cm-1)   Eexp        DE      Ewgt     f(/10^-6))        Iexp       DI       Iwgt")') 
          IF (EDunits.EQ.2) WRITE(IO,'(23X," E(cm-1)   Eexp        DE      Ewgt  D(/10^-6 Debye^2))   Iexp       DI       Iwgt")') 
          IF (EDunits.EQ.3) WRITE(IO,'(23X," E(cm-1)   Eexp        DE      Ewgt   I(/10^-24 cm^2))    Iexp       DI       Iwgt")') 
          IF (EDunits.EQ.4) WRITE(IO,'(23X," E(cm-1)   Eexp        DE      Ewgt      e(Area))         Iexp       DI       Iwgt")') 
        else
          WRITE(IO,'(" Judd-Ofelt Intensities:  Dielectric constant=",F6.4,/,  &
                     " Omega(2)=",E10.4," Omega(4)=",E10.4," Omega(6)=",E10.4," x10^-20 cm^2")') Dielectric,(JO_omega(i),i=1,3)
          iF (EDunits.eq.3 .or. EDunits.eq.4) WRITE(IO,'(" e(Area) is the calculated integrated area of the molar extinction", &
                                                 " coefficient plotted again cm-1. Units: mol-1 L cm-1 x cm-1")') 
          if (EngInt.eq.0) then                                       
            IF (EDunits.EQ.1) WRITE(IO,'(25X," E(cm-1)  f(/10^-6)")') 
            IF (EDunits.EQ.2) WRITE(IO,'(25X," E(cm-1)  D(/10^-6 Debye^2)")') 
            IF (EDunits.EQ.3) WRITE(IO,'(25X," E(cm-1)  I(/10^-24 cm^2)")') 
            IF (EDunits.EQ.4) WRITE(IO,'(25X," E(cm-1)  e(Area)")') 
            IF (EDunits.EQ.5) WRITE(IO,'(25X," E(cm-1)  f(/10^-6)   D(/10^-6 Debye^2)  e(Area)")') 
          else
            IF (EDunits.EQ.1) WRITE(IO,'(23X," E(cm-1)   Eexp        DE      Ewgt     f(/10^-6))        Iexp       DI       Iwgt")') 
            IF (EDunits.EQ.2) WRITE(IO,'(23X," E(cm-1)   Eexp        DE      Ewgt  D(/10^-6 Debye^2))   Iexp       DI       Iwgt")') 
            IF (EDunits.EQ.3) WRITE(IO,'(23X," E(cm-1)   Eexp        DE      Ewgt   I(/10^-24 cm^2))    Iexp       DI       Iwgt")') 
            IF (EDunits.EQ.4) WRITE(IO,'(23X," E(cm-1)   Eexp        DE      Ewgt      e(Area))         Iexp       DI       Iwgt")') 
          endif
        endif  
!
        DelSqE=0.0d+00; DelSqI=0.0d+00
        DO IG=ED(1,1),ED(2,1)   
          TwoJ=jBasis(3,wfg(ig)); LS=jBasis(4,wfg(ig))  !  2J, Place in the tables of Nielson&Koster
          if (mod(TwoJ,2).eq.0) then
            write(label,'(I4,1X,A3,"(",I2,")")') IG,TS_labels(LS,nelectrons),TwoJ/2  
          else
            write(label,'(I2,1X,A3,"(",I2,"/2)")') IG,TS_labels(LS,nelectrons),TwoJ 
          endif
          WRITE(IO,'(115("-"),/,A12,2X,F9.2," -->")') label,engs(IG)-engs(1)
          SumInt=0.0d+00; 
          DO IE=ED(1,2),ED(2,2)
            if (PrangeL.ne.PrangeH .and. (EngJO(IG,IE).lt.PrangeL .or. EngJO(IG,IE).gt.PrangeH)) goto 110
            TwoJ=jBasis(3,wfe(ie)); LS=jBasis(4,wfe(ie))  !  2J, Place in the tables of Nielson&Koster
            if (mod(TwoJ,2).eq.0) then
              write(label,'(I5,1X,A3,"(",I2,")")') IE,TS_labels(LS,nelectrons),TwoJ/2 
            else
              write(label,'(I3,1X,A3,"(",I2,"/2)")') IE,TS_labels(LS,nelectrons),TwoJ 
            endif
!  D(Debye^2) = 2.127d+6 * f(oscillator strenth) / dE
            if (EDunits.EQ.1) write(CalcInt,'(f12.4)') REDJO(IG,IE,1)*1.0d+6                            ! f
            if (EDunits.EQ.2) write(CalcInt,'(F12.2)') 2.127d+12*REDJO(IG,IE,1)/EngJO(IG,IE)            ! D
            if (EDunits.EQ.3) write(CalcInt,'(F12.2)') 2.127d+12*REDJO(IG,IE,1)/EngJO(IG,IE)/2.307d-07  ! I
            if (EDunits.EQ.4) write(CalcInt,'(F12.2)') IntJO(IG,IE)                                     ! Area
            if (engJO(IG,IE).gt.1.0d-10) then
              if (EngInt.eq.0) then
                if (EDunits.EQ.1) WRITE(IO,'(4X,A13,2X,F12.2,f13.4)') label,EngJO(IG,IE),REDJO(IG,IE,1)*1.0d+6
                if (EDunits.EQ.2) WRITE(IO,'(4X,A13,2X,F12.2,f13.2)') label,EngJO(IG,IE),2.127d+12*REDJO(IG,IE,1)/EngJO(IG,IE)
                if (EDunits.EQ.3) WRITE(IO,'(4X,A13,2X,F12.2,f13.2)') label,EngJO(IG,IE),2.127d+12*REDJO(IG,IE,1)/EngJO(IG,IE) &
                                                                                                                 /2.307d-07
                if (EDunits.EQ.4) WRITE(IO,'(4X,A13,2X,F12.2,f13.2)') label,EngJO(IG,IE),IntJO(IG,IE)
                if (EDunits.EQ.5) WRITE(IO,'(4X,A13,2X,F12.2,f13.4,2f13.2)') label,EngJO(IG,IE),  & 
                                                      REDJO(IG,IE,1)*1d+6, 2.127d+12*REDJO(IG,IE,1)/EngJO(IG,IE), IntJO(IG,IE)
              endif               
              if (EngInt.ne.0) then
                ExpEng=""; DelEng=""; WgtEng=""; ExpInt=""; DelInt=""; WgtInt=""
                bit="  "; last=.false.; notLast=.true.
                do k=1,nexp
                  do k1=1,Nass(k)
                    if (NAssign(k,k1).eq.ie) then
!                  write(io,'("k=",i2,"; Nassign(k)=",i2,"; EngJO(1,ie)=",F10.2,"; Exp_E(k)=",F10.2)') &
!                              k, Nassign(k),EngJO(1,ie),Exp_E(k)
                      if (                  Nass(k).eq.1) then; bit="- "; last=.true.;  notLast=.false.; endif
                      if (k1.ne.Nass(k).and.Nass(k).ne.1) then; bit=" |"; last=.false.; notLast=.true.;  endif
                      if (k1.eq.Nass(k).and.Nass(k).ne.1) then; bit="-|"; last=.true.;  notLast=.false.; endif
                      write(ExpEng,'(F10.2)') EXP_E(k)
                      write(DelEng,'(F10.2)') EngJO(IG,IE) - EXP_E(k)
                      DelSqE=DelSqE+WGT_E(k)*(EngJO(IG,IE) - EXP_E(k))**2
                      write(WgtEng,'(F8.2)') WGT_E(k)
                      write(ExpInt,'(F12.4)') EXP_I(k)
                      if (EDunits.EQ.1) SumInt = SumInt + REDJO(IG,IE,1)*1d+6
                      if (EDunits.EQ.2) SumInt = SumInt + 2.127d+12*REDJO(IG,IE,1)/EngJO(IG,IE)
                      if (EDunits.EQ.3) SumInt = SumInt + 2.127d+12*REDJO(IG,IE,1)/EngJO(IG,IE)/2.307d-07
                      if (EDunits.EQ.4) SumInt = SumInt + IntJO(IG,IE)
                      write(WgtInt,'(F8.2)') WGT_I(k)
                      if (last) then
                        if (EDunits.eq.1) write(DelInt,'(F12.4)') SumInt - EXP_I(k)
                        if (EDunits.eq.3) write(DelInt,'(F12.4)') SumInt - EXP_I(k)
                        if (EDunits.eq.2 .or. EDunits.eq.4) write(DelInt,'(F12.4)') SumInt - EXP_I(k)
                        DelSqI=DelSqI+WGT_I(k)*(SumInt - EXP_I(k))**2
                        SumInt=0.0d+00
                      endif
                    endif
                  enddo ! k1
                enddo ! k
                if (notLast) then; ExpInt=""; DelInt=""; endif
                WRITE(IO,'(5X,A13,2X,F12.2,2A10,A8,A12,A2,2A12,A8)') label,EngJO(ig,ie),ExpEng,DelEng,WgtEng,  & 
                                                                      CalcInt,bit,ExpInt,DelInt,WgtInt
              endif ! EngInt.ne.0 
            endif ! engJO(IG,IE).gt.1.0d-10
 110        continue
          enddo ! IE
        enddo ! IG
        WRITE(IO,'(115("-"))')
        if (EngInt.ne.0) then 
          WRITE(IO,'("Sqrt(Sdel2/N)=",F10.1,5X,"Sqrt(SdelE/N)=",F10.1,18X,"Sqrt(SdelI/N)=",F10.1)')  &
                                                          SQRT((DelSqE+DelSqI)/nexp), SQRT(DelSqE/nexp), SQRT(DelSqI/nexp)
          WRITE(IO,'(115("-"))')
        endif  
        goto 500  ! return
      endif !  "J-O" 
!
!  Write the characters for symmetry operations for each wavefunction.
      IF (IGP.EQ.1) THEN
        WRITE(IO,5420) GROUP
        WRITE(IO,5421) (GROUPO(I),I=1,NCLASS)
        WRITE(IO,'(65("-"))')
        DO 300 I=1,neng
          J=n_matrix-I+1
          IF (neng.NE.n_matrix) J=neng-I+1
          if (  PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                            &
           .or. PrangeL.ne.PrangeH .and. (engs(J).gt.PrangeL .and. engs(J).lt.PrangeH) &
           .or. NrangeL.ne.PrangeH .and. (      J.ge.NrangeL .and.       J.le.NrangeH)) then
            WRITE(IO,'(I4,F10.2,1X,2A4,24(1X,2F4.1))') J,engs(J),IRREP(J,1),IRREP(J,2),(CSYM(J,K),K=1,NCLASS)
          endif
 300    CONTINUE
        WRITE(IO,'(65("-"))')
!
!  Write the fractional irreps for each wavefunction.
        WRITE(IO,5430) GROUP
        WRITE(IO,5431) (GRPR1(I),I=1,NCLASS)
        WRITE(IO,'(24X,20(A4,1X))') (GRPR2(I),I=1,NCLASS)
        WRITE(IO,'(65("-"))')
        DO 350 I=1,neng
          J=n_matrix-I+1
          IF (neng.NE.n_matrix) J=neng-I+1
          if (  PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                            &
           .or. PrangeL.ne.PrangeH .and. (engs(J).gt.PrangeL .and. engs(J).lt.PrangeH) &
           .or. NrangeL.ne.PrangeH .and. (      J.ge.NrangeL .and.       J.le.NrangeH)) then
            IF (REPS(J,1).EQ.-1.0D+00) THEN
              WRITE(IO,'(I4,F10.2,1X,2A4,24F5.1)') J,engs(J),IRREP(J,1),IRREP(J,2)
             ELSE
              WRITE(IO,'(I4,F10.2,1X,2A4,24F5.1)') J,engs(J),IRREP(J,1),IRREP(J,2),(REPS(J,K),K=1,NCLASS)
            endif  
          ENDIF
 350    CONTINUE
        WRITE(IO,'(65("-"))')
      ENDIF
!
!  Print eigenvectors. 
      IF (vectLow.ne.0) THEN
        WRITE(IO,'(/,"Eigenvectors:",i4,"-",i4)') vectLow,vectHigh
        IF (.not.MatComplex) then
          n3=(vectHigh-vectLow)/16        
          n2=vectLow-1
          if (n3.gt.0) then
            do i3=1,n3
              n1=n2+1; n2=n1-1+16
              WRITE(IO,'(15X,15(2X,I4,3X))') (J,j=n1,n2)
              WRITE(IO,'("|2S,L,2J,2MJ,i>",15F9.1)') (engs(J),j=n1,n2)
              DO i=1,n_matrix; WRITE(IO,'((5I3,15F9.5))') (FullBasis(k,i),k=1,5),(MAT(i,J),j=n1,n2); enddo  
            enddo 
          endif
          n1=n2+1; n2=vectHigh
          WRITE(IO,'(15X,15(2X,I4,3X))') (J,j=n1,n2)
          WRITE(IO,'("|2S,L,2J,2MJ,i>",15F9.1)') (engs(J),j=n1,n2)
          DO i=1,n_matrix; WRITE(IO,'((5I3,15F9.5))') (FullBasis(k,i),k=1,5),(MAT(i,J),j=n1,n2); enddo 
        else ! MatComplex
          n3=(vectHigh-vectLow)/8        
          n2=vectLow-1
          if (n3.gt.0) then
            do i3=1,n3
              n1=n2+1; n2=n1-1+8
              WRITE(IO,'(15X,15(8X,I4,8X))') (J,j=n1,n2)
              WRITE(IO,'("|2S,L,2J,2MJ,i>",15(6X,F9.1,5X))') (engs(J),j=n1,n2)
              DO i=1,n_matrix; WRITE(IO,'(5I3,15(2X,2F9.5))')  (FullBasis(k,i),k=1,5),(CMAT(K,J),j=n1,n2); enddo  
            enddo 
          endif
          n1=n2+1; n2=vectHigh
          WRITE(IO,'(15X,15(8X,I4,8X))') (J,j=n1,n2)
          WRITE(IO,'("|2S,L,2J,2MJ,i>",15(6X,F9.1,5X))') (engs(J),j=n1,n2)
          DO i=1,n_matrix; WRITE(IO,'(5I3,15(2X,2F9.5))')  (FullBasis(k,i),k=1,5),(CMAT(K,J),j=n1,n2); enddo  
        endif

        WRITE(IO,'(115("-"))')
      ENDIF
     IF (n_expect.ne.0) call Output_Expect(Neng)
!
 500  call flush(IO)  ! Allows output to be seen during multiple calculations.
      RETURN
! 
 5420 FORMAT(/,12X,'Characters of the symmetry operations of the ',A4,' point group',/,65('-'))
 5421 FORMAT(3X,'Level',17X,'Transformations <X:RX>',/,24X,24(3X,A4,2X))
 5430 FORMAT(/,12X,'Projections of representations in the ',A4,' point group',/,65('-'))
 5431 FORMAT(25X,'Irreducible representation',/,24X,24(A4,1X))
!
      END subroutine OUTPUT
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OUTPUT_EXP1()
!
!  Writes results of fitting parameters to 1-electron real f-orbital LF matrix.
!      
      IMPLICIT none
      INTEGER*4 i,j,k
      real*8 sum1,sum2,M(7,7),fv1(7),EC(7),EG(7)       
      
      if (fitType.ne.6) then; write(io,'("OUTPUT_EXP1 should not be called with fitType=",i2)') fitType; stop; endif
      k=0
      do i=1,2*Lvalue+1
        do j=1,i
          k=k+1
          M(i,j)=EXP_E(k)
          M(j,i)=M(i,j)
        enddo
      enddo  
      if (Lvalue.eq.2) then  ! d-orbital LF1e
      sum1=(EXP_E(1)+EXP_E(3)+EXP_E(6)+EXP_E(10)+EXP_E(15))/5.0d+00
      write(io,'(" The experimental 1-electron d-orbital LF MEs given in EXP1 (x indicates zero weighting) Diag(ave)=",F10.2)') sum1 
      write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>   ")') 
      write(io,'(2X," <sig|",7(F11.2,A1))')  EXP_E( 1),  WGTchar( 1)
      write(io,'(2X," <piS|",7(F11.2,A1))') (EXP_E( 1+i),WGTchar( 1+i),i=1,2)
      write(io,'(2X," <piC|",7(F11.2,A1))') (EXP_E( 3+i),WGTchar( 3+i),i=1,3)
      write(io,'(2X,"<delS|",7(F11.2,A1))') (EXP_E( 6+i),WGTchar( 6+i),i=1,4)
      write(io,'(2X,"<delC|",7(F11.2,A1))') (EXP_E(10+i),WGTchar(10+i),i=1,5)
      WRITE(IO,'(115("-"))'); 
      write(io,'(" The final calculated 1-electron d-orbital LF MEs with final parameters (diagonal adjusted)")') 
      write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>  ")') 
      write(io,'(2X," <sig|",7(F12.2))')  LF1eMat(1,1) 
      write(io,'(2X," <piS|",7(F12.2))') (LF1eMat(2,j),j=1,2)
      write(io,'(2X," <piC|",7(F12.2))') (LF1eMat(3,j),j=1,3)
      write(io,'(2X,"<delS|",7(F12.2))') (LF1eMat(4,j),j=1,4)
      write(io,'(2X,"<delC|",7(F12.2))') (LF1eMat(5,j),j=1,5)
      WRITE(IO,'(115("-"))'); 
      write(io,'(" The 1-electron d-orbital LF MEs differences (Calc-Given)*wgt (x indicates zero weighting)")') 
      write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>     ")') 
      write(io,'(2X," <sig|",7(F11.2,A1))')   (LF1eMat(1,1)-EXP_E( 1))*WGT_E( 1),  WGTchar( 1)
      write(io,'(2X," <piS|",7(F11.2,A1))') ((LF1eMat(2,i)-EXP_E( 1+i))*WGT_E( 1+i),WGTchar( 1+i),i=1,2)
      write(io,'(2X," <piC|",7(F11.2,A1))') ((LF1eMat(3,i)-EXP_E( 3+i))*WGT_E( 3+i),WGTchar( 3+i),i=1,3)
      write(io,'(2X,"<delS|",7(F11.2,A1))') ((LF1eMat(4,i)-EXP_E( 6+i))*WGT_E( 6+i),WGTchar( 6+i),i=1,4)
      write(io,'(2X,"<delC|",7(F11.2,A1))') ((LF1eMat(5,i)-EXP_E(10+i))*WGT_E(10+i),WGTchar(10+i),i=1,5)
      else if (Lvalue.eq.3) then ! f-electron LF1e 
!        write(io,'(/,"M")') 
!        do i=1,7; write(io,'(7F12.2)') (M(i,j),j=1,7); enddo
        call diagrs(M,7,7,EG,fv1,1,io)  ! calc eng only
        call zeroE(EG,7)
      sum1=(EXP_E(1)+EXP_E(3)+EXP_E(6)+EXP_E(10)+EXP_E(15)+EXP_E(21)+EXP_E(28))/7.0d+00
      write(io,'(" The experimental 1-electron f-orbital LF MEs given in EXP1 (x indicates zero weighting) Diag(ave)=",F10.2)') sum1 
      write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>      |psiS>      |psiC>     Eng(given)")') 
      write(io,'(2X," <sig|",1(F11.2,A1),6(12X),F15.2)')  EXP_E( 1),  WGTchar( 1),          EG(1)
      write(io,'(2X," <piS|",2(F11.2,A1),5(12X),F15.2)') (EXP_E( 1+i),WGTchar( 1+i),i=1,2), EG(2)
      write(io,'(2X," <piC|",3(F11.2,A1),4(12X),F15.2)') (EXP_E( 3+i),WGTchar( 3+i),i=1,3), EG(3)
      write(io,'(2X,"<delS|",4(F11.2,A1),3(12X),F15.2)') (EXP_E( 6+i),WGTchar( 6+i),i=1,4), EG(4)
      write(io,'(2X,"<delC|",5(F11.2,A1),2(12X),F15.2)') (EXP_E(10+i),WGTchar(10+i),i=1,5), EG(5)
      write(io,'(2X,"<psiS|",6(F11.2,A1),1(12X),F15.2)') (EXP_E(15+i),WGTchar(15+i),i=1,6), EG(6)
      write(io,'(2X,"<psiC|",7(F11.2,A1),       F15.2)') (EXP_E(21+i),WGTchar(21+i),i=1,7), EG(7)
      WRITE(IO,'(115("-"))'); 
!        write(io,'(/,"LF1eMat")') 
!        do i=1,7; write(io,'(7F12.2)') (LF1eMat(i,j),j=1,7); enddo
        do i=1,7; do j=1,7; M(i,j)=LF1eMat(i,j); enddo; enddo  ! 
        call diagrs(M,7,7,EC,fv1,1,io)  ! calc eng only
        call zeroE(EC,7)
      write(io,'(" The final calculated 1-electron f-orbital LF MEs with final parameters (diagonal adjusted)")') 
      write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>      |psiS>      |psiC>     Eng(calc)")') 
      write(io,'(2X," <sig|",1(F12.2),6(12X),F15.2)')  LF1eMat(1,1),        EC(1)
      write(io,'(2X," <piS|",2(F12.2),5(12X),F15.2)') (LF1eMat(2,j),j=1,2), EC(2)
      write(io,'(2X," <piC|",3(F12.2),4(12X),F15.2)') (LF1eMat(3,j),j=1,3), EC(3)
      write(io,'(2X,"<delS|",4(F12.2),3(12X),F15.2)') (LF1eMat(4,j),j=1,4), EC(4)
      write(io,'(2X,"<delC|",5(F12.2),2(12X),F15.2)') (LF1eMat(5,j),j=1,5), EC(5)
      write(io,'(2X,"<psiS|",6(F12.2),1(12X),F15.2)') (LF1eMat(6,j),j=1,6), EC(6)
      write(io,'(2X,"<psiC|",7(F12.2),       F15.2)') (LF1eMat(7,j),j=1,7), EC(7)
      WRITE(IO,'(115("-"))');                               
      write(io,'(" The 1-electron f-orbital LF MEs differences (Calc-Given)*wgt (x indicates zero weighting)")') 
      write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>      |psiS>      |psiC>     E(calc)-E(given)")') 
      write(io,'(2X," <sig|",1(F11.2,A1),6(12X),F15.2)')  (LF1eMat(1,1)-EXP_E( 1))*WGT_E( 1),    WGTchar( 1),          EC(1)-EG(1)
      write(io,'(2X," <piS|",2(F11.2,A1),5(12X),F15.2)') ((LF1eMat(2,i)-EXP_E( 1+i))*WGT_E( 1+i),WGTchar( 1+i),i=1,2), EC(2)-EG(2)
      write(io,'(2X," <piC|",3(F11.2,A1),4(12X),F15.2)') ((LF1eMat(3,i)-EXP_E( 3+i))*WGT_E( 3+i),WGTchar( 3+i),i=1,3), EC(3)-EG(3)
      write(io,'(2X,"<delS|",4(F11.2,A1),3(12X),F15.2)') ((LF1eMat(4,i)-EXP_E( 6+i))*WGT_E( 6+i),WGTchar( 6+i),i=1,4), EC(4)-EG(4)
      write(io,'(2X,"<delC|",5(F11.2,A1),2(12X),F15.2)') ((LF1eMat(5,i)-EXP_E(10+i))*WGT_E(10+i),WGTchar(10+i),i=1,5), EC(5)-EG(5)
      write(io,'(2X,"<psiS|",6(F11.2,A1),1(12X),F15.2)') ((LF1eMat(6,i)-EXP_E(15+i))*WGT_E(15+i),WGTchar(15+i),i=1,6), EC(6)-EG(6)
      write(io,'(2X,"<psiC|",7(F11.2,A1),       F15.2)') ((LF1eMat(7,i)-EXP_E(21+i))*WGT_E(21+i),WGTchar(21+i),i=1,7), EC(7)-EG(7)
      endif                           
!      
      WRITE(IO,'(115("-"))'); 
      sum1=0.0d+00; k=0
      do i=1,2*Lvalue+1 ! 5 or 7
        do j=1,i
          k=k+1
          sum1=sum1+(LF1eMat(i,j)-EXP_E(k))**2*WGT_E(k)
        enddo
      enddo
      sum1=sqrt(sum1/(dble(k)))
      WRITE(IO,'("sqrt(sum(wgt*(MEcalc-MEexp)^2))=",E12.5)') SUM1 
      WRITE(IO,'(115("-"))'); 
!
      RETURN
      END SUBROUTINE OUTPUT_EXP1
!
!234567890123456789012345678901234567890123456789012345678901234567890**
!-----------------------------------------------------------------------
!
      SUBROUTINE OUTPUT_Expect(neng)
!
!  Write the results of an EXPT calculation.
! 
      IMPLICIT none
      REAL*8 E0,EE 
      INTEGER*4 I,J,K,II,IDEG,IWRITE,JHIGH,neng,in_e
      CHARACTER LINE1*150,LINE2*120,bit*81, Format1a*40,Format1b*40 
      
      line1=" "
      line2="  Level    Rel.Eng. (deg)|"
      Do i=1,N_expect
        line2=trim(line2)//"  <"//Plab(P_expect(i))//">"
      enddo  
      write(io,'(A)') trim(line1)
      write(io,'(A)') trim(line2)
      WRITE(IO,'(115("-"))')
!      
      IDEG=1 
      write(Format1a,'("(I4,5X,F11.",i1,",5X,""|"")")') Ndecpts                    ! single energy levels
      write(format1b,'("(I4,""-"",I4,F11.",i1,",1X,""("",I2,"")|"")")') Ndecpts    ! summed over deg energy levels
!
      DO I=1,neng
        IWRITE=1
        J=n_matrix-I+1
        IF (neng.NE.n_matrix) J=neng-I+1
        IF (IDEG.EQ.1) JHIGH=J
        EE=engs(J)
        if  ( PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                    &
         .or. PrangeL.ne.PrangeH .and. (EE.gt.PrangeL .and. EE.lt.PrangeH)   & 
         .or. NrangeL.ne.NrangeH .and. ( J.ge.NrangeL .and.  J.le.NrangeH)  ) then
          IF (J.GT.1) THEN
            IF (.not.OUTP(6) .AND. ABS(EE-engs(J-1)).LT.Edegen) THEN
              IDEG=IDEG+1
              IWRITE=0
            ENDIF
          ENDIF
          IF (IWRITE.EQ.1) THEN
            IF (IDEG.EQ.1) WRITE(line1,Format1a) J,      engs(J)
            IF (IDEG.GT.1) WRITE(line1,Format1b) J,JHIGH,engs(J),IDEG
            Do in_e=1,N_expect
               write(bit,'(4X,F8.4,4X)') Expect(j,in_e)
               line1=trim(line1)//bit
            enddo
            write(io,'(A)') trim(line1)
!
            IDEG=1
          endif  ! IF (IWRITE.EQ.1) 
        ENDIF  ! if (EE.gt.PrangeL .and. EE.lt.PrangeH) or (i.gt.PrangeL .and. i.lt.PrangeH) 
      enddo ! I=1,neng
      WRITE(IO,'(115("-"))')
!
 500  call flush(IO)  ! Allows output to be seen during multiple calculations.
      RETURN
!
      END subroutine OUTPUT_Expect
!
!-----------------------------------------------------------------------
!
      SUBROUTINE PlotfileOUT(neng)
!
! Writes output to plotfile.
! 
      IMPLICIT none
      REAL*8 AMV(max_N)
      INTEGER*4 NCLASS1,I,J,IJ,K,II,TMJ,neng,IG,IE,IS, indexAMV(max_N),iblock
      character symlab*6 
!
      if (option(8)) goto 100  ! skip main calc; could be JO
!            Eng            S        spin-allowed     MJ
      IF (PLOTout1(1).or.PLOTout1(2).or.PLOTout1(3).or.PLOTout1(4)) THEN
        IF (IDG.EQ.0) THEN ! No symmetry
          if (.not.PrintCQN) then       
            II=0
            DO J=1,neng
              if (  PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                            &
               .or. PrangeL.ne.PrangeH .and. (engs(J).gt.PrangeL .and. engs(J).lt.PrangeH) &
               .or. NrangeL.ne.PrangeH .and. (      J.ge.NrangeL .and.       J.le.NrangeH)) then
                II=II+1
                AMV(II)=engs(J)
                indexAMV(II)=J 
              endif
            enddo 
            WRITE(IP1,'(2I4,"  ! Nblocks, ID")') 1,0  
            WRITE(IP1,'(I4,2X,"Eng     ! N, lab ")') II 
            WRITE(IP1,'((5E16.8))') (AMV(J),J=1,II)
            if (PLOTout1(2)) then
              WRITE(IP1,'(I4)') NSPIN
              do I=1,NSPIN   ! S_mult(I)=2S+1
                WRITE(IP1,'(2I4,"  ! I, 2S+1",/,(5E16.8))') I,S_mult(I),(SPIN(indexAMV(J),I),J=1,II)
              enddo
            endif
            if (PLOTout1(4)) then
              TMJ=TJp1(Nelectrons)-1
              if (N_odd) then
                WRITE(IP1,'(I4)') TJp1(Nelectrons)
                do IJ=1,TJp1(Nelectrons)
                  WRITE(IP1,'(I4,I3,"/2  ! I, MJ",/,(5E16.8))') IJ,-TMJ+2*(IJ-1),    &
                                                                     (binMJ(indexAMV(J),-TMJ+2*(IJ-1)),J=1,II)
                enddo
              else 
                WRITE(IP1,'(I4)') TJp1(Nelectrons)
                do IJ=1,TJp1(Nelectrons)
                  WRITE(IP1,'(I4,I5,"  ! I, MJ",/,(5E16.8))') IJ,-TMJ+1*(IJ-1),      &
                                                                     (binMJ(indexAMV(J),-TMJ+2*(IJ-1)),J=1,II)
                  enddo
              endif  ! N_odd
            endif  ! PLOTout1(4)
          else ! use Crystal Quantum numbers
            WRITE(IP1,'(I4)') nblocks
            DO Iblock=1,nblocks
              II=0
              DO J=1,neng
                IF (allCQN(J).EQ.CQN(Iblock)) then
                  if (  PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                            &
                   .or. PrangeL.ne.PrangeH .and. (engs(J).gt.PrangeL .and. engs(J).lt.PrangeH) &
                   .or. NrangeL.ne.PrangeH .and. (      J.ge.NrangeL .and.       J.le.NrangeH)) then
                    II=II+1
                    AMV(II)=engs(J)
                    indexAMV(II)=J 
                  endif
                ENDIF
              enddo ! J  
!   Filenames cannot be +1/2, etc, replace "/" by "_"      
              symlab="M"//ADJUSTL(CQN(Iblock)); i=scan(symlab,"/"); if (i.gt.0 .and. i.lt.6) symlab(i:i)="_" 
              IF (II.EQ.0) WRITE(IP1,'(I4,2X,A6)') II,symlab
              IF (II.NE.0) WRITE(IP1,'(I4,2X,A6,/,(5E16.8))') II,symlab,(AMV(J),J=1,II)
              if (PLOTout1(2)) then
                WRITE(IP1,'(I4)') NSPIN
                do I=1,NSPIN   ! S_mult(I)=2S+1
                  WRITE(IP1,'(2I4,"  ! I, 2S+1",/,(5E16.8))') I,S_mult(I),(SPIN(indexAMV(J),I),J=1,II)
                enddo
              endif
              if (PLOTout1(4)) then
                TMJ=TJp1(Nelectrons)-1
                if (N_odd) then
                  WRITE(IP1,'(I4)') TJp1(Nelectrons)
                  do IJ=1,TJp1(Nelectrons)
                    WRITE(IP1,'(I4,I3,"/2  ! I, MJ",/,(5E16.8))') IJ,-TMJ+2*(IJ-1),  &
                                                                (binMJ(indexAMV(J),-TMJ+2*(IJ-1)),J=1,II)
                  enddo
                else 
                  WRITE(IP1,'(I4)') TJp1(Nelectrons)
                  do IJ=1,TJp1(Nelectrons)
                    WRITE(IP1,'(I4,I5,"  ! I, MJ",/,(5E16.8))') IJ,-TMJ+1*(IJ-1),   &
                                                                (binMJ(indexAMV(J),-TMJ+2*(IJ-1)),J=1,II)
                    enddo
                endif  ! N_odd
              endif  ! PLOTout1(4)
            enddo  ! Iblock
          endif          
         ELSE  ! use IRREPs
          IF (ISO.EQ.0) THEN   ! spin-orbit coupling=0, Mulliken labels.
            NCLASS1=0
            DO I=1,NCLASS
              IF (GRPR2(I).NE."++++".AND.GRPR2(I).NE." +++") NCLASS1=NCLASS1+1
            enddo
            WRITE(IP1,'(3I4)') NCLASS1*NSPIN,IDG,ISO
            DO I=1,NCLASS1
              DO K=1,NSPIN
                II=0
                DO J=1,neng
                  IF (IRREP(J,1).EQ.GRPR2(I).AND.ABS(SPIN(J,K)-1.0).LT.0.01 .and.          &
                 (  PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                            &
               .or. PrangeL.ne.PrangeH .and. (engs(J).gt.PrangeL .and. engs(J).lt.PrangeH) &
               .or. NrangeL.ne.PrangeH .and. (      J.ge.NrangeL .and.       J.le.NrangeH))) THEN
                    II=II+1
                    AMV(II)=engs(J)
                    indexAMV(II)=J 
                  ENDIF
                enddo
                write(symlab,'(I2,A4)') S_mult(K),adjustl(GRPR2(I))
                WRITE(IP1,'(I4,2X,A6)') II,symlab 
                IF (II.NE.0) then
                  if (PLOTout1(1)) WRITE(IP1,'((5E16.8))') (AMV(J),J=1,II)
                  if (PLOTout1(2)) then  ! Doesn't make sense to do this, but...
                    WRITE(ip1,'(I4)') NSPIN
                    do IS=1,NSPIN   ! S_mult(I)=2S+1
                      WRITE(IP1,'(2I4,"  ! I, 2S+1",/,(5E16.8))') I,S_mult(IS),(SPIN(indexAMV(J),IS),J=1,II)
                    enddo ! IS
                  endif ! PLOTout1(2)
                  if (PLOTout1(4)) then
                    TMJ=TJp1(Nelectrons)-1
                    if (N_odd) then
                      WRITE(IP1,'(I4)') TJp1(Nelectrons)
                      do IJ=1,TJp1(Nelectrons)
                        WRITE(IP1,'(I4,I3,"/2  ! I, MJ",/,(5E16.8))') IJ,-TMJ+2*(IJ-1),  &
                                                                            (binMJ(indexAMV(J),-TMJ+2*(IJ-1)),J=1,II)
                      enddo
                    else 
                      WRITE(IP1,'(I4)') TJp1(Nelectrons)
                      do IJ=1,TJp1(Nelectrons)
                        WRITE(IP1,'(I4,I5,"  ! I, MJ",/,(5E16.8))') IJ,-TMJ+1*(IJ-1),  &
                                                                          (binMJ(indexAMV(J),-TMJ+2*(IJ-1)),J=1,II)
                        enddo
                    endif  ! N_odd
                  endif  ! PLOTout1(4)
                endif 
              enddo  ! K=1,NSPIN
            enddo ! I=1,NCLASS1              
          ELSE IF (ISO.EQ.1) THEN  ! s.o. #0; Bethe notation
            WRITE(IP1,'(3I4)') NCLASS,IDG,ISO
            DO I=1,NCLASS
              II=0
              DO J=1,neng
                IF (IRREP(J,1).EQ.GRPR1(I).and.                                            &
                 (  PrangeL.eq.PrangeH .and. NrangeL.eq.NrangeH                            &
               .or. PrangeL.ne.PrangeH .and. (engs(J).gt.PrangeL .and. engs(J).lt.PrangeH) &
               .or. NrangeL.ne.PrangeH .and. (      J.ge.NrangeL .and.       J.le.NrangeH))) then
                  II=II+1
                  AMV(II)=engs(J)
                  indexAMV(II)=J 
                ENDIF
              enddo ! J  
              WRITE(IP1,'(I4,4X,A4)') II,GRPR1(I) 
              IF (II.NE.0) then
                if (PLOTout1(1)) WRITE(IP1,'((5E16.8))') (AMV(J),J=1,II)
                if (PLOTout1(2)) then
                  WRITE(IP1,'(I4)') NSPIN
                  do IS=1,NSPIN   ! S_mult(I)=2S+1
                    WRITE(IP1,'(2I4,"  ! I, 2S+1",/,(5E16.8))') I,S_mult(IS),(SPIN(indexAMV(J),IS),J=1,II)
                  enddo ! IS
                endif ! PLOTout1(2)
                if (PLOTout1(3)) then
                  write(io,'(2I4,"  ! SpAllow(1), SpAllow(2)",/,(5E16.8))') SpAllow(1),SpAllow(2),(SpinAllowed(j), j=1,n_matrix)
                endif  ! PLOTout1(3)
                if (PLOTout1(4)) then
                  TMJ=TJp1(Nelectrons)-1
                  if (N_odd) then
                    WRITE(IP1,'(I4)') TJp1(Nelectrons)
                    do IJ=1,TJp1(Nelectrons)
                      WRITE(IP1,'(I4,I3,"/2  ! I, MJ",/,(5E16.8))') IJ,-TMJ+2*(IJ-1),  &
                                                                          (binMJ(indexAMV(J),-TMJ+2*(IJ-1)),J=1,II)
                    enddo
                  else 
                    WRITE(IP1,'(I4)') TJp1(Nelectrons)
                    do IJ=1,TJp1(Nelectrons)
                      WRITE(IP1,'(I4,I5,"  ! I, MJ",/,(5E16.8))') IJ,-TMJ+1*(IJ-1),  &
                                                                        (binMJ(indexAMV(J),-TMJ+2*(IJ-1)),J=1,II)
                    enddo
                  endif  ! N_odd
                endif  ! PLOTout1(4)
              endif ! II.ne.0   
            enddo  ! I
          ENDIF  ! ISO
        ENDIF  !
      ENDIF   !  PLOTout1(1).or.PLOTout1(2).or.PLOTout1(3).or.PLOTout1(4) = Eng, S, Spin-allowed, MJ
!
!  Magnetic Dipole Transitions
      if (PLOTout1(5) .and. IMD.NE.0) then
        WRITE(IP1,'(2I4)') MDG,MDE
        do I=1,MDG
          WRITE(IP1,'(1X,F10.2," (",I2,") -->")') RENGG(I),IDEGG(I)
          DO J=1,MDE
            IF(IMD.EQ.1) WRITE(IP1,'(F12.2," (",I2,")",7F10.6)') RENGE(J),IDEGE(J),(RMD(I,J,K),K=1,3)
            IF(IMD.EQ.2) WRITE(IP1,'(F12.2," (",I2,")",7F10.6)') RENGE(J),IDEGE(J),(RMD(I,J,K),K=1,5), &
                                                                 RMD(I,J,5)-RMD(I,J,4),RMD(I,J,1)-RMD(I,J,2)
            IF(IMD.EQ.3) WRITE(IP1,'(F12.2," (",I2,")",3(2X,2F8.5))')RENGE(J),IDEGE(J),(CMD(I,J,K),K=1,3)
          enddo ! J
        enddo ! I
      endif ! PLOTout1(5)
!
!  Electric Dipole Transitions: Judd-Ofelt
 100  continue
      if (PLOTout1(6) .and. JuddOfelt) then
        WRITE(IP1,'(2I4)') ED(2,1)-ED(1,1),ED(2,2)-ED(1,2)
        DO IG=ED(1,1),ED(2,1)   
          WRITE(IP1,'(1X,F10.2," -->")') engs(IG)-engs(1)
          DO IE=ED(1,2),ED(2,2)
            if (PrangeL.ne.PrangeH .and. (EngJO(IG,IE).lt.PrangeL .or. EngJO(IG,IE).gt.PrangeH)) goto 110
!  D(Debye^2) = 2.127d+6 * f(oscillator strenth) / dE
            if (engJO(IG,IE).gt.1.0d-10) then
              WRITE(IP1,'(F12.2,5x,3E12.6)') EngJO(ig,ie),REDJO(IG,IE,1)*1d+6,2.127d+12*REDJO(IG,IE,1)/EngJO(ig,ie),IntJO(ig,ie)
            endif
 110        continue
          enddo ! IE
        enddo ! IG
        goto 500  ! return
      endif !  PLOTout1(6) & "J-O" 
!
!  Write fit parameters
     ! if (PLOTout1(7)) ! this done in subroutine fcn(m,n,x,fvec,iflag) of FIT
!
!  Write Bkq parameters
!                 1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16   
!DATA BkqR_label/"B00 ","B20 ","B21 ","B22 ","B40 ","B41 ","B42 ","B43 ","B44 ","B60 ","B61 ","B62 ","B63 ","B64 ","B65 ","B66 "/
!DATA BkqI_label/"****","****","B21'","B22'","****","B41'","B42'","B43'","B44'","****","B61'","B62'","B63'","B64'","B65'","B66'"/
      IF (PLOTout1(8)) THEN
        if (Lvalue.eq.2) then
          write(ip1,'(14F12.5)') (BKQR(i),i=2,9),BKQI(3),BKQI(4),BKQI(6),BKQI(7),BKQI(8),BKQI(9)
        else
          write(ip1,'(15F12.5)') (BKQR(i),i=2,16)
          write(ip1,'(12F12.5)') BKQI(3),BKQI(4),BKQI(6),BKQI(7),BKQI(8),BKQI(9),(BKQI(i),i=11,16)
        endif
      endif
!
!  write g-values
      if (PLOTout1(9) .and. GEXS(1).ne.0) then
        do i=1,ngexs  ! (GEXS(2)-GEXS(1)+1)/2
          write(ip1,'(i3,"-",i3,F9.2,4X,3F8.4,3(2X,3F7.3))') gexs(2*i-1),gexs(2*i),engs(gexs(2*i-1)),   &
                                           (gval(j,i),j=1,3),(g_prin(1,j,i),g_prin(2,j,i),g_prin(3,j,i),j=1,3)
        enddo ! ! i
      endif
!
!  write d/f-orbital energies
      if (PLOTout1(10)) then
        write(ip1,'(7E16.8)') (OrbEngs(i),i=1,2*Lvalue+1)
      endif
!
 500  call flush(ip1)  ! Allows output to be seen during multiple calculations.
      RETURN
! 
      END subroutine PlotfileOUT
!
!-----------------------------------------------------------------------
!
      subroutine TO_LOWER(KEY)
      IMPLICIT none
      integer i,j,KEYlen
      character(len=*), intent(inout)  :: KEY
      CHARACTER lower*26, upper*26
      PARAMETER(lower="abcdefghijklmnopqrstuvwxyz")
      PARAMETER(upper="ABCDEFGHIJKLMNOPQRSTUVWXYZ")
!      
      KEYlen = len_trim(KEY)
      DO I=1,KEYlen
        DO J=1,26
          IF (key(i:i).EQ.UPPER(j:j)) then
            KEY(I:I)=lower(J:J)
            goto 20
          endif
        enddo  
 20     continue
      enddo
! 
      return
      end subroutine TO_LOWER
!
!-----------------------------------------------------------------------
!
      subroutine TO_UPPER(KEY)
      IMPLICIT none
      integer i,j,KEYlen
      character(len=*), intent(inout)  :: KEY
      CHARACTER lower*26, upper*26
      PARAMETER(lower="abcdefghijklmnopqrstuvwxyz")
      PARAMETER(upper="ABCDEFGHIJKLMNOPQRSTUVWXYZ")
!      
      KEYlen = len_trim(KEY)
      DO I=1,KEYlen
        DO J=1,26
          IF (key(i:i).EQ.lower(j:j)) then
            KEY(I:I)=UPPER(J:J)
            goto 20
          endif
        enddo  
 20     continue
      enddo
! 
      return
      end subroutine TO_UPPER
!
!-----------------------------------------------------------------------
END MODULE f_e_input
!    3240