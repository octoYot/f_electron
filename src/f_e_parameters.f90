MODULE f_e_parameters

USE f_e_data
USE f_e_interpreter

IMPLICIT NONE

PUBLIC 

character*4 LFtype/"    "/

!  AOM stuff
INTEGER, PARAMETER :: max_ligands=20
integer NLIGANDS
! The only way for AOMchanged to be set to true is through subroutines setpEq() or makeLinks()
! flags that a call to AOMmatrixD() or  AOMmatrixF() is required at which AOMchanged set back to .false.
LOGICAL AOMchanged/.false./, AOM_Angles/.false./
character*4 Lname(max_ligands)
real*8 eSigma(max_ligands),ePiX(max_ligands),ePiY(max_ligands),Theta(max_ligands),Phi(max_ligands),Chi(max_ligands)
real*8 xLig(max_ligands),yLig(max_ligands),zLig(max_ligands)
! Extended AOM
Logical AOMX/.false./
real*8 eDelS(max_ligands),eDelC(max_ligands),ePsiS(max_ligands),ePsiC(max_ligands)

!  Extended Stevens Operators
INTEGER*4 ESO2J/0/,ESOkmax/0/,ESONorm/0/,ESOskip/0/   !  Extended Stevens Operators
REAL*8 ESOBkq(20,41)
complex*16 ESOcmat(18,18)

!  Crystal field stuff
REAL*8 BKQR(16),BKQI(16),savedBKQR(16),savedBKQI(16),sumCLF,RotLF(3)/3*0.0d0/
integer Altp_lookup(3,-1:1,0:7)
character*4 Para_label(20),BkqR_label(16),BkqI_label(16),Altp_label(45)
DATA Para_label/"EAVE","F2  ","F4  ","F6  ","ZETA","ALPH","BETA","GAMM","M0  ","M2  ","M4  ","P2  ","P4  ","P6  ","T2  ","T3  ",  &
                "T4  ","T6  ","T7  ","T8  "/
!                1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16   
DATA BkqR_label/"B00 ","B20 ","B21 ","B22 ","B40 ","B41 ","B42 ","B43 ","B44 ","B60 ","B61 ","B62 ","B63 ","B64 ","B65 ","B66 "/
DATA BkqI_label/"****","****","B21'","B22'","****","B41'","B42'","B43'","B44'","****","B61'","B62'","B63'","B64'","B65'","B66'"/

!  CCF labels
character CCF_label(12)*3
DATA CCF_label/"1  ","2  ","3  ","4  ","5  ","6  ","7  ","8  ","9  ","10A","10B","11 "/

!  Intensity parameters
DATA Altp_label/"A210","A211","A220","A221","A222","A230","A231","A232","A233","A430","A431","A432","A433","A440","A441","A442", &
                "A443","A444","A450","A451","A452","A453","A454","A455","A650","A651","A652","A653","A654","A655","A660","A661", &
                "A662","A663","A664","A665","A666","A670","A671","A672","A673","A674","A675","A676","A677"/
Data Altp_lookup/ 1,10,25,  3,14,31,  6,19,38,  &
                  2,11,26,  4,15,32,  7,20,39,  &
                 -1,12,27,  5,16,33,  8,21,40,  &
                 -1,13,28, -1,17,34,  9,22,41,  &
                 -1,-1,29, -1,18,35, -1,23,42,  &
                 -1,-1,30, -1,-1,36, -1,24,43,  &
                 -1,-1,-1, -1,-1,37, -1,-1,44,  &
                 -1,-1,-1, -1,-1,-1, -1,-1,45/
REAL*8 Tli(max_ligands,9),Cli(max_ligands,7)

! one electron LF matrix elements
REAL*8 LF1emat(7,7)

! Atomic parameters
REAL*8, PARAMETER :: Z=0.0d+00
REAL*8 EAVE/z/,F2/z/,F4/z/,F6/z/,Bracah/z/,Cracah/z/,ALPHA/z/,BETA/z/,GAMMA/z/,ZETA/z/
REAL*8 M0/z/,M2/z/,M4/z/, P2/z/,P4/z/,P6/z/, T2/z/,T3/z/,T4/z/,T6/z/,T7/z/,T8/z/
LOGICAL Racah/.false./, xref/.false./, doAtomic/.false./
REAL*8 Atomic(20)/20*z/

! Expectation values
INTEGER N_expect/0/,P_expect(20)
Real*8 Expect(max_N,20)

INTEGER, PARAMETER :: maxTotPar=450
logical FitOpt(10)/10*.false./
character*1 V1(maxTotPar)/maxTotPar*" "/, V2(maxTotPar)/maxTotPar*" "/
real*8 P(maxTotPar), rlink(2,50)

integer nlink(2,50), nlinks/0/, nRef(3)/2,1,3/, FitType/0/

! Offset parameters (set as parameters 350-400):
integer offsetType,offsetN, iOffset(49,2)  ! offsetType=1/2 for diagonalised |LSJ> multiplets/full diagonalised basis.
real*8 ROffset(49)

! Fitting data:
INTEGER, PARAMETER :: maxFit=30   !  maximum number of parameters that can be varied in a fit at one time.
INTEGER, PARAMETER :: maxNF= 630000  !  maximum size of FitMatR,FitMatI = maxFit*maxNF  ! f6 B61 requires 205290 for maxNF, 
                                     !  but AOM parameters are a lot more, since they are linear combinations of Bkqs.
real*8 FitP(maxFit),FitL(maxTotPar),FitH(maxTotPar),cov(maxFit,maxFit),FitTol,fnorm_init
!REAL*8 ARfit(6,(max_N+24)/6,(max_N+24)/6),AIfit(6,(max_N+24)/6,(max_N+24)/6)
REAL*8, dimension(:,:,:), allocatable :: ARfit,AIfit,FitMatR,FitMatI
integer, dimension(:,:,:), allocatable :: iIndexR,jIndexR,iIndexI,jIndexI  ! indices for FitMatR,FitMatI MEs
Integer FitN(maxFit),NFit/0/  ! FitN(Nfit) contains the NFit parameters "standard order" numbers 
INTEGER Its/0/,maxIts,nMatFR(maxFit,6),nMatFI(maxFit,6)  
logical nMatC(maxFit) 
logical linearFit/.true./ ! true if the parameters being varied can be factored outside a matrix. 
!  This will be true for the atomic,CF,Esig/Epi parameters, but not for the Ligand angular positions. 
!  It is also possible to assign a parameter to be a nonlinear function of a variable.

! Constants: Evaluated once in the program
INTEGER, PARAMETER :: maxCon=30       ! maximum number of Constants
INTEGER      NCon/0/
character*10 Constants(maxCon)
character*100      cEq(maxCon) 
real*8       ConstantV(maxCon)

! Variables used by the evaluator to build and to determine the value of the expressions
INTEGER, PARAMETER :: maxTotVar=30       ! total number of variables
INTEGER, PARAMETER :: maxTotVarVals=501  ! total number of values each variable can have
logical grid
INTEGER      Nvar/0/, nVarVals(maxTotVar)/maxTotVar*0/, varNum(2)/1,2/
character*10 variables(maxTotVar)
real*8       variablesvalues(maxTotVar),allVarVals(maxTotVar,maxTotVarVals)
character*100 pEq(maxTotPar)     !  Equation for the parameters. (Can be a number or zero).
logical       pEqF(maxTotPar)    !  The parameter is an formula that cannot be immediately evaluated to 
                                 !  a number, but is a formula of a variable.
logical       pEqFF(maxTotPar)   !  The parameter is an formula of a variable that is to be fitted.
LOGICAL       atLeastOneVarForm 

CONTAINS
!
!-----------------------------------------------------------------------
!  Determines which parameter/variables being varied in a fit are linear 
!  of non-linear functions of a matrix.
!  Atomic, CF & AOM parameters linear. AOM angles non-linear.
!  However a linear parameter can be a non-linear function of a variable.
!     VARI  x
!     L1    x^2 ePi theta phi chi
!  makes the eSigma parameter of L1 equal to x^2
      SUBROUTINE parLinearity()
      IMPLICIT none
      integer i,j
      real*8 S1,P1,P2,P3
      character func*100, statusflag*5
!
      linearFit=.true.
      do i=1,maxTotPar
        if (V1(i).eq."*") then
          V2(i)="l"
! CF matrices are non-linear fucntions of the AOM angles. 
          if (LFtype.eq."AOM " .and. i.gt.80  .and. i.le.140) V2(i)="n"  
          if (AOMX             .and. i.gt.140 .and. i.le.180) V2(i)="l"   ! AOMX  
          if (i.gt.400) then  ! this is a variable being fitted.

            do j=1,399  ! Find the parameters that are a function of this variable.
              if (pEqF(j)) then
                if (LFtype.eq."AOM " .and. j.gt.80 .and. j.le.140) then    ! AOM angle
                  V2(i)="n" 
                  linearFit=.false.
                else
                  func=pEq(j)
                  call init(func, variables, statusflag)
                  if(statusflag == "ok") Then
                    S1 = variablesvalues(i-400)
                    P2 = evaluate(variablesvalues)
                    variablesvalues(i-400)=S1-1.0d+00 
                    P1 = evaluate(variablesvalues)
                    variablesvalues(i-400)=S1+1.0d+00 
                    P3 = evaluate(variablesvalues)
                    variablesvalues(i-400)=S1
                    if (abs((P2-P1)-(P3-P2)).gt.0.00001) then
                      V2(i)="n" 
                      linearFit=.false.
!                      write(IO,'("Parameter:",I3,"(",A9,") is a non-linear function of the variable ",A10)') j,Plab(j), &
!                                                                                                             variables(i-400)
!                      write(IO,'("non-linearity: value - (value-1)=",F10.4,"; (value+1)-value=",F10.4)') P2-P1,P3-P2
                    else
!                      write(IO,'("Parameter:",I3," is a linear function of the variable ",A10)') j,variables(i-400)
!                      write(IO,'("linearity: value - (value-1)=",F10.4,"; (value+1)-value=",F10.4)') P2-P1,P3-P2
                    endif
                  else
                    write(IO,'("***FATAL: Parameter:",i3," ",A9," contains unknown variables: ",A)') i,pLab(i),trim(func)
                    stop              
                  endif  
                  call destroyfunc()
                endif ! LFtype.eq."AOM " .and. j.gt.80 .and. j.le.140
              endif  ! pEqF(j)           
            enddo ! j=1,300

          endif ! i.gt.400
        endif
      enddo
!
      RETURN
      END SUBROUTINE parLinearity
!
!-----------------------------------------------------------------------
!
      SUBROUTINE setpEq(ipar)
!
!  Sets the parameters from a given set of variablevalues.
!  Assumes that matrix P() is up to date.
!  Must load variablevalues from the latest values in P(400-450), which may have changed.
!  If ipar=0; set all variable values.
!  If ipar>0 & <=400; then all variable values set to zero.
!  If ipar>400 & <451; then all variable values set to zero except for variable ipar. 
!
!  Note: 
!     ipar will be non-zero only when called from PrepCalc()<prepareFit()
!     ipar will never refer to a constant. Constants always stay the same
!
      IMPLICIT none
      integer i,n,ipar
      character variableString*200, vv*12, statusflag*5, func*100
      
      if (.not.atLeastOneVarForm) return ! There much be at least one parameter dependent on a variable
      if (ipar.eq.0) then 
        do i=1,nVar; variablesvalues(i)=P(400+i); enddo ! update variablesvalues from P().
      else if (ipar.gt.0 .and. ipar.le.400) then 
        do i=1,nVar; variablesvalues(i)=0.0d+00; enddo ! set all variablesvalues to 0. 
      else if (ipar.gt.400 .and. ipar.le.450) then 
        do i=1,nVar; variablesvalues(i)=0.0d+00; enddo; ! set all variablesvalues to 0, except for ipar. Called from FIT
       variablesvalues(ipar-400)=P(ipar); ! (ipar should never refer to a constnt)
      endif
!
      do i=1,maxTotPar
        if (pEqF(i)) then
!          write(io,'("Equation for parameter:",i3,2X,A9,"=",A100)') i,pLab(i),pEq(i)
          if (LFtype.eq."AOM " .and. i.ge.21 .and. i.le.180) AOMchanged=.true.  ! AOMX changed 140 to 180
          func=pEq(i)
          call init (func, variables, statusflag)
          if(statusflag == "ok") Then
            P(i) = evaluate(variablesvalues)   ! variablesvalues(i) 
           else
            write(io,'("***FATAL: Error is setting equation for parameter",A9,"; statusflag=",A5)') pLab(i),statusflag
          endif            
!          write(IO,'("Parameter:",i2,", pEqF(i)=",L1,",  P(i)=",F12.5,", func=",A100)') i,pEqF(i),P(i),func
          call destroyfunc()
        endif
      enddo
!
!      write(*,'("Return from  setpEq;  P(121)=",F12.2)') P(121)
      RETURN
      END SUBROUTINE setpEq
!
!-----------------------------------------------------------------------
!
      SUBROUTINE checkLFcomplex()
!  Checks to see if there are any imaginary terms in the ligand field.
!  Assumes that matrix P() is up to date.
!  Assumes MFcomplex is up to date.
!
      IMPLICIT none
      integer i
      logical fC/.true./  ! firstCall
      save fC
!      
      sumCLF=0.0d+00
      do i=1,16; if (abs(BKQI(i)).lt.CFbit) BKQI(i)=0.0d+00; sumCLF=sumCLF+abs(BKQI(i)); enddo
      if (sumCLF.ge.CFbit) then
        if (.not.CFcomplex) then
          if (.not.fC) write(io,'("***Warning: The ligand field has been changed from real to complex in checkLFcomplex().",/)') 
          CFcomplex=.true.
          MatComplex=CFcomplex.or.MFcomplex
        endif 
      else ! sumCLF < CFbit
        if (CFcomplex) then
          if (.not.fC) write(io,'("***Warning: The ligand field has been changed from complex to real in checkLFcomplex().")') 
!          write(io,'("***BKQI(15)=",F12.6," BKQI(16)=",F12.6," sumCLF=",F12.6,/ )') BKQI(15),BKQI(16),sumCLF
          CFcomplex=.false.
          MatComplex=CFcomplex.or.MFcomplex
        endif
      endif  
      fC=.false. 
!
      RETURN
      END SUBROUTINE checkLFcomplex
!
!-----------------------------------------------------------------------
!
!  Assumes that matrix P() is up to date.
!
      SUBROUTINE makeLinks()
      IMPLICIT none
      integer i
!      
      if (nlinks.eq.0) return
      do i=1,nlinks
        if (nlink(1,i).gt.20 .and. nlink(1,i).le.140 .and. LFtype.eq."AOM ") AOMchanged=.true.
        p(nlink(1,i)) = p(nlink(2,i))*rlink(1,i) + rlink(2,i)
      enddo
!
      RETURN
      END SUBROUTINE makeLinks
!
!-----------------------------------------------------------------------
!
      SUBROUTINE loadP()
! Loads the parameters into the array P
      IMPLICIT none
      integer i,j,k,t,icount,lamda
!      
      P(1)=EAVE
      if (Racah) then
        P(2)=Bracah;P(3)=Cracah;P(4)=0.0d+00
      else
        P(2)=F2;    P(3)=F4;    P(4)=F6
      endif
      P(5)=ZETA;   P(6)=ALPHA  
      P(7)=BETA;  P(8)=GAMMA;P(9)=M0;   P(10)=M2;  P(11)=M4;    P(12)=P2
      P(13)=P4;   P(14)=P6;  P(15)=T2;  P(16)=T3;  P(17)=T4;    P(18)=T6
      P(19)=T7;   P(20)=T8
      if (LFtype.eq."CF  ") then
        P(21)=BKQR(2); P(22)=BKQR(3); P(23)=BKQR(4); P(24)=BKQR(5); P(25)=BKQR(6)
        P(26)=BKQR(7); P(27)=BKQR(8); P(28)=BKQR(9); P(29)=BKQR(10);P(30)=BKQR(11)
        P(31)=BKQR(12);P(32)=BKQR(13);P(33)=BKQR(14);P(34)=BKQR(15);P(35)=BKQR(16)
        P(42)=BKQI(3); P(43)=BKQI(4); P(45)=BKQI(6); P(46)=BKQI(7); P(47)=BKQI(8)
        P(48)=BKQI(9); P(50)=BKQI(11);P(51)=BKQI(12);P(52)=BKQI(13);P(53)=BKQI(14)
        P(54)=BKQI(15);P(55)=BKQI(16)
      else if (LFtype.eq."AOM ") then
        do i=1,NLIGANDS
          P(20+i)=eSigma(i); P(40+i)=ePiX(i); P(60+i)=ePiY(i); P(80+i)=Theta(i); P(100+i)=Phi(i); P(120+i)=Chi(i)
        enddo
        if (AOMX) then
          do i=1,NLIGANDS
            P(140+i)=eDelS(i); P(150+i)=eDelC(i); P(160+i)=ePsiS(i); P(170+i)=ePsiC(i)
          enddo
        endif
      else if (LFtype.eq."LF1E") then
        k=0
        do i=1,2*Lvalue+1
          do j=1,i
            k=k+1
            P(20+k)=LF1eMat(i,j) 
          enddo
        enddo
      endif  ! LFtype
      if (CCFtype.ne.0) then
        do i=1,N_CCF
          P(140+i)=CCF(i)
        enddo
      endif
!      
      if (ED_IntMech.eq."Full") then
        if (IntN1.eq.1) then
          icount=0
          do i=1,3
            lamda=2*i
            do j=1,3
              t=lamda-2+j  !  1,2,3 or 3,4,5 or 5,6,7
              do k=1,t+1
                icount=icount+1
                P(200+icount*2-1)=AltpR(i,j,k)
                P(200+icount*2  )=AltpI(i,j,k)
              enddo
            enddo
          enddo
          if (icount.ne.45) then; write(io,'("Something wrong in loadP; icount=",i4)') icount; Stop; endIF 
        else if (IntN1.eq.2) then
          do i=1,15
            do j=1,3
              P(200+(i-1)*6+j*2-1)=BlkiR(i,j)
              P(200+(i-1)*6+j*2  )=BlkiI(i,j)
            enddo  ! j
          enddo ! i
        else if (IntN1.eq.3) then
          do i=1,NLIGANDS
            do j=1,9; P(200+(i-1)*9+j)=Tli(i,j);  enddo  ! j
          enddo ! i
        else if (IntN1.eq.4) then
          do i=1,NLIGANDS
            do j=1,7; P(200+(i-1)*7+j)=Cli(i,j);  enddo  ! j
          enddo ! i
        endif ! IntN1  
      endif
      if (JuddOfelt) then
        do i=1,3; P(300+i)=JO_omega(i); enddo
      endif
      P(304)=MAGF(1);  P(305)=MAGF(2);  P(306)=MAGF(3)
      P(307)=ROTLF(1); P(308)=ROTLF(2); P(309)=ROTLF(3)
      P(310)=RK(1);    P(311)=RK(2);    P(312)=RK(3)
      if (offsetN.gt.0) then
        do i=1,offsetN; P(350+i)=ROffset(i); enddo
      endif

      if (nVar.gt.0) then
        do i=1,nVar
          P(400+i)=variablesvalues(i)
        enddo
      endif 
      if (NCon.gt.0) then
        do i=1,NCon
          P(400+nVar+i)=ConstantV(i)
        enddo
      endif
!
      RETURN
      END SUBROUTINE loadP
!
!-----------------------------------------------------------------------
!
      SUBROUTINE unloadP1()
! Defines atomic parameters for a PRED calculation where the atomic parameters are different.
      IMPLICIT none
      integer i,j,k,t,icount,lamda
!      
      EAVE=Atomic(1);   F2=Atomic(2);   F4=Atomic(3);   F6=Atomic(4);   
      ZETA=Atomic(5);ALPHA=Atomic(6); BETA=Atomic(7);GAMMA=Atomic(8)
      M0=  Atomic(9);   M2=Atomic(10);  M4=Atomic(11);  P2=Atomic(12);  P4=Atomic(13);  P6=Atomic(14)  
      T2=  Atomic(15);  T3=Atomic(16);  T4=Atomic(17);  T6=Atomic(18)
      T7=  Atomic(19);  T8=Atomic(20)
!
      RETURN
      END SUBROUTINE unloadP1
!
!-----------------------------------------------------------------------
!
      SUBROUTINE unloadP()
! Unloads the parameters from the array P
!
! An actual calculation uses the named parameters in BuildMatrix
! In a fit the parameters in P are changed, then unloadP() called.
!
      IMPLICIT none
      integer i,j,k,t,icount,lamda
      character func*100, statusflag*5
      real*8 dtr
      dtr=4.0d+00*atan(1.0d+00)/180.0d+00
 !     write(*, '("P(I)=")') 
 !     write(*, '((10F8.2))') (P(i),i=1,maxTotPar) 
 !     write(io,'("P(I)=")') 
 !     write(io,'((10F8.2))') (P(i),i=1,maxTotPar) 
            
      EAVE=P(1)
      if (Racah) then
        Bracah=P(2); Cracah=P(3)
      else
        F2=P(2);   F4=P(3);   F6=P(4);   
      endif
      ZETA=P(5);ALPHA=P(6)  
      BETA=P(7);GAMMA=P(8);  M0=P(9);   M2=P(10);  M4=P(11);    P2=P(12)
      P4=P(13);   P6=P(14);  T2=P(15);  T3=P(16);  T4=P(17);    T6=P(18)
      T7=P(19);   T8=P(20)
      if (LFtype.eq."CF  ") then
        BKQR(2) =P(21); BKQR(3) =P(22); BKQR(4) =P(23); BKQR(5) =P(24); BKQR(6)=P(25)
        BKQR(7) =P(26); BKQR(8) =P(27); BKQR(9) =P(28); BKQR(10)=P(29); BKQR(11)=P(30)
        BKQR(12)=P(31); BKQR(13)=P(32); BKQR(14)=P(33); BKQR(15)=P(34); BKQR(16)=P(35)
        BKQI(3) =P(42); BKQI(4) =P(43); BKQI(6) =P(45); BKQI(7) =P(46); BKQI(8) =P(47)
        BKQI(9) =P(48); BKQI(11)=P(50); BKQI(12)=P(51); BKQI(13)=P(52); BKQI(14)=P(53)
        BKQI(15)=P(54); BKQI(16)=P(55)                                      
      else if (LFtype.eq."AOM ") then
        do i=1,NLIGANDS
          if (eSigma(i).ne.P(20+i) .or. ePiX(i).ne.P(40+i) .or. ePiY(i).ne.P(60+i).or. Theta(i).ne.P(80+i) .or.  &
                Phi(i).ne.P(100+i)  .or. Chi(i).ne.P(120+i) ) AOMchanged=.true.
          eSigma(i)=P(20+i); ePiX(i)=P(40+i); ePiY(i)=P(60+i)
          if (AOM_angles) then
            Theta(i)=P(80+i); Phi(i)=P(100+i); Chi(i)=P(120+i)
            xLig(i)=sin(theta(i)*dtr)*cos(phi(i)*dtr)
            yLig(i)=sin(theta(i)*dtr)*sin(phi(i)*dtr) 
            zLig(i)=cos(theta(i)*dtr) 
          else
            xLig(i)=P(80+i); yLig(i)=P(100+i); zLig(i)=P(120+i)
            phi(i)=atan2(yLig(i),xLig(i))/dtr
            theta(i)=atan2(sqrt(xLig(i)**2+yLig(i)**2),zLig(i))/dtr
          endif
        enddo
        if (AOMX) then
          do i=1,NLIGANDS
            eDelS(i)=P(140+i); eDelC(i)=P(150+i); ePsiS(i)=P(160+i); ePsiC(i)=P(170+i)
          enddo
        endif
      else if (LFtype.eq."LF1E") then
        k=0
        do i=1,2*Lvalue+1
          do j=1,i
            k=k+1
            LF1eMat(i,j) =  P(20+k)
          enddo
        enddo
      endif  !  LFtype
      if (CCFtype.ne.0) then
        do i=1,N_CCF
          CCF(i)=P(140+i)
        enddo
      endif
!      
      if (ED_IntMech.eq."Full") then
        if (IntN1.eq.1) then     !  Altp 
          icount=0
          do i=1,3
            lamda=2*i
            do j=1,3
              t=lamda-2+j  !  1,2,3 or 3,4,5 or 5,6,7
              do k=1,t+1
                icount=icount+1
                AltpR(i,j,k) = P(200+icount*2-1)
                AltpI(i,j,k) = P(200+icount*2  )
              enddo
            enddo
          enddo
          if (icount.ne.45) then; write(io,'("Something wrong in unloadP; icount=",i4)') icount; Stop; endIF 
        elseif (IntN1.eq.2) then    !  Blki
          do i=1,15
            do j=1,3
              BlkiR(i,j) = P(200+(i-1)*6+j*2-1)
              BlkiI(i,j) = P(200+(i-1)*6+j*2  )
            enddo
          enddo
        else if (IntN1.eq.3) then
          do i=1,NLIGANDS
            do j=1,9; Tli(i,j)=P(200+(i-1)*9+j);  enddo  ! j
          enddo ! i
        else if (IntN1.eq.4) then
          do i=1,NLIGANDS
            do j=1,7; Cli(i,j)=P(200+(i-1)*7+j);  enddo  ! j
          enddo ! i
        endif ! IntN1  
      endif
      if (JuddOfelt) then
        do i=1,3; JO_omega(i)=P(300+i); enddo
      endif
      MAGF(1)=P(304);  MAGF(2)=P(305);  MAGF(3)=P(306)
      ROTLF(1)=P(307); ROTLF(2)=P(308); ROTLF(3)=P(309)
      RK(1)=P(310);    RK(2)=P(311);    RK(3)=P(312)
      if (offsetN.gt.0) then
        do i=1,offsetN; ROffset(i)=P(350+i); enddo
      endif
      if (nVar.gt.0) then
        do i=1,nVar
          variablesvalues(i)=P(400+i)
        enddo
      endif
      if (nCon.gt.0) then
        do i=1,nCon
          ConstantV(i)=P(400+nVar+i)
        enddo
      endif
!
      RETURN
      END SUBROUTINE unloadP
!
!-----------------------------------------------------------------------
!
      character*9 function PLab(n)
! Returns the parameter label of parameter number n      
      IMPLICIT none
      integer n,i,j,k
      character lab*9, lab2*3
!      
!      write(io,'("function PLab(n),  n=:",I4)') n
      if (n.le.0 .or. n.gt.maxTotPar) then
        write(io,'("***FATAL: Invalid parameter number:",I4)') n
        stop
      endif
      lab="****"      
      if (n.gt.0  .and. n.le.20) then
        lab=Para_label(n)
      endif
      if (LFtype.eq."CF  ") then
        if      (n.gt.20 .and. n.le.36) then; lab=BkqR_label(n-20+1) 
        else if (n.gt.36 .and. n.le.55) then; lab=BkqI_label(n-40+1)
        endif     
      Else if (LFtype.eq."AOM ") then
        if      (n.ge.21 .and. n.le.29) then; write(lab,'("Esig(",I1,")")') n-20 
        else if (n.ge.30 .and. n.le.40) then; write(lab,'("Esig(",I2,")")') n-20 
        else if (n.ge.41 .and. n.le.49) then; write(lab,'("EpiX(",I1,")")') n-40 
        else if (n.ge.50 .and. n.le.60) then; write(lab,'("EpiX(",I2,")")') n-40 
        else if (n.ge.61 .and. n.le.69) then; write(lab,'("EpiY(",I1,")")') n-60 
        else if (n.ge.70 .and. n.le.80) then; write(lab,'("EpiY(",I2,")")') n-70 
        else if (n.ge.81 .and. n.le.89) then; write(lab,'("theta(",I1,")")') n-80 
        else if (n.ge.90 .and. n.le.100)  then; write(lab,'("theta(",I2,")")') n-80 
        else if (n.ge.101 .and. n.le.109) then; write(lab,'("phi(",I1,")")') n-100 
        else if (n.ge.110 .and. n.le.120) then; write(lab,'("phi(",I2,")")') n-100 
        else if (n.ge.121 .and. n.le.129) then; write(lab,'("chi(",I1,")")') n-120 
        else if (n.ge.130 .and. n.le.140) then; write(lab,'("chi(",I2,")")') n-120 
        endif
      else if (LFtype.eq."LF1E") then
        k=0
        do i=1,2*Lvalue+1
          do j=1,i
            k=k+1
            if (n.eq.k+20 .and. Lvalue.eq.2) write(lab,'("d(",2I1,")")') i,j 
            if (n.eq.k+20 .and. Lvalue.eq.3) write(lab,'("f(",2I1,")")') i,j 
          enddo
        enddo  
      endif      
      if (ED_IntMech.eq."Full".and. n.gt.200 .and. n.le.290) then
        if (IntN1.eq.1) lab = Altp_label(Int((n-200+1)/2))
        if (IntN1.eq.2) lab = "undefined"
        if (IntN1.eq.3) lab = "undefined"
        if (IntN1.eq.4) lab = "undefined"
      endif
      IF (CCFtype.eq.1 .and. n.ge.141 .and. n.le.200) then
        write(lab,'("D_",i2,"_",i2)') (CCFindex(n-140,j),j=2,3)
      endif
      IF (CCFtype.eq.3 .and. n.ge.141 .and. n.le.200) then
        if (CCFindex(n-140,1).le.9) then; write(lab,'("G",i1,"_",i2,"_",i2)') (CCFindex(n-140,j),j=1,3)
        elseif (CCFindex(n-140,1).eq.10) then; write(lab,'("10A_",i2,"_",i2)') (CCFindex(n-140,j),j=2,3)
        elseif (CCFindex(n-140,1).eq.11) then; write(lab,'("10B_",i2,"_",i2)') (CCFindex(n-140,j),j=2,3)
        elseif (CCFindex(n-140,1).eq.12) then; write(lab,'("G11_",i2,"_",i2)') (CCFindex(n-140,j),j=2,3)
        endif        
      endif
      
      if (n.eq.301) lab="Omega(2) "
      if (n.eq.302) lab="Omega(4) "
      if (n.eq.303) lab="Omega(6) "
      if (n.eq.304) lab=" Hx      "
      if (n.eq.305) lab=" Hy      "
      if (n.eq.306) lab=" Hz      "
      if (n.eq.307) lab=" alpha   "
      if (n.eq.308) lab=" beta    "
      if (n.eq.309) lab=" gamma   "
      if (n.eq.310) lab=" kx      "
      if (n.eq.311) lab=" ky      "
      if (n.eq.312) lab=" kz      "
      
      if (n.ge.351 .and. n.le.359) write(lab,'("offset",I1)') n-350
      if (n.ge.360 .and. n.le.400) write(lab,'("offset",I2)') n-350
 
      if (n.gt.400 .and. n.le.450) lab=variables(n-400)
      PLab = lab
!      
      RETURN
      END function PLab
!
!-----------------------------------------------------------------------
!
      real*8 function CCFvalue(i,k,q)
! Returns the CCF parameter  with the values i,k,q     
! i=1-12, k=2(2)12, q=0,k
      IMPLICIT none
      integer i1,i,k,q
      REAL*8 x
!      
      X=0.0d+00
      do i1 = 1, N_CCF
        if (i.eq.CCFindex(i1,1).and.k.eq.CCFindex(i1,2).and.q.eq.CCFindex(i1,3)) X=CCF(i1)
      enddo
      CCFvalue=X
      RETURN
      END function CCFvalue
!
!-----------------------------------------------------------------------

END MODULE f_e_parameters
!  634