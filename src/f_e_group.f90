MODULE f_e_group
USE f_e_data
USE f_e_parameters
USE f_e_wigner
!
!  This module contains the point group calculations.
!

IMPLICIT NONE

! PRIVATE

PUBLIC 

CONTAINS
!-----------------------------------------------------------------------
!
      SUBROUTINE DJMPM(D,I2J,BETA)
!
!  Returns the Wigner rotation matrix d(j,m',m,beta).
!  In this program beta = -pi/2, thetaT=cos-1[1/sqrt(3)]
!
!  Each d(j,m',m,pi/2) submatrix is of the size ((2J+1)x(2J+1)).
!  It is given by:
!  eq. (2.7)  Silver, "Irreducible Tensor Methods", (1976).
!  eq. (3.57) Zare, "Angular Momentum", (1989).
!
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      REAL*8 D(max_N,max_N)
      PARAMETER(Z=0.0D+00, H=0.5D+00)
!
      C = COS(H*BETA)
      S = SIN(H*BETA)
      DO I=1,I2J+1
        DO J=1,I2J+1
          D(I,J) = Z
        enddo
      enddo
        
!
      AJ = DBLE(I2J)*H
      DO I1=1,I2J+1
        AMP = AJ - (I1-1)
        JMMP= INT(AJ-AMP + 0.001)
        JPMP= INT(AJ+AMP + 0.001)
        DO I2=1,I2J+1
          AM = AJ - (I2-1)
          JPM = INT(AJ+AM + 0.001)
          JMM = INT(AJ-AM + 0.001)
          X = AM-AMP
          IF (X.LT.Z) MMMP=INT(X-0.001)
          IF (X.GE.Z) MMMP=INT(X+0.001)
          ISTART = MAX(0,  MMMP)
          ISTOP  = MIN(JPM,JMMP)
          X = Z
          IF (ISTART.GT.ISTOP) GOTO 100
          DO I=ISTART,ISTOP
            X = X + (-1)**I*C**(I2J+MMMP-2*I)*(-S)**(2*I-MMMP) /   &
                (FACT(JPM-I)*FACT(JMMP-I)*FACT(I)*FACT(I-MMMP))
          enddo
          X = SQRT(FACT(JPM)*FACT(JMM)*FACT(JPMP)*FACT(JMMP))*X
          D(I1,I2) = X
 100      continue
        enddo  
      enddo
!
      RETURN
      END SUBROUTINE DJMPM
!
!-----------------------------------------------------------------------
!
      REAL*8 FUNCTION FACT(N)
!  Returns n!
      REAL*8 X
      INTEGER*4 N,I
      FACT=1.0D+00
      IF (N.LT.0) THEN
        WRITE(*,'(A)') '***FATAL: N is negative in function FACT.'
        STOP
       ELSE IF (N.EQ.0) THEN
        RETURN
       ELSE IF (N.GT.30) THEN
        WRITE(*,'(A)') '***FATAL: N is too big in function FACT.'
        STOP
       ELSE 
        X = 1.0D+00
        DO I=1,N
          X = X*DBLE(I)
        enddo
        FACT = X
        RETURN
      ENDIF
!
      END FUNCTION FACT
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GETGRP()
!
!  Returns information about the group given by IDG.
!  by reading the file: GROUP.DAT. If the file is not 
!  found then the subroutine returns with IDG=ISO=IGP=0.
!
!  IDG is an integer between 1 and 37 that identifies the group.
!  The first 32 groups are defined according to the tables of 
!  Koster, Dimmock, Wheeler & Statz; "Properties of the 32 crystallographic 
!          point groups",MIT Press,(1963).
!  The additional groups TdT,OT,OhT are the Td,O,Oh with trigonal axes.
!  The additional groups C2Y is the C2 point group with the C2(y) symmetry element.
!  The additional groups CsY is the Cs point group with the s(y) symmetry element.
!
!   1. C1         9. C4         17. C3i        25. C6v      33. TdT  
!   2. Ci        10. S4         18. D3         26. D3h      34. OT 
!   3. C2        11. C4h        19. C3v        27. D6h      35. OhT 
!   4. Cs        12. D4         20. D3d        28. T        36. C2Y 
!   5. C2h       13. C4v        21. C6         29. Th       37. CsY
!   6. D2        14. D2d        22. C3h        30. O
!   7. C2v       15. D4h        23. C6h        31. Td
!   8. D2h       16. C3         24. D6         32. Oh
!
!  ISO=(0/1) The orbital/complete wavefunctions are used.
!  IGP=(0/1) The transformation properties of the wavefunctions 
!            under the symmetry operations (are not/are) printed.
!
!  GROUP is the label of group IDG above which is returned.
!
!  NCLASS is the number of classes in the group (= total number  
!         of irreducible representations.)
!         The maximum number of classes is given by NGRP=24,
!         which is the number of classes in the C6h double group.
!  NNUM: Defines a constant w = exp(2*Pi/NNUM) 
!  NROT:   is the number of rotational angle needed to describe
!         the symmetry operation.
!  NSINGLE: Number of irreps in the single group.
!
!  GROUPT(NCLASS,NCLASS) contains the character table.
!
!  GROUPO(NCLASS) contains the characters of the symmetry operations.
!
!  GRPR1(NCLASS), GRPR2(NCLASS) contain the characters which 
!        denote the irreducible representations of the group.
!  GRPR1: Contains Bethe's notation.
!  GRPR2: Contains Schoenflies (Mulliken's) notation.
!
!  ISYMEL describes the symmetry operations of the group.
!         All operations can be described in terms of rotations,
!         inversion and time reversal.
!  ISYMEL(K,1) = 1 Rotation;  2 Inversion.Rotation
!  ISYMEL(K,2) = The number of elements in each class.
!  ISYMEL(K,3) = The degeneracy of the irreps.
!  ISYMEL(K,4) = 1,2,3 = a,b,c for behaviour under Time Reversal: KDWS pg 11
!
!  RSYMEL describes the rotation axis of the symmetry element.
!  RSYMEL(K,1) = alpha
!  RSYMEL(K,2) = beta
!  RSYMEL(K,3) = gamma    [in units of pi]
!  where alpha, beta are the angles in the Wigner rotational
!  matrices D(alpha,beta,gamma).
!
!  These angles are given in Table 2.1 and figures 1.1, 1.2
!  and 1.3 of Bradley & Cracknell, "Mathematical theory of 
!  Symmetry in Solids", Claredon, Oxford, (1972).
! 
!  NSINGLE The number of irreps in the single group. (will be <NCLASS)
!
!  NGRP The maximum number of classes for arrays. 
!       This has been set to NGRP=24 in main program, which is the number 
!       of classes in the C6h double group.
!
!  The Character Table is printed if global check is true.
!
!--------1---------2---------3---------4---------5---------6---------7--
      IMPLICIT REAL*8(A-H,O-Z), INTEGER*4(I-N)
      COMPLEX*16 Ortho(24,24),II,W
      real*8 Z,D2
      PARAMETER(Z=0.0D+00, D2=2.0d+00, D3=3.0d+00)
      INTEGER*4 I1(24),I2(24),I3(24)
      CHARACTER Ugroup*4, ASTER*4, DATE*10, TimeRev(3)*1/"a","b","c"/, FMT*40
      CHARACTER lower*26, upper*26
      PARAMETER(lower="abcdefghijklmnopqrstuvwxyz")
      PARAMETER(upper="ABCDEFGHIJKLMNOPQRSTUVWXYZ")

!
      II=DCMPLX(0.0D+00,1.0D+00)
      PI=4.0D+00*ATAN(1.0D+00)
      TRIG=ACOS(1.0D+00/SQRT(3.0D+00))   ! trigonal angle
!
!  Set pointer to the beginning of the group IDG in file GROUP.DAT:
      REWIND(IG)

      IDG=0
      do i=1,37
        Ugroup=Groups(i)
        DO K=1,4
          DO J=1,26
            IF (Ugroup(K:K).EQ.lower(J:J)) then
              Ugroup(K:K)=UPPER(J:J)
              goto 5
            endif
          enddo
 5       continue
        enddo      
        if (TRIM(GROUP).eq.TRIM(Ugroup)) IDG=i
      enddo
      if (IDG.eq.0) then
        write(io,'("***FATAL: Could not find point group:",A4," in file GROUP.DAT")') GROUP
        STOP
      endif
      
      I=0
 10   I=I+1; READ(IG,'(A4)',ERR=60,END=30) ASTER
      IF (ASTER.EQ.'****') goto 10
      IF (ASTER.EQ.'DATE') then
        I=I+1;  READ(IG,'(A10)') DATE
        if (check1(10)) WRITE(iDebug,'("The GROUP.DAT is date stamped:",A10)') DATE
        goto 10
      endif
      I=I+1; READ(ASTER,'(I4)') ID
!      WRITE(IO,'("ASTER,ID,I=",A4,2I4)') ASTER,ID,I
      IF (ID.EQ.IDG) GOTO 100
      do j=1,100
        I=I+1;  READ(IG,'(A4)',ERR=60,END=30) ASTER
        IF (ASTER.EQ.'****') goto 10
      enddo
!
 30   CONTINUE
      WRITE(IO,'(" Point group =",I4," has not been found in the file:GROUP.DAT",/,I3," lines read.")') IDG,I
      GOTO 500
 60   WRITE(IO,'(" ERROR: Reading line",I4," of file GROUP.DAT")') I
      GOTO 500
!
!  Read the information for the appropriate group.
 100  READ(IG,110) GROUP,NCLASS,NNUM,NROT,NSINGLE
 110  FORMAT(A4,4I4)
      IF (NNUM.NE.0) W = EXP(2.0D+00*PI*II/ABS(NNUM))
      DO J=1,3
        DO I=1,NCLASS
          RSYMEL(I,J) = Z
        enddo
      enddo  
      DO J=1,4
        READ(IG,'(20I4)') (ISYMEL(I,J),I=1,NCLASS)
      enddo
      DO J=1,NROT
        READ(IG,'(40I2)') (I1(I),I2(I),I=1,NCLASS)
        DO I=1,NCLASS
          IF (I2(I).EQ.0) THEN
            WRITE(IO,'("***FATAL: RSYMEL(",I2,",",I2,") infinite.")') I,J
            WRITE(IO,'(" File: GROUP.DAT is corrupted.")')
            STOP
          ENDIF
!          IF (ABS(I1(I)).EQ.9.AND.ABS(I2(I)).EQ.9) THEN
!            RSYMEL(I,J) = DBLE(SIGN(I1(I),1))*TRIG  ! special angle thetaT = 54.735..
          IF (I1(I).EQ.11 .and. I2(I).EQ.11) THEN
            RSYMEL(I,J) = D2*TRIG       ! special angle 2*thetaT used for trigonal axes TdT,OT,OhT
           ELSEIF (I1(I).EQ.12 .and. I2(I).EQ.12) THEN
            RSYMEL(I,J) = PI - D2*TRIG  ! special angle pi - 2*thetaT
           ELSEIF (I1(I).EQ.13 .and. I2(I).EQ.13) THEN
            RSYMEL(I,J) = TRIG  ! special angle thetaT
           ELSE 
            RSYMEL(I,J) = DBLE(I1(I))/DBLE(I2(I))*PI   ! This changed!!!
          ENDIF
        enddo
      enddo  

!  The data in file GROUP.DAT is from Koster, Dimmock, Wheeler & Statz, "Properties of the 32 Point Groups".
!  They use passive rotations, so you would expect the -ve of these angles to be taken, to convert them into 
!  active rotations. However they are used in a -ve sense down the line: Exp(+i.m.alpha)djm'm(-beta)Exp(+i.m.gamma)
!  instead of Exp(-i.m.alpha)djm'm(beta)Exp(-i.m.gamma), so positive angles are kept here.   
! 160      RSYMEL(I,J) = -DBLE(I1(I))/DBLE(I2(I))*PI  ! original 
!      write(io,161)
! 161  format("***Warning the sign of RSYMEL has been change to",
!     +       " +ve in subroutine GETGRP***")
      READ(IG,'(20A4)') (GROUPO(I),I=1,NCLASS)   ! Operators
      READ(IG,'(20A4)') (GRPR1(I),I=1,NCLASS)    ! Double PointGroup irreps
      READ(IG,'(20A4)') (GRPR2(I),I=1,NCLASS)    ! Single PointGroup irreps
!  Three ways of reading the character table, depending on NNUM.
      IF (NNUM.LT.0) THEN  
        DO I=1,NCLASS
          READ(IG,'(40I2)') (I1(J),I2(J),I3(J),J=1,NCLASS)
          DO J=1,NCLASS
            IF (I2(J).LE.0) THEN
              WRITE(IO,'("***FATAL: GROUPT(",I2,",",I2,") SQRT(-).")') I,J
              WRITE(IO,'(" File: GROUP.DAT is corrupted.")')
              STOP
            ENDIF
            GROUPT(I,J) = DBLE(I1(J))*SQRT(DBLE(I2(J)))*W**I3(J)
          enddo
        enddo  
       ELSE IF (NNUM.EQ.0) THEN
        DO I=1,NCLASS
          READ(IG,'(40I2)') (I1(J),I2(J),J=1,NCLASS)
          DO J=1,NCLASS
            IF (I2(J).LE.0) THEN
              WRITE(IO,'("***FATAL: GROUPT(",I2,",",I2,") SQRT(-).")') I,J
              WRITE(IO,'(" File: GROUP.DAT is corrupted.")')
              STOP
            ENDIF
            GROUPT(I,J) = DCMPLX(DBLE(I1(J))*SQRT(DBLE(I2(J))),Z)
          enddo
        enddo  
       ELSE IF (NNUM.GT.0) THEN
        DO I=1,NCLASS
          READ(IG,'(40I2)') (I1(J),J=1,NCLASS)
          DO J=1,NCLASS
            GROUPT(I,J) = W**I1(J)
          enddo
        enddo  
      ENDIF
!
!  Check for consistency: 
!  a) The number of symmetry classes should equal the number
!  of irreducible representations.
!  b) The sum over the number of symmetry elements in each 
!  class should equal the sum of the squares of the degeneracies
!  of the irreducible representations.
!
      ISYM = 0
      IDEG = 0
      DO I=1,NCLASS
        ISYM = ISYM + ISYMEL(I,2)  
        IDEG = IDEG + ISYMEL(I,3)**2  
      enddo  
      IF (ISYM.NE.IDEG) THEN
        WRITE(IO,310)
 310    FORMAT('***WARNING Point group:',I3,' is inconsistent.',/,              &
               ' The sum of the symmetry elements does not equal',/,  &
               ' the sum of the squares of the degeneracies.')
      ENDIF
!
! Check symmetry operations, no.of elements/class, degeneracy of irreps:
      DO I=1,NCLASS
        IF ((ISYMEL(I,1).NE.1 .AND. ISYMEL(I,1).NE.2) .OR.  &   
            (ISYMEL(I,2).LE.0 .OR.  ISYMEL(I,2).GT.13) .OR. &  
            (ISYMEL(I,3).LE.0 .OR.  ISYMEL(I,3).GT.4) )THEN  
          WRITE(IO,'("***FATAL: ISYMEL=",3I2)') (ISYMEL(I,J),J=1,3)
          WRITE(IO,'(" File: GROUP.DAT is corrupted.")')
          STOP
        ENDIF
! Check that beta is one of: 0, pi/2, pi, ThetaT, 2*ThetaT, pi-2*ThetaT, 2pi, 3pi.
        IF ((ABS(RSYMEL(I,2))           .GT.1.0D-10).AND.  &
            (ABS(RSYMEL(I,2)-0.5D+00*PI).GT.1.0D-10).AND.  &
            (ABS(RSYMEL(I,2)-TRIG)      .GT.1.0D-10).AND.  &
            (ABS(RSYMEL(I,2)-D2*TRIG)   .GT.1.0D-10).AND.  &
            (ABS(RSYMEL(I,2)-PI+D2*TRIG).GT.1.0D-10).AND.  &
            (ABS(RSYMEL(I,2)-PI)        .GT.1.0D-10).AND.  &
            (ABS(RSYMEL(I,2)-D2*PI)     .GT.1.0D-10).AND.  &
            (ABS(RSYMEL(I,2)-D3*PI)     .GT.1.0D-10) ) THEN
          WRITE(IO,'("***FATAL: Invalid beta=",F10.4)') RSYMEL(I,2)
          WRITE(IO,'(" File: GROUP.DAT is corrupted.")')
          STOP
        ENDIF
      enddo  
!
!  Print character table
      if (check1(10)) then
        write(idebug,'(/,"Character Table for:",A4,4X,   &
                     "GROUP.DAT date stamped:",A10)') GROUP,DATE
        write(idebug,'(3X,A4,3X,"|",20(3X,A4,3X))') GROUP,(GROUPO(i),i=1,NCLASS)
        write(idebug,'("-------",20A10)') ("----------",i=1,NCLASS+1)
        write(FMT,'("(1X,A4,1X,A4,""|"",",i2,"(2F5.2)"" | "",A1,"" |"")")') NCLASS
!        write(idebug,'("FMT=",A40)') FMT
        do i=1,NCLASS
          write(idebug,FMT) GRPR1(i),GRPR2(i),(GROUPT(i,j),j=1,NCLASS),TimeRev(ISYMEL(i,4))
        enddo
        write(idebug,'("-------",20A10)') ("----------",i=1,NCLASS+1)

        write(idebug,'(//,"Description of the symmetry elements:")')
        write(idebug,'(3X,A4,3X,"|",20(3X,A4,3X))') GROUP,(GROUPO(i),i=1,NCLASS)
        write(idebug,'("-",20A10)') ("----------",i=1,NCLASS+1)
        write(idebug,456) "ISYMEL(1)",(ISYMEL(i,1),i=1,NCLASS)
        write(idebug,456) "ISYMEL(2)",(ISYMEL(i,2),i=1,NCLASS)
        write(idebug,456) "ISYMEL(3)",(ISYMEL(i,3),i=1,NCLASS)
        write(idebug,456) "ISYMEL(4)",(ISYMEL(i,4),i=1,NCLASS)
        write(idebug,457) " ALPHA   ",(180.0d+00/PI*RSYMEL(i,1),i=1,NCLASS)
        write(idebug,457) " BETA    ",(180.0d+00/PI*RSYMEL(i,2),i=1,NCLASS)
        write(idebug,457) " GAMMA   ",(180.0d+00/PI*RSYMEL(i,3),i=1,NCLASS)
        write(idebug,'("-",20A10)') ("----------",i=1,NCLASS+1)
 456    format(1X,A9,"|",20(3X,I4,3X)) 
 457    format(1X,A9,"|",20(F10.2))    

        do i=1, NCLASS
          do j=1, NCLASS
            Ortho(i,j)=DCMPLX(0.0D+00,0.0D+00)
            do k=1, NCLASS
              Ortho(i,j)=Ortho(i,j)+ISYMEL(k,2)*DCONJG(GROUPT(i,k))*GROUPT(j,k) 
            enddo
            Ortho(i,j)=Ortho(i,j)/ISYM
          enddo
        enddo  
        write(idebug,'(//,"Orthogonality Table for:",A4)') GROUP
        write(idebug,'(10X,"|",20(3X,A4,3X))') (GRPR1(j),j=1,NCLASS)
        write(idebug,'(10X,"|",20(3X,A4,3X))') (GRPR2(j),j=1,NCLASS)
        write(idebug,'("-",20A10)') ("----------",i=1,NCLASS+1)
        do i=1, NCLASS
          write(idebug,469) GRPR1(i),GRPR2(i),(Ortho(i,j),j=1,NCLASS)
 469      format(A4,"*",A4,"*|",20(2F5.2)) 
        enddo

      endif
      RETURN
!
 500  WRITE(IO,510) 
 510  FORMAT('***WARNING: The following values are reset: IDG=0, ISO=0, IGP=0.')
      IDG=0
      ISO=0
      IGP=0
      RETURN
!
      END SUBROUTINE GETGRP
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GROUPL(NPRT)
!
!  Returns how the wavefunctions transforms under the different
!  symmetry operations of the point group that is being used to 
!  model the calculation. These characters are returned in the 
!  complex array: CSYM(II,I) for the II wavefunctions and the I
!  symmetry classes.
!
!  If ISO = 0, then only the angular part of the wavefunctions
!  is considered. If ISO = 1, then the total angular momentum of
!  the wavefunctions is considered. These situations are
!  appropriate when the spin-orbit coupling is zero and non-zero
!  respectively.
!
!  It also uses the transformations to find the projection onto the 
!  irreducible representations of the point group. These are returned 
!  in REAL array REPS(II,I) for the II wavefunctions and the Ith irrep. 
! 
!  The irreducible representations with the largest component
!  in each of the NPRT levels is the assignment of that level.
!  This is returned in IRREP(NPRT,1). For perfect symmetry, REPS 
!  will be 1 for the Ith irrep and 0 for all others.
!
!  If an assignment is dubious then it will be marked by asterisks.
!  The criteria for being acceptable is that the projection of the 
!  largest irrep is greater than 0.95 and the sum of the others is
!  less than 0.05.
!
!  If check is true, then the Wigner beta matrix is printed for beta=Pi/2.
!
      IMPLICIT none 
      COMPLEX*16 ANS,W(max_N),ROTW(max_N),VEC,X(24)
      COMPLEX*16 EXPI,FAC,CX(24),PMJ,C1,CI,CZ  
!      REAL*8 T(max_N,max_N),R(max_N,max_N)  T->VR  R->VI
      REAL*8 FACINV,Pi, TRIG, SMALL, AL,AS,AJ,AMJ,BMJ,AMS,AML,AJP,AMJP,PJMJ
      REAL*8 ALPHA,BETA,GAMMA, Z,H,D1,D2,D3,BIT,PHASE,PHASE2,XX,EM1,EP1,SDEG,GMAX,GMAX2,GMIN
      logical check
      INTEGER*4 NDEG(max_N),NPRT,IORDER,IT,I,II,J,K,ISTART,ILD,ISD,ITD,L,IDEG,ISP,IX,iphase
      INTEGER*4 ISEL,ITR,IS1,IS2,IL1,IL2,IJMJ,IMAX,IMAX2,IMMJ,I1,I2,N1,N2,M,MM,MLP,LL,I2MJ,ML
      integer*4 quickBeta(24)
      PARAMETER(SMALL=1.0D-10)
      PARAMETER(Z=0.0D+00, H=0.5D+00, D1=1.0D+00, D2=2.0D+00, D3=3.0D+00, BIT=0.001D+00)
!
      check=.true.
      CZ=DCMPLX(0.0D+00,0.0D+00)
      C1=DCMPLX(1.0D+00,0.0D+00)
      CI=DCMPLX(0.0D+00,1.0D+00)
!
      PI=4.0D+00*ATAN(D1)
      TRIG=ACOS(1.0D+00/SQRT(3.0D+00))   ! trigonal angle
      EXPI=EXP(CI)
      IORDER=0
      DO I=1,NCLASS
        IORDER = IORDER + ISYMEL(I,2)
      enddo 
      IF (check1(10)) WRITE(idebug,'(/,"Order of group=",I4,";  Number of classes=",I4)') IORDER,NCLASS
!
! Pre-calculate the Rotational matrix for beta if a group contains a beta (pi/2, 2*ThetaT or pi-2*ThetaT).
      IF (IDG.GE.28 .AND. IDG.LE.35) THEN
        if (IDG.GE.28 .and. IDG.LE.32) BETA=-PI/2
        if (IDG.GE.33 .and. IDG.LE.35) BETA=D2*TRIG       
        CALL ROTB(VI,VR,BETA)  !  VI is the rotation matrix VR is workspace
        IF (check.and.check1(10)) call WriteBmat(BETA,1)
        if (IDG.GE.33 .and. IDG.LE.35) then
          BETA=PI-D2*TRIG
          CALL ROTB(AI,AR,BETA)  !  AI is the rotation matrix AR is workspace
          IF (check.and.check1(10)) call WriteBmat(BETA,2)
        endif     
      ENDIF
!      
! Determine which symmetry elements this applies to.
      DO IX=1,NCLASS
        quickBeta(IX)=0
        BETA=RSYMEL(IX,2)
        IF (ABS(BETA-PI/2).lt.SMALL) quickBeta(IX)=1
        IF (ABS(BETA-D2*TRIG).LT.SMALL) quickBeta(IX)=1
        IF (ABS(BETA-(PI-D2*TRIG)).LT.SMALL) quickBeta(IX)=2
      enddo
      IF (check1(10)) WRITE(idebug,'("quickBeta=",24I4)') (quickBeta(IX),IX=1,NCLASS) ! 24 for C6h
!
! Prepare transformation matrix that decouples the |J,MJ> basis into |L,Ml>|S,Ms> basis if the symmetry 
! properties of only orbital parts are wanted.
! One can then operate on S & L states independently.
!
      IF (ISO.EQ.0) THEN
        DO I=1,n_matrix
          DO J=1,n_matrix
            VR(I,J)=Z
          enddo
        enddo  
        ISTART=0
        DO IT=1,n_states
          ILD=2*TS_Bases(2,IT,nelectrons) + 1  !  2L + 1
          ISD=  TS_Bases(1,IT,nelectrons) + 1  !  2S + 1 
          ITD=ILD*ISD
          IF (IT.NE.1) ISTART=ISTART+(2*TS_Bases(2,IT-1,nelectrons)+1)*(TS_Bases(1,IT-1,nelectrons)+1)
          DO I=1,ITD
            II=I+ISTART
            AS= DBLE(FullBasis(1,II))/D2
            AL= DBLE(FullBasis(2,II))
            AJ= DBLE(FullBasis(3,II))/D2
            AMJ=DBLE(FullBasis(4,II))/D2
            BMJ=-AMJ
            PHASE=SQRT(2.0D+00*AJ+D1)
            XX=ABS(AL-AS+AMJ)+BIT
            XX=MOD(XX,2.0D+00)
            IF (XX.GT.0.02D+00) PHASE=-PHASE
            PHASE2=D1; iphase=int(AJ-AMJ); if (mod(iphase,2).ne.0) PHASE2=-D1 
            L=ISTART
            DO K=1,ISD
              AMS=AS-K+1
              DO J=1,ILD
                L=L+1
                AML=AL-J+1
                VR(II,L)=PHASE*PHASE2*WIG3J(AL,AS,AJ,AML,AMS,BMJ)
              enddo
            enddo   
            IF (check1(10)) then
              IF (ii.eq.1) then
                write(idebug,'("matrix to decouple |J,Mj> into |L,Ml>|S,Ms> states")') 
                write(idebug,'("(Note an extra phase2=(-1)^(J-MJ) is used in this matrix to compensate",  &
                                    /," for the WFs here differing to those of CAMMAGS in sign.)")') 
                write(idebug,'(12X,10("|",F4.1,",",F3.1,"> "))') (((AL-J+D1),(AS-K+D1),J=1,ILD),K=1,ISD) 
              endif                         
              write(idebug,'("|",F3.1,",",F4.1,">",14F11.5)') AJ,AMJ,(VR(ii,j),j=1,n_matrix)
            endif  
          enddo  ! i
        enddo  !  it
        IF (check1(10)) write(idebug,'("Real part of WF expressed as uncoupled |L,Ml>|S,Ms> states")') 
      ENDIF
!
!  Preparation complete,
!    ISO=0:        VR contains the decoupling matrix.
!    Groups 28-32: VI contains the rotation matrix for beta = pi/2
!    Groups 33-35: VI contains the rotation matrix for beta = 2*thetaT
!    Groups 33-35: AI contains the rotation matrix for beta = pi-2*thetaT
!  TODO: the matrix for pi-2*thetaT should be related to 2*thetaT, could save time by not having 
!        to actually calculate AI. 
! 
!  Start calculation....
!  Loop for every eigenvector:
      DO 500 II=1,NPRT
!        WRITE(IO,'("Calculating symmetry properties of eigenvector:",I4)') II
!
!  Determine the degeneracy:
!  If state is nondegenerate:                     IDEG=0.
!  If degenerate & the first member:              IDEG=1.
!  If degenerate & neither first nor last member: IDEG=2.
!  If degenerate & the last member:               IDEG=3.
!
        IDEG=0
        EM1=D1
        EP1=D1
        IF (II.NE.1)    EM1=ABS(ENGS(II)-ENGS(II-1))
        IF (II.NE.NPRT) EP1=ABS(ENGS(II)-ENGS(II+1))
        IF (EM1.GT.Edegen .AND. EP1.LT.Edegen) IDEG=1
        IF (EM1.LT.Edegen .AND. EP1.LT.Edegen) IDEG=2
        IF (EM1.LT.Edegen .AND. EP1.GT.Edegen) IDEG=3
!
!  Spin degeneracy.
        SDEG=D1
        IF (ISO.EQ.0) THEN  ! SOC=0; Orbital degeneracy only, split into |L,Ml>|S,Ms> basis 
          SDEG = Z
!          WRITE(Idebug,'("NSPIN:",I4)') NSPIN
          DO ISP=1,Nspin
!            WRITE(Idebug,'("SPIN(II,ISP),S_mult(ISP)",F8.3,I4)') SPIN(II,ISP),S_mult(ISP)
            SDEG = SDEG + SPIN(II,ISP)*S_mult(ISP)
          enddo
        ENDIF
!        WRITE(Idebug,'("Spin degeneracy:",F8.3)') SDEG
!
        X(1)=C1  ! Identity symmetry element
        DO 400 IX=2,NCLASS
          X(IX)=CZ
          ISEL= ISYMEL(IX,1)
          ITR = ISYMEL(IX,2)
          ALPHA=RSYMEL(IX,1)
          BETA =RSYMEL(IX,2)
          GAMMA=RSYMEL(IX,3)
!   Inversion:
          FAC=D1
          FACINV=D1
          IF (ISEL.EQ.2 .and. Lvalue.eq.3 .and. n_odd) FACINV = -D1  ! d-electrons even, odd number of f-electrons odd wrt inversion
!
          IF (ISO.EQ.0) THEN
!  Uncouple into |L,Ml>|S,Ms> states
            DO K=1,n_matrix
              W(K)=CZ
              DO J=1,n_matrix
                IF (.not.MatComplex) W(K)=W(K)+DCMPLX(MAT(J,II),Z)*VR(J,K)
                IF (     MatComplex) W(K)=W(K) +        CMAT(J,II)*VR(J,K)
              enddo
            enddo            
            if (check1(10) .and. ix.eq.2) write(idebug,'("W(",I2,")=",14F9.5)') II,(dble(W(K)),k=1,10)
!
            IF (quickBeta(IX).eq. 0) THEN      ! beta is simple.. 0,+/-pi, 
              DO M=1,n_matrix
                IF (.not.MatComplex) VEC=MAT(M,II)
                IF (     MatComplex) VEC=DCONJG(CMAT(M,II))
                IF (ABS(VEC).LT.SMALL) GOTO 250
                ISTART=0
                DO IT=1,n_states
                  ILD=2*TS_Bases(2,IT,nelectrons) + 1  !  2L + 1
                  ISD=  TS_Bases(1,IT,nelectrons) + 1  !  2S + 1 
                  IF (IT.NE.1) ISTART=ISTART+(2*TS_Bases(2,IT-1,nelectrons)+1)*(TS_Bases(1,IT-1,nelectrons)+1)
!  LL can only take integer values....
                  LL=TS_Bases(2,IT,nelectrons)
                  I=ISTART
                  DO I1=1,ISD
                    DO I2=1,ILD
                      I=I+1
                      L=ILD-I2+ISTART+1+(I1-1)*ILD   ! L
                      IF (ABS(VR(M,I)).LT.SMALL) GO TO 200
!  ML can only take integer values...
                      ML=-LL-1+I2  ! negative to CAMMAG
                      MM=I
!   Rotations:
                      IF (ABS(BETA+PI).LT.SMALL .or. ABS(BETA-PI).LT.SMALL .or. ABS(BETA-D3*PI).LT.SMALL) THEN   !  beta = +/-pi
                        FAC=(-D1)**(LL-ML)*EXPI**(ML*(ALPHA+GAMMA))
                        MM = L
                       ELSE IF (ABS(BETA).LT.SMALL .or. ABS(BETA-D2*PI).LT.SMALL) THEN  !  beta = 0, 2Pi 
                        FAC=EXPI**(ML*(ALPHA+GAMMA))
                       ELSE
                        WRITE(IO,'("***FATAL: Invalid value of beta = ",F10.4, " in subroutine: GROUPL")') BETA
                        STOP
                      ENDIF
                      X(IX)=X(IX)+FAC*VEC*VR(M,I)*W(MM)
 200                  CONTINUE 
                    enddo  !  I2=1,ILD
                  enddo  !  I1=1,ISD
                enddo  ! IT=1,n_states
 250            CONTINUE
              enddo  ! M=1,n_matrix
            ELSE  !  a more complicated beta
              ISTART=0
              DO IT=1,n_states
                ILD=2*TS_Bases(2,IT,nelectrons) + 1  !  2L + 1
                ISD=  TS_Bases(1,IT,nelectrons) + 1  !  2S + 1 
                ITD=ILD*ISD
                IF (IT.NE.1) ISTART=ISTART+(2*TS_Bases(2,IT-1,nelectrons)+1)*(TS_Bases(1,IT-1,nelectrons)+1)
                I=ISTART
                LL=TS_Bases(2,IT,nelectrons)
                DO IS1=1,ISD
                  DO IL1=1,ILD
                    I=I+1
                    ML =-LL - 1 + IL1 ! negative to CAMMAG
                    J=ISTART
                    ROTW(I) = CZ
                    DO IS2=1,ISD
                      DO IL2=1,ILD
                        J=J+1
                        MLP =-LL - 1 + IL2 ! negative to CAMMAG
!  VI/AI are the rotation matrix (collection of d(j,m',m,beta) ) for precalculated values of beta
                        if (quickBeta(IX).eq.1) ROTW(I) = ROTW(I) + EXPI**(MLP*GAMMA)*VI(J,I)*EXPI**(ML*ALPHA)*W(J)
                        if (quickBeta(IX).eq.2) ROTW(I) = ROTW(I) + EXPI**(MLP*GAMMA)*AI(J,I)*EXPI**(ML*ALPHA)*W(J)
                      enddo 
                    enddo
                  enddo
                enddo
                DO J=ISTART+1,ISTART+ITD
                  X(IX) = X(IX) + DCONJG(W(J))*ROTW(J)
                enddo
              enddo ! IT=1,n_states
            ENDIF
!
          ELSE IF (ISO.EQ.1) THEN ! SOC#0; Find symmetry propoertis of |J,MJ> basis 
!
            ISTART=0
            DO IT=1,n_states
              ITD=(TS_Bases(1,IT,nelectrons)+1)*(2*TS_Bases(2,IT,nelectrons) + 1)  !  (2S+1)(2L+1) 
!  ITD is the total degeneracy of the term = 2J+1.
              IF (IT.NE.1) ISTART=ISTART+(2*TS_Bases(2,IT-1,nelectrons)+1)*(TS_Bases(1,IT-1,nelectrons)+1)
              DO I=1,ITD
                IF (.not.MatComplex) W(I) = DCMPLX(MAT(ISTART+I,II),Z)
                IF (     MatComplex) W(I) =        CMAT(ISTART+I,II)
              enddo
              DO I=1,ITD
                ROTW(I) = Z
                AJ  = DBLE(FullBasis(3,ISTART+I))/D2
                AMJ = DBLE(FullBasis(4,ISTART+I))/D2
                BMJ=-AMJ
                IF (quickBeta(IX) .eq. 0) THEN     ! beta simple
                  IF (.not.MatComplex) VEC = DCMPLX(MAT(ISTART+I,II),Z)
                  IF (     MatComplex) VEC = DCONJG(CMAT(ISTART+I,II))
                  IF (ABS(VEC).LT.SMALL) GOTO 380
!  (-1)**(J-MJ) note that this is real since (J-MJ) is always integral.
                  IJMJ = INT(AJ-AMJ+BIT)
                  PJMJ = D1
                  IF (MOD(IJMJ,2).GT.0) PJMJ=-D1
!  (-1)**MJ can be complex, depending whether MJ is integral or half
!  integral. MJ can also be positive or negative.
!not used                  I2MJ = INT(2*AMJ+SIGN(BIT,AMJ))
!not used                  PMJ =  CI**(MOD(I2MJ,4))
                  DO I1=1,ITD
                    IMMJ=I1
                    I2=ISTART+I1
         IF (DBLE(FullBasis(3,I2))/D2.EQ.AJ .AND. DBLE(FullBasis(4,I2))/D2.EQ.BMJ) GOTO 330
                  enddo
 330              CONTINUE
                  M=I
!   Rotations:
                  IF (ABS(BETA+PI).LT.SMALL .or. ABS(BETA-PI).LT.SMALL) THEN       !  beta=+/-pi
                    FAC = PJMJ*EXPI**(AMJ*(ALPHA+GAMMA))
                    M = IMMJ
                   ELSE IF (ABS(BETA).LT.SMALL) THEN    !  beta=0
                    FAC = EXPI**(AMJ*(ALPHA+GAMMA))
                  ENDIF
!
                  X(IX) = X(IX) + FAC*VEC*W(M)
!
                ELSE IF (quickBeta(IX) .ne. 0) THEN  !  beta= pi/2, 2*thetaT, pi-2*thetaT
                  DO J=1,ITD
                    AJP= DBLE(FullBasis(3,ISTART+J))/D2
                    AMJP=DBLE(FullBasis(4,ISTART+J))/D2
!  VI/AI are the rotation matrix (collection of d(j,m',m,beta) ) for precalculated values of beta
                    if (quickBeta(IX).eq.1) ROTW(I) = ROTW(I) + EXPI**(AMJP*GAMMA)*VI(ISTART+J,ISTART+I)*EXPI**(AMJ*ALPHA)*W(J)
                    if (quickBeta(IX).eq.2) ROTW(I) = ROTW(I) + EXPI**(AMJP*GAMMA)*AI(ISTART+J,ISTART+I)*EXPI**(AMJ*ALPHA)*W(J)  ! changed VI to AI 16/8/18
                  enddo
                ENDIF
!
 380            CONTINUE 
              enddo ! I=1,ITD loop over degeneracy of each state

              IF (quickBeta(IX) .ne. 0) THEN         !  beta= pi/2, 2*thetaT, pi-2*thetaT
                DO J=1,ITD
                  X(IX) = X(IX) + DCONJG(W(J))*ROTW(J)
                enddo
              ENDIF
            enddo  ! IT=1,n_states
!
          ENDIF    ! End ISO IF.
          X(IX) = X(IX)*FACINV
 400    CONTINUE ! IX=2,NCLASS
!
        DO I=1,NCLASS
          REPS(II,I)=Z     ! fraction of each irrep
          CSYM(II,I)=X(I)  ! characters
        enddo  
!
!  Here we sum over the irreps all above or all below the double group line.
!        if (check1(10)) write(idebug,'("W(",I2,"), SDEG=",F5.2,"; IORDER=",I4)') II,SDEG,IORDER
        N1=1
        N2=NSINGLE
        IF (ISO.eq.1 .and. n_odd) then
          N1=NSINGLE+1
          N2=NCLASS
        endif
        N1=1           ! sum over all irreps
        N2=NCLASS  
        IF (IDEG.EQ.0) THEN   
          NDEG(II)=1                !  (NDEG counts the degeneracy.)
          DO I=N1,N2
            ANS=CZ
            DO J=1,NCLASS
!              ANS = ANS + GROUPT(I,J)*DCONJG(X(J))*ISYMEL(J,3) !  original
              ANS = ANS + DCONJG(GROUPT(I,J))*X(J)*ISYMEL(J,2)  ! ISYMEL(J,2) is the number of elements in each class
            enddo
            REPS(II,I) = ABS(ANS)/IORDER
          enddo
        ELSE IF (IDEG.EQ.1) THEN   !  Degenerate & 1st member 
          NDEG(II)=1
          REPS(II,1)=-D1
          DO I=N1,N2
            CX(I) = X(I)
          enddo
        ELSE IF (IDEG.EQ.2) THEN   !  Degenerate & neither 1st or last
          NDEG(II)=NDEG(II-1)+1
          REPS(II,1)=-D1
          DO I=N1,N2
            CX(I) = CX(I) + X(I)
          enddo 
        ELSE IF (IDEG.EQ.3) THEN   !  Degenerate & last member
          NDEG(II)=NDEG(II-1)+1
          DO I=1,NDEG(II)-1
            NDEG(II-I)=(NDEG(II)+BIT)/SDEG   ! SDEG spin degeneracy if ISO=0, otherwise SDEG=1
          enddo
          NDEG(II)=(NDEG(II)+BIT)/SDEG
          DO I=1,NCLASS
            CX(I) = CX(I) + X(I)
          enddo
          DO I=N1,N2
            ANS=CZ
            DO J=1,NCLASS
              ANS = ANS + ISYMEL(J,2)*GROUPT(I,J)*DCONJG(CX(J)) ! Nclass(I)*character(irrep J, R I) * R.V
            enddo
            REPS(II,I) = ABS(ANS)/(IORDER*SDEG) 
!            REPS(II,I) = ABS(ANS)/(IORDER*SDEG*NDEG(II))   ! NDEG(II) added 14/9/11
          enddo 
        ENDIF
!
 500  CONTINUE  !  II=1,NPRT  loop over eigenvectors
!
      DO II=1,NPRT
        IRREP(II,1) = '    '
        IRREP(II,2) = '    '
        IRREP(II,3) = '    '
        IF (REPS(II,1).EQ.-D1) GOTO 550   ! marked as not the last member of a degenerate set.
        GMAX=REPS(II,1); IMAX=1
        DO I=2,NCLASS
          IF (REPS(II,I).GT.GMAX) THEN
            GMAX = REPS(II,I)
            IMAX = I
          ENDIF
        enddo 
        GMAX2 = Z; IMAX2=0
        DO I=1,NCLASS
          IF (I.EQ.IMAX) GOTO 520
          IF (REPS(II,I).GT.GMAX2) THEN
            GMAX2 = REPS(II,I)
            IMAX2 = I
          endif  
 520      CONTINUE
        enddo
        GMIN = Z
        DO I=1,NCLASS
          IF (I.EQ.IMAX .or. I.eq.IMAX2) GOTO 530
          GMIN = GMIN+REPS(II,I)
 530      CONTINUE
        enddo
        IF (ISO.EQ.0) IRREP(II,1) = GRPR2(IMAX)
        IF (ISO.EQ.1) IRREP(II,1) = GRPR1(IMAX)
        IF (IMAX2.gt.0 .and. GMAX2.GE.0.95) then
          IF (ISO.EQ.0) IRREP(II,2) = GRPR2(IMAX2)
          IF (ISO.EQ.1) IRREP(II,2) = GRPR1(IMAX2)
        endif
! Degenerate states will have two irreps if they are type "b" OR have an odd number of electrons 
! and are type "a". The only degenerate states that have a single irrep (majority) are those of 
! type "a" (even) or "c" (odd).
 !       write(io,'("II=",I3,", GMAX,GMAX2,GMIN=",3F8.3,"NDEG(II),ISYMEL(IMAX,3)=",2I3)')   &
 !                                             II,GMAX,GMAX2,GMIN,NDEG(II),ISYMEL(IMAX,3)
        IF (GMAX.GE.0.95 .AND. GMIN.LT.0.05) IRREP(II,3) = '    '
        IF (GMAX.LT.0.95 .OR.  GMIN.GT.0.05) IRREP(II,3) = '****'
  !      IF (NDEG(II).NE.ISYMEL(IMAX,3)) IRREP(II,3) = '****' 
 550    continue
      enddo ! II=1,NPRT
!
      RETURN
      END SUBROUTINE GROUPL
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ROTB(R,D,beta)
!
!  Called by subroutine: GROUPL
!
!  Returns the matrix R, which contains the Wigner rotation matrix
!  d(j,m',m,-pi/2) on the diagonals for different j values of the
!  basis functions being used. This is then used to perform symmetry
!  rotations on the wavefunctions in subroutine GROUPL.
!
!  For ISO = 0, the L of the basis is used for j.
!      ISO = 1, the J of the basis is used for j.
!
!  Each d(j,m',m,-pi/2) submatrix is of the size ((2j+1)x(2j+1)).
!
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      REAL*8 R(max_N,max_N),D(max_N,max_N)
      PARAMETER(Z=0.0D+00, H=0.5D+00, D1=1.0D+00, D2=2.0D+00)
!
      PI = 4.0D+00*ATAN(D1)
      BETA = -beta   ! (it works)    
      DO I=1,n_matrix
        DO J=1,n_matrix
          R(I,J) = Z
        enddo
      enddo  
!
      IOFF = 0
      DO IT=1,n_states
        ISD=  TS_Bases(1,IT,nelectrons) + 1  !  2S + 1 
        ILD=2*TS_Bases(2,IT,nelectrons) + 1  !  2L + 1
        ITD=ILD*ISD ! total degeneracy
        I2L=ILD-1   ! 2L
        AL= DBLE(TS_Bases(2,IT,nelectrons))     ! (DBLE(ILD)-D1+0.001)/D2
        AS= DBLE(TS_Bases(1,IT,nelectrons))/D2  ! (DBLE(ISD)-D1+0.001)/D2
        IF (ISO.EQ.0) THEN
          CALL DJMPM(D,I2L,BETA)
          DO IS = 1,ISD
            DO IL1 = 1,ILD
              DO IL2 = 1,ILD
                R(IOFF+ILD*(IS-1)+IL1,IOFF+ILD*(IS-1)+IL2)=D(IL1,IL2)
              enddo
            enddo
          enddo  
          IOFF = IOFF + ITD
         ELSE IF (ISO.EQ.1) THEN
          STARTJ = ABS(AL-AS)
          STOPJ  =    (AL+AS)
          JSUM=0
          JJ=INT(STOPJ-STARTJ+0.001)
          DO J1=0,JJ
            JSUM=JSUM + 2*(STARTJ+J1+0.001) + 1
          enddo 
          IF (ITD.NE.JSUM) THEN
            WRITE(IO,'("***FATAL: Something wrong in subroutine ROTB.")')
            STOP
          ENDIF
          DO J1 = 0,JJ      ! JJ,0,-1
            I2J=INT(D2*(STARTJ+J1+0.001))
            CALL DJMPM(D,I2J,BETA)
            DO MJ = 1,I2J+1
              DO MJP = 1,I2J+1
                R(IOFF+MJP,IOFF+MJ) = D(MJP,MJ)
              enddo
            enddo  
            IOFF = IOFF + I2J + 1
          enddo
        ENDIF
      enddo
!
      RETURN
      END SUBROUTINE ROTB
!
!-----------------------------------------------------------------------
!
      SUBROUTINE WriteBmat(BETA,icall)
!
!  Called by subroutine: ROTB
!
!  Prints the transformational matrix R calculated in ROTB for the BETA rotation. 
!  If IDG=28-32 then beta = pi/2, icall =1, R is in VI.
!
!  If IDG=33-35 
!               beta = 2*thetaT;      icall=1, R is in VI.
!               beta = pi - 2*thetaT; icall=2, R is in RI.
!
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)
      REAL*8 beta,PI
      integer icall,i,j,it,isd,ild,itd,istart
!
      PI = 4.0D+00*ATAN(1.0d+0)
      IF (ISO.EQ.0) WRITE(idebug,30) BETA/PI,'L,S '
      IF (ISO.EQ.1) WRITE(idebug,30) BETA/PI,'J,MJ'
 30   FORMAT(//," Wigner rotation matrices for beta=",F7.3,"*pi for wavefunctions in a ",A5," basis.",/,41('-'))
      ISTART=0
      DO IT=1,n_states
        ISD=  TS_Bases(1,IT,nelectrons) + 1  ! NSD(IT)  !  2S + 1 
        ILD=2*TS_Bases(2,IT,nelectrons) + 1  ! NLD(IT)  !  2L + 1
        ITD=ISD*ILD
        IF (IT.NE.1) ISTART=ISTART+(2*TS_Bases(2,IT-1,nelectrons)+1)*(TS_Bases(1,IT-1,nelectrons)+1)
        IF (ISO.EQ.0) WRITE(idebug,35) IT,ISD,ILD,(ISTART+I,I=1,ITD)
 35     FORMAT(" Term:",I4,", 2S+1, 2L+1 = ",2I4,/,21I5)
        IF (ISO.EQ.1) then
          if (.not.N_odd) WRITE(idebug,40) IT,(ILD+ISD)/2-1,ABS(ILD-ISD)/2,(ISTART+I,I=1,ITD)
          if (     N_odd) WRITE(idebug,41) IT,(ILD+ISD)-2,ABS(ILD-ISD),(ISTART+I,I=1,ITD)
 40       FORMAT(" Term: ",I2,",  J=",I2," to ",I2,/,21I5)
 41       FORMAT(" Term: ",I2,",  J=",I2,"/2 to ",I2,"/2",/,21I5)
        endif
        DO I=1,ITD
          if (icall.eq.1) WRITE(idebug,'(1X,21F5.2)') (VI(ISTART+I,ISTART+J),J=1,ITD)
          if (icall.eq.2) WRITE(idebug,'(1X,21F5.2)') (AI(ISTART+I,ISTART+J),J=1,ITD)
        enddo 
      enddo  ! IT=1,n_states
!
      RETURN
      END SUBROUTINE WriteBmat
!
!-----------------------------------------------------------------------
END MODULE f_e_group
!  973 lines