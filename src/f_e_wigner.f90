MODULE f_e_wigner

CONTAINS

!----------------------------------------------------------------------
!
      REAL*8 FUNCTION DEL(U,V,W)
!  DEL is required in the calculation of the Wigner 6-j symbol.
      IMPLICIT REAL*8(A-H,O-Z), INTEGER*4(I-N)
      INTEGER*4 J(3)
      SP = 1.0D+00
      PP = 1.0D+00
      J(1) = U+V-W+0.02D+00
      J(2) = U-V+W+0.02D+00
      J(3) = V+W-U+0.02D+00
      J4 = U+V+W+1.02D+00
      DO 1 M = 1,3
        N = J(M)
        IF (N.EQ.0) GO TO 2
        DO 3 I = 1,N
          SP=SP*DBLE(I)
    3   CONTINUE
    2   CONTINUE
    1 CONTINUE
      DO 4 K = 1,J4
        PP = PP*DBLE(K)
    4 CONTINUE
      DEL = SQRT(SP/PP)
      RETURN
      END FUNCTION DEL
!
!----------------------------------------------------------------------
!
      REAL*8 FUNCTION WIG3J (A,B,C,X,Y,Z)
!   Selection rules are built into WIG3J.
!   Determined from corresponding Clebsch-Gordon coefficient using
!   Racah's Equation. Ref: G.Racah, Phys.Rev.,62,440,(1942).
      IMPLICIT REAL*8(A-H,O-Z), INTEGER*4(I-N)
      REAL*8 D(6)
      INTEGER*4 K(9)
      Q = 0.001D+00
      WIG3J = 0.0D+00
      HA = ABS(X+Y+Z)
      HB = A+B-C+0.02D+00
      HC = B+C-A+0.02D+00
      HD = C+A-B+0.02D+00
      HE=A+X+0.02D+00
      HF=A-X+0.02D+00
      HG=B+Y+0.02D+00
      HH=B-Y+0.02D+00
      HI=C+Z+0.02D+00
      HJ=C-Z+0.02D+00
      HK=A+B+C+1.02D+00
      IF ((HA.GT.Q).OR.(HB.LT.Q).OR.(HC.LT.Q).OR.(HD.LT.Q).OR.(HE.LT.Q)  &
      .OR.(HF.LT.Q).OR.(HG.LT.Q).OR.(HH.LT.Q).OR.(HI.LT.Q).OR.(HJ.LT.Q)  &
      .OR.(HK.LT.Q))  GO TO 13
      SN = 1.0D+00
      RN = 1.0D+00
      SUM = 0.0D+00
      V = 0.0D+00
      SM = A + B + C
      GO TO 2
    1 CONTINUE
      V = V+1.0D+00
      IF (V.GT.SM) GO TO 10
    2 CONTINUE
      PN = 1.0D+00
      D(1) = A+B-C-V+0.02D+00
      D(2) = A-X-V+0.02D+00
      D(3) = B+Y-V+0.02D+00
      D(4) = C-B+X+V+0.02D+00
      D(5) = C-A-Y+V+0.02D+00
      D(6) = V +0.02
      IF ((D(1).LT.Q).OR.(D(2).LT.Q).OR.(D(3).LT.Q).OR.(D(4).LT.Q).OR.  &
          (D(5).LT.Q).OR.(D(6).LT.Q))  GO TO 1
      DO 3 M = 1,6
        N = INT(D(M))
        IF (N.EQ.0) GO TO 14
        DO 4 I=1,N
          PN = PN*DBLE(I)
    4   CONTINUE
   14   CONTINUE
    3 CONTINUE
      PX=0.5D+00*(V+2.02D+00)
      FACT = -1.0D+00
      IF((PX-AINT(PX)-0.3D+00).LT.0.0D+00)   FACT =1.0D+00
      SUM = SUM + FACT/PN
      GO TO 1
   10 CONTINUE
      K(1)=HD
      K(2)=HC
      K(3)=HB
      K(4)=HE
      K(5)=HF
      K(6)=HG
      K(7)=HH
      K(8)=HI
      K(9)=HJ
      K10=HK
      DO 5 M=1,9
        N=K(M)
        IF (N.EQ.0) GO TO 16
        DO 6 I=1,N
          RN = RN*DBLE(I)
    6   CONTINUE
   16 CONTINUE
    5 CONTINUE
      DO 7 I = 1,K10
        SN = SN*DBLE(I)
    7 CONTINUE
      QX = 0.5D+00*(A-B-Z+100.02D+00)
      PHF = -1.0D+00
      IF ((QX-AINT(QX)-0.3D+00).LT.0.0D+00)   PHF = 1.0D+00
      WIG3J = PHF*SUM*SQRT(RN/SN)
   13 CONTINUE
      RETURN
      END FUNCTION WIG3J
!
!----------------------------------------------------------------------
!
      REAL*8 FUNCTION WIG6J (A,B,C,D,E,F)
!  Selection rules are built into WIG6J.
!  Calculated from an equation given by Judd.
!  Ref: B.R.Judd, "Operator Techniques in Atomic Spectroscopy",
!       McGraw Hill, (1963); pg. 57.
      IMPLICIT REAL*8(A-H,O-Z), INTEGER*4(I-N)
      REAL*8 G(7)
      Q = 0.001D+00
      WIG6J = 0.0D+00
      
      HA = A+B-C+0.02D+00
      HB = B+C-A+0.02D+00
      HC = C+A-B+0.02D+00
      HD = C+D-E+0.02D+00
      HE = D+E-C+0.02D+00
      HF = E+C-D+0.02D+00
      HG = A+E-F+0.02D+00
      HI = E+F-A+0.02D+00
      HJ = F+A-E+0.02D+00
      HK = D+B-F+0.02D+00
      HL = B+F-D+0.02D+00
      HM = F+D-B+0.02D+00
      IF ((HA.LT.Q).OR.(HB.LT.Q).OR.(HC.LT.Q).OR.(HD.LT.Q).OR.(HE.LT.Q)  &
      .OR.(HF.LT.Q).OR.(HG.LT.Q).OR.(HI.LT.Q).OR.(HJ.LT.Q).OR.(HK.LT.Q)  &
      .OR.(HL.LT.Q).OR.(HM.LT.Q))  GO TO 13
      Z = 0.0D+00
      SM = A+B+C+D+E+F
      SUM = 0.0D+00
      GO TO 2
    1 CONTINUE
      Z = Z + 1.0D+00
      IF(Z.GT.SM) GO TO 10
    2 CONTINUE
      PN = 1.0D+00
      QN = 1.0D+00
      G(1) = Z-A-B-C+0.02D+00
      G(2) = Z-A-E-F+0.02D+00
      G(3) = Z-D-B-F+0.02D+00
      G(4) = Z-D-E-C+0.02D+00
      G(5) = A+B+D+E-Z+0.02D+00
      G(6) = B+C+E+F-Z+0.02D+00
      G(7) = C+A+F+D-Z+0.02D+00
      IF ((G(1).LT.Q).OR.(G(2).LT.Q).OR.(G(3).LT.Q).OR.(G(4).LT.Q).OR.  &
          (G(5).LT.Q).OR.(G(6).LT.Q).OR.(G(7).LT.Q))  GO TO 1
      DO 3 M = 1,7
        N = INT(G(M))
        IF (N.EQ.0) GO TO 14
        DO 4 I = 1,N
          PN = PN*DBLE(I)
    4   CONTINUE
   14   CONTINUE
    3 CONTINUE
      JR = Z+1.02D+00
      DO 17 J = 1,JR
        QN = QN*DBLE(J)
   17 CONTINUE
      PX = 0.5D+00*(Z+2.02D+00)
      FACT = -1.0D+00
      IF ((PX-AINT(PX)-0.3D+00).LT.0.0D+00)   FACT =1.0D+00
      SUM = SUM+FACT*QN/PN
      GO TO 1
   10 CONTINUE
      WIG6J = DEL(A,B,C)*DEL(A,E,F)*DEL(D,B,F)*DEL(D,E,C)*SUM
   13 CONTINUE
      RETURN
      END 
!
!----------------------------------------------------------------------
END MODULE f_e_wigner
!  190 lines
