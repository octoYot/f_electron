C----------------------------------------------------------------------
C
      SUBROUTINE DIAGCH(AR,IAR,AI,IAI,N,WR,VR,IVR,VI,IVI,
     +                  WK1,WK2,WK3,IDIAG,IO)
      IMPLICIT REAL*8(A-H,O-Z), INTEGER(I-N)
      INTEGER IFAIL, I, N, J, IVR, IVI, IAR, IAI
      REAL*8 AR(IAR,N),AI(IAI,N),WR(N),VR(IVR,N),VI(IVI,N)
      REAL*8 WK1(N),WK2(N),WK3(N)
      PARAMETER(ZERO=0.0D+00, ONE=1.0D+00)
C
C  EIGENVALUES AND EIGENVECTORS OF A COMPLEX*16 HERMITIAN MATRIX
C
C  Note: VR and VI are used for workspace for the calculation of 
C        the eigenvectors. The actual eigenvectors are returned in 
C        AR and AI.
      CALL HTRIDI(IAR,N,AR,IAI,AI,WR,WK1,WK1,WK2,WK3)
      IF (IDIAG.EQ.1) THEN
        CALL TQL1(N,WR,WK1,IFAIL)
       ELSE
        DO 40 I=1,N
          DO 20 J=1,N
   20       VR(I,J) = ZERO
   40     VR(I,I) = ONE
        CALL TQL2(IVR,N,WR,WK1,VR,IFAIL)
        CALL HTRIBK(IAR,N,AR,IAI,AI,WK2,WK3,IVR,VR,IVI,VI)
C  Wavefunctions are normalised so that the sum of the squares of 
C  the moduli of the elements is equal to 1, and the element with
C  the largest moduli is real.
        DO 100 I=1,N
          SUM=ZERO
          CMAX=ZERO
          DO 80 J=1,N
            SQ=VR(J,I)*VR(J,I) + VI(J,I)*VI(J,I)
            SUM=SUM+SQ
            IF (SQ.LE.CMAX) GOTO 80
            CMAX=SQ
            A=VR(J,I)
            B=VI(J,I)
 80       CONTINUE
          IF (SUM.EQ.ZERO) GOTO 100
          SUM = ONE/SQRT(SUM*CMAX)
          DO 90 J=1,N
            SQ = SUM*(VR(J,I)*A+VI(J,I)*B)
            AI(J,I) = SUM*(VI(J,I)*A-VR(J,I)*B)
            AR(J,I) = SQ
 90       CONTINUE
 100    CONTINUE
      ENDIF
      IF (IFAIL.NE.0) THEN
        WRITE(IO,50) IFAIL
        STOP
      ENDIF
C
 50   FORMAT(//,'***FATAL: More then 30 iterations for eigenvalue ',I3,/,
     +          ' Program terminated in subroutine TQL1/TQL2.',/)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DIAGRS(MATRIX,NM,N,VALUES,FV1,IDIAG,IO)
C
C   THIS SUBROUTINE DIAGONALISES A REAL SYMMETRIC NxN MATRIX, FINDING
C  ALL THE EIGENVALUES AND EIGENVECTORS DEPENDING ON THE OPTIONS TAKEN:
C
C            IDIAG
C              1            ALL EIGENVALUES.
C              2            ALL EIGENVALUES AND EIGENVECTORS.
C
C    NM      DIMENSION OF MATRIX SPECIFIED IN THE DIMENSION
C            STATEMENT OF THE CALLING PROGRAM.
C    N       ORDER OF INPUT MATRIX. N MUST NOT BE GREATER THAN NM.
C    MATRIX  WORKING PRECISION REAL 2-D MATRIX WITH ROW DIMENSION
C            NM AND COLUMN DIMENSION AT LEAST N.
C            ON INPUT ONLY THE FULL LOWER TRIANGLE NEED BE GIVEN.
C            ON OUTPUT THE COLUMNS WILL CONTAIN THE EIGENVECTORS(if IDIAG=2).
C            (IF IDIAG=1, the original matrix will have been changed.)        
C    VALUES  EIGENVALUES IN ASCENDING ORDER.
C
      IMPLICIT REAL*8(A-H,O-Z), INTEGER(I-N)
      INTEGER NM,N,IFAIL,IDIAG
      REAL*8 MATRIX(NM,N),VALUES(N),FV1(N)
C
      if (N.lt.1 .or. n.gt.NM) then
        write(io,'("Invalid dimension of matrix in DIAGRS; N=",i6)') N
        STOP
      endif
      IF (IDIAG.EQ.1) THEN
        CALL TRED1(NM,N,MATRIX,VALUES,FV1,FV1)
        CALL TQL1(N,VALUES,FV1,IFAIL)
       ELSE IF (IDIAG.NE.1) THEN
        CALL TRED2(NM,N,MATRIX,VALUES,FV1)
        CALL TQL2(NM,N,VALUES,FV1,MATRIX,IFAIL)
      ENDIF
      IF (IFAIL.NE.0) THEN
        WRITE(IO,50) IFAIL
        STOP
      ENDIF
C
 50   FORMAT(//,'***FATAL: More then 30 iterations for eigenvalue ',I3,/,
     +          ' Program terminated in subroutine TQL1/TQL2.',/)
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE HTRIDI(IAR,N,AR,IAI,AI,D,E,E2,TAU1,TAU2)
C
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)
      DIMENSION AR(IAR,N),AI(IAI,N),D(N),E(N),E2(N),TAU1(N),TAU2(N)
      PARAMETER(ONE=1.0D+00, ZERO=0.0D+00)
C
      TAU1(N)=ONE
      TAU2(N)=ZERO
C
      DO 100 I = 1, N
 100  D(I) = AR(I,I)
C
C    ********* FOR I=N STEP -1 UNTIL 1 DO --  ******************
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GOTO 130
C   *********  SCALE ROW (ALGOL TOL THEN NOT NEEDED)  **********
         DO 120 K = 1, L
 120     SCALE = SCALE + ABS(AR(I,K)) + ABS(AI(I,K))
C
         IF (SCALE .NE. ZERO) GOTO 140
         TAU1(L) = ZERO  ! ONE  correction MJR 13/12/15
         TAU2(L) = ZERO
 130     E(I) = ZERO
         E2(I) = ZERO
         GOTO 290
C
 140     DO 150 K = 1, L
            AR(I,K) = AR(I,K) / SCALE
            AI(I,K) = AI(I,K) / SCALE
            H = H + AR(I,K) * AR(I,K) + AI(I,K) * AI(I,K)
 150     CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G = SQRT(H)
         E(I) = SCALE * G
         F = CDABS(DCMPLX(AR(I,L),AI(I,L)))    ! trap MJR 13/12/15  must use DCMPLX, not CMPLX
                                               !                             CDABS, not CABS   
C     ********  FORM NEXT DIAGONAL ELEMENT OF MATRIX T  **********
         IF (F .EQ. ZERO) GOTO 160
         TAU1(L) = (AI(I,L) * TAU2(I) - AR(I,L) * TAU1(I)) / F
         SI = (AR(I,L) * TAU2(I) + AI(I,L) * TAU1(I)) / F
         H = H + F * G
         G = ONE + G / F
         AR(I,L) = G * AR(I,L)
         AI(I,L) = G * AI(I,L)
         IF (L .EQ. 1) GOTO 270
         GOTO 170
 160     TAU1(L) = -TAU1(I)
         SI = TAU2(I)
         AR(I,L) = G
 170     F = ZERO
C
         DO 240 J = 1, L
            G = ZERO
            GI = ZERO
C      *************  FORM ELEMENT OF A*U  *************
            DO 180 K = 1, J
               G = G + AR(J,K) * AR(I,K) + AI(J,K) * AI(I,K)
               GI = GI - AR(J,K) * AI(I,K) + AI(J,K) * AR(I,K)
 180        CONTINUE
C
            JP1 = J + 1
            IF (L .LT. JP1) GOTO 220
C
            DO 200 K = JP1, L
               G = G + AR(K,J) * AR(I,K) - AI(K,J) * AI(I,K)
               GI = GI - AR(K,J) * AI(I,K) - AI(K,J) * AR(I,K)
 200        CONTINUE
C      ***********  FORM ELEMENT OF P  **************
 220        E(J) = G / H
            TAU2(J) = GI / H
            F = F + E(J) * AR(I,J) - TAU2(J) * AI(I,J)
 240     CONTINUE
C
         HH = F / (H + H)
C     ************  FORM REDUCED A  ******************
         DO 260 J = 1, L
            F = AR(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -AI(I,J)
            GI = TAU2(J) - HH * FI
            TAU2(J) = -GI
C
            DO 260 K = 1, J
              AR(J,K)=AR(J,K)-F*E(K)-G*AR(I,K)+FI*TAU2(K)+GI*AI(I,K)
              AI(J,K)=AI(J,K)-F*TAU2(K)-G*AI(I,K)-FI*E(K)-GI*AR(I,K)
 260     CONTINUE
C
 270     DO 280 K = 1, L
            AR(I,K) = SCALE * AR(I,K)
            AI(I,K) = SCALE * AI(I,K)
 280     CONTINUE
C
         TAU2(L) = -SI
 290     HH = D(I)
         D(I) = AR(I,I)
         AR(I,I) = HH
         AI(I,I) = SCALE * SCALE * H
 300  CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE HTRIBK(IAR,N,AR,IAI,AI,TAU1,TAU2,IZR,ZR,IZI,ZI)
C
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)
      DIMENSION AR(IAR,N),AI(IAI,N),TAU1(N),TAU2(N),ZR(IZR,N),ZI(IZI,N)
      PARAMETER(ZERO=0.0D+00)
C
C   *******  TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
C            TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
C            TRIDIAGONAL MATRIX.  *****************
      DO 50 K = 1, N
C
         DO 50 J = 1, N
            ZI(K,J) = -ZR(K,J) * TAU2(K)
            ZR(K,J) =  ZR(K,J) * TAU1(K)
 50   CONTINUE
C
      IF (N .EQ. 1) GOTO 200
C   ********  RECOVER AND APPLY THE HOUSEHOLDER MATRICES.  ***********
      DO 140 I = 2, N
         L = I - 1
         H = AI(I,I)
         IF (H .EQ. ZERO) GOTO 140
C
         DO 130 J = 1, N
            S = ZERO
            SI = ZERO
C
            DO 110 K = 1, L
               S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
               SI = SI + AR(I,K) *ZI(K,J) + AI(I,K) * ZR(K,J)
 110        CONTINUE
C
            S = S/ H
            SI = SI / H
C
            DO 120 K = 1, L
               ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
               ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
 120        CONTINUE
C
 130     CONTINUE
C
 140  CONTINUE
C
 200  RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
      IMPLICIT REAL*8(A-H,O-Z), INTEGER(I-N)
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 A(NM,N),D(N),E(N),E2(N),F,G,H,SCALE,SQRT,ABS,SIGN
      PARAMETER(ZERO=0.0D+00)
C
      DO 100 I=1,N
 100  D(I)=A(I,I)
C*************** FOR I=N STEP -1 UNTIL 1 DO *************************
      DO 300 II=1,N
         I=N+1-II
         L=I-1
         H=ZERO
         SCALE=ZERO
         IF (L.LT.1) GOTO 130
C************** SCALE ROW (ALGOL TOL THEN NOT NEEDED) ****************
         DO 120 K=1,L
 120     SCALE=SCALE+ABS(A(I,K))
C
         IF (SCALE.NE.ZERO) GOTO 140
 130     E(I)=ZERO
         E2(I)=ZERO
         GOTO 290
C
 140     DO 150 K=1,L
            A(I,K)=A(I,K)/SCALE
            H=H+A(I,K)*A(I,K)
 150     CONTINUE
C
         E2(I)=SCALE*SCALE*H
         F=A(I,L)
         G=-SIGN(SQRT(H),F)
         E(I)=SCALE*G
         H=H-F*G
         A(I,L)=F-G
         IF (L.EQ.1) GOTO 270
         F=ZERO
C
         DO 240 J=1,L
            G=ZERO
C***************** FORM ELEMENT OF A*U *******************************
            DO 180 K=1,J
 180        G=G+A(J,K)*A(I,K)
C
            JP1=J+1
            IF (L.LT.JP1) GOTO 220
C
            DO 200 K=JP1,L
 200        G=G+A(K,J)*A(I,K)
C**************** FORM ELEMENT OF P **********************************
 220        E(J)=G/H
            F=F+E(J)*A(I,J)
 240     CONTINUE
C
         H=F/(H+H)
C**************** FORM REDUCED A *************************************
         DO 260 J=1,L
            F=A(I,J)
            G=E(J)-H*F
            E(J)=G
C
            DO 260 K=1,J
               A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
 260     CONTINUE
C
 270     DO 280 K=1,L
 280     A(I,K)=SCALE*A(I,K)
C
 290     H=D(I)
         D(I)=A(I,I)
         A(I,I)=H
 300  CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE TRED2(NM,N,Z,D,E)
C
C   This subroutine was adapted from: B.T.SMITH ET.AL., "MATRIX EIGENSYSTEM
C  ROUTINES- EISPACK GUIDE", SPRINGER-VERLAG, (1974).
C  Original call:
C      SUBROUTINE TRED2(NM,N,A,D,E,Z,IZ)
C      REAL*8 A(NM,N),Z(IZ,N)
C  Modified to always overwrite the original matrix.
C
      IMPLICIT REAL*8(A-H,O-Z), INTEGER(I-N)
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 Z(NM,N),D(N),E(N),F,G,H,HH,SCALE,SQRT,ABS,SIGN
      PARAMETER(ZERO=0.0D+00, ONE=1.0D+00)
C
C     DO 100 I=1,N
C       DO 100 J=1,I
C         Z(I,J)=A(I,J)
C100  CONTINUE
C
      IF (N.EQ.1) GOTO 320
C ********** FOR I=N STEP -1 UNTIL 2 DO __ **************************
      DO 300 II=2,N
         I=N+2-II
         L=I-1
         H=ZERO
         SCALE=ZERO
         IF (L.LT.2) GOTO 130
C ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) ********************
         DO 120 K=1,L
 120     SCALE=SCALE+ABS(Z(I,K))
C
         IF (SCALE.NE.ZERO) GOTO 140
 130     E(I)=Z(I,L)
         GOTO 290
C
 140     DO 150 K=1,L
            Z(I,K)=Z(I,K)/SCALE
            H=H+Z(I,K)*Z(I,K)
 150     CONTINUE
C
         F=Z(I,L)
         G=-SIGN(SQRT(H),F)
         E(I)=SCALE*G
         H=H-F*G
         Z(I,L)=F-G
         F=ZERO
C
         DO 240 J=1,L
            Z(J,I)=Z(I,J)/(SCALE*H)
            G=ZERO
C ********** FORM ELEMENT OF A*U ***************************************
            DO 180 K=1,J
 180        G=G+Z(J,K)*Z(I,K)
C
            JP1=J+1
            IF (L.LT.JP1) GOTO 220
C
            DO 200 K=JP1,L
 200        G=G+Z(K,J)*Z(I,K)
C ********** FORM ELEMENT OF P *****************************************
 220        E(J)=G/H
            F=F+E(J)*Z(I,J)
 240     CONTINUE
C
         HH=F/(H+H)
C ********** FORM REDUCED A ********************************************
         DO 260 J=1,L
            F=Z(I,J)
            G=E(J)-HH*F
            E(J)=G
C
            DO 260 K=1,J
               Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
 260     CONTINUE
C
         DO 280 K=1,L
 280     Z(I,K)=SCALE*Z(I,K)
C
 290     D(I)=H
 300  CONTINUE
C
 320  D(1)=ZERO
      E(1)=ZERO
C ********** ACCUMULATION OF TRANSFORMED MATRICES **********************
      DO 500 I=1,N
         L=I-1
         IF (D(I).EQ.ZERO) GOTO 380
C
         DO 360 J=1,L
            G=ZERO
C
            DO 340 K=1,L
 340        G=G+Z(I,K)*Z(K,J)
C
            DO 360 K=1,L
               Z(K,J)=Z(K,J)-G*Z(K,I)
 360     CONTINUE
C
 380     D(I)=Z(I,I)
         Z(I,I)=ONE
         IF (L.LT.1) GOTO 500
C
         DO 400 J=1,L
            Z(I,J)=ZERO
            Z(J,I)=ZERO
 400     CONTINUE
C
 500  CONTINUE
C
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE TQL1(N,D,E,IERR)
C
      IMPLICIT REAL*8(A-H,O-Z), INTEGER(I-N)
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL*8 D(N),E(N),B,C,F,G,P,R,S,MACHEP,SQRT,ABS,SIGN,X02AAF
      PARAMETER(ZERO=0.0D+00, ONE=1.0D+00)
C************** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C             THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC. *****
      MACHEP = X02AAF(ONE)
C
      IERR=0
      IF (N.EQ.1) GOTO 1001
C
      DO 100 I=2,N
 100  E(I-1)=E(I)
C
      F=ZERO
      B=ZERO
      E(N)=ZERO
C
      DO 290 L=1,N
         J=0
         H=MACHEP*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B=H
C************** LOOK FOR SMALL SUB-DIAGONAL ELEMENT ******************
         DO 110 M=L,N
            IF (ABS(E(M)).LE.B) GOTO 120
C************** E(N) IS ALWAYS ZERO SO THERE IS NO EXIT
C           THROUGH THE BOTTEM OF THE LOOP ***************************
 110     CONTINUE
C
 120     IF (M.EQ.L) GOTO 210
 130     IF (J.EQ.30) GOTO 1000
         J=J+1
C************** FORM SHIFT *******************************************
         L1=L+1
         G=D(L)
         P=(D(L1)-G)/(2.0*E(L))
         R=SQRT(P*P+ONE)
         D(L)=E(L)/(P+SIGN(R,P))
         H=G-D(L)
C
         DO 140 I=L1,N
 140     D(I)=D(I)-H
C
         F=F+H
C************* QL TRANSFORM ******************************************
         P=D(M)
         C=ONE
         S=ZERO
         MML=M-L
C************* FOR I=M-1 STEP -1 UNTIL L DO __ ***********************
         DO 200 II=1,MML
            I=M-II
            G=C*E(I)
            H=C*P
            IF (ABS(P).LT.ABS(E(I))) GOTO 150
            C=E(I)/P
            R=SQRT(C*C+ONE)
            E(I+1)=S*P*R
            S=C/R
            C=ONE/R
            GOTO 160
 150        C=P/E(I)
            R=SQRT(C*C+ONE)
            E(I+1)=S*E(I)*R
            S=ONE/R
            C=C*S
 160        P=C*D(I)-S*G
            D(I+1)=H+S*(C*G+S*D(I))
 200     CONTINUE
C
         E(L)=S*P
         D(L)=C*P
         IF (ABS(E(L)).GT.B) GOTO 130
 210     P=D(L)+F
C***************** ORDER EIGENVALUES *********************************
         IF (L.EQ.1) GOTO 250
C***************** FOR I=L STEP -1 UNTIL 2 DO __ *********************
         DO 230 II=2,L
            I=L+2-II
            IF (P.GE.D(I-1)) GOTO 270
            D(I)=D(I-1)
 230     CONTINUE
C
 250     I=1
 270     D(I)=P
 290  CONTINUE
C
      GOTO 1001
C**************** SET ERROR IF NO CONVERGENCE TO AN EIGENVALUE
C                      AFTER 30 ITERATIONS *****************************
 1000 IERR=L
 1001 RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      IMPLICIT REAL*8(A-H,O-Z), INTEGER(I-N)
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N),B,C,F,G,H,P,R,S
      REAL*8 MACHEP,SQRT,ABS,SIGN,X02AAF
      PARAMETER(ZERO=0.0D+00, ONE=1.0D+00, TWO=2.0D+00)
C
C ************* MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C        THE RELATIVE PRECISION OF FLOATING POINT ARITHMATIC. **********
      MACHEP = X02AAF(ONE)
C
      IERR=0
      IF (N.EQ.1) GOTO 1001
C
      DO 100 I=2,N
 100  E(I-1)=E(I)
C
      F=ZERO
      B=ZERO
      E(N)=ZERO
C
      DO 240 L=1,N
         J=0
         H=MACHEP*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B=H
C ************* LOOK FOR SMALL SUB-DIAGONAL ELEMENT ********************
         DO 110 M=L,N
            IF (ABS(E(M)).LE.B) GOTO 120
C ************* E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C             THROUGH THE BOTTEM OF THE LOOP. **************************
 110  CONTINUE
C
 120  IF (M.EQ.L) GOTO 220
 130  IF (J.EQ.30) GOTO 1000
      J=J+1
C ************* FORM SHIFT *********************************************
      L1=L+1
      G=D(L)
      P=(D(L1)-G)/(TWO*E(L))
      R=SQRT(P*P+ONE)
      D(L)=E(L)/(P+SIGN(R,P))
      H=G-D(L)
C
      DO 140 I=L1,N
 140  D(I)=D(I)-H
C
      F=F+H
C ************ QL TRANSFORMATIONS **************************************
      P=D(M)
      C=ONE
      S=ZERO
      MML=M-L
C ************* FOR I=M-1 STEP -1 UNTIL L DO __ ************************
      DO 200 II=1,MML
         I=M-II
         G=C*E(I)
         H=C*P
         IF (ABS(P).LT.ABS(E(I))) GOTO 150
         C=E(I)/P
         R=SQRT(C*C+ONE)
         E(I+1)=S*P*R
         S=C/R
         C=ONE/R
         GOTO 160
 150     C=P/E(I)
         R=SQRT(C*C+ONE)
         E(I+1)=S*E(I)*R
         S=ONE/R
         C=C*S
 160     P=C*D(I)-S*G
         D(I+1)=H+S*(C*G+S*D(I))
C ************** FORM VECTOR *******************************************
         DO 180 K=1,N
            H=Z(K,I+1)
            Z(K,I+1)=S*Z(K,I)+C*H
            Z(K,I)=C*Z(K,I)-S*H
 180     CONTINUE
C
 200  CONTINUE
C
      E(L)=S*P
      D(L)=C*P
      IF (ABS(E(L)).GT.B) GOTO 130
 220  D(L)=D(L)+F
 240  CONTINUE
C ************** ORDER EIGENVALUES AND EIGENVECTORS ********************
      DO 300 II=2,N
         I=II-1
         K=I
         P=D(I)
C
         DO 260 J=II,N
            IF (D(J).GE.P) GOTO 260
            K=J
            P=D(J)
 260     CONTINUE
C
         IF (K.EQ.I) GOTO 300
         D(K)=D(I)
         D(I)=P
C
         DO 280 J=1,N
            P=Z(J,I)
            Z(J,I)=Z(J,K)
            Z(J,K)=P
 280     CONTINUE
C
 300  CONTINUE
C
      GOTO 1001
C ****************** SET ERROR __ NO CONVERGENCE TO AN EIGENVALUE
C                    AFTER 30 ITERATIONS. ******************************
 1000 IERR=L
 1001 RETURN
      END
C
C----------------------------------------------------------------------
C
      REAL*8 FUNCTION X02AAF(X)
      REAL*8 XX,X,Z
C     NAG COPYRIGHT 1975
C     IBM DOUBLE PRECISION
C     MARK 4.5 RELEASE
C     * MACHEPS *
C     RETURNS THE VALUE MACHEPS WHERE MACHEPS IS THE SMALLEST POSITIVE
C     NUMBER SUCH THAT 1.0 + EPS > 1.0
C     THE X PARAMETER IS NOT USED
C     FOR ICL 1900
C     X02AAF = 2.0**(-37.0)
C     FOR IBM 370
C     X02AAF = 2.0D0**(-52.0D0)
C     SET IN HEX FOR ACCURACY
C     FOR BURROUGHS DBL PREC
C     X02AAF = 8**(-25)
C     FOR RSC VAX/750 
C     X02AAF = 2**(-55)
C     SET IN HEX FOR ACCURACY:
C     DATA Z/'0000000000002500'X/
C
C  For SUN Workstation, IEEE Standard 754, X02AAF= 2**(-52)
C  Set in HEX:
C     DATA Z/'3CB0000000000000'X/

      XX = X
      Z = 2.220447D-16  ! = 1/(2^52)

      X02AAF=Z
      RETURN
      END
C-----------------------------------------------------------------------
!  707 lines