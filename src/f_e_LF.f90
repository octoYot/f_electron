MODULE f_e_LF 

USE f_e_data
USE f_e_parameters
!
!  This module does ligand field manipulations.
!

IMPLICIT NONE

! PRIVATE

PUBLIC 

CONTAINS
!
!-----------------------------------------------------------------------
!
      INTEGER FUNCTION Bkq_index(k,q)
 
!  Returns the index of BkqR and BkqI array appropriate to particular values of k & q.
!  k | 0  2  2  2  4  4  4  4  4  6  6  6  6  6  6  6  
!  q | 0  0  1  2  0  1  2  3  4  0  1  2  3  4  5  6
!  i | 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16
!
      IMPLICIT none
      integer k,q,i
      if (k.eq.0) i=1
      if (k.eq.2) i=2
      if (k.eq.4) i=5
      if (k.eq.6) i=10
      Bkq_index=i+q
      if (Bkq_index.lt.1.or.Bkq_index.gt.16) then      
        write(*,'("***FATAL: Invalid index (",I4,") found for k,q=",2I4)') Bkq_index,k,q
        stop
      endif
      Return
      end function Bkq_index
!      
!-----------------------------------------------------------------------
!
      SUBROUTINE zeroE(e,n)
!
!  Makes the eignevalues relative to the lowest value (if option(1) true)
!  EAVE added here rather than diagonal elements of matrix, to avoid the matrix diagonalisation.
!  Makes it very quick to fit a value of EAVE, if this is the only parameter being fitted.
!
      IMPLICIT NONE
      integer i,n
      real*8 e(n),e0
      
      do i=1,n; e(i)=e(i)+EAVE; enddo
      if (Option(1)) then
        e0=e(1)
        do i=2,n
          if (e(i).lt.e0) e0=e(i)
        enddo
        do i=1,n
          e(i)=e(i)-e0
        enddo
      endif  

      return
!
      end subroutine zeroE
!
!-----------------------------------------------------------------------
!
      SUBROUTINE AOMmatrixD()
! Calculates the AOM matrix for d-orbitals, and the equivalent crystal field parameters.
      IMPLICIT none
      INTEGER*4 iL,i,j,k,q
      REAL*8 T,P,Ps,Pi,fv1(5),M(5,5),M1(5,5),F(5,5),D0,D1,D2,D3,D4,D5,D6,D7,D10,R2,R3,E1,sum1,sum2
      PARAMETER(D0=0.0D+00,D1=1.0D+00,D2=2.0D+00,D3=3.0D+00,D4=4.0D+00,D5=5.0D+00,D6=6.0D+00,D7=7.0D+00)
      PARAMETER(D10=10.0D+00)
      
      Pi=D4*ATAN(D1)
      R2=Sqrt(D2); R3=Sqrt(D3)
      
      if (LFtype.eq."AOM ") then

        M=D0; M1=D0; F=D0

        do iL=1,NLIGANDS
          T=Theta(iL)*Pi/180.d+00
          P=Phi(iL)*Pi/180.d+00
          Ps=Chi(iL)*Pi/180.d+00
!  The <d|V|L> elements from table 2b Schaeffer,Struct&Bond,5,68,(1968). 
!  Note: "Standard" ordering of orbitals = z2, yz, xz, xy, x2-y2.
          F(1,1)=+(D1 + D3*Cos(D2*T))/D4
          F(1,2)=-R3/D2*Sin(D2*T)*Cos(Ps)  ! PiX
          F(1,3)=+R3/D2*Sin(D2*T)*Sin(Ps)  ! PiY
          F(2,1)= R3/D2*Sin(D2*T)*Sin(P)
          F(2,2)= Cos(T)*Cos(P)*Sin(Ps) + Cos(D2*T)*Sin(P)*Cos(Ps)
          F(2,3)= Cos(T)*Cos(P)*Cos(Ps) - Cos(D2*T)*Sin(P)*Sin(Ps)
          F(3,1)= R3/D2*Sin(D2*T)*Cos(P)
          F(3,2)=-Cos(T)*Sin(P)*Sin(Ps) + Cos(D2*T)*Cos(P)*Cos(Ps)
          F(3,3)=-Cos(T)*Sin(P)*Cos(Ps) - Cos(D2*T)*Cos(P)*Sin(Ps)
          F(4,1)= R3/D4*(D1-Cos(D2*T))*Sin(D2*P)
          F(4,2)= Sin(T)*Cos(D2*P)*Sin(Ps) + Sin(D2*T)*Sin(D2*P)*Cos(Ps)/D2
          F(4,3)= Sin(T)*Cos(D2*P)*Cos(Ps) - Sin(D2*T)*Sin(D2*P)*Sin(Ps)/D2
          F(5,1)= R3/D4*(D1-Cos(D2*T))*Cos(D2*P)
          F(5,2)=-Sin(T)*Sin(D2*P)*Sin(Ps) + Sin(D2*T)*Cos(D2*P)*Cos(Ps)/D2
          F(5,3)=-Sin(T)*Sin(D2*P)*Cos(Ps) - Sin(D2*T)*Cos(D2*P)*Sin(Ps)/D2
!        
          do i=1,5
            do j=1,5
              M(i,j) = M(i,j) + eSigma(iL)*F(i,1)*F(j,1) + ePiX(iL)*F(i,2)*F(j,2) + ePiY(iL)*F(i,3)*F(j,3)
              M1(i,j)= M(i,j) ! M1 gets destroyed below
!     write(io,'("eS,ePx,ePy,T,P,Ps=",2(3F8.1,5x,3F8.1))') eSigma(iL),ePiX(iL),ePiY(iL), Theta(iL),Phi(iL),Chi(iL)   
            enddo
          enddo
        enddo  ! iL=1,NLIGANDS       
 
        if (outp(3).or.outp(9).or.check1(8)) then
          OrbEngs=D0; fv1=D0
          call diagrs(M1,5,5,OrbEngs,fv1,1,io)
          E1=EAVE; EAVE=D0; call zeroE(OrbEngs,5); EAVE=E1
          if (outp(3)) then
            write(io,'(/," AOM matrix in terms of the 5x5 real 5-orbital matrix")') 
            write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>     orbital energies")') 
            write(io,'(2X," <sig|",5F12.2,5X,F12.2)') (M(1,j),j=1,5),OrbEngs(5)
            write(io,'(2X," <piS|",5F12.2,5X,F12.2)') (M(2,j),j=1,5),OrbEngs(4)
            write(io,'(2X," <piC|",5F12.2,5X,F12.2)') (M(3,j),j=1,5),OrbEngs(3)
            write(io,'(2X,"<delS|",5F12.2,5X,F12.2)') (M(4,j),j=1,5),OrbEngs(2)
            write(io,'(2X,"<delC|",5F12.2,5X,F12.2)') (M(5,j),j=1,5),OrbEngs(1)
          endif
          if (outp(9)) then
            write(io,'(/,2X,"Energies:",5F12.2)') (OrbEngs(j),j=1,5)
            write(io,'(2X,"  |z^2> :",5F12.2)') (M1(1,j),j=1,5)
            write(io,'(2X,"  |yz>  :",5F12.2)') (M1(2,j),j=1,5)
            write(io,'(2X,"  |xz>  :",5F12.2)') (M1(3,j),j=1,5)
            write(io,'(2X,"  |xy>  :",5F12.2)') (M1(4,j),j=1,5)
            write(io,'(2X,"|x2-y2)>:",5F12.2)') (M1(5,j),j=1,5)
          endif
          if (check1(8)) then
            write(idebug,'(/," AOM matrix in terms of the 5x5 real 5-orbital matrix")') 
            write(idebug,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>     orbital energies")') 
            write(idebug,'(2X," <sig|",5F12.2,5X,F12.2)') (M(1,j),j=1,5),OrbEngs(5)
            write(idebug,'(2X," <piS|",5F12.2,5X,F12.2)') (M(2,j),j=1,5),OrbEngs(4)
            write(idebug,'(2X," <piC|",5F12.2,5X,F12.2)') (M(3,j),j=1,5),OrbEngs(3)
            write(idebug,'(2X,"<delS|",5F12.2,5X,F12.2)') (M(4,j),j=1,5),OrbEngs(2)
            write(idebug,'(2X,"<delC|",5F12.2,5X,F12.2)') (M(5,j),j=1,5),OrbEngs(1)
          endif
        
          if (FitType.eq.6) then
            do i=1,5; do j=1,i; LF1eMat(i,j)=M(i,j); enddo; enddo
            sum1=(EXP_E(1)+EXP_E(3)+EXP_E(6)+EXP_E(10)+EXP_E(15))/D5  ! Adjust traces of matrices to be the same
            sum2=(LF1eMat(1,1)+LF1eMat(2,2)+LF1eMat(3,3)+LF1eMat(4,4)+LF1eMat(5,5))/D5
            do i=1,5; LF1eMat(i,i)=LF1eMat(i,i)-sum2+sum1; enddo
          endif
        Endif

      else if (LFtype.eq."LF1E") then
        do i=1,5; do j=1,i;  M(i,j)=LF1eMat(i,j); enddo; enddo  ! Fill lower triangle only
      endif
!
!  Relationships below are given by Gerloch & McMeeking, Dalton, 2443, (1975). Table 2.
!  These are the Ckq coefficents of the ligand field expanded as spherical harmonics.
! B00 
      BKQR(1)=D2/D5*Sqrt(Pi)*(M(1,1) + M(2,2) + M(3,3) + M(4,4) + M(5,5))
      BKQI(1)=D0
! B20     
      BKQR(2)= Sqrt(Pi/D5)*(D2*M(1,1) + M(2,2) + M(3,3) - D2*M(4,4) -D2*M(5,5))
      BKQI(2)=D0
!  Note: "Standard" ordering of orbitals = z2, yz, xz, xy, x2-y2.
! B21      
      BKQR(3)=-Sqrt(D4*Pi/D5)*( R3/R2*(M(4,2) + M(5,3)) + M(3,1)/R2)  ! (M(2,4) + M(3,5)) + M(3,1)/R2)
      BKQI(3)= Sqrt(D4*Pi/D5)*(-R3/R2*(M(5,2) - M(4,3)) + M(2,1)/R2)  ! (M(2,5) - M(3,4)) + M(2,1)/R2)
! B22      
      BKQR(4)=-Sqrt(D4*Pi/D5)*( R2*M(5,1) + R3/(D2*R2)*(M(2,2) - M(3,3)))
      BKQI(4)=-Sqrt(D4*Pi/D5)*(-R2*M(4,1) + R3/R2*M(3,2))
! B40     
      BKQR(5)= Sqrt(Pi)/D5*(D6*M(1,1) - D4*M(2,2) - D4*M(3,3) + M(4,4) + M(5,5))
      BKQI(5)=D0
! B41      
      BKQR(6)=D2*Sqrt(D2*Pi/D5)*(-R3/R2*M(3,1) + D1/(D2*R2)*(M(4,2) + M(5,3)))   !  (M(2,4) + M(3,5))
      BKQI(6)=D2*Sqrt(D2*Pi/D5)*( R3/R2*M(2,1) + D1/(D2*R2)*(M(5,2) - M(4,3)))   !  (M(2,5) - M(3,4))
! B42      
      BKQR(7)=D2*Sqrt(D2*Pi/D5)*( R3/D2*M(5,1) + (M(3,3) - M(2,2))/D2)
      BKQI(7)=D2*Sqrt(D2*Pi/D5)*(-R3/D2*M(4,1) - M(3,2))
! B43      
      BKQR(8)=Sqrt(D7*Pi/D5)*(M(4,2) - M(5,3))   !  (M(2,4) - M(3,5))
      BKQI(8)=Sqrt(D7*Pi/D5)*(M(4,3) + M(5,2))   !  (M(3,4) + M(2,5))
! B44      
      BKQR(9)= Sqrt(D7*Pi/D10)*(M(5,5) - M(4,4))
      BKQI(9)=-Sqrt(D7*Pi/D5)*R2*M(5,4)  ! changed D10 to D5 
      
      do i=10,16
        BKQR(i)= D0
        BKQI(i)= D0
      enddo
!   Adjust phases so the Bkq and Bkq' correspond to those of Goeller-Walrand.
      Do k=0,4,2
        do q=0,k
          i=Bkq_index(k,q)
          BKQR(i)= sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQR(i)
          BKQI(i)=-sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQI(i)  ! minus sign added 6/2/15
        enddo
      enddo
      AOMchanged=.false.
!
      return
      end subroutine AOMmatrixD
!
!-----------------------------------------------------------------------
!
      SUBROUTINE AOMintD()
! Calculates the AOM intensity matrix for d-orbitals, 
! and the equivalent crystal field parameters.
      IMPLICIT none
      INTEGER*4 iL,i,j,k,q
      REAL*8 T,P,Ps,Pi,fv1(5),M(5,5),M1(5,5),F(5,5),D0,D1,D2,D3,D4,D5,D6,D7,D10,R2,R3
      PARAMETER(D0=0.0D+00,D1=1.0D+00,D2=2.0D+00,D3=3.0D+00,D4=4.0D+00,D5=5.0D+00,D6=6.0D+00,D7=7.0D+00)
      PARAMETER(D10=10.0D+00)
      Pi=D4*ATAN(D1)
      R2=Sqrt(D2); R3=Sqrt(D3)
      do iL=1,NLIGANDS
        do i=1,5; do j=1,5; M(i,j) = D0; Enddo; Enddo
! Need to do: write in M matrix for each ligand.
        T=Theta(iL)*Pi/180.d+00
        P=Phi(iL)*Pi/180.d+00
        Ps=Chi(iL)*Pi/180.d+00
!  The <d|V|L> elements from table 2b Schaeffer,Struct&Bond,5,68,(1968). 
!  Note: "Standard" ordering of orbitals = z2, yz, xz, xy, x2-y2.
!  Note that in this case we need all 5 columns, not just the first 3.
!  This is because the AOM transition moment matrix is non-diagonal.
        F(1,1)=+(D1 + D3*Cos(D2*T))/D4
        F(1,2)=+R3/D2*Sin(D2*T)*Sin(Ps)  ! PiY
        F(1,3)=-R3/D2*Sin(D2*T)*Cos(Ps)  ! PiX
        F(1,4)=+R3/D4*(D1-Cos(D2*T))*(-Sin(D2*Ps)) 
        F(1,5)=+R3/D4*(D1-Cos(D2*T))*( Cos(D2*Ps)) 
        F(2,1)= R3/D2*Sin(D2*T)*Sin(P)
        F(2,2)= Cos(T)*Cos(P)*Cos(Ps) - Cos(D2*T)*Sin(P)*Sin(Ps)
        F(2,3)= Cos(T)*Cos(P)*Sin(Ps) + Cos(D2*T)*Sin(P)*Cos(Ps)
        F(2,4)= Cos(P)*(-Sin(T))*Cos(D2*Ps) + Sin(P)*(-Sin(D2*T)/D2)*(-Sin(D2*Ps))
        F(2,5)= Cos(P)*(-Sin(T))*Sin(D2*Ps) + Sin(P)*(-Sin(D2*T)/D2)*( Cos(D2*Ps))
        F(3,1)= R3/D2*Sin(D2*T)*Cos(P)
        F(3,2)=-Cos(T)*Sin(P)*Cos(Ps) - Cos(D2*T)*Cos(P)*Sin(Ps)
        F(3,3)=-Cos(T)*Sin(P)*Sin(Ps) + Cos(D2*T)*Cos(P)*Cos(Ps)
        F(3,4)=(-Sin(P))*(-Sin(T))*Cos(D2*Ps) + Cos(P)*(-Sin(D2*T)/D2)*(-Sin(D2*Ps))
        F(3,5)=(-Sin(P))*(-Sin(T))*Sin(D2*Ps) + Cos(P)*(-Sin(D2*T)/D2)*( Cos(D2*Ps))
        F(4,1)= R3/D4*(D1-Cos(D2*T))*Sin(D2*P)
        F(4,2)= Sin(T)*Cos(D2*P)*Cos(Ps) - Sin(D2*T)*Sin(D2*P)*Sin(Ps)/D2
        F(4,3)= Sin(T)*Cos(D2*P)*Sin(Ps) + Sin(D2*T)*Sin(D2*P)*Cos(Ps)/D2
        F(4,4)=(Cos(D2*P))*(Cos(T))*Cos(D2*Ps) + Sin(D2*P)*(D3+Cos(D2*T))/D4*(-Sin(D2*Ps))
        F(4,5)=(Cos(D2*P))*(Cos(T))*Sin(D2*Ps) + Sin(D2*P)*(D3+Cos(D2*T))/D4*( Cos(D2*Ps))
        F(5,1)= R3/D4*(D1-Cos(D2*T))*Cos(D2*P)
        F(5,2)=-Sin(T)*Sin(D2*P)*Cos(Ps) - Sin(D2*T)*Cos(D2*P)*Sin(Ps)/D2
        F(5,3)=-Sin(T)*Sin(D2*P)*Sin(Ps) + Sin(D2*T)*Cos(D2*P)*Cos(Ps)/D2
        F(5,4)=(-Sin(D2*P))*(Cos(T))*Cos(D2*Ps) + Cos(D2*P)*(D3+Cos(D2*T))/D4*(-Sin(D2*Ps))
        F(5,5)=(-Sin(D2*P))*(Cos(T))*Sin(D2*Ps) + Cos(D2*P)*(D3+Cos(D2*T))/D4*( Cos(D2*Ps))
!        
        do i=1,5
          do j=1,5
!transform each M matrix:
            M(i,j) = M(i,j) + eSigma(iL)*F(i,1)*F(j,1) + ePiX(iL)*F(i,2)*F(j,2) + ePiY(iL)*F(i,3)*F(j,3)
! Add each transformed matrix:
            M1(i,j)= M1(i,j)+ M(i,j) ! M1 gets destroyed below
          enddo
        enddo
      enddo  ! iL=1,NLIGANDS  
 
      if (outp(3).or.check1(8)) then
        if (outp(3)) then
          write(io,'(/," AOM transition matrix in terms of the 5x5 real 5-orbitals")') 
          write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>  ")') 
          write(io,'(2X," <sig|",5F12.2)') (M(1,j),j=1,5)
          write(io,'(2X," <piS|",5F12.2)') (M(2,j),j=1,5)
          write(io,'(2X," <piC|",5F12.2)') (M(3,j),j=1,5)
          write(io,'(2X,"<delS|",5F12.2)') (M(4,j),j=1,5)
          write(io,'(2X,"<delC|",5F12.2)') (M(5,j),j=1,5)
        endif
      Endif
!
!  Relationships below are given by Gerloch & McMeeking, Dalton, 2443, (1975). Table 2.
!  These are the Ckq coefficents of the ligand field expanded as spherical harmonics.
! B00 
      BKQR(1)=D2/D5*Sqrt(Pi)*(M(1,1) + M(2,2) + M(3,3) + M(4,4) + M(5,5))
      BKQI(1)=D0
! B20     
      BKQR(2)= Sqrt(Pi/D5)*(D2*M(1,1) + M(2,2) + M(3,3) - D2*M(4,4) -D2*M(5,5))
      BKQI(2)=D0
!  Note: "Standard" ordering of orbitals = z2, yz, xz, xy, x2-y2.
! B21      
      BKQR(3)=-Sqrt(D4*Pi/D5)*( R3/R2*(M(2,4) + M(3,5)) + M(3,1)/R2)
      BKQI(3)= Sqrt(D4*Pi/D5)*(-R3/R2*(M(2,5) - M(3,4)) + M(2,1)/R2)
! B22      
      BKQR(4)=-Sqrt(D4*Pi/D5)*( R2*M(5,1) + R3/(D2*R2)*(M(2,2) - M(3,3)))
      BKQI(4)=-Sqrt(D4*Pi/D5)*(-R2*M(4,1) + R3/R2*M(3,2))
! B40     
      BKQR(5)= Sqrt(Pi)/D5*(D6*M(1,1) - D4*M(2,2) - D4*M(3,3) + M(4,4) + M(5,5))
      BKQI(5)=D0
! B41      
      BKQR(6)=D2*Sqrt(D2*Pi/D5)*(-R3/R2*M(3,1) + D1/(D2*R2)*(M(2,4) + M(3,5)))
      BKQI(6)=D2*Sqrt(D2*Pi/D5)*( R3/R2*M(2,1) + D1/(D2*R2)*(M(2,5) - M(3,4)))
! B42      
      BKQR(7)=D2*Sqrt(D2*Pi/D5)*( R3/D2*M(5,1) + (M(3,3) - M(2,2))/D2)
      BKQI(7)=D2*Sqrt(D2*Pi/D5)*(-R3/D2*M(4,1) - M(3,2))
! B43      
      BKQR(8)=Sqrt(D7*Pi/D5)*(M(2,4) - M(3,5))
      BKQI(8)=Sqrt(D7*Pi/D5)*(M(3,4) + M(2,5))
! B44      
      BKQR(9)= Sqrt(D7*Pi/D10)*(M(5,5) - M(4,4))
      BKQI(9)=-Sqrt(D7*Pi/D5)*R2*M(5,4)  ! changed D10 to D5 
      do i=10,16
        BKQR(i)= D0
        BKQI(i)= D0
      enddo
!   Adjust phases so the Bkq and Bkq' correspond to those of Goeller-Walrand.
      Do k=0,4,2
        do q=0,k
          i=Bkq_index(k,q)
          BKQR(i)= sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQR(i)
          BKQI(i)=-sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQI(i)  ! minus sign added 6/2/15
        enddo
      enddo
      AOMchanged=.false.
!
      return
      end subroutine AOMintD
!
!-----------------------------------------------------------------------
!
      SUBROUTINE AOMmatrixFX()
! ** TODO: not activated yet.      
! Extended AOM parameters.   
! New Subroutine for non-zero Chi values (anisotropic pi, delta, psi bonding).
! Calculates the AOM matrix for f-orbitals, and the equivalent crystal field parameters.
      IMPLICIT none
      INTEGER*4 iL,i,j,k,q
      REAL*8 T,P,Ps,Pi,fv1(7), E1, sum1,sum2
      REAL*8 M(7,7),M1(7,7),F(7,7),D0,D1,D2,D3,D4,D5,D6,D7,D8,D10,D11,D13,D14,D15,D16
      PARAMETER(D0=0.0D+00,D1=1.0D+00,D2=2.0D+00,D3=3.0D+00,D4=4.0D+00,D5=5.0D+00,D6=6.0D+00,D7=7.0D+00)
      PARAMETER(D8=8.0D+00,D10=10.0D+00,D11=11.0D+00,D13=13.0D+00,D14=14.0D+00,D15=15.0D+00,D16=16.0D+00)
!
      Pi=D4*ATAN(D1)
      
      if (LFtype.eq."AOM ") then
        do i=1,7
          do j=1,7
            M(i,j) = D0
          enddo
        enddo  
        do iL=1,NLIGANDS
          T=Theta(iL)*Pi/180.d+00
          P=Phi(iL)*Pi/180.d+00
          Ps=Chi(iL)*Pi/180.d+00
!
!  Note: "Standard" ordering of real f-orbitals |sigma>, |piS>, |piC>, |deltaS>, |deltaC>, |phiS>, |phiC>.
!  For a definition of these orbitals: see Harnung & Schaffer, Struct&Bond,12,201,(1972).
      F(1,1)=(D3*Cos(T) +D5*Cos(D3*T))/D8
      F(1,2)= (Sqrt(D3/D2)*(Sin(T)+D5*Sin(D3*T)))*Sin(Ps)/D8
      F(1,3)=-(Sqrt(D3/D2)*(Sin(T)+D5*Sin(D3*T)))*Cos(Ps)/D8
      F(1,4)=-(Sqrt(D15)*Cos(T)*Sin(T)**2)*Sin(D2*Ps)/D2
      F(1,5)= (Sqrt(D15)*Cos(T)*Sin(T)**2)*Cos(D2*Ps)/D2
      F(1,6)= (Sqrt(D5/D2)*Sin(T)**3)*Sin(D3*Ps)/D2
      F(1,7)=-(Sqrt(D5/D2)*Sin(T)**3)*Cos(D3*Ps)/D2
      
      F(2,1)=(Sqrt(D3/D2)*Sin(P)*(Sin(T) +D5*Sin(D3*T)))/D8
      F(2,2)=(Cos(Ps)*Cos(P)*(D3+D5*Cos(D2*T))+Cos(T)*( D7-D15*Cos(D2*T))*Sin(Ps)*Sin(P))/D8
      F(2,3)=(Cos(P)*(D3+D5*Cos(D2*T))*Sin(Ps)+Cos(Ps)*Cos(T)*(-D7+D15*Cos(D2*T))*Sin(P))/D8
      F(2,4)=(Sqrt(D5/D2)*(-D4*Cos(D2*Ps)*Cos(P)*Cos(T) + (D1+D3*Cos(D2*T))*Sin(D2*Ps)*Sin(P))*Sin(T))/D4
      F(2,5)=(Sqrt(D5/D2)*(-D4*Cos(P)*Sin(D2*Ps)*Sin(D2*T) + Cos(D2*Ps)*Sin(P)*(Sin(T)-D3*Sin(D3*T))))/D8
      F(2,6)=(Sqrt(D15)*(Cos(D3*Ps)*Cos(P) - Cos(T)*Sin(D3*Ps)*Sin(P))*Sin(T)**2)/D4
      F(2,7)=(Sqrt(D15)*(Cos(P)*Sin(D3*Ps) + Cos(D3*Ps)*Cos(T)*Sin(P))*Sin(T)**2)/D4
      
      F(3,1)=(Sqrt(D3/D2)*Cos(P)*(Sin(T) +D5*Sin(D3*T)))/D8
      F(3,2)=(-Cos(P)*(Cos(T)+D15*Cos(D3*T))*Sin(Ps) - D2*Cos(Ps)*(D3+D5*Cos(D2*T))*Sin(P))/D16
      F(3,3)=( Cos(Ps)*Cos(P)*(Cos(T)+D15*Cos(D3*T)) - D2*(D3+D5*Cos(D2*T))*Sin(Ps)*Sin(P))/D16
      F(3,4)=(Sqrt(D5/D2)*(Cos(P)*(D1+D3*Cos(D2*T))*Sin(D2*Ps) +D4*Cos(D2*Ps)*Cos(T)*Sin(P))*Sin(T))/D4
      F(3,5)=(Sqrt(D5/D2)*(D4*Sin(D2*Ps)*Sin(P)*Sin(D2*T) + Cos(D2*Ps)*Cos(P)*(Sin(T)-D3*Sin(D3*T))))/D8
      F(3,6)=-(Sqrt(D15)*(Cos(P)*Cos(T)*Sin(D3*Ps) + Cos(D3*Ps)*Sin(P))*Sin(T)**2)/D4
      F(3,7)= (Sqrt(D15)*(Cos(D3*Ps)*Cos(P)*Cos(T) - Sin(D3*Ps)*Sin(P))*Sin(T)**2)/D4
      
      F(4,1)=Sqrt(D15)*Cos(P)*Cos(T)*Sin(P)*Sin(T)**2
      F(4,2)=(Sqrt(D5/D2)*(D4*Cos(Ps)*Cos(D2*P)*Cos(T) - (D1+D3*Cos(D2*T))*Sin(Ps)*Sin(D2*P))*Sin(T))/D4
      F(4,3)=(Sqrt(D5/D2)*(D4*Cos(D2*P)*Sin(Ps)*Sin(D2*T) + Cos(Ps)*Sin(D2*P)*(-Sin(T)+D3*Sin(D3*T))))/D8
      F(4,4)=Cos(D2*Ps)*Cos(D2*P)*Cos(D2*T) - ((D5*Cos(T)+D3*Cos(D3*T))*Sin(D2*Ps)*Sin(D2*P))/D8
      F(4,5)=Cos(D2*P)*Cos(D2*T)*Sin(D2*Ps) + (Cos(D2*Ps)*(D5*Cos(T)+D3*Cos(D3*T))*Sin(D2*P))/D8
      F(4,6)=(Sqrt(D3/D2)*(-D4*Cos(D3*Ps)*Cos(D2*P)*Sin(D2*T) + Sin(D3*Ps)*Sin(D2*P)*(D5*Sin(T)+Sin(D3*T))))/D8
      F(4,7)=(Sqrt(D3/D2)*(-D4*Cos(D2*P)*Sin(D3*Ps)*Sin(D2*T) - Cos(D3*Ps)*Sin(D2*P)*(D5*Sin(T)+Sin(D3*T))))/D8
      
      F(5,1)=(Sqrt(D15)*Cos(D2*P)*Cos(T)*Sin(T)**2)/D2 
      F(5,2)=(Sqrt(D5/D2)*(-D4*Cos(Ps)*Sin(D2*P)*Sin(D2*T) + Cos(D2*P)*Sin(Ps)*(Sin(T)-D3*Sin(D3*T))))/D8
      F(5,3)=(Sqrt(D5/D2)*(Cos(Ps)*Cos(D2*P)*(D1+D3*Cos(D2*T)) - D8*Cos(P)*Cos(T)*Sin(Ps)*Sin(P))*Sin(T))/D4
      F(5,4)=-(Cos(D2*P)*(D5*Cos(T) +D3*Cos(D3*T))*Sin(D2*Ps))/D8 - Cos(D2*Ps)*Cos(D2*T)*Sin(D2*P)
      F(5,5)= (Cos(D2*Ps)*Cos(D2*P)*(D5*Cos(T) +D3*Cos(D3*T)))/D8 - Cos(D2*T)*Sin(D2*Ps)*Sin(D2*P)
      F(5,6)= (Sqrt(D3/D2)*( D4*Cos(D3*Ps)*Sin(D2*P)*Sin(D2*T) + Cos(D2*P)*Sin(D3*Ps)*(D5*Sin(T)+Sin(D3*T))))/D8
      F(5,7)=-(Sqrt(D3/D2)*(-D4*Sin(D3*Ps)*Sin(D2*P)*Sin(D2*T) + Cos(D3*Ps)*Cos(D2*P)*(D5*Sin(T)+Sin(D3*T))))/D8
      
      F(6,1)=(Sqrt(D5/D2)*(D1+D2*Cos(D2*P))*Sin(P)*Sin(T)**3)/D2
      F(6,2)=(Sqrt(D15)*(Cos(Ps)*Cos(D3*P)-Cos(T)*Sin(Ps)*Sin(D3*P))*Sin(T)**2)/D4
      F(6,3)=(Sqrt(D15)*(Cos(D3*P)*Sin(Ps)+Cos(Ps)*Cos(T)*Sin(D3*P))*Sin(T)**2)/D4
      F(6,4)=(Sqrt(D3/D2)*(D4*Cos(D2*Ps)*Cos(D3*P)*Sin(D2*T) - Sin(D2*Ps)*Sin(D3*P)*(D5*Sin(T)+Sin(D3*T))))/D8
      F(6,5)=(Sqrt(D3/D2)*(D4*Cos(D3*P)*Sin(D2*Ps)*Sin(D2*T) + Cos(D2*Ps)*Sin(D3*P)*(D5*Sin(T)+Sin(D3*T))))/D8
      F(6,6)=(Cos(D3*Ps)*Cos(D3*P)*(D5+D3*Cos(D2*T)))/D8 - ((D15*Cos(T)+Cos(D3*T))*Sin(D3*Ps)*Sin(D3*P))/D16
      F(6,7)=(Cos(D3*P)*(D5+D3*Cos(D2*T))*Sin(D3*Ps))/D8 + (Cos(D3*Ps)*(D15*Cos(T)+Cos(D3*T))*Sin(D3*P))/D16
      
      F(7,1)=(Sqrt(D5/D2)*Cos(P)*(-D1+D2*Cos(D2*P))*Sin(T)**3)/D2
      F(7,2)=-(Sqrt(D15)*(Cos(D3*P)*Cos(T)*Sin(Ps) + Cos(Ps)*Sin(D3*P))*Sin(T)**2)/D4
      F(7,3)= (Sqrt(D15)*(Cos(Ps)*Cos(D3*P)*Cos(T) - Sin(Ps)*Sin(D3*P))*Sin(T)**2)/D4
      F(7,4)=-(Sqrt(D3/D2)*( D4*Cos(D2*Ps)*Sin(D3*P)*Sin(D2*T) + Cos(D3*P)*Sin(D2*Ps)*(D5*Sin(T)+Sin(D3*T))))/D8
      F(7,5)= (Sqrt(D3/D2)*(-D4*Sin(D2*Ps)*Sin(D3*P)*Sin(D2*T) + Cos(D2*Ps)*Cos(D3*P)*(D5*Sin(T)+Sin(D3*T))))/D8
      F(7,6)=-(Cos(D3*P)*(D15*Cos(T)+Cos(D3*T))*Sin(D3*Ps))/D16 - (Cos(D3*Ps)*(D5+D3*Cos(D2*T))*Sin(D3*P))/D8
      F(7,7)= (Cos(D3*Ps)*Cos(D3*P)*(D15*Cos(T)+Cos(D3*T)))/D16 - ((D5+D3*Cos(D2*T))*Sin(D3*Ps)*Sin(D3*P))/D8
              
          do i=1,7
            do j=1,7
              M(i,j) = M(i,j) + eSigma(iL)*F(i,1)*F(j,1) + ePiY(iL)*F(i,2)*F(j,2) + ePiX(iL)*F(i,3)*F(j,3)    &
                                                         + eDelS(iL)*F(i,4)*F(j,4) + eDelC(iL)*F(i,5)*F(j,5)  &
                                                         + ePsiS(iL)*F(i,6)*F(j,6) + ePsiC(iL)*F(i,7)*F(j,7)
              M1(i,j)= M(i,j) ! M1 gets destroyed below
            enddo
          enddo
        enddo  !  iL=1,NLIGANDS
! 
        if (outp(3) .or. outp(9) .or.check1(8)) then
          call diagrs(M1,7,7,OrbEngs,fv1,2,io)  ! calc vects
          E1=EAVE; EAVE=D0; call zeroE(OrbEngs,7); EAVE=E1
          sum1=(M(1,1)+M(2,2)+M(3,3)+M(4,4)+M(5,5)+M(6,6)+M(7,7))/D7
          do i=1,7; M(i,i)=M(i,i)-sum1; enddo  ! subtract baricentre
          write(io,'(/," AOM matrix in terms of the 7x7 real f-orbital matrix (baricentre=",F12.2," subtracted)")') sum1
          write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>      |psiS>      |psiC>     orbital energies")') 
          write(io,'(2X," <sig|",7F12.2,5X,F12.2)') (M(1,j),j=1,7),OrbEngs(7)
          write(io,'(2X," <piS|",7F12.2,5X,F12.2)') (M(2,j),j=1,7),OrbEngs(6)
          write(io,'(2X," <piC|",7F12.2,5X,F12.2)') (M(3,j),j=1,7),OrbEngs(5)
          write(io,'(2X,"<delS|",7F12.2,5X,F12.2)') (M(4,j),j=1,7),OrbEngs(4)
          write(io,'(2X,"<delC|",7F12.2,5X,F12.2)') (M(5,j),j=1,7),OrbEngs(3)
          write(io,'(2X,"<psiS|",7F12.2,5X,F12.2)') (M(6,j),j=1,7),OrbEngs(2)
          write(io,'(2X,"<psiC|",7F12.2,5X,F12.2)') (M(7,j),j=1,7),OrbEngs(1)
          if (outp(9)) then  ! print energies and eigenvectors
            write(io,'(/,2X,"Energies:",7F12.2)') (OrbEngs(j),j=1,7)
            write(io,'(2X,"   |z^3>   :",7F12.2)') (M1(1,j),j=1,7)
            write(io,'(2X,"   |yz^2>  :",7F12.2)') (M1(2,j),j=1,7)
            write(io,'(2X,"   |xz^2>  :",7F12.2)') (M1(3,j),j=1,7)
            write(io,'(2X," |z(xy)>   :",7F12.2)') (M1(4,j),j=1,7)
            write(io,'(2X," |z(x2-y2)>:",7F12.2)') (M1(5,j),j=1,7)
            write(io,'(2X,"|y(3x2-y2)>:",7F12.2)') (M1(6,j),j=1,7)
            write(io,'(2X,"|x(x2-3y2)>:",7F12.2)') (M1(7,j),j=1,7)
          endif
        Endif
        if (FitType.eq.6) then
          do i=1,7; do j=1,i; LF1eMat(i,j)=M(i,j); enddo; enddo   !?????  changed =M(1,j) to =M(i,j)
          sum1=(EXP_E(1)+EXP_E(3)+EXP_E(6)+EXP_E(10)+EXP_E(15)+EXP_E(21)+EXP_E(28))/D7
          sum2=(LF1eMat(1,1)+LF1eMat(2,2)+LF1eMat(3,3)+LF1eMat(4,4)+LF1eMat(5,5)+LF1eMat(6,6)+LF1eMat(7,7))/D7
          do i=1,7; LF1eMat(i,i)=LF1eMat(i,i)-sum2+sum1; enddo
        endif
        
      elseif (LFtype.eq."LF1E") then
        do i=1,7
          do j=1,i;  M(i,j)=LF1eMat(i,j); if (i.ne.j) M(j,i)=LF1eMat(i,j); enddo
        enddo
      endif
!      do i=1,7; write(io,'(7f12.4)') (M(i,j),j=1,7); enddo

!
!  Relationships below are given by Urland, Chem.Phys.14, 393,(1976). Table 3.
!  These are the Ckq coefficents of the ligand field expanded as spherical harmonics.
! B00 
      BKQR(1)=(D2*Sqrt(Pi)*(M(1,1) + M(2,2) + M(3,3) + M(4,4) + M(5,5) + M(6,6) + M(7,7)))/D7
      BKQI(1)=D0
! B20     
      BKQR(2)=(D2*Sqrt(D5*Pi)*M(1,1))/D7 + (D3*Sqrt(D5*Pi)*(M(2,2) + M(3,3)))/D14 - (D5*Sqrt(D5*Pi)*(M(6,6) + M(7,7)))/D14
      BKQI(2)=D0
!      CF(3)=(Sqrt(D5*Pi)*((0,1)*M(2,1) - M(3,1)))/D7 + (D5*Sqrt(D3*Pi)*(-M(4,2) + (0,1)*M(4,3) - (0,1)*M(5,2) - M(5,3)))/D14 + (D5*Sqrt(D5*Pi)*(-M(6,4) + (0,1)*M(6,5) - (0,1)*M(7,4) - M(7,5)))/D14
      BKQR(3)=(Sqrt(D5*Pi)*(-M(3,1)))/D7 + (D5*Sqrt(D3*Pi)*(-M(4,2) - M(5,3)))/D14 + (D5*Sqrt(D5*Pi)*(-M(6,4) - M(7,5)))/D14
      BKQI(3)=(Sqrt(D5*Pi)*( M(2,1)))/D7 + (D5*Sqrt(D3*Pi)*( M(4,3) - M(5,2)))/D14 + (D5*Sqrt(D5*Pi)*( M(6,5) - M(7,4)))/D14
!      CF(4)=(Sqrt(D5*D6*Pi)*(-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2))/D7 + (D5*Sqrt(D2*Pi)*((0,1)*M(4,1) - M(5,1)))/D7 + (D5*Sqrt(Pi/D2)*(-M(6,2) + (0,1)*M(6,3) - (0,1)*M(7,2) - M(7,3)))/D7
      BKQR(4)=Sqrt(D5*D6*Pi)*(-M(2,2)/D2 + M(3,3)/D2)/D7 + D5*Sqrt(D2*Pi)*(-M(5,1))/D7 + D5*Sqrt(Pi/D2)*(-M(6,2) - M(7,3))/D7
      BKQI(4)=Sqrt(D5*D6*Pi)*(-M(2,3))/D7 + D5*Sqrt(D2*Pi)*(M(4,1))/D7 + D5*Sqrt(Pi/D2)*(M(6,3) - M(7,2))/D7
! B40
      BKQR(5)=-(Sqrt(Pi)*(M(4,4) + M(5,5))) + (Sqrt(Pi)*(D6*M(1,1) + M(2,2) + M(3,3) + D3*M(6,6) + D3*M(7,7)))/D7
      BKQI(5)=D0      
!      CF(6)=(Sqrt(D5*D6*Pi)*((0,1)*M(2,1) - M(3,1)))/D7 + (D4*Sqrt(D2*Pi)*(-M(4,2) + (0,1)*M(4,3) - (0,1)*M(5,2) - M(5,3)))/D7 + (Sqrt(D5*D6*Pi)*(M(6,4) - (0,1)*M(6,5) + (0,1)*M(7,4) + M(7,5)))/D7
      BKQR(6)=(Sqrt(D5*D6*Pi)*(-M(3,1)))/D7 + (D4*Sqrt(D2*Pi)*(-M(4,2) - M(5,3)))/D7 + (Sqrt(D5*D6*Pi)*( M(6,4) + M(7,5)))/D7
      BKQI(6)=(Sqrt(D5*D6*Pi)*( M(2,1)))/D7 + (D4*Sqrt(D2*Pi)*( M(4,3) - M(5,2)))/D7 + (Sqrt(D5*D6*Pi)*(-M(6,5) + M(7,4)))/D7
!      CF(7)=(D2*Sqrt(D10*Pi)*(-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2))/D7 + (Sqrt(D6*Pi)*((0,1)*M(4,1) - M(5,1)))/D7 + (3*Sqrt(D6*Pi)*(M(6,2) - (0,1)*M(6,3) + (0,1)*M(7,2) + M(7,3)))/D7
      BKQR(7)=(D2*Sqrt(D10*Pi)*(-M(2,2)/D2 + M(3,3)/D2))/D7 + (Sqrt(D6*Pi)*(-M(5,1)))/D7 + (D3*Sqrt(D6*Pi)*(M(6,2) + M(7,3)))/D7
      BKQI(7)=(D2*Sqrt(D10*Pi)*(-M(2,3)))/D7 + Sqrt(D6*Pi)/D7*(M(4,1) + D3*(-M(6,3) + M(7,2)))
!      CF(8)=Sqrt((D2*Pi)/D7)*(M(4,2) + (0,1)*M(4,3) + (0,1)*M(5,2) - M(5,3)) + D3*Sqrt((D2*Pi)/D7)*((0,-1)*M(6,1) + M(7,1))
      BKQR(8)=Sqrt((D2*Pi)/D7)*(M(4,2) - M(5,3)) + D3*Sqrt((D2*Pi)/D7)*(M(7,1))
      BKQI(8)=Sqrt((D2*Pi)/D7)*(M(4,3) + M(5,2)) - D3*Sqrt((D2*Pi)/D7)*(M(6,1))
!      CF(9)=Sqrt((D10*Pi)/D7)* (-M(4,4)/D2 - (0,1)*M(4,5) + M(5,5)/D2) + Sqrt((D6*Pi)/D7)*(M(6,2) + (0,1)*M(6,3) + (0,1)*M(7,2) - M(7,3))
      BKQR(9)=Sqrt((D10*Pi)/D7)*(-M(4,4)/D2 + M(5,5)/D2) + Sqrt((D6*Pi)/D7)*(M(6,2) - M(7,3))
      BKQI(9)=Sqrt((D10*Pi)/D7)*(-M(5,4)) + Sqrt((D6*Pi)/D7)*(M(6,3) + M(7,2))
! B60
      BKQR(10)=Sqrt(D13*Pi)/D7*(D2*M(1,1) - D3*(M(2,2) + M(3,3))/D2 + D3*(M(4,4) + M(5,5))/D5 - (M(6,6) + M(7,7))/(D2*D5))
      BKQI(10)=D0
!      CF(11)=Sqrt((D13*Pi)/D7)*((0,1)*M(2,1) - M(3,1)) + (Sqrt((39*Pi)/35.)*(M(4,2) - (0,1)*M(4,3) + (0,1)*M(5,2) + M(5,3)))/D2 + (Sqrt((D13*Pi)/D7)*(-M(6,4) + (0,1)*M(6,5) - (0,1)*M(7,4) - M(7,5)))/(D2*D5)
      BKQR(11)=Sqrt((D13*Pi)/D7)*(-M(3,1) + Sqrt(D3/D5)*( M(4,2) + M(5,3))/D2 + (-M(6,4) - M(7,5))/(D2*D5))
      BKQI(11)=Sqrt((D13*Pi)/D7)*( M(2,1) + Sqrt(D3/D5)*(-M(4,3) + M(5,2))/D2 + ( M(6,5) - M(7,4))/(D2*D5))      
!      CF(12)=Sqrt((D3*D13*Pi)/(D5*D7))* (-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2) + (D4*Sqrt((D13*Pi)/D7)*((0,-1)*M(4,1) + M(5,1)))/D5 + (Sqrt((D13*Pi)/D7)*(-M(6,2) + (0,1)*M(6,3) - (0,1)*M(7,2) - M(7,3)))/D5
      BKQR(12)=Sqrt((D13*Pi)/D7)*(-Sqrt(D3/D5)*(M(2,2) - M(3,3))/D2 + D4*(M(5,1))/D5 - (M(6,2) + M(7,3))/D5)
      BKQI(12)=Sqrt((D13*Pi)/D7)*(-Sqrt(D3/D5)*M(3,2) - D4*(M(4,1))/D5 + (M(6,3) - M(7,2))/D5)
!      CF(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,2) + (0,1)*M(4,3) + (0,1)*M(5,2) - M(5,3)))/D5 + (Sqrt((D6*D13*Pi)/D7)*((0,1)*M(6,1) - M(7,1)))/D5
      BKQR(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,2) - M(5,3)))/D5 + (Sqrt((D6*D13*Pi)/D7)*(-M(7,1)))/D5
      BKQI(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,3) + M(5,2)))/D5 + (Sqrt((D6*D13*Pi)/D7)*( M(6,1)))/D5      
!      CF(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(4,4)/D2 - (0,1)*M(4,5) + M(5,5)/D2))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,2) - (0,1)*M(6,3) - (0,1)*M(7,2) + M(7,3))
      BKQR(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(4,4)/D2 + M(5,5)/D2))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,2) + M(7,3))
      BKQI(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(5,4)))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,3) - M(7,2))
!      CF(15)=(Sqrt((D3*D11*D13*Pi)/D14)* (M(6,4) + (0,1)*M(6,5) + (0,1)*M(7,4) - M(7,5)))/D5
      BKQR(15)=(Sqrt((D3*D11*D13*Pi)/D14)*(M(6,4) - M(7,5)))/D5
      BKQI(15)=(Sqrt((D3*D11*D13*Pi)/D14)*(M(6,5) + M(7,4)))/D5    
!      CF(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(6,6)/D2 - (0,1)*M(6,7) + M(7,7)/D2))/D5
      BKQR(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(6,6)/D2 + M(7,7)/D2))/D5
      BKQI(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(7,6)))/D5
      
!   Adjust phases so the Bkq and Bkq' correspond to those of Goeller-Walrand.
      Do k=0,6,2
        do q=0,k
          i=Bkq_index(k,q)
          BKQR(i)= sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQR(i)
          BKQI(i)=-sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQI(i)  ! minus sign added 6/2/15
        enddo
      enddo
      AOMchanged=.false.
!      write(io,'("AOMmatrixF: eSigma(1),eSigma(2),BKQR(2)=",3F10.4)') eSigma(1),eSigma(2),BKQR(2)
!
      return
      end subroutine AOMmatrixFX
!
!-----------------------------------------------------------------------
!
      SUBROUTINE AOMmatrixF()
! New Subroutine for non-zero Chi values (anisotropic pi bonding).
! Calculates the AOM matrix for f-orbitals, and the equivalent crystal field parameters.
      IMPLICIT none
      INTEGER*4 iL,i,j,k,q
      REAL*8 T,P,Ps,Pi,fv1(7), A1,A2,A3,B1,B2,B3,G1,G2,G3,E1, sum1,sum2
      REAL*8 M(7,7),M1(7,7),F(7,7),D0,D1,D2,D3,D4,D5,D6,D7,D8,D10,D11,D13,D14,D15
      PARAMETER(D0=0.0D+00,D1=1.0D+00,D2=2.0D+00,D3=3.0D+00,D4=4.0D+00,D5=5.0D+00,D6=6.0D+00,D7=7.0D+00)
      PARAMETER(D8=8.0D+00,D10=10.0D+00,D11=11.0D+00,D13=13.0D+00,D14=14.0D+00,D15=15.0D+00)
!
      Pi=D4*ATAN(D1)
      
      if (LFtype.eq."AOM ") then
        do i=1,7
          do j=1,7
            M(i,j) = D0
          enddo
        enddo  
        do iL=1,NLIGANDS
          T=Theta(iL)*Pi/180.d+00
          P=Phi(iL)*Pi/180.d+00
          Ps=Chi(iL)*Pi/180.d+00
!
! Directional cosines: Use matrix given in Zare: (3.37)
          A1 = Cos(P)*Cos(T)*Cos(Ps)-Sin(P)*Sin(Ps);  B1 = -Cos(P)*Cos(T)*Sin(Ps)-Sin(P)*Cos(Ps);  G1 = Cos(P)*Sin(T) 
          A2 = Sin(P)*Cos(T)*Cos(Ps)+Cos(P)*Sin(Ps);  B2 = -Sin(P)*Cos(T)*Sin(Ps)+Cos(P)*Cos(Ps);  G2 = Sin(P)*Sin(T) 
          A3 =       -Sin(T)*Cos(Ps);                 B3 =         Sin(T)*Sin(Ps);                 G3 =        Cos(T)

      
!  The <f|V|L> elements from Table 1: Urland,Chem.Phys.14,393,(1976).
!  Note: "Standard" ordering of real f-orbitals |sigma>, |piS>, |piC>, |deltaS>, |deltaC>, |phiS>, |phiC>.
!  For a definition of these orbitals: see Harnung & Schaffer, Struct&Bond,12,201,(1972).
          F(1,1)= G3*(D5*G3**2-D3)/D2
          F(1,2)= Sqrt(D3/D2)*B3*(D5*G3**2-D1)/D2
          F(1,3)= Sqrt(D3/D2)*A3*(D5*G3**2-D1)/D2

          F(2,1)= Sqrt(D6)*G2*(D5*G3**2-D1)/D4
          F(2,2)= B2*(D5*G3**2-D1)/D4+D5/D2*G2*G3*B3
          F(2,3)= A2*(D5*G3**2-D1)/D4+D5/D2*G2*G3*A3
        
          F(3,1)= Sqrt(D6)*G1*(D5*G3**2-D1)/D4
          F(3,2)= B1*(D5*G3**2-D1)/D4+D5/D2*G1*G3*B3
          F(3,3)= A1*(D5*G3**2-D1)/D4+D5/D2*G1*G3*A3

          F(4,1)= Sqrt(D15)*G1*G2*G3
          F(4,2)= Sqrt(D5/D2)*(B1*G2*G3+B2*G1*G3+B3*G1*G2)
          F(4,3)= Sqrt(D5/D2)*(A1*G2*G3+A2*G1*G3+A3*G1*G2)

          F(5,1)= Sqrt(D15)*G3*(G1**2-G2**2)/D2
          F(5,2)= Sqrt(D5/D8)*(G1**2*B3-G2**2*B3+D2*B1*G1*G3-D2*B2*G2*G3)
          F(5,3)= Sqrt(D5/D8)*(G1**2*A3-G2**2*A3+D2*A1*G1*G3-D2*A2*G2*G3)
        
          F(6,1)= Sqrt(D10)*G2*(D3*G1**2-G2**2)/D4
          F(6,2)= Sqrt(D15)*(B2*(G1**2-G2**2)/D4+G1*G2*B1/D2)
          F(6,3)= Sqrt(D15)*(A2*(G1**2-G2**2)/D4+G1*G2*A1/D2)
        
          F(7,1)= Sqrt(D10)*G1*(G1**2-D3*G2**2)/D4
          F(7,2)= Sqrt(D15)*(B1*(G1**2-G2**2)/D4-G1*G2*B2/D2)
          F(7,3)= Sqrt(D15)*(A1*(G1**2-G2**2)/D4-G1*G2*A2/D2)
        
          do i=1,7
             do j=1,7
              M(i,j) = M(i,j) + eSigma(iL)*F(i,1)*F(j,1) + ePiY(iL)*F(i,2)*F(j,2) + ePiX(iL)*F(i,3)*F(j,3)
              M1(i,j)= M(i,j) ! M1 gets destroyed below
            enddo
          enddo 
        enddo  ! iL=1,NLIGANDS
! 
        if (outp(3) .or. outp(9) .or.check1(8)) then
          call diagrs(M1,7,7,OrbEngs,fv1,2,io)  ! calc vects
          E1=EAVE; EAVE=D0; call zeroE(OrbEngs,7); EAVE=E1
          sum1=(M(1,1)+M(2,2)+M(3,3)+M(4,4)+M(5,5)+M(6,6)+M(7,7))/D7
          do i=1,7; M(i,i)=M(i,i)-sum1; enddo  ! subtract baricentre
          write(io,'(/," AOM matrix in terms of the 7x7 real f-orbital matrix (baricentre=",F12.2," subtracted)")') sum1
          write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>      |psiS>      |psiC>     orbital energies")') 
          write(io,'(2X," <sig|",7F12.2,5X,F12.2)') (M(1,j),j=1,7),OrbEngs(7)
          write(io,'(2X," <piS|",7F12.2,5X,F12.2)') (M(2,j),j=1,7),OrbEngs(6)
          write(io,'(2X," <piC|",7F12.2,5X,F12.2)') (M(3,j),j=1,7),OrbEngs(5)
          write(io,'(2X,"<delS|",7F12.2,5X,F12.2)') (M(4,j),j=1,7),OrbEngs(4)
          write(io,'(2X,"<delC|",7F12.2,5X,F12.2)') (M(5,j),j=1,7),OrbEngs(3)
          write(io,'(2X,"<psiS|",7F12.2,5X,F12.2)') (M(6,j),j=1,7),OrbEngs(2)
          write(io,'(2X,"<psiC|",7F12.2,5X,F12.2)') (M(7,j),j=1,7),OrbEngs(1)
          if (outp(9)) then
            write(io,'(/,2X,"Energies:",7F12.2)') (OrbEngs(j),j=1,7)
            write(io,'(2X,"   |z^3>   :",7F12.2)') (M1(1,j),j=1,7)
            write(io,'(2X,"   |yz^2>  :",7F12.2)') (M1(2,j),j=1,7)
            write(io,'(2X,"   |xz^2>  :",7F12.2)') (M1(3,j),j=1,7)
            write(io,'(2X," |z(xy)>   :",7F12.2)') (M1(4,j),j=1,7)
            write(io,'(2X," |z(x2-y2)>:",7F12.2)') (M1(5,j),j=1,7)
            write(io,'(2X,"|y(3x2-y2)>:",7F12.2)') (M1(6,j),j=1,7)
            write(io,'(2X,"|x(x2-3y2)>:",7F12.2)') (M1(7,j),j=1,7)
          endif
        Endif
        if (FitType.eq.6) then
          do i=1,7; do j=1,i; LF1eMat(i,j)=M(i,j); enddo; enddo   !?????  changed =M(1,j) to =M(i,j)
          sum1=(EXP_E(1)+EXP_E(3)+EXP_E(6)+EXP_E(10)+EXP_E(15)+EXP_E(21)+EXP_E(28))/D7
          sum2=(LF1eMat(1,1)+LF1eMat(2,2)+LF1eMat(3,3)+LF1eMat(4,4)+LF1eMat(5,5)+LF1eMat(6,6)+LF1eMat(7,7))/D7
          do i=1,7; LF1eMat(i,i)=LF1eMat(i,i)-sum2+sum1; enddo
        endif
        
      elseif (LFtype.eq."LF1E") then
        do i=1,7
          do j=1,i;  M(i,j)=LF1eMat(i,j); if (i.ne.j) M(j,i)=LF1eMat(i,j); enddo
        enddo
      endif
!      do i=1,7; write(io,'(7f12.4)') (M(i,j),j=1,7); enddo

!
!  Relationships below are given by Urland, Chem.Phys.14, 393,(1976). Table 3.
!  These are the Ckq coefficents of the ligand field expanded as spherical harmonics.
! B00 
      BKQR(1)=(D2*Sqrt(Pi)*(M(1,1) + M(2,2) + M(3,3) + M(4,4) + M(5,5) + M(6,6) + M(7,7)))/D7
      BKQI(1)=D0
! B20     
      BKQR(2)=(D2*Sqrt(D5*Pi)*M(1,1))/D7 + (D3*Sqrt(D5*Pi)*(M(2,2) + M(3,3)))/D14 - (D5*Sqrt(D5*Pi)*(M(6,6) + M(7,7)))/D14
      BKQI(2)=D0
!      CF(3)=(Sqrt(D5*Pi)*((0,1)*M(2,1) - M(3,1)))/D7 + (D5*Sqrt(D3*Pi)*(-M(4,2) + (0,1)*M(4,3) - (0,1)*M(5,2) - M(5,3)))/D14 + (D5*Sqrt(D5*Pi)*(-M(6,4) + (0,1)*M(6,5) - (0,1)*M(7,4) - M(7,5)))/D14
      BKQR(3)=(Sqrt(D5*Pi)*(-M(3,1)))/D7 + (D5*Sqrt(D3*Pi)*(-M(4,2) - M(5,3)))/D14 + (D5*Sqrt(D5*Pi)*(-M(6,4) - M(7,5)))/D14
      BKQI(3)=(Sqrt(D5*Pi)*( M(2,1)))/D7 + (D5*Sqrt(D3*Pi)*( M(4,3) - M(5,2)))/D14 + (D5*Sqrt(D5*Pi)*( M(6,5) - M(7,4)))/D14
!      CF(4)=(Sqrt(D5*D6*Pi)*(-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2))/D7 + (D5*Sqrt(D2*Pi)*((0,1)*M(4,1) - M(5,1)))/D7 + (D5*Sqrt(Pi/D2)*(-M(6,2) + (0,1)*M(6,3) - (0,1)*M(7,2) - M(7,3)))/D7
      BKQR(4)=Sqrt(D5*D6*Pi)*(-M(2,2)/D2 + M(3,3)/D2)/D7 + D5*Sqrt(D2*Pi)*(-M(5,1))/D7 + D5*Sqrt(Pi/D2)*(-M(6,2) - M(7,3))/D7
      BKQI(4)=Sqrt(D5*D6*Pi)*(-M(2,3))/D7 + D5*Sqrt(D2*Pi)*(M(4,1))/D7 + D5*Sqrt(Pi/D2)*(M(6,3) - M(7,2))/D7
! B40
      BKQR(5)=-(Sqrt(Pi)*(M(4,4) + M(5,5))) + (Sqrt(Pi)*(D6*M(1,1) + M(2,2) + M(3,3) + D3*M(6,6) + D3*M(7,7)))/D7
      BKQI(5)=D0      
!      CF(6)=(Sqrt(D5*D6*Pi)*((0,1)*M(2,1) - M(3,1)))/D7 + (D4*Sqrt(D2*Pi)*(-M(4,2) + (0,1)*M(4,3) - (0,1)*M(5,2) - M(5,3)))/D7 + (Sqrt(D5*D6*Pi)*(M(6,4) - (0,1)*M(6,5) + (0,1)*M(7,4) + M(7,5)))/D7
      BKQR(6)=(Sqrt(D5*D6*Pi)*(-M(3,1)))/D7 + (D4*Sqrt(D2*Pi)*(-M(4,2) - M(5,3)))/D7 + (Sqrt(D5*D6*Pi)*( M(6,4) + M(7,5)))/D7
      BKQI(6)=(Sqrt(D5*D6*Pi)*( M(2,1)))/D7 + (D4*Sqrt(D2*Pi)*( M(4,3) - M(5,2)))/D7 + (Sqrt(D5*D6*Pi)*(-M(6,5) + M(7,4)))/D7
!      CF(7)=(D2*Sqrt(D10*Pi)*(-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2))/D7 + (Sqrt(D6*Pi)*((0,1)*M(4,1) - M(5,1)))/D7 + (3*Sqrt(D6*Pi)*(M(6,2) - (0,1)*M(6,3) + (0,1)*M(7,2) + M(7,3)))/D7
      BKQR(7)=(D2*Sqrt(D10*Pi)*(-M(2,2)/D2 + M(3,3)/D2))/D7 + (Sqrt(D6*Pi)*(-M(5,1)))/D7 + (D3*Sqrt(D6*Pi)*(M(6,2) + M(7,3)))/D7
      BKQI(7)=(D2*Sqrt(D10*Pi)*(-M(2,3)))/D7 + Sqrt(D6*Pi)/D7*(M(4,1) + D3*(-M(6,3) + M(7,2)))
!      CF(8)=Sqrt((D2*Pi)/D7)*(M(4,2) + (0,1)*M(4,3) + (0,1)*M(5,2) - M(5,3)) + D3*Sqrt((D2*Pi)/D7)*((0,-1)*M(6,1) + M(7,1))
      BKQR(8)=Sqrt((D2*Pi)/D7)*(M(4,2) - M(5,3)) + D3*Sqrt((D2*Pi)/D7)*(M(7,1))
      BKQI(8)=Sqrt((D2*Pi)/D7)*(M(4,3) + M(5,2)) - D3*Sqrt((D2*Pi)/D7)*(M(6,1))
!      CF(9)=Sqrt((D10*Pi)/D7)* (-M(4,4)/D2 - (0,1)*M(4,5) + M(5,5)/D2) + Sqrt((D6*Pi)/D7)*(M(6,2) + (0,1)*M(6,3) + (0,1)*M(7,2) - M(7,3))
      BKQR(9)=Sqrt((D10*Pi)/D7)*(-M(4,4)/D2 + M(5,5)/D2) + Sqrt((D6*Pi)/D7)*(M(6,2) - M(7,3))
      BKQI(9)=Sqrt((D10*Pi)/D7)*(-M(5,4)) + Sqrt((D6*Pi)/D7)*(M(6,3) + M(7,2))
! B60
      BKQR(10)=Sqrt(D13*Pi)/D7*(D2*M(1,1) - D3*(M(2,2) + M(3,3))/D2 + D3*(M(4,4) + M(5,5))/D5 - (M(6,6) + M(7,7))/(D2*D5))
      BKQI(10)=D0
!      CF(11)=Sqrt((D13*Pi)/D7)*((0,1)*M(2,1) - M(3,1)) + (Sqrt((39*Pi)/35.)*(M(4,2) - (0,1)*M(4,3) + (0,1)*M(5,2) + M(5,3)))/D2 + (Sqrt((D13*Pi)/D7)*(-M(6,4) + (0,1)*M(6,5) - (0,1)*M(7,4) - M(7,5)))/(D2*D5)
      BKQR(11)=Sqrt((D13*Pi)/D7)*(-M(3,1) + Sqrt(D3/D5)*( M(4,2) + M(5,3))/D2 + (-M(6,4) - M(7,5))/(D2*D5))
      BKQI(11)=Sqrt((D13*Pi)/D7)*( M(2,1) + Sqrt(D3/D5)*(-M(4,3) + M(5,2))/D2 + ( M(6,5) - M(7,4))/(D2*D5))      
!      CF(12)=Sqrt((D3*D13*Pi)/(D5*D7))* (-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2) + (D4*Sqrt((D13*Pi)/D7)*((0,-1)*M(4,1) + M(5,1)))/D5 + (Sqrt((D13*Pi)/D7)*(-M(6,2) + (0,1)*M(6,3) - (0,1)*M(7,2) - M(7,3)))/D5
      BKQR(12)=Sqrt((D13*Pi)/D7)*(-Sqrt(D3/D5)*(M(2,2) - M(3,3))/D2 + D4*(M(5,1))/D5 - (M(6,2) + M(7,3))/D5)
      BKQI(12)=Sqrt((D13*Pi)/D7)*(-Sqrt(D3/D5)*M(3,2) - D4*(M(4,1))/D5 + (M(6,3) - M(7,2))/D5)
!      CF(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,2) + (0,1)*M(4,3) + (0,1)*M(5,2) - M(5,3)))/D5 + (Sqrt((D6*D13*Pi)/D7)*((0,1)*M(6,1) - M(7,1)))/D5
      BKQR(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,2) - M(5,3)))/D5 + (Sqrt((D6*D13*Pi)/D7)*(-M(7,1)))/D5
      BKQI(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,3) + M(5,2)))/D5 + (Sqrt((D6*D13*Pi)/D7)*( M(6,1)))/D5      
!      CF(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(4,4)/D2 - (0,1)*M(4,5) + M(5,5)/D2))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,2) - (0,1)*M(6,3) - (0,1)*M(7,2) + M(7,3))
      BKQR(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(4,4)/D2 + M(5,5)/D2))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,2) + M(7,3))
      BKQI(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(5,4)))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,3) - M(7,2))
!      CF(15)=(Sqrt((D3*D11*D13*Pi)/D14)* (M(6,4) + (0,1)*M(6,5) + (0,1)*M(7,4) - M(7,5)))/D5
      BKQR(15)=(Sqrt((D3*D11*D13*Pi)/D14)*(M(6,4) - M(7,5)))/D5
      BKQI(15)=(Sqrt((D3*D11*D13*Pi)/D14)*(M(6,5) + M(7,4)))/D5    
!      CF(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(6,6)/D2 - (0,1)*M(6,7) + M(7,7)/D2))/D5
      BKQR(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(6,6)/D2 + M(7,7)/D2))/D5
      BKQI(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(7,6)))/D5
!   Adjust phases so the Bkq and Bkq' correspond to those of Goeller-Walrand.
      Do k=0,6,2
        do q=0,k
          i=Bkq_index(k,q)
          BKQR(i)= sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQR(i)
          BKQI(i)=-sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQI(i)  ! minus sign added 6/2/15
        enddo
      enddo
      AOMchanged=.false.
!      write(io,'("AOMmatrixF: eSigma(1),eSigma(2),BKQR(2)=",3F10.4)') eSigma(1),eSigma(2),BKQR(2)
!
      return
      end subroutine AOMmatrixF
!
!-----------------------------------------------------------------------
!
      SUBROUTINE AOMmatrixF1()
!****old replaced by: AOMmatrixF *****
! Calculates the AOM matrix for f-orbitals, and the equivalent crystal field parameters.
!  TODO add expression for non-zero Chi values (anisotropic pi bonding).
      IMPLICIT none
      INTEGER*4 iL,i,j,k,q
      REAL*8 T,P,Ps,Pi,fv1(7),E1,sum1
      REAL*8 M(7,7),M1(7,7),F(7,7),D0,D1,D2,D3,D4,D5,D6,D7,D8,D10,D11,D13,D14,D15
      PARAMETER(D0=0.0D+00,D1=1.0D+00,D2=2.0D+00,D3=3.0D+00,D4=4.0D+00,D5=5.0D+00,D6=6.0D+00,D7=7.0D+00)
      PARAMETER(D8=8.0D+00,D10=10.0D+00,D11=11.0D+00,D13=13.0D+00,D14=14.0D+00,D15=15.0D+00)
!
      Pi=D4*ATAN(D1)
      do i=1,7
        do j=1,7
           M(i,j) = D0
        enddo
      enddo  
      do iL=1,NLIGANDS
        T=Theta(iL)*Pi/180.d+00
        P=Phi(iL)*Pi/180.d+00
        Ps=Chi(iL)*Pi/180.d+00
!  The <f|V|L> elements from Table 2: Urland,Chem.Phys.14,393,(1976).
!  Note: isotropic pi bonding assumed. Only 2 AOM angles (theta,phi) are used.
!  Note: "Standard" ordering of real f-orbitals |sigma>, |piS>, |piC>, |deltaS>, |deltaC>, |phiS>, |phiC>.
!  For a definition of these orbitals: see Harnung & Schaffer, Struct&Bond,12,201,(1972).
        F(1,1)=+(Cos(T)*(-D3 + D5*Cos(T)**2))/D2
        F(1,2)= D0
        F(1,3)=+Sqrt(D6)*Sin(T)*(-D1 + (D5*Sin(T)**2)/D4)
!        F(1,4)= D0
!        F(1,5)=+(Sqrt(D15)*Sin(T)*Sin(D2*T))/D4
!        F(1,6)= D0
!        F(1,7)=-(Sqrt(D5/D2)*Sin(T)**3)/D2
        F(2,1)=+Sqrt(D6)*Sin(P)*Sin(T)*(D1 - (D5*Sin(T)**2)/D4)
        F(2,2)=+(Cos(P)*(-D1 + D5*Cos(T)**2))/D4
        F(2,3)=+(Cos(T)*(-D11 + D15*Cos(T)**2)*Sin(P))/D4
!        F(2,4)=-(Sqrt(D5/D2)*Cos(P)*Sin(D2*T))/D2
!        F(2,5)=+Sqrt(D5/D2)*Sin(P)*Sin(T)*(-D1 + (D3*Sin(T)**2)/D2)
!        F(2,6)=+(Sqrt(D15)*Cos(P)*Sin(T)**2)/D4
!        F(2,7)=+(Sqrt(D15)*Sin(P)*Sin(T)*Sin(D2*T))/D8
        F(3,1)=+Sqrt(D6)*Cos(P)*Sin(T)*(D1 - (D5*Sin(T)**2)/D4)
        F(3,2)=-((-D1 + D5*Cos(T)**2)*Sin(P))/D4
        F(3,3)=+(Cos(P)*Cos(T)*(-D11 + D15*Cos(T)**2))/D4
!        F(3,4)=+(Sqrt(D5/D2)*Sin(P)*Sin(D2*T))/D2
!        F(3,5)=+Sqrt(D5/D2)*Cos(P)*Sin(T)*(-D1 + (D3*Sin(T)**2)/D2)
!        F(3,6)=-(Sqrt(D15)*Sin(P)*Sin(T)**2)/D4
!        F(3,7)=+(Sqrt(D15)*Cos(P)*Sin(T)*Sin(D2*T))/D8
        F(4,1)=+(Sqrt(D15)*Sin(D2*P)*Sin(T)*Sin(D2*T))/D4
        F(4,2)=+(Sqrt(D5/D2)*Cos(D2*P)*Sin(D2*T))/D2
        F(4,3)=+Sqrt(D5/D2)*Sin(D2*P)*Sin(T)*(Cos(T)**2 - Sin(T)**2/D2)
!        F(4,4)=+Cos(D2*P)*Cos(D2*T)
!        F(4,5)=+(Cos(T)*(-D1 + D3*Cos(T)**2)*Sin(D2*P))/D2
!        F(4,6)=-(Sqrt(D3/D2)*Cos(D2*P)*Sin(D2*T))/D2
!        F(4,7)=+Sqrt(D3/D2)*Sin(D2*P)*Sin(T)*(-D1 + Sin(T)**2/D2)
        F(5,1)=+(Sqrt(D15)*Cos(D2*P)*Sin(T)*Sin(D2*T))/D4
        F(5,2)=-(Sqrt(D5/D2)*Sin(D2*P)*Sin(D2*T))/D2
        F(5,3)=+Sqrt(D5/D2)*Cos(D2*P)*Sin(T)*(Cos(T)**2 - Sin(T)**2/D2)
!        F(5,4)=-(Cos(D2*T)*Sin(D2*P))
!        F(5,5)=+(Cos(D2*P)*Cos(T)*(-D1 + D3*Cos(T)**2))/D2
!        F(5,6)=+(Sqrt(D3/D2)*Sin(D2*P)*Sin(D2*T))/D2
!        F(5,7)=+Sqrt(D3/D2)*Cos(D2*P)*Sin(T)*(-D1 + Sin(T)**2/D2)
        F(6,1)=+(Sqrt(D5/D2)*Sin(D3*P)*Sin(T)**3)/D2
        F(6,2)=+(Sqrt(D15)*Cos(D3*P)*Sin(T)**2)/D4
        F(6,3)=+(Sqrt(D15)*Sin(D3*P)*Sin(T)*Sin(D2*T))/D8
!        F(6,4)=+(Sqrt(D3/D2)*Cos(D3*P)*Sin(D2*T))/D2
!        F(6,5)=+Sqrt(D3/D2)*Sin(D3*P)*(Sin(T) - Sin(T)**3/D2)
!        F(6,6)=+Cos(D3*P)*(D1 - (D3*Sin(T)**2)/D4)
!        F(6,7)=+((D3*Cos(T) + Cos(T)**3)*Sin(D3*P))/D4
        F(7,1)=+(Sqrt(D5/D2)*Cos(D3*P)*Sin(T)**3)/D2
        F(7,2)=-(Sqrt(D15)*Sin(D3*P)*Sin(T)**2)/D4
        F(7,3)=+(Sqrt(D15)*Cos(D3*P)*Sin(T)*Sin(D2*T))/D8
!        F(7,4)=-(Sqrt(D3/D2)*Sin(D3*P)*Sin(D2*T))/D2
!        F(7,5)=+Sqrt(D3/D2)*Cos(D3*P)*(Sin(T) - Sin(T)**3/D2)
!        F(7,6)=-(Sin(D3*P)*(D1 - (D3*Sin(T)**2)/D4))
!        F(7,7)=+(Cos(D3*P)*(D3*Cos(T) + Cos(T)**3))/D4
        do i=1,7
          do j=1,7
            M(i,j) = M(i,j) + eSigma(iL)*F(i,1)*F(j,1) + ePiX(iL)*F(i,2)*F(j,2) + ePiY(iL)*F(i,3)*F(j,3)
            M1(i,j)= M(i,j) ! M1 gets destroyed below
          enddo
        enddo 
      enddo  !  iL=1,NLIGANDS
 
      if (check1(8)) then
        call diagrs(M1,7,7,engs,fv1,1,io)
        E1=EAVE; EAVE=D0; call zeroE(OrbEngs,7); EAVE=E1
        sum1=(M(1,1)+M(2,2)+M(3,3)+M(4,4)+M(5,5)+M(6,6)+M(7,7))/D7
        do i=1,7; M(i,i)=M(i,i)-sum1; enddo  ! subtract baricentre
        write(io,'(/," AOM matrix in terms of the 7x7 real f-orbital matrix (baricentre=",F12.2," subtracted)")') sum1
        write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>      |psiS>      |psiC>     orbital energies")') 
        write(io,'(2X," <sig|",7F12.2,5X,F12.2)') (M(1,j),j=1,7),OrbEngs(7)
        write(io,'(2X," <piS|",7F12.2,5X,F12.2)') (M(2,j),j=1,7),OrbEngs(6)
        write(io,'(2X," <piC|",7F12.2,5X,F12.2)') (M(3,j),j=1,7),OrbEngs(5)
        write(io,'(2X,"<delS|",7F12.2,5X,F12.2)') (M(4,j),j=1,7),OrbEngs(4)
        write(io,'(2X,"<delC|",7F12.2,5X,F12.2)') (M(5,j),j=1,7),OrbEngs(3)
        write(io,'(2X,"<psiS|",7F12.2,5X,F12.2)') (M(6,j),j=1,7),OrbEngs(2)
        write(io,'(2X,"<psiS|",7F12.2,5X,F12.2)') (M(7,j),j=1,7),OrbEngs(1)
      Endif
!
!  Relationships below are given by Urland, Chem.Phys.14, 393,(1976). Table 3.
!  These are the Ckq coefficents of the ligand field expanded as spherical harmonics.
! B00 
      BKQR(1)=(D2*Sqrt(Pi)*(M(1,1) + M(2,2) + M(3,3) + M(4,4) + M(5,5) + M(6,6) + M(7,7)))/D7
      BKQI(1)=D0
! B20     
      BKQR(2)=(D2*Sqrt(D5*Pi)*M(1,1))/D7 + (D3*Sqrt(D5*Pi)*(M(2,2) + M(3,3)))/D14 - (D5*Sqrt(D5*Pi)*(M(6,6) + M(7,7)))/D14
      BKQI(2)=D0
!      CF(3)=(Sqrt(D5*Pi)*((0,1)*M(2,1) - M(3,1)))/D7 + (D5*Sqrt(D3*Pi)*(-M(4,2) + (0,1)*M(4,3) - (0,1)*M(5,2) - M(5,3)))/D14 + (D5*Sqrt(D5*Pi)*(-M(6,4) + (0,1)*M(6,5) - (0,1)*M(7,4) - M(7,5)))/D14
      BKQR(3)=(Sqrt(D5*Pi)*(-M(3,1)))/D7 + (D5*Sqrt(D3*Pi)*(-M(4,2) - M(5,3)))/D14 + (D5*Sqrt(D5*Pi)*(-M(6,4) - M(7,5)))/D14
      BKQI(3)=(Sqrt(D5*Pi)*( M(2,1)))/D7 + (D5*Sqrt(D3*Pi)*( M(4,3) - M(5,2)))/D14 + (D5*Sqrt(D5*Pi)*( M(6,5) - M(7,4)))/D14
!      CF(4)=(Sqrt(D5*D6*Pi)*(-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2))/D7 + (D5*Sqrt(D2*Pi)*((0,1)*M(4,1) - M(5,1)))/D7 + (D5*Sqrt(Pi/D2)*(-M(6,2) + (0,1)*M(6,3) - (0,1)*M(7,2) - M(7,3)))/D7
      BKQR(4)=Sqrt(D5*D6*Pi)*(-M(2,2)/D2 + M(3,3)/D2)/D7 + D5*Sqrt(D2*Pi)*(-M(5,1))/D7 + D5*Sqrt(Pi/D2)*(-M(6,2) - M(7,3))/D7
      BKQI(4)=Sqrt(D5*D6*Pi)*(-M(2,3))/D7 + D5*Sqrt(D2*Pi)*(M(4,1))/D7 + D5*Sqrt(Pi/D2)*(M(6,3) - M(7,2))/D7
! B40
      BKQR(5)=-(Sqrt(Pi)*(M(4,4) + M(5,5))) + (Sqrt(Pi)*(D6*M(1,1) + M(2,2) + M(3,3) + D3*M(6,6) + D3*M(7,7)))/D7
      BKQI(5)=D0      
!      CF(6)=(Sqrt(D5*D6*Pi)*((0,1)*M(2,1) - M(3,1)))/D7 + (D4*Sqrt(D2*Pi)*(-M(4,2) + (0,1)*M(4,3) - (0,1)*M(5,2) - M(5,3)))/D7 + (Sqrt(D5*D6*Pi)*(M(6,4) - (0,1)*M(6,5) + (0,1)*M(7,4) + M(7,5)))/D7
      BKQR(6)=(Sqrt(D5*D6*Pi)*(-M(3,1)))/D7 + (D4*Sqrt(D2*Pi)*(-M(4,2) - M(5,3)))/D7 + (Sqrt(D5*D6*Pi)*( M(6,4) + M(7,5)))/D7
      BKQI(6)=(Sqrt(D5*D6*Pi)*( M(2,1)))/D7 + (D4*Sqrt(D2*Pi)*( M(4,3) - M(5,2)))/D7 + (Sqrt(D5*D6*Pi)*(-M(6,5) + M(7,4)))/D7
!      CF(7)=(D2*Sqrt(D10*Pi)*(-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2))/D7 + (Sqrt(D6*Pi)*((0,1)*M(4,1) - M(5,1)))/D7 + (3*Sqrt(D6*Pi)*(M(6,2) - (0,1)*M(6,3) + (0,1)*M(7,2) + M(7,3)))/D7
      BKQR(7)=(D2*Sqrt(D10*Pi)*(-M(2,2)/D2 + M(3,3)/D2))/D7 + (Sqrt(D6*Pi)*(-M(5,1)))/D7 + (D3*Sqrt(D6*Pi)*(M(6,2) + M(7,3)))/D7
      BKQI(7)=(D2*Sqrt(D10*Pi)*(-M(2,3)))/D7 + Sqrt(D6*Pi)/D7*(M(4,1) + D3*(-M(6,3) + M(7,2)))
!      CF(8)=Sqrt((D2*Pi)/D7)*(M(4,2) + (0,1)*M(4,3) + (0,1)*M(5,2) - M(5,3)) + D3*Sqrt((D2*Pi)/D7)*((0,-1)*M(6,1) + M(7,1))
      BKQR(8)=Sqrt((D2*Pi)/D7)*(M(4,2) - M(5,3)) + D3*Sqrt((D2*Pi)/D7)*(M(7,1))
      BKQI(8)=Sqrt((D2*Pi)/D7)*(M(4,3) + M(5,2)) - D3*Sqrt((D2*Pi)/D7)*(M(6,1))
!      CF(9)=Sqrt((D10*Pi)/D7)* (-M(4,4)/D2 - (0,1)*M(4,5) + M(5,5)/D2) + Sqrt((D6*Pi)/D7)*(M(6,2) + (0,1)*M(6,3) + (0,1)*M(7,2) - M(7,3))
      BKQR(9)=Sqrt((D10*Pi)/D7)*(-M(4,4)/D2 + M(5,5)/D2) + Sqrt((D6*Pi)/D7)*(M(6,2) - M(7,3))
      BKQI(9)=Sqrt((D10*Pi)/D7)*(-M(4,5)) + Sqrt((D6*Pi)/D7)*(M(6,3) + M(7,2))
! B60
      BKQR(10)=Sqrt(D13*Pi)/D7*(D2*M(1,1) - D3*(M(2,2) + M(3,3))/D2 + D3*(M(4,4) + M(5,5))/D5 - (M(6,6) + M(7,7))/(D2*D5))
      BKQI(10)=D0
!      CF(11)=Sqrt((D13*Pi)/D7)*((0,1)*M(2,1) - M(3,1)) + (Sqrt((39*Pi)/35.)*(M(4,2) - (0,1)*M(4,3) + (0,1)*M(5,2) + M(5,3)))/D2 + (Sqrt((D13*Pi)/D7)*(-M(6,4) + (0,1)*M(6,5) - (0,1)*M(7,4) - M(7,5)))/(D2*D5)
      BKQR(11)=Sqrt((D13*Pi)/D7)*(-M(3,1) + Sqrt(D3/D5)*( M(4,2) + M(5,3))/D2 + (-M(6,4) - M(7,5))/(D2*D5))
      BKQI(11)=Sqrt((D13*Pi)/D7)*( M(2,1) + Sqrt(D3/D5)*(-M(4,3) + M(5,2))/D2 + ( M(6,5) - M(7,4))/(D2*D5))      
!      CF(12)=Sqrt((D3*D13*Pi)/(D5*D7))* (-M(2,2)/D2 - (0,1)*M(2,3) + M(3,3)/D2) + (D4*Sqrt((D13*Pi)/D7)*((0,-1)*M(4,1) + M(5,1)))/D5 + (Sqrt((D13*Pi)/D7)*(-M(6,2) + (0,1)*M(6,3) - (0,1)*M(7,2) - M(7,3)))/D5
      BKQR(12)=Sqrt((D13*Pi)/D7)*(-Sqrt(D3/D5)*(M(2,2) - M(3,3))/D2 + D4*(M(5,1))/D5 - (M(6,2) + M(7,3))/D5)
      BKQI(12)=Sqrt((D13*Pi)/D7)*(-Sqrt(D3/D5)*M(2,3) - D4*(M(4,1))/D5 + (M(6,3) - M(7,2))/D5)
!      CF(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,2) + (0,1)*M(4,3) + (0,1)*M(5,2) - M(5,3)))/D5 + (Sqrt((D6*D13*Pi)/D7)*((0,1)*M(6,1) - M(7,1)))/D5
      BKQR(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,2) - M(5,3)))/D5 + (Sqrt((D6*D13*Pi)/D7)*(-M(7,1)))/D5
      BKQI(13)=(D3*Sqrt((D3*D13*Pi)/D14)*(M(4,3) + M(5,2)))/D5 + (Sqrt((D6*D13*Pi)/D7)*( M(6,1)))/D5      
!      CF(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(4,4)/D2 - (0,1)*M(4,5) + M(5,5)/D2))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,2) - (0,1)*M(6,3) - (0,1)*M(7,2) + M(7,3))
      BKQR(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(4,4)/D2 + M(5,5)/D2))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,2) + M(7,3))
      BKQI(14)=(D3*Sqrt((D2*D13*Pi)/D7)*(-M(4,5)))/D5 + Sqrt((D3*D13*Pi)/(D2*D5*D7))*(-M(6,3) - M(7,2))
!      CF(15)=(Sqrt((D3*D11*D13*Pi)/D14)* (M(6,4) + (0,1)*M(6,5) + (0,1)*M(7,4) - M(7,5)))/D5
      BKQR(15)=(Sqrt((D3*D11*D13*Pi)/D14)*(M(6,4) - M(7,5)))/D5
      BKQI(15)=(Sqrt((D3*D11*D13*Pi)/D14)*(M(6,5) + M(7,4)))/D5    
!      CF(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(6,6)/D2 - (0,1)*M(6,7) + M(7,7)/D2))/D5
      BKQR(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(6,6)/D2 + M(7,7)/D2))/D5
      BKQI(16)=(Sqrt((D3*D11*D13*Pi)/D7)*(-M(6,7)))/D5
      
!   Adjust phases so the Bkq and Bkq' correspond to those of Goeller-Walrand.
      Do k=0,6,2
        do q=0,k
          i=Bkq_index(k,q)
          BKQR(i)= sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQR(i)
          BKQI(i)=-sqrt((D2*k+D1)/(D4*Pi))*(-D1)**q*BKQI(i)  ! minus sign added 6/2/15
        enddo
      enddo
      AOMchanged=.false.
!      write(io,'("AOMmatrixF: eSigma(1),eSigma(2),BKQR(2)=",3F10.4)') eSigma(1),eSigma(2),BKQR(2)
!
      return
      end subroutine AOMmatrixF1
!
!-----------------------------------------------------------------------
!
      subroutine convertToAOM(DistMat,AngMat,DistLab,Distbonds)
! 
!  Converts the AOM angles theta & phi into x, y, z and 
!  the x, y, z into theta, phi and bondlength.
!
      IMPLICIT none
      real*8 DistMat(max_ligands,max_ligands), AngMat(max_ligands,max_ligands), Distbonds(max_ligands)
      real*8 y(3,3), w(3,3),s, DTR, Xnew,Ynew,Znew,xx
      character*1 DistLab(max_ligands,max_ligands)
      integer i,j,L,iL
      DTR=(4.0d+00*ATAN(1.0d+00))/180.0d+00  ! converts degrees to radians
!
      if (xref) then
!        if (AOM_angles) then; write(io,'("****FATAL: cannot use XREF with AOM specified by angles.")'); STOP; endif
!   (Should now be possible, xLig,yLig,zLig defined from angles)
        if (nref(1).lt.1 .or. nref(2).lt.1) then
          write(io,'("****FATAL: Ligand numbers too small in XREF command: N1,N2=",2i4)') nref(1),nref(2); STOP  
        endif
        if (nref(1).gt.NLIGANDS .or. nref(2).gt.NLIGANDS) then
          write(io,'("****FATAL: Ligand numbers too large in XREF command: N1,N2=",2i4)') nref(1),nref(2); STOP  
        endif
        Y(1,1)=xLig(nref(1)); Y(1,2)=yLig(nref(1)); Y(1,3)=zLig(nref(1));
        Y(2,1)=xLig(nref(2)); Y(2,2)=yLig(nref(2)); Y(2,3)=zLig(nref(2));
        do i=1,2
          S=SQRT(Y(I,1)**2+Y(I,2)**2+Y(I,3)**2)
          IF (S.LE.1.0D-10) THEN; WRITE(IO,*) 'Axes of coordinate system are colinear !'; STOP; ENDIF
          DO J=1,3
            W(I,J)=Y(I,J)/S
          enddo 
        enddo
        W(3,1)=W(1,2)*W(2,3)-W(1,3)*W(2,2); W(3,2)=W(1,3)*W(2,1)-W(1,1)*W(2,3); W(3,3)=W(1,1)*W(2,2)-W(1,2)*W(2,1)
        S=SQRT(W(3,1)**2+W(3,2)**2+W(3,3)**2)
        IF (S.LE.1.0D-10) THEN; WRITE(IO,*) 'Axes of coordinate system are colinear !'; STOP; ENDIF
        DO j=1,3; W(3,j)=W(3,j)/S; enddo
        W(2,1)=W(1,3)*W(3,2)-W(1,2)*W(3,3); W(2,2)=W(1,1)*W(3,3)-W(1,3)*W(3,1); W(2,3)=W(1,2)*W(3,1)-W(1,1)*W(3,2)
        DO I=1,3; DO J=1,3; L=MOD(I,3)+1; Y(I,J)=W(L,J); enddo; enddo
        write(IO,'("Transformation matrix from XREF",2I3)') Nref(1),Nref(2) 
        do i=1,3; write(IO,'(3f12.6)') (Y(i,j),j=1,3); enddo 
        write(IO,'("Ligand    Transformed Coords             Original Coords",/,  &
                   "          X        Y        Z            X        Y        Z" )')
        do il=1,Nligands
          Xnew=Y(1,1)*XLig(iL)+Y(1,2)*YLig(iL)+Y(1,3)*ZLig(iL)        
          Ynew=Y(2,1)*XLig(iL)+Y(2,2)*YLig(iL)+Y(2,3)*ZLig(iL)        
          Znew=Y(3,1)*XLig(iL)+Y(3,2)*YLig(iL)+Y(3,3)*ZLig(iL)
          write(IO,'(1X,A4,1X,3F9.5,4X,3F9.5)') Lname(iL),Xnew,Ynew,Znew,XLig(iL),YLig(iL),ZLig(iL)
          XLig(iL)=Xnew; YLig(iL)=Ynew; ZLig(iL)=Znew
          phi(i)=atan2(yLig(i),xLig(i))/DTR
          theta(i)=atan2(sqrt(xLig(i)**2+yLig(i)**2),zLig(i))/DTR
        ENDDO  
      endif  ! xref
!      
      do i=1,NLIGANDS
        Distbonds(i)=sqrt(xLig(i)**2+yLig(i)**2+zLig(i)**2)
        if (Distbonds(i).lt.0.1) write(io,'("****Warning: Bond length to atom:",A4," too small=",F5.3)') Lname(i),Distbonds(i)
        do j=1,i
          if (i.eq.j) then
            DistMat(i,j)=0.0d+00;  AngMat(i,j)=0.0d+00
            Distlab(i,j)=" "
          elseif (i.ne.j) then
            DistMat(i,j)=sqrt((xLig(i)-xLig(j))**2+(yLig(i)-yLig(j))**2+(zLig(i)-zLig(j))**2)
            Distlab(i,j)=" "
            xx=(Distbonds(i)**2+Distbonds(j)**2-DistMat(i,j)**2)/(2.0d+00*Distbonds(i)*Distbonds(j))
            if (xx.gt. 1.0d+00) then
              if (xx-1.0d+00.gt.1.0d-6) write(io,'("Geometry Warning, angle near 0, atoms",i2,"and",i2, &
                                                                                  "; xx >1; xx-1=",E16.8)') i,j,xx-1.0d+00
              xx= 1.0d+00;
            endif
            if (xx.lt.-1.0d+00) then
              if (xx+1.0d+00.lt.-1.0d-6) write(io,'("Geometry Warning, angle near 180, atoms",i2,"and",i2,  &
                                                                                  "; xx<-1; xx+1=",E16.8)') i,j,xx+1.0d+00
              xx=-1.0d+00
            endif
            AngMat(i,j)=ACOS(xx)/DTR
            if (DistMat(i,j).lt.0.25) then
              write(io,'("****Warning: distance between atoms:",A4," and ",A4," too close=",F5.3)') Lname(i),Lname(j),DistMat(i,j)
              Distlab(i,j)="*"
            endif
          endif ! i.eq.j   
        enddo ! j
      enddo ! i
!
      Return
      end subroutine convertToAOM
!      
!-----------------------------------------------------------------------
!
      SUBROUTINE RotateLF()
!
!  This subroutine rotates the ligand field. Since all other atomic parameters 
!  are symmetric with respect to orientation, rotating the Ligand Field is  
!  equivalent to rotating the crystal in fixed magnetic field and polarisation 
!  orientations.  
!
!  Ligand field in terms of spherical harmonics:
!     k   q    Bkq
!     0   0   BKQR(1)
!     2   0   BKQR(2)
!     2   1   BKQR(3)  BKQI(3)
!     2   2   BKQR(4)  BKQI(4)
!     4   0   BKQR(5)
!     4   1   BKQR(6)  BKQI(6)
!     4   2   BKQR(7)  BKQI(7)
!     4   3   BKQR(8)  BKQI(8)
!     4   4   BKQR(9)  BKQI(9)
!     6   0   BKQR(10)
!     6   1   BKQR(11)  BKQI(11)
!     6   2   BKQR(12)  BKQI(12)
!     6   3   BKQR(13)  BKQI(13)
!     6   4   BKQR(14)  BKQI(14)
!     6   5   BKQR(15)  BKQI(15)
!     6   6   BKQR(16)  BKQI(16)
!
      IMPLICIT none
      INTEGER*4 i,j,k,q, mp,mm
      REAL*8 TR(9),TI(9),a,b,g,PI,DTR
      real*8 Z,D1,D2,D3,D4,D5,D6,D7,D8,D9, D11,D14,D15,D16,D20,D21,D28,D30,D32,D35,D64,  R2,R3,R5,R7,R14,R35
      real*8 CB2,SB2,CB,SB,C2B,C3B,C4B
      COMPLEX*16 V(13),VN(13),M(13,13),CZ,CII 
      LOGICAL DEBUG
      PARAMETER(  Z=0.0d+00,  D1=1.0d+00,  D2=2.0d+00,   D3=3.0d+00,  D4=4.0d+00,  D5=5.0d+00,  D6=6.0d+00,   D7=7.0d+00)
      PARAMETER( D8=8.0d+00,  D9=9.0d+00, D11=11.0d+00, D14=14.0d+00,D15=15.0d+00,D16=16.0d+00,D20=20.0d+00,D21=21.0d+00)
      PARAMETER(D28=28.0d+00,D30=30.0d+00,D32=32.0d+00, D35=35.0d+00,D64=64.0d+00)
!
      a=RotLF(1); b=RotLF(2); g=RotLF(3)  ! alpha, beta, gamma
      DEBUG=.false.
      CZ=DCMPLX(Z, Z);     CII=DCMPLX(Z, 1.0d+00)
      PI=4.0d+00*ATAN(1.0d+00);  DTR = PI/180.0D+00
      if (debug) write(io,'("Ligand field rotation:  alpha,beta,gamma=",3F14.4)') a,b,g   
      a=-a*DTR;   b=-b*DTR; g=-g*DTR
!   
      if (OUTP(5).or.debug) then 
        WRITE(IO,'(" Before rotation, Crystal field parameters ")')
        do k=2,2*Lvalue,2
          do q=0,k
            i=Bkq_index(k,q)
            if (q.eq.0) then 
              WRITE(IO,'(7X,A4,F14.6)')    BkqR_label(i),BkqR(i)
            else if (q.ne.0) then
              WRITE(IO,'(2(7X,A4,F14.6))') BkqR_label(i),BkqR(i),BkqI_label(i),BkqI(i) 
            endif
          enddo
        enddo
      endif
!
      R2=SQRT(D2);  R3=SQRT(D3);  R5=SQRT(5.0d+00);  R7=SQRT(7.0d+00);  R14=SQRT(14.0D+00);  R35=SQRT(35.0d+00)
      IF (ABS(B).gt.1.0d-6 .and. ABS(B-PI).gt.1.0d-6) then 
        CB2=Cos(B/D2); SB2=Sin(B/D2); CB=Cos(B); SB=Sin(B); C2B=Cos(D2*B);  C3B=Cos(D3*B); C4B=Cos(D4*B)
        M(1,1)=CB2**4;                 M(1,2)=D2*CB2**3*SB2;        M(1,3)=(Sqrt(D3/D2)*SB**2)/D2 
        M(1,4)=SB2**2*SB;              M(1,5)=SB2**4
        M(2,1)=-D2*CB2**3*SB2;         M(2,2)=CB2**2*(-D1+D2*CB);   M(2,3)=Sqrt(D3/D2)*CB*SB      
        M(2,4)=(D1+D2*CB)*SB2**2;      M(2,5)=SB2**2*SB
        M(3,1)=(Sqrt(D3/D2)*SB**2)/D2; M(3,2)=-(Sqrt(D3/D2)*CB*SB); M(3,3)=(D1+D3*C2B)/D4 
        M(3,4)=Sqrt(D3/D2)*CB*SB;      M(3,5)=(Sqrt(D3/D2)*SB**2)/D2
        M(4,1)=-(SB2**2*SB);           M(4,2)=(D1 + D2*CB)*SB2**2;  M(4,3)=-(Sqrt(D3/D2)*CB*SB)   
        M(4,4)=CB2**2*(-D1+D2*CB);     M(4,5)=D2*CB2**3*SB2
        M(5,1)=SB2**4;                 M(5,2)=-(SB2**2*SB);         M(5,3)=(Sqrt(D3/D2)*SB**2)/D2  
        M(5,4)=-D2*CB2**3*SB2;         M(5,5)=CB2**4
      ELSE
        do i=1,5
          mp=-2+i-1
          do j=1,5
            M(i,j)=CZ
          enddo
          if (ABS(B).le.1.0d-6) M(i,i)=DCMPLX(D1,Z)
          if (ABS(B-PI).le.1.0d-6) M(i,5-i+1)=(-D1)**(2-mp)*DCMPLX(D1,Z)
        enddo
      ENDIF

      if (debug) then 
        write(io,'(/,"This is as Zare (3.54) d2mpm   b=",F8.2)') b/dtr
        write(io,'(8X,5("  |",I2,">  "))') ((-2+j-1),j=1,5)
        do i=1,5
          mp=-2+i-1
          write(io,'("  <",I2,"|  ",5F8.4)') mp,(DBLE(M(i,j)),j=1,5)
        enddo
      endif
!
!  Two things are happening here:
!  1. Use the transpose of M so you can post-multiply 
!     a column vector rather then pre-multipy a row vector for the kets. 
!  2. Use the inverse of M (since unitary, this is the conjugate transpose).
!     Inverse used because want to describe the rotated functions (Ykq) in terms 
!     of the original coordinate system, rather than in terms of the original functions.   

      do I=1,5
        do J=1,5
          mp=-2+I-1
          mm=-2+J-1
          M(I,J)=CDExp(-CII*mp*a)*M(I,J)*CDExp(-CII*mm*g)
        enddo
      enddo  
!
      V(1)=DCMPLX( BKQR(4), BKQI(4)); V(2)=DCMPLX( BKQR(3), BKQI(3)); V(3)=DCMPLX( BKQR(2), Z)
      V(4)=DCMPLX(-BKQR(3), BKQI(3)); V(5)=DCMPLX( BKQR(4),-BKQI(4))
      do i=1,5
        VN(i)=CZ
        do k=1,5
!          VN(i) = VN(i) + M(k,i)*v(k)
          VN(i) = VN(i) + DCONJG(M(i,k))*v(k)
        enddo
      enddo
      BKQR(2)=DBLE(VN(3));   BKQR(3)=DBLE(VN(2));    BKQR(4)=DBLE(VN(1));   BKQI(3)=DIMAG(VN(2));    BKQI(4)=DIMAG(VN(1))
!
      IF (ABS(B).gt.1.0d-6 .and. ABS(B-PI).gt.1.0d-6) then 
        M(1,1)=CB2**8;                                   M(1,2)=D2*R2*CB2**7*SB2
        M(1,3)=(R7*(D1/SB2**4)*SB**6)/D32;               M(1,4)=(Sqrt(D7/D2)*(D1/SB2**2)*SB**5)/D8
        M(1,5)=(Sqrt(D35/D2)*SB**4)/D8;                  M(1,6)=(Sqrt(D7/D2)*SB2**2*SB**3)/D2
        M(1,7)=(R7*SB2**4*SB**2)/D2;                     M(1,8)=R2*SB2**6*SB
        M(1,9)=SB2**8
        M(2,1)=-D2*R2*CB2**7*SB2;                        M(2,2)=CB2**6*(-D3 + D4*CB)
        M(2,3)=R14*CB2**5*(-D1 + D2*CB)*SB2;             M(2,4)=R7*CB2**4*(-D1 + D4*CB)*SB2**2
        M(2,5)=(R35*CB*SB**3)/D4;                        M(2,6)=R7*CB2**2*(D1 + D4*CB)*SB2**4
        M(2,7)=R14*CB2*(D1 + D2*CB)*SB2**5;              M(2,8)=(D3 + D4*CB)*SB2**6
        M(2,9)=R2*SB2**6*SB
        M(3,1)=(R7*(D1/SB2**4)*SB**6)/D32;               M(3,2)=-(R14*CB2**5*(-D1+D2*CB)*SB2)
        M(3,3)=-(CB2**4*(-D9 + D14*CB - D7*C2B))/D2;     M(3,4)=-((CB2**3*(-D6+D7*CB-D7*C2B)*SB2)/R2)
        M(3,5)=(Sqrt(D5/D2)*(D5+D7*C2B)*SB**2)/D8;       M(3,6)=(CB2*(D6 + D7*CB+D7*C2B)*SB2**3)/R2
        M(3,7)=(D1 + D7*CB + D7*CB**2)*SB2**4;           M(3,8)=R14*CB2*(D1 + D2*CB)*SB2**5
        M(3,9)=(R7*SB2**4*SB**2)/D2
        M(4,1)=-(Sqrt(D7/D2)*(D1/SB2**2)*SB**5)/D8;      M(4,2)=R7*CB2**4*(-D1 + D4*CB)*SB2**2
        M(4,3)=(CB2**3*(-D6+D7*CB-D7*C2B)*SB2)/R2;       M(4,4)=(CB2**2*(-D15+D30*CB-D21*C2B+D14*C3B))/D8
        M(4,5)=(R5*(D9*CB + D7*C3B)*SB)/D16;             M(4,6)=((D15 + D30*CB+D21*C2B+D14*C3B)*SB2**2)/D8
        M(4,7)=(CB2*(D6+D7*CB+D7*C2B)*SB2**3)/R2;        M(4,8)=R7*CB2**2*(D1 + D4*CB)*SB2**4
        M(4,9)=(Sqrt(D7/D2)*SB2**2*SB**3)/D2
        M(5,1)=(Sqrt(D35/D2)*SB**4)/D8;                  M(5,2)=-(R35*CB*SB**3)/D4
        M(5,3)=(Sqrt(D5/D2)*(D5 + D7*C2B)*SB**2)/D8;     M(5,4)=-(R5*(D9*CB + D7*C3B)*SB)/D16
        M(5,5)=(D9+D20*C2B+D35*C4B)/D64;                 M(5,6)=(R5*(D9*CB + D7*C3B)*SB)/D16
        M(5,7)=(Sqrt(D5/D2)*(D5 + D7*C2B)*SB**2)/D8;     M(5,8)=(R35*CB*SB**3)/D4
        M(5,9)=(Sqrt(D35/D2)*SB**4)/D8
        M(6,1)=-(Sqrt(D7/D2)*SB2**2*SB**3)/D2;           M(6,2)=R7*CB2**2*(D1 + D4*CB)*SB2**4
        M(6,3)=-((CB2*(D6+D7*CB+D7*C2B)*SB2**3)/R2);     M(6,4)=((D15+D30*CB+D21*C2B+D14*C3B)*SB2**2)/D8
        M(6,5)=-(R5*(D9*CB + D7*C3B)*SB)/D16;            M(6,6)=(CB2**2*(D3-D6*CB-D21*CB**2+D28*CB**3))/D4
        M(6,7)=-((CB2**3*(-D6+D7*CB-D7*C2B)*SB2)/R2);    M(6,8)=R7*CB2**4*(-D1 + D4*CB)*SB2**2
        M(6,9)=(Sqrt(D7/D2)*(D1/SB2**2)*SB**5)/D8
        M(7,1)=(R7*SB2**4*SB**2)/D2;                     M(7,2)=-(Sqrt(D7/D2)*(D1+D2*CB)*SB2**4*SB)
        M(7,3)=((D9 + D14*CB + D7*C2B)*SB2**4)/D2;       M(7,4)=-((CB2*(D6+D7*CB+D7*C2B)*SB2**3)/R2)
        M(7,5)=(Sqrt(D5/D2)*(D5+D7*C2B)*SB**2)/D8;       M(7,6)=(CB2**3*(-D6+D7*CB-D7*C2B)*SB2)/R2
        M(7,7)=CB2**4*(D1 - D7*CB + D7*CB**2);           M(7,8)=R14*CB2**5*(-D1 + D2*CB)*SB2
        M(7,9)=(R7*(D1/SB2**4)*SB**6)/D32
        M(8,1)=-(R2*SB2**6*SB);                          M(8,2)=(D3 + D4*CB)*SB2**6
        M(8,3)=-(Sqrt(D7/D2)*(D1 + D2*CB)*SB2**4*SB);    M(8,4)=(R7*(D1 + D4*CB)*SB2**2*SB**2)/D4
        M(8,5)=-(R35*CB*SB**3)/D4;                       M(8,6)=R7*CB2**4*(-D1 + D4*CB)*SB2**2
        M(8,7)=-(R14*CB2**5*(-D1 + D2*CB)*SB2);          M(8,8)=CB2**6*(-D3 + D4*CB)
        M(8,9)=D2*R2*CB2**7*SB2
        M(9,1)=SB2**8;                                   M(9,2)=-(R2*SB2**6*SB)
        M(9,3)=(R7*SB2**4*SB**2)/D2;                     M(9,4)=-(Sqrt(D7/D2)*SB2**2*SB**3)/D2
        M(9,5)=(Sqrt(D35/D2)*SB**4)/D8;                  M(9,6)=-(Sqrt(D7/D2)*(D1/SB2**2)*SB**5)/D8
        M(9,7)=(R7*(D1/SB2**4)*SB**6)/D32;               M(9,8)=-D2*R2*CB2**7*SB2
        M(9,9)=CB2**8
      ELSE
        do i=1,9
          mp=-4+i-1
          do j=1,9
            M(i,j)=CZ
          enddo
          if (ABS(B).le.1.0d-6) M(i,i)=DCMPLX(D1,Z)
          if (ABS(B-PI).le.1.0d-6) M(i,9-i+1)=(-D1)**(4-mp)*DCMPLX(D1,Z)
        enddo
      ENDIF

      if (debug) then 
        write(io,'(/,"This is as  Zare (3.54) d4mpm   b=",F8.2)') b/dtr
        write(io,'(8X,9("  |",I2,">  "))') ((-4+j-1),j=1,9)
        do i=1,9
          mp=-4+i-1
          write(io,'("  <",I2,"|  ",9F8.4)') mp,(DBLE(M(i,j)),j=1,9)
        enddo
      endif

      do i=1,9
        do j=1,9
          mp=-4+i-1
          mm=-4+j-1
          M(i,j)=CDExp(-CII*mp*a)*M(i,j)*CDExp(-CII*mm*g)
        enddo
      enddo         
!
      V(1)=DCMPLX(BKQR(9), BKQI(9)); V(2)=DCMPLX( BKQR(8),BKQI(8)); V(3)=DCMPLX(BKQR(7), BKQI(7)); V(4)=DCMPLX( BKQR(6),BKQI(6))
      V(5)=DCMPLX(BKQR(5), Z);       V(6)=DCMPLX(-BKQR(6),BKQI(6)); V(7)=DCMPLX(BKQR(7),-BKQI(7)); V(8)=DCMPLX(-BKQR(8),BKQI(8))
      V(9)=DCMPLX(BKQR(9),-BKQI(9))
      do i=1,9
        VN(i)=CZ
        do k=1,9
!          VN(i) = VN(i) + M(k,i)*v(k)
          VN(i) = VN(i) + DCONJG(M(i,k))*v(k)
        enddo
      enddo        
      BKQR(5)=DBLE(VN(5)); BKQR(6)=DBLE(VN(4));  BKQR(7)=DBLE(VN(3));  BKQR(8)=DBLE(VN(2));  BKQR(9)=DBLE(VN(1))
                           BKQI(6)=DIMAG(VN(4)); BKQI(7)=DIMAG(VN(3)); BKQI(8)=DIMAG(VN(2)); BKQI(9)=DIMAG(VN(1))

      IF (ABS(B).gt.1.0d-6 .and. ABS(B-PI).gt.1.0d-6) then 
        M(1,1)=CB2**12
        M(1,2)=D2*R3*CB2**11*SB2
        M(1,3)=Sqrt(66.d+0)*CB2**10*SB2**2
        M(1,4)=D2*Sqrt(55.d+0)*CB2**9*SB2**3
        M(1,5)=D3*Sqrt(55.d+0)*CB2**8*SB2**4
        M(1,6)=D6*Sqrt(22.d+0)*CB2**7*SB2**5
        M(1,7)=(Sqrt(231.d+0)*SB**6)/D32
        M(1,8)=D6*Sqrt(22.d+0)*CB2**5*SB2**7
        M(1,9)=(D3*Sqrt(55.d+0)*SB2**4*SB**4)/D16
        M(1,10)=(Sqrt(55.d+0)*SB2**6*SB**3)/D4
        M(1,11)=Sqrt(66.d+0)*CB2**2*SB2**10
        M(1,12)=R3*SB2**10*SB
        M(1,13)=SB2**12

        M(2,1)=-D2*R3*CB2**11*SB2
        M(2,2)=CB2**10*(-D5 + D6*CB)
        M(2,3)=Sqrt(22.d+0)*CB2**9*(-D2 + D3*CB)*SB2
        M(2,4)=Sqrt(165.d+0)*CB2**8*(-D1 + D2*CB)*SB2**2
        M(2,5)=Sqrt(165.d+0)*CB2**7*(-D1 + D3*CB)*SB2**3
        M(2,6)=Sqrt(66.d+0)*CB2**6*(-D1 + D6*CB)*SB2**4
        M(2,7)=(D3*Sqrt(77.d+0)*CB*SB**5)/D16
        M(2,8)=Sqrt(66.d+0)*CB2**4*(D1 + D6*CB)*SB2**6
        M(2,9)=Sqrt(165.d+0)*CB2**3*(D1 + D3*CB)*SB2**7
        M(2,10)=Sqrt(165.d+0)*CB2**2*(D1 + D2*CB)*SB2**8
        M(2,11)=Sqrt(22.d+0)*CB2*(D2 + D3*CB)*SB2**9
        M(2,12)=(D5 + D6*CB)*SB2**10
        M(2,13)=R3*SB2**10*SB
      
        M(3,1)=Sqrt(66.d+0)*CB2**10*SB2**2
        M(3,2)=-(Sqrt(22.d+0)*CB2**9*(-D2 + D3*CB)*SB2)
        M(3,3)=-(CB2**8*(-59.d+0 + 88.d+0*CB - 33.d+0*C2B))/D4
        M(3,4)=-(Sqrt(7.5d+0)*CB2**7*(-D15 + 22.d+0*CB - D11*C2B)*SB2)/D2
        M(3,5)=-(Sqrt(7.5d+0)*CB2**6*(-D35 + 44.d+0*CB - 33.d+0*C2B)*SB2**2)/D4
        M(3,6)=-(R3*CB2**5*(-29.d+0 + 22.d+0*CB - 33.d+0*C2B)*SB2**3)/D2
        M(3,7)=(D3*Sqrt(3.5d+0)*(D9 + D11*C2B)*SB**4)/D32
        M(3,8)=(R3*CB2**3*(29.d+0 + 22.d+0*CB + 33.d+0*C2B)*SB2**5)/D2
        M(3,9)=(Sqrt(7.5d+0)*CB2**2*(D35 + 44.d+0*CB + 33.d+0*C2B)*SB2**6)/D4
        M(3,10)=(Sqrt(7.5d+0)*CB2*(D15 + 22.d+0*CB + D11*C2B)*SB2**7)/D2
        M(3,11)=((59.d+0 + 88.d+0*CB + 33.d+0*C2B)*SB2**8)/D4
        M(3,12)=Sqrt(22.d+0)*CB2*(D2 + D3*CB)*SB2**9
        M(3,13)=Sqrt(66.d+0)*CB2**2*SB2**10

        M(4,1)=-D2*Sqrt(55.d+0)*CB2**9*SB2**3
        M(4,2)=Sqrt(165.d+0)*CB2**8*(-D1 + D2*CB)*SB2**2
        M(4,3)=(Sqrt(7.5d+0)*CB2**7*(-D15 + 22.d+0*CB - D11*C2B)*SB2)/D2
        M(4,4)=(CB2**6*(-167.d+0 + 285.d+0*CB - 165.d+0*C2B + 55.d+0*C3B))/D8
        M(4,5)=(D3*CB2**5*(-98.d+0 + 185.d+0*CB - 110.d+0*C2B + 55.d+0*C3B)*SB2)/D16
        M(4,6)=(D3*Sqrt(2.5d+0)*CB2**4*(-D9 + 25.d+0*CB - D11*C2B + D11*C3B)*SB2**2)/D4
        M(4,7)=(Sqrt(105.d+0)*(D21*CB + D11*C3B)*SB**3)/D64
        M(4,8)=(D3*Sqrt(2.5d+0)*CB2**2*(D9 + 25*CB + D11*C2B + D11*C3B)*SB2**4)/D4
        M(4,9)=(D3*CB2*(98.d+0 + 185.d+0*CB + 110.d+0*C2B + 55.d+0*C3B)*SB2**5)/D16
        M(4,10)=((D1 + 60.d+0*CB + 165.d+0*CB**2 + 110.d+0*CB**3)*SB2**6)/D4
        M(4,11)=Sqrt(7.5d+0)*CB2*(D2 + D11*CB + D11*CB**2)*SB2**7
        M(4,12)=Sqrt(165.d+0)*CB2**2*(D1 + D2*CB)*SB2**8
        M(4,13)=(Sqrt(55.d+0)*SB2**6*SB**3)/D4

        M(5,1)=D3*Sqrt(55.d+0)*CB2**8*SB2**4
        M(5,2)=-(Sqrt(165.d+0)*CB2**7*(-D1 + D3*CB)*SB2**3)
        M(5,3)=-(Sqrt(7.5d+0)*CB2**6*(-D35 + 44.d+0*CB - 33.d+0*C2B)*SB2**2)/D4
        M(5,4)=(-3*CB2**5*(-98.d+0 + 185.d+0*CB - 110.d+0*C2B + 55.d+0*C3B)*SB2)/D16
        M(5,5)=-(CB2**4*(-1709.d+0 + 3096.d+0*CB - 2340.d+0*C2B + 1320.d+0*C3B - 495.d+0*C4B))/128.
        M(5,6)=-(Sqrt(2.5d+0)*CB2**3*(-161.d+0 + 252.d+0*CB - 252.d+0*C2B + 132.d+0*C3B - 99.d+0*C4B)*SB2)/D32
        M(5,7)=(Sqrt(105.d+0)*(D35 + 60.d+0*C2B + 33.d+0*C4B)*SB**2)/256.d+0
        M(5,8)=(Sqrt(2.5d+0)*CB2*(161.d+0 + 252.d+0*CB + 252.d+0*C2B + 132.d+0*C3B + 99.d+0*C4B)*SB2**3)/D32
        M(5,9)=((-17.d+0 - 108.d+0*CB + 90.d+0*CB**2 + 660.d+0*CB**3 + 495.d+0*CB**4)*SB2**4)/D16
        M(5,10)=(D3*CB2*(-D3 + D5*CB + 55.d+0*CB**2 + 55.d+0*CB**3)*SB2**5)/D4
        M(5,11)=(Sqrt(7.5d+0)*CB2**2*(D1 + 22.d+0*CB + 33.d+0*CB**2)*SB2**6)/D2
        M(5,12)=Sqrt(165.d+0)*CB2**3*(D1 + D3*CB)*SB2**7
        M(5,13)=(D3*Sqrt(55.d+0)*SB2**4*SB**4)/D16

        M(6,1)=-D6*Sqrt(22.d+0)*CB2**7*SB2**5
        M(6,2)=Sqrt(66.d+0)*CB2**6*(-D1 + D6*CB)*SB2**4
        M(6,3)=(R3*CB2**5*(-29.d+0 + 22.d+0*CB - 33.d+0*C2B)*SB2**3)/D2
        M(6,4)=(D3*Sqrt(2.5d+0)*CB2**4*(-D9 + 25.d+0*CB - D11*C2B + D11*C3B)*SB2**2)/D4
        M(6,5)=(Sqrt(2.5d+0)*CB2**3*(-161.d+0 + 252.d+0*CB - 252.d+0*C2B + 132.d+0*C3B - 99.d+0*C4B)*SB2)/D32
        M(6,6)=(CB2**2*(-175.d+0 + 350.d+0*CB - 300.d+0*C2B + 255.d+0*C3B - 165.d+0*C4B + 99.d+0*Cos(5*B)))/D64
        M(6,7)=(Sqrt(10.5d+0)*(50.d+0*CB + 45.d+0*C3B + 33.d+0*Cos(5*B))*SB)/128.d+0
        M(6,8)=((D5 + 10.d+0*CB - 90.d+0*CB**2 - 120.d+0*CB**3 + 165.d+0*CB**4 + 198.d+0*CB**5)*SB2**2)/D8
        M(6,9)=(Sqrt(2.5d+0)*CB2*(D1 - 18.d+0*CB - 36.d+0*CB**2 + 66.d+0*CB**3 + 99.d+0*CB**4)*SB2**3)/D4
        M(6,10)=(3*Sqrt(2.5d+0)*CB2**2*(D9 + 25.d+0*CB + D11*C2B + D11*C3B)*SB2**4)/D4
        M(6,11)=R3*CB2**3*(-D2 + D11*CB + 33.d+0*CB**2)*SB2**5
        M(6,12)=Sqrt(66.d+0)*CB2**4*(D1 + D6*CB)*SB2**6
        M(6,13)=D6*Sqrt(22.d+0)*CB2**5*SB2**7

        M(7,1)=(Sqrt(231.d+0)*SB**6)/D32
        M(7,2)=(-D3*Sqrt(77.d+0)*CB*SB**5)/D16
        M(7,3)=(D3*Sqrt(3.5d+0)*(D9 + D11*C2B)*SB**4)/D32
        M(7,4)=-(Sqrt(105.d+0)*(D21*CB + D11*C3B)*SB**3)/D64
        M(7,5)=(Sqrt(105.d+0)*(D35 + 60.d+0*C2B + 33.d+0*C4B)*SB**2)/256.d+0
        M(7,6)=-(Sqrt(10.5d+0)*(50.d+0*CB + 45.d+0*C3B + 33.d+0*Cos(5*B))*SB)/128.d+0
        M(7,7)=(50.d+0 + 105.d+0*C2B + 126.d+0*C4B + 231.d+0*Cos(6*B))/512.d+0
        M(7,8)=(Sqrt(10.5)*(50.d+0*CB + 45.d+0*C3B + 33.d+0*Cos(5*B))*SB)/128.d+0
        M(7,9)=(Sqrt(105.d+0)*(D1 - 18.d+0*CB**2 + 33.d+0*CB**4)*SB**2)/D32
        M(7,10)=(Sqrt(105.d+0)*(D21*CB + D11*C3B)*SB**3)/D64
        M(7,11)=(D3*Sqrt(3.5d+0)*(D9 + D11*C2B)*SB**4)/D32
        M(7,12)=(D3*Sqrt(77.d+0)*CB*SB**5)/D16
        M(7,13)=(Sqrt(231.d+0)*SB**6)/D32

        M(8,1)=-D6*Sqrt(22.d+0)*CB2**5*SB2**7
        M(8,2)=Sqrt(66.d+0)*CB2**4*(D1 + D6*CB)*SB2**6
        M(8,3)=-(R3*CB2**3*(29.d+0 + 22.d+0*CB + 33.d+0*C2B)*SB2**5)/D2
        M(8,4)=(D3*Sqrt(2.5d+0)*CB2**2*(D9 + 25.d+0*CB + D11*C2B + D11*C3B)*SB2**4)/D4
        M(8,5)=-(Sqrt(2.5d+0)*CB2*(161.d+0 + 252.d+0*CB + 252.d+0*C2B + 132.d+0*C3B + 99.d+0*C4B)*SB2**3)/D32
        M(8,6)=((175.d+0 + 350.d+0*CB + 300.d+0*C2B + 255.d+0*C3B + 165.d+0*C4B + 99.d+0*Cos(5*B))*SB2**2)/D64
        M(8,7)=-(Sqrt(10.5d+0)*(50.d+0*CB + 45.d+0*C3B + 33.d+0*Cos(5*B))*SB)/128.d+0
        M(8,8)=(CB2**2*(-D5 + 10.d+0*CB + 90.d+0*CB**2 - 120.d+0*CB**3 - 165.d+0*CB**4 + 198.d+0*CB**5))/D8
        M(8,9)=(Sqrt(2.5d+0)*CB2**3*(D1 + 18.d+0*CB - 36.d+0*CB**2 - 66.d+0*CB**3 + 99.d+0*CB**4)*SB2)/D4
        M(8,10)=(D3*Sqrt(2.5d+0)*CB2**4*(D1 - D4*CB - D11*CB**2 + 22.d+0*CB**3)*SB2**2)/D2
        M(8,11)=R3*CB2**5*(-D2 - D11*CB + 33.d+0*CB**2)*SB2**3
        M(8,12)=Sqrt(66.d+0)*CB2**6*(-D1 + D6*CB)*SB2**4
        M(8,13)=D6*Sqrt(22.d+0)*CB2**7*SB2**5

        M(9,1)=(D3*Sqrt(55.d+0)*SB2**4*SB**4)/D16
        M(9,2)=-(Sqrt(165.d+0)*CB2**3*(D1 + D3*CB)*SB2**7)
        M(9,3)=(Sqrt(7.5d+0)*CB2**2*(D35 + 44*CB + 33*C2B)*SB2**6)/D4
        M(9,4)=(-D3*CB2*(98.d+0 + 185.d+0*CB + 110.d+0*C2B + 55.d+0*C3B)*SB2**5)/D16
        M(9,5)=((1709.d+0 + 3096.d+0*CB + 2340.d+0*C2B + 1320.d+0*C3B + 495.d+0*C4B)*SB2**4)/128.d+0
        M(9,6)=-(Sqrt(2.5)*CB2*(161.d+0 + 252.d+0*CB + 252.d+0*C2B + 132.d+0*C3B + 99.d+0*C4B)*SB2**3)/D32
        M(9,7)=(Sqrt(105.d+0)*(D35 + 60.d+0*C2B + 33.d+0*C4B)*SB**2)/256.d+0
        M(9,8)=(Sqrt(2.5d+0)*CB2**3*(-161.d+0 + 252.d+0*CB - 252.d+0*C2B + 132.d+0*C3B - 99.d+0*C4B)*SB2)/D32
        M(9,9)=(CB2**4*(-17.d+0 + 108.d+0*CB + 90.d+0*CB**2 - 660.d+0*CB**3 + 495.d+0*CB**4))/D16
        M(9,10)=(D3*CB2**5*(D3 + D5*CB - 55.d+0*CB**2 + 55.d+0*CB**3)*SB2)/D4
        M(9,11)=(Sqrt(7.5d+0)*CB2**6*(D1 - 22.d+0*CB + 33.d+0*CB**2)*SB2**2)/D2
        M(9,12)=Sqrt(165.d+0)*CB2**7*(-D1 + D3*CB)*SB2**3
        M(9,13)=D3*Sqrt(55.d+0)*CB2**8*SB2**4

        M(10,1)=-(Sqrt(55.d+0)*SB2**6*SB**3)/D4
        M(10,2)=(Sqrt(165.d+0)*SB2**6*SB*(SB + Sin(D2*B)))/D4
        M(10,3)=-(Sqrt(7.5d+0)*CB2*(D15 + 22.d+0*CB + D11*C2B)*SB2**7)/D2
        M(10,4)=((167.d+0 + 285.d+0*CB + 165.d+0*C2B + 55.d+0*C3B)*SB2**6)/D8
        M(10,5)=(-D3*CB2*(98.d+0 + 185.d+0*CB + 110.d+0*C2B + 55.d+0*C3B)*SB2**5)/D16
        M(10,6)=(D3*Sqrt(2.5d+0)*CB2**2*(D9 + 25.d+0*CB + D11*C2B + D11*C3B)*SB2**4)/D4
        M(10,7)=-(Sqrt(105.d+0)*(D21*CB + D11*C3B)*SB**3)/D64
        M(10,8)=(D3*Sqrt(2.5d+0)*CB2**4*(-D9 + 25.d+0*CB - D11*C2B + D11*C3B)*SB2**2)/D4
        M(10,9)=(-D3*CB2**5*(-98.d+0 + 185.d+0*CB - 110.d+0*C2B + 55.d+0*C3B)*SB2)/D16
        M(10,10)=(CB2**6*(-D1 + 60.d+0*CB - 165.d+0*CB**2 + 110.d+0*CB**3))/D4
        M(10,11)=Sqrt(7.5d+0)*CB2**7*(D2 - D11*CB + D11*CB**2)*SB2
        M(10,12)=Sqrt(16.d+05)*CB2**8*(-D1 + D2*CB)*SB2**2
        M(10,13)=D2*Sqrt(55.d+0)*CB2**9*SB2**3

        M(11,1)=Sqrt(66.d+0)*CB2**2*SB2**10
        M(11,2)=-(Sqrt(22.d+0)*CB2*(D2 + D3*CB)*SB2**9)
        M(11,3)=((59.d+0 + 88.d+0*CB + 33.d+0*C2B)*SB2**8)/D4
        M(11,4)=-(Sqrt(7.5d+0)*CB2*(D15 + 22.d+0*CB + D11*C2B)*SB2**7)/D2
        M(11,5)=(Sqrt(7.5d+0)*CB2**2*(D35 + 44.d+0*CB + 33.d+0*C2B)*SB2**6)/D4
        M(11,6)=-(R3*CB2**3*(29.d+0 + 22.d+0*CB + 33.d+0*C2B)*SB2**5)/D2
        M(11,7)=(D3*Sqrt(3.5d+0)*(D9 + D11*C2B)*SB**4)/D32
        M(11,8)=(R3*CB2**5*(-29.d+0 + 22.d+0*CB - 33.d+0*C2B)*SB2**3)/D2
        M(11,9)=-(Sqrt(7.5d+0)*CB2**6*(-D35 + 44.d+0*CB - 33.d+0*C2B)*SB2**2)/D4
        M(11,10)=(Sqrt(7.5d+0)*CB2**7*(-D15 + 22.d+0*CB - D11*C2B)*SB2)/D2
        M(11,11)=-(CB2**8*(-59.d+0 + 88.d+0*CB - 33.d+0*C2B))/D4
        M(11,12)=Sqrt(22.d+0)*CB2**9*(-D2 + D3*CB)*SB2
        M(11,13)=Sqrt(66.d+0)*CB2**10*SB2**2

        M(12,1)=-(R3*SB2**10*SB)
        M(12,2)=(D5 + D6*CB)*SB2**10
        M(12,3)=-(Sqrt(22.d+0)*CB2*(D2 + D3*CB)*SB2**9)
        M(12,4)=(Sqrt(165.d+0)*SB2**6*SB*(SB + Sin(D2*B)))/D4
        M(12,5)=-(Sqrt(165.d+0)*CB2**3*(D1 + D3*CB)*SB2**7)
        M(12,6)=Sqrt(66.d+0)*CB2**4*(D1 + D6*CB)*SB2**6
        M(12,7)=(-D3*Sqrt(77.d+0)*CB*SB**5)/D16
        M(12,8)=Sqrt(66.d+0)*CB2**6*(-D1 + D6*CB)*SB2**4
        M(12,9)=-(Sqrt(165.d+0)*CB2**7*(-D1 + D3*CB)*SB2**3)
        M(12,10)=Sqrt(165.d+0)*CB2**8*(-D1 + D2*CB)*SB2**2
        M(12,11)=-(Sqrt(22.d+0)*CB2**9*(-D2 + D3*CB)*SB2)
        M(12,12)=CB2**10*(-D5 + D6*CB)
        M(12,13)=D2*R3*CB2**11*SB2

        M(13,1)=SB2**12
        M(13,2)=-(R3*SB2**10*SB)
        M(13,3)=Sqrt(66.d+0)*CB2**2*SB2**10
        M(13,4)=-(Sqrt(55.d+0)*SB2**6*SB**3)/D4
        M(13,5)=(D3*Sqrt(55.d+0)*SB2**4*SB**4)/D16
        M(13,6)=-D6*Sqrt(22.d+0)*CB2**5*SB2**7
        M(13,7)=(Sqrt(231.d+0)*SB**6)/D32
        M(13,8)=-D6*Sqrt(22.d+0)*CB2**7*SB2**5
        M(13,9)=D3*Sqrt(55.d+0)*CB2**8*SB2**4
        M(13,10)=-D2*Sqrt(55.d+0)*CB2**9*SB2**3
        M(13,11)=Sqrt(66.d+0)*CB2**10*SB2**2
        M(13,12)=-D2*R3*CB2**11*SB2
        M(13,13)=CB2**12
      ELSE
        do i=1,13
          mp=-6+i-1
          do j=1,13
            M(i,j)=CZ
          enddo
          if (ABS(B).le.1.0d-6) M(i,i)=DCMPLX(D1,Z)
          if (ABS(B-PI).le.1.0d-6) M(i,13-i+1)=(-D1)**(6-mp)*DCMPLX(D1,Z)
        enddo
      ENDIF

      if (debug) then 
        write(io,'(/,"This is as  Zare (3.54) d6mpm   b=",F8.2)') b/dtr
        write(io,'(8X,13("  |",I2,">  "))') ((-6+j-1),j=1,13)
        do i=1,13
          mp=-6+i-1
          write(io,'("  <",I2,"|  ",13F8.4)') mp,(DBLE(M(i,j)),j=1,13)
        enddo
      endif

      do i=1,13
        do j=1,13
          mp=-6+i-1
          mm=-6+j-1
          M(i,j)=CDExp(-CII*mp*a)*M(i,j)*CDExp(-CII*mm*g)
        enddo
      enddo         
!
       V(1)=DCMPLX( BKQR(16), BKQI(16));  V(2)=DCMPLX( BKQR(15), BKQI(15));  V(3)=DCMPLX( BKQR(14), BKQI(14)); 
       V(4)=DCMPLX( BKQR(13), BKQI(13));  V(5)=DCMPLX( BKQR(12), BKQI(12));  V(6)=DCMPLX( BKQR(11), BKQI(11));
       V(7)=DCMPLX( BKQR(10), Z);         V(8)=DCMPLX(-BKQR(11), BKQI(11));  V(9)=DCMPLX( BKQR(12),-BKQI(12)); 
      V(10)=DCMPLX(-BKQR(13), BKQI(13)); V(11)=DCMPLX( BKQR(14),-BKQI(14)); V(12)=DCMPLX(-BKQR(15), BKQI(15));
      V(13)=DCMPLX( BKQR(16),-BKQI(16))

      do i=1,13
        VN(i)=CZ
        do k=1,13
!          VN(i) = VN(i) + M(k,i)*v(k)
          VN(i) = VN(i) + DCONJG(M(i,k))*v(k)
        enddo
      enddo        
      BKQR(10)=DBLE(VN(7)); BKQR(11)= DBLE(VN(6));  BKQR(12)= DBLE(VN(5));  BKQR(13)= DBLE(VN(4)) 
                            BKQR(14)= DBLE(VN(3));  BKQR(15)= DBLE(VN(2));  BKQR(16)= DBLE(VN(1))
                            BKQI(11)=DIMAG(VN(6));  BKQI(12)=DIMAG(VN(5));  BKQI(13)=DIMAG(VN(4)) 
                            BKQI(14)=DIMAG(VN(3));  BKQI(15)=DIMAG(VN(2));  BKQI(16)=DIMAG(VN(1))
!      
      if (OUTP(5).or.debug) then 
        WRITE(IO,'(" After rotation, Crystal field parameters ")')
        do k=2,2*Lvalue,2
          do q=0,k
            i=Bkq_index(k,q)
            if (q.eq.0) then 
              WRITE(IO,'(7X,A4,F14.6)')    BkqR_label(i),BkqR(i)
            else if (q.ne.0) then
              WRITE(IO,'(2(7X,A4,F14.6))') BkqR_label(i),BkqR(i),BkqI_label(i),BkqI(i) 
            endif
          enddo
        enddo
      endif
!
      RETURN
      END SUBROUTINE RotateLF
!
!-----------------------------------------------------------------------
!
      SUBROUTINE A_B_Mat()
!
!  Reads the ALtp or BLli intensity parameters and converts them to BLli or ALtp.
!  l=2,4,6; l=0,L; i=x,y,z; ALtp & BLli complex 
!      
      IMPLICIT none
      integer i,j,k, lamda,t,ii
      COMPLEX*16 AV2(9),AV4(15),AV6(21),BV2(9),BV4(15),BV6(21)
      real*8 B2(9,9),B4(15,15),B6(21,21), A2(9,9),A4(15,15),A6(21,21)  
      real*8 A2_1(4,5), A2_2(5,4), A4_1(7,8), A4_2(8,7), A6_1(10,11), A6_2(11,10) 
      real*8 B2_1(5,4), B2_2(4,5), B4_1(8,7), B4_2(7,8), B6_1(11,10), B6_2(10,11) 
      ReAL*8, PARAMETER::  Z=0.0d+00,  D1= 1.0d+00, D2= 2.0d+00, D3= 3.0d+00, D4= 4.0d+00, D5= 5.0d+00
      real*8, Parameter:: D6=6.0d+00,  D7= 7.0d+00, D8= 8.0d+00, D9= 9.0d+00,D10=10.0d+00,D11=11.0d+00
      real*8, Parameter::D12=12.0d+00,D13=13.0d+00,D14=14.0d+00,D15=15.0d+00,D16=16.0d+00
      real*8, Parameter::D18=18.0d+00,D20=20.0d+00,D21=21.0d+00,D24=24.0d+00,D25=25.0d+00,D26=26.0d+00
      real*8, Parameter::D28=28.0d+00,D30=30.0d+00,D33=33.0d+00,D35=35.0d+00
      real*8, Parameter::D36=36.0d+00,D39=39.0d+00,D40=40.0d+00,D42=42.0d+00,D45=45.0d+00,D48=48.0d+00
      real*8, Parameter::D52=52.0d+00,D55=55.0d+00,D72=72.0d+00,D78=78.0d+00,D84=84.0d+00,D90=90.0d+00
      real*8, Parameter::D91=91.0d+00,D182=182.0d+00,D156=156.0d+00
      COMPLEX*16, Parameter:: CZ=(Z,Z), C1=(D1,Z), CI=(Z,D1)
      integer d21r(5),d21c(4),d22r(4),d22c(5), d41r(8),d41c(7),d42r(7),d42c(8), d61r(11),d61c(10),d62r(10),d62c(11)
      data d21r/1,2,6,7,8/, d21c/2,4,7,9/, d22r/3,4,5,9/, d22c/1,3,5,6,8/
      data d41r/1,2,6,7,8,12,13,14/, d41c/2,4,6,8,11,13,15/, d42r/3,4,5,9,10,11,15/, d42c/1,3,5,7,9,10,12,14/
      data d61r/1,2,6,7,8,12,13,14,18,19,20/, d61c/2,4,6,8,10,12,15,17,19,21/ 
      data d62r/3,4,5,9,10,11,15,16,17,21/,   d62c/1,3,5,7,9,11,13,14,16,18,20/
      
      if (IntN1.eq.1) then  ! Altp convert to Blki
      
        B2_1=z; B2_2=z; B2=0
        B2_1(1,1) =-sqrt(d1/d5);  B2_1(1,2) = d1; B2_1(1,3) =-sqrt(d4/d5)
        B2_1(2,1) =-sqrt(d1/d5);  B2_1(2,2) = d1; B2_1(2,3) =-sqrt(d4/d5)
        B2_1(3,1) =-sqrt(d3/d10); B2_1(3,2) = sqrt(d1/d6); B2_1(3,3) = sqrt(d8/d15)
        B2_1(4,1) = sqrt(d3/d10); B2_1(4,2) = sqrt(d1/d6); B2_1(4,3) = sqrt(d1/d30); B2_1(4,4) =-sqrt(d1/d2)
        B2_1(5,1) =-sqrt(d3/d10); B2_1(5,2) =-sqrt(d1/d6); B2_1(5,3) =-sqrt(d1/d30); B2_1(5,4) =-sqrt(d1/d2)
      
        B2_2(1,1) =-sqrt(d2/d5);  B2_2(1,4) = sqrt(d3/d5)
        B2_2(2,1) = sqrt(d3/d20); B2_2(2,2) = sqrt(d1/d4); B2_2(2,3) = sqrt(d1/d6); B2_2(2,4) = sqrt(d1/d10) 
        B2_2(2,5) =-sqrt(d1/d3)  
        B2_2(3,1) =-sqrt(d3/d20); B2_2(3,2) =-sqrt(d1/d4); B2_2(3,3) = sqrt(d1/d6); B2_2(3,4) =-sqrt(d1/d10)
        B2_2(3,5) =-sqrt(d1/d3)  
        B2_2(4,3) = sqrt(d2/d3);  B2_2(4,5) = sqrt(d1/d3)  
        do i=1,5; do j=1,4; B2(d21r(i),d21c(j))=B2_1(i,j); enddo; enddo
        do i=1,4; do j=1,5; B2(d22r(i),d22c(j))=B2_2(i,j); enddo; enddo
!
        B4_1=z; B4_2=z; B4=0 
        B4_1(1,1) =-sqrt(d1/d3);  B4_1(1,3) = d1;           B4_1(1,5) =-sqrt(d2/d3)
        B4_1(2,1) =-sqrt(d1/d3);  B4_1(2,3) = d1;           B4_1(2,5) =-sqrt(d2/d3)
        B4_1(3,1) =-sqrt(d5/d12); B4_1(3,3) = sqrt(d1/d20); B4_1(3,5) = sqrt(d8/d15)
        B4_1(4,1) = sqrt(d5/d24); B4_1(4,2) =-sqrt(d1/d72); B4_1(4,3) = sqrt(d9/d40); B4_1(4,4) = sqrt(d7/d40)
        B4_1(4,5) = sqrt(d1/d15); B4_1(4,6) =-sqrt(d14/d45)
        B4_1(5,1) =-sqrt(d5/d24); B4_1(5,2) =-sqrt(d1/d72); B4_1(5,3) =-sqrt(d9/d40); B4_1(5,4) = sqrt(d7/d40)
        B4_1(5,5) =-sqrt(d1/d15); B4_1(5,6) =-sqrt(d14/d45)
        B4_1(6,2) =-sqrt(d7/d36); B4_1(6,4) = sqrt(d9/d20); B4_1(6,6) = sqrt(d16/d45)
        B4_1(7,2) = sqrt(d7/d18); B4_1(7,4) = sqrt(d1/d10); B4_1(7,6) = sqrt(d1/d90); B4_1(7,7) =-sqrt(d1/d2)
        B4_1(8,2) =-sqrt(d7/d18); B4_1(8,4) =-sqrt(d1/d10); B4_1(8,6) =-sqrt(d1/d90); B4_1(8,7) =-sqrt(d1/d2)
        
        B4_2(1,1) =-sqrt(d4/d9);  B4_2(1,6) = sqrt(d5/d9)
        B4_2(2,1) = sqrt(d5/d36); B4_2(2,2) =-sqrt(d1/d24); B4_2(2,3) = sqrt(d1/d4);  B4_2(2,4) = sqrt(d9/d40)
        B4_2(2,6) = sqrt(d1/d9);  B4_2(2,7) =-sqrt(d7/d30)
        B4_2(3,1) =-sqrt(d5/d36); B4_2(3,2) =-sqrt(d1/d24); B4_2(3,3) =-sqrt(d1/d4);  B4_2(3,4) = sqrt(d9/d40)
        B4_2(3,6) =-sqrt(d1/d9);  B4_2(3,7) =-sqrt(d7/d30)
        B4_2(4,2) =-sqrt(d1/d3);  B4_2(4,4) = sqrt(d1/d5);  B4_2(4,7) = sqrt(d7/d15)
        B4_2(5,2) = sqrt(d7/d24); B4_2(5,4) = sqrt(d7/d40); B4_2(5,5) = sqrt(d1/d10); B4_2(5,7) = sqrt(d1/d30)
        B4_2(5,8) =-sqrt(d2/d5)
        B4_2(6,2) =-sqrt(d7/d24); B4_2(6,4) =-sqrt(d7/d40); B4_2(6,5) = sqrt(d1/d10); B4_2(6,7) =-sqrt(d1/d30)
        B4_2(6,8) =-sqrt(d2/d5)
        B4_2(7,5) = sqrt(d4/d5);  B4_2(7,8) = sqrt(d1/d5); 
        do i=1,8; do j=1,7; B4(d41r(i),d41c(j))=B4_1(i,j); enddo; enddo
        do i=1,7; do j=1,8; B4(d42r(i),d42c(j))=B4_2(i,j); enddo; enddo
     
        B6_1=z; B6_2=z; B6=0
        B6_1(1,1)=-sqrt(d5/d13);   B6_1(1,4)= d1;            B6_1(1,7)=-sqrt(d8/d13)
        B6_1(2,1)=-sqrt(d5/d13);   B6_1(2,4)= d1;            B6_1(2,7)=-sqrt(d8/d13)
        B6_1(3,1)=-sqrt(d35/d78);  B6_1(3,4)= sqrt(d1/d42);  B6_1(3,7)= sqrt(d48/d91)
        B6_1(4,1)= sqrt(d7/d39);   B6_1(4,2)=-sqrt(d1/d26);  B6_1(4,4)= sqrt(d5/d21); B6_1(4,5)= sqrt(d3/d14)
        B6_1(4,7)= sqrt(d15/d182); B6_1(4,8)=-sqrt(d45/d182)
        B6_1(5,1)=-sqrt(d7/d39);   B6_1(5,2)=-sqrt(d1/d26);  B6_1(5,4)=-sqrt(d5/d21); B6_1(5,5)= sqrt(d3/d14)
        B6_1(5,7)=-sqrt(d15/d182); B6_1(5,8)=-sqrt(d45/d182)
        B6_1(6,2)=-sqrt(d9/d26);   B6_1(6,5)= sqrt(d3/d14);  B6_1(6,8)= sqrt(d40/d91)
        B6_1(7,2)= sqrt(d15/d52);  B6_1(7,3)=-sqrt(d1/d156); B6_1(7,5)= sqrt(d5/d28); B6_1(7,6)= sqrt(d11/d84)
        B6_1(7,8)= sqrt(d3/d91);   B6_1(7,9)=-sqrt(d33/d91)
        B6_1(8,2)=-sqrt(d15/d52);  B6_1(8,3)=-sqrt(d1/d156); B6_1(8,5)=-sqrt(d5/d28); B6_1(8,6)= sqrt(d11/d84)
        B6_1(8,8)=-sqrt(d3/d91);   B6_1(8,9)=-sqrt(d33/d91)
        B6_1(9,3)=-sqrt(d11/d78);  B6_1(9,6)= sqrt(d25/d42); B6_1(9,9)= sqrt(d24/d91)
        B6_1(10,3)= sqrt(d11/d26); B6_1(10,6)= sqrt(d1/d14); B6_1(10,9)= sqrt(d1/d182); B6_1(10,10)=-sqrt(d1/d2)
        B6_1(11,3)=-sqrt(d11/d26); B6_1(11,6)=-sqrt(d1/d14); B6_1(11,9)=-sqrt(d1/d182); B6_1(11,10)=-sqrt(d1/d2)
        
        B6_2(1,1)=-sqrt(d6/d13);   B6_2(1,8)= sqrt(d7/d13)
        B6_2(2,1)= sqrt(d7/d52);   B6_2(2,2)=-sqrt(d5/d78);  B6_2(2,4)= d1/d2; B6_2(2,5)= sqrt(d5/d21);
        B6_2(2,8)= sqrt(d3/d26);   B6_2(2,9)=-sqrt(d18/d91)
        B6_2(3,1)=-sqrt(d7/d52);   B6_2(3,2)=-sqrt(d5/d78);  B6_2(3,4)=-d1/d2; B6_2(3,5)= sqrt(d5/d21); 
        B6_2(3,8)=-sqrt(d3/d26);   B6_2(3,9)=-sqrt(d18/d91)
        B6_2(4,2)=-sqrt(d16/d39);  B6_2(4,5)= sqrt(d2/d21);  B6_2(4,9)= sqrt(d45/d91)
        B6_2(5,2)= sqrt(d3/d13);   B6_2(5,3)=-sqrt(d1/d52);  B6_2(5,5)= sqrt(d3/d14);  B6_2(5,6)= sqrt(d5/d28)
        B6_2(5,9)= sqrt(d5/d91);   B6_2(5,10)=-sqrt(d55/d182)
        B6_2(6,2)=-sqrt(d3/d13);   B6_2(6,3)=-sqrt(d1/d52);  B6_2(6,5)=-sqrt(d3/d14);  B6_2(6,6)= sqrt(d5/d28)
        B6_2(6,9)=-sqrt(d5/d91);   B6_2(6,10)=-sqrt(d55/d182)
        B6_2(7,3)=-sqrt(d10/d39);  B6_2(7,6)= sqrt(d8/d21);  B6_2(7,10)= sqrt(d33/d91)
        B6_2(8,3)= sqrt(d55/d156); B6_2(8,6)= sqrt(d11/d84); B6_2(8,7)= sqrt(d1/d14);  B6_2(8,10)= sqrt(d3/d182)
        B6_2(8,11)=-sqrt(d3/d7)
        B6_2(9,3)=-sqrt(d55/d156); B6_2(9,6)=-sqrt(d11/d84); B6_2(9,7)= sqrt(d1/d14);  B6_2(9,10)=-sqrt(d3/d182)
        B6_2(9,11)=-sqrt(d3/d7)
        B6_2(10,7)= sqrt(d6/d7);   B6_2(10,11)= sqrt(d1/d7)
        do i=1,11; do j=1,10; B6(d61r(i),d61c(j))=B6_1(i,j); enddo; enddo
        do i=1,10; do j=1,11; B6(d62r(i),d62c(j))=B6_2(i,j); enddo; enddo
!  
        do i=1,3
          lamda=2*i
          ii=0
          do j=1,3
            t=lamda-2+j  ! t=lamda-1,lamda,lamda+1 
            do k=1,t+1
              ii=ii+1
              if (i.eq.1) AV2(ii)=DCMPLX(AltpR(i,j,k),AltpI(i,j,k))
              if (i.eq.2) AV4(ii)=DCMPLX(AltpR(i,j,k),AltpI(i,j,k))
              if (i.eq.3) AV6(ii)=DCMPLX(AltpR(i,j,k),AltpI(i,j,k))
            enddo
          enddo
        enddo

        if (check2(4)) then
          write(Idebug,'("A2tp p=",5X,"0",14X,"+1",13X,"+2",13X,"+3")')
          write(Idebug,'("t=1   ",2(2F6.1,3x))')(AV2(I),I=1,2);
          write(Idebug,'("t=2   ",3(2F6.1,3x))')(AV2(I),I=3,5)
          write(Idebug,'("t=3   ",4(2F6.1,3x))')(AV2(I),I=6,9)
          write(Idebug,'("A4tp p=",5X,"0",14X,"+1",13X,"+2",13X,"+3",13X,"+4",13X,"+5")')
          write(Idebug,'("t=3   ",4(2F6.1,3x))')(AV4(I),I=1,4);
          write(Idebug,'("t=4   ",5(2F6.1,3x))')(AV4(I),I=5,9)
          write(Idebug,'("t=5   ",6(2F6.1,3x))')(AV4(I),I=10,15)
          write(Idebug,'("A6tp p=",5X,"0",14X,"+1",13X,"+2",13X,"+3",13X,"+4",13X,"+5",13X,"+6",13X,"+7")')
          write(Idebug,'("t=5   ",6(2F6.1,3x))')(AV6(I),I=1,6)
          write(Idebug,'("t=6   ",7(2F6.1,3x))')(AV6(I),I=7,13)
          write(Idebug,'("t=7   ",8(2F6.1,3x))')(AV6(I),I=14,21)
          write(Idebug,'(/,"B2")'); Do i=1,9;  write(Idebug,'(9F6.3, 5X,2F6.1)')(B2(i,j),j=1,9), AV2(i);enddo   
          write(Idebug,'(/,"B4")'); Do i=1,15; write(Idebug,'(15F6.3,5X,2F6.1)')(B4(i,j),j=1,15),AV4(i);enddo   
          write(Idebug,'(/,"B6")'); Do i=1,21; write(Idebug,'(21F6.3,5X,2F6.1)')(B6(i,j),j=1,21),AV6(i);enddo   
        endif

        BV2=CZ
        do j=1,9;  BV2(1)=BV2(1)+DCMPLX(B2(1,j))*DREAL(AV2(j)); enddo  ! 1st 2 rows must be treated explicitly
        do j=1,9;  BV2(2)=BV2(2)+DCMPLX(B2(2,j))*DIMAG(AV2(j))*CI; enddo
        Do i=3,9;  do j=1,9;  BV2(i)=BV2(i)+DCMPLX(B2(i,j))*AV2(j); enddo; enddo   
        do i=2,8,3; BV2(i)=BV2(i)*CI; enddo  ! all Blky are multiplied by i
        BV4=CZ 
        do j=1,15; BV4(1)=BV4(1)+DCMPLX(B4(1,j))*DREAL(AV4(j)); enddo
        do j=1,15; BV4(2)=BV4(2)+DCMPLX(B4(2,j))*DIMAG(AV4(j))*CI; enddo
        Do i=3,15; do j=1,15; BV4(i)=BV4(i)+DCMPLX(B4(i,j))*AV4(j); enddo; enddo   
        do i=2,14,3; BV4(i)=BV4(i)*CI; enddo  
        BV6=CZ 
        do j=1,21; BV6(1)=BV6(1)+DCMPLX(B6(1,j))*DREAL(AV6(j)); enddo 
        do j=1,21; BV6(2)=BV6(2)+DCMPLX(B6(2,j))*DIMAG(AV6(j))*CI; enddo
        Do i=3,21; do j=1,21; BV6(i)=BV6(i)+DCMPLX(B6(i,j))*AV6(j); enddo; enddo   
        do i=2,20,3; BV6(i)=BV6(i)*CI; enddo  

        ii=0; do i=1,3; do j=1,3; ii=ii+1; BlkiR(i,j)=DREAL(BV2(ii));   BlkiI(i,j)=DIMAG(BV2(ii));   endDO; endDO
        ii=0; do i=1,5; do j=1,3; ii=ii+1; BlkiR(3+i,j)=DREAL(BV4(ii)); BlkiI(3+i,j)=DIMAG(BV4(ii)); enddo; enddo
        ii=0; do i=1,7; do j=1,3; ii=ii+1; BlkiR(8+i,j)=DREAL(BV6(ii)); BlkiI(8+i,j)=DIMAG(BV6(ii)); enddo; enddo

        if (check2(4)) then
          write(Idebug,'(/,12x,"x",15X,"y",15X,"z")')
          ii=0
          do lamda=2,6,2
            do i= 0,lamda
              ii=ii+1
              write(Idebug,'("B",2i1,"i=",3(2F7.2,3x))') lamda,i,(BlkiR(ii,j),BlkiI(ii,j),j=1,3)
            enddo
          enddo
        endif ! check2(4)
        
      endif ! IntN1.eq.1
      
      if (IntN1.eq.2) then  ! Blki convert to Altp 
      
        ii=0; do i=1,3; do j=1,3; ii=ii+1; BV2(ii)=DCMPLX(BlkiR(i,j),BlkiI(i,j)); enddo; enddo
        ii=0; do i=1,5; do j=1,3; ii=ii+1; BV4(ii)=DCMPLX(BlkiR(3+i,j),BlkiI(3+i,j)); enddo; enddo
        ii=0; do i=1,7; do j=1,3; ii=ii+1; BV6(ii)=DCMPLX(BlkiR(8+i,j),BlkiI(8+i,j)); enddo; enddo
        
        if (check2(4)) then
          write(Idebug,'(/,12x,"x",15X,"y",15X,"z")')
          ii=0
          do lamda=2,6,2
            do i= 0,lamda
              ii=ii+1
              write(Idebug,'("B",2i1,"i=",3(2F7.2,3x))') lamda,i,(BlkiR(ii,j),BlkiI(ii,j),j=1,3)
            enddo
          enddo
        endif  ! check2(4)       

        A2_1=z; A2_2=z; A2=0
        A2_1(1,1)=-sqrt(d1/d20);A2_1(1,2)=-sqrt(d1/d20);A2_1(1,3)=-sqrt(d3/d10);A2_1(1,4)=sqrt(d3/d10);A2_1(1,5)=-sqrt(d3/d10)
        A2_1(2,1)= sqrt(d1/d4); A2_1(2,2)= sqrt(d1/d4); A2_1(2,3)= sqrt(d1/d6); A2_1(2,4)=sqrt(d1/d6); A2_1(2,5)=-sqrt(d1/d6)
        A2_1(3,1)=-sqrt(d1/d5); A2_1(3,2)=-sqrt(d1/d5); A2_1(3,3)= sqrt(d8/d15);A2_1(3,4)=sqrt(d1/d30);A2_1(3,5)=-sqrt(d1/d30)
        A2_1(4,4)=-sqrt(d1/d2); A2_1(4,5)=-sqrt(d1/d2)
        A2_2(1,1) =-sqrt(d2/d5); A2_2(1,2) = sqrt(d3/d5);  A2_2(1,3) =-sqrt(d3/d5)
        A2_2(2,2) = d1; A2_2(2,3) =-d1
        A2_2(3,2) = sqrt(d1/d6); A2_2(3,3) = sqrt(d1/d6); A2_2(3,4) = sqrt(d2/d3)
        A2_2(4,1) = sqrt(d3/d5); A2_2(4,2) = sqrt(d2/d5); A2_2(4,3) =-sqrt(d2/d5)  
        A2_2(5,2) =-sqrt(d1/d3); A2_2(5,3) =-sqrt(d1/d3); A2_2(5,4) = sqrt(d1/d3)  
        do i=1,4; do j=1,5; A2(d21c(i),d21r(j))=A2_1(i,j); enddo; enddo
        do i=1,5; do j=1,4; A2(d22c(i),d22r(j))=A2_2(i,j); enddo; enddo
!
        A4_1=z; A4_2=z; A4=0 
        A4_1(1,1)=-sqrt(d1/d12); A4_1(1,2)=-sqrt(d1/d12); A4_1(1,3)=-sqrt(d5/d12); A4_1(1,4)=sqrt(d5/d24);A4_1(1,5)=-sqrt(d5/d24)
        A4_1(2,4)=-sqrt(d1/d72); A4_1(2,5)=-sqrt(d1/d72); A4_1(2,6)=-sqrt(d7/d36); A4_1(2,7)=sqrt(d7/d18);A4_1(2,8)=-sqrt(d7/d18) 
        A4_1(3,1)= sqrt(d1/d4);  A4_1(3,2)=sqrt(d1/d4);  A4_1(3,3)=sqrt(d1/d20); A4_1(3,4)=sqrt(d9/d40); A4_1(3,5)=-sqrt(d9/d40)
        A4_1(4,4)= sqrt(d7/d40); A4_1(4,5)=sqrt(d7/d40); A4_1(4,6)=sqrt(d9/d20); A4_1(4,7)=sqrt(d1/d10); A4_1(4,8)=-sqrt(d1/d10)
        A4_1(5,1)=-sqrt(d1/d6); A4_1(5,2)=-sqrt(d1/d6); A4_1(5,3)=sqrt(d8/d15); A4_1(5,4)=sqrt(d1/d15); A4_1(5,5)=-sqrt(d1/d15)
        A4_1(6,4)=-sqrt(d14/d45);A4_1(6,5)=-sqrt(d14/d45);A4_1(6,6)=sqrt(d16/d45);A4_1(6,7)=sqrt(d1/d90);A4_1(6,8)=-sqrt(d1/d90)
        A4_1(7,7)=-sqrt(d1/d2);  A4_1(7,8)=-sqrt(d1/d2)
        A4_2(1,1)=-sqrt(d4/d9);  A4_2(1,2)= sqrt(d5/d9);  A4_2(1,3)=-sqrt(d5/d9)
        A4_2(2,2)=-sqrt(d1/d24); A4_2(2,3)=-sqrt(d1/d24); A4_2(2,4)=-sqrt(d1/d3); A4_2(2,5)= sqrt(d7/d24); A4_2(2,6)=-sqrt(d7/d24)
        A4_2(3,2)= d1;           A4_2(3,3)=-d1
        A4_2(4,2)= sqrt(d9/d40); A4_2(4,3)= sqrt(d9/d40); A4_2(4,4)= sqrt(d1/d5); A4_2(4,5)= sqrt(d7/d40); A4_2(4,6)=-sqrt(d7/d40)
        A4_2(5,5)= sqrt(d1/d10); A4_2(5,6)= sqrt(d1/d10); A4_2(5,7)= sqrt(d4/d5);
        A4_2(6,1)= sqrt(d5/d9);  A4_2(6,2)= sqrt(d4/d9);  A4_2(6,3)=-sqrt(d4/d9)
        A4_2(7,2)=-sqrt(d7/d30); A4_2(7,3)=-sqrt(d7/d30); A4_2(7,4)= sqrt(d7/d15); A4_2(7,5)= sqrt(d1/d30); A4_2(7,6)=-sqrt(d1/d30)
        A4_2(8,5)=-sqrt(d2/d5);  A4_2(8,6)=-sqrt(d2/d5);  A4_2(8,7)= sqrt(d1/d5)
        do i=1,7; do j=1,8; A4(d41c(i),d41r(j))=A4_1(i,j); enddo; enddo
        do i=1,8; do j=1,7; A4(d42c(i),d42r(j))=A4_2(i,j); enddo; enddo

        A6_1=z; A6_2=z; A6=0
        A6_1(1,1)=-sqrt(d5/d52); A6_1(1,2)=-sqrt(d5/d52); A6_1(1,3)=-sqrt(d35/d78);A6_1(1,4)= sqrt(d7/d39);  A6_1(1,5)=-sqrt(d7/d39)
        A6_1(2,4)=-sqrt(d1/d26); A6_1(2,5)=-sqrt(d1/d26); A6_1(2,6)=-sqrt(d9/d26); A6_1(2,7)= sqrt(d15/d52);A6_1(2,8)=-sqrt(d15/d52)
       A6_1(3,7)=-sqrt(d1/d156);A6_1(3,8)=-sqrt(d1/d156);A6_1(3,9)=-sqrt(d11/d78);A6_1(3,10)=sqrt(d11/d26);A6_1(3,11)=-sqrt(d11/d26)
        A6_1(4,1)= sqrt(d1/d4);  A6_1(4,2)= sqrt(d1/d4);  A6_1(4,3)= sqrt(d1/d42); A6_1(4,4)= sqrt(d5/d21);  A6_1(4,5)=-sqrt(d5/d21)
        A6_1(5,4)= sqrt(d3/d14); A6_1(5,5)=sqrt(d3/d14);  A6_1(5,6)=sqrt(d3/d14);  A6_1(5,7)= sqrt(d5/d28);  A6_1(5,8)=-sqrt(d5/d28)
        A6_1(6,7)= sqrt(d11/d84);A6_1(6,8)=sqrt(d11/d84);A6_1(6,9)=sqrt(d25/d42); A6_1(6,10)= sqrt(d1/d14); A6_1(6,11)=-sqrt(d1/d14)
        A6_1(7,1)=-sqrt(d2/d13); A6_1(7,2)=-sqrt(d2/d13);A6_1(7,3)=sqrt(d48/d91);A6_1(7,4)= sqrt(d15/d182);A6_1(7,5)=-sqrt(d15/d182)
        A6_1(8,4)=-sqrt(d45/d182);A6_1(8,5)=-sqrt(d45/d182);A6_1(8,6)=sqrt(d40/d91);A6_1(8,7)= sqrt(d3/d91); A6_1(8,8)=-sqrt(d3/d91)
        A6_1(9,7)=-sqrt(d33/d91);A6_1(9,8)=-sqrt(d33/d91);A6_1(9,9)=sqrt(d24/d91);A6_1(9,10)=sqrt(d1/d182);A6_1(9,11)=-sqrt(d1/d182)
        A6_1(10,10)=-sqrt(d1/d2);A6_1(10,11)=-sqrt(d1/d2)
        A6_2(1,1)=-sqrt(d6/d13); A6_2(1,2)= sqrt(d7/d13); A6_2(1,3)=-sqrt(d7/d13)
        A6_2(2,2)=-sqrt(d5/d78); A6_2(2,3)=-sqrt(d5/d78); A6_2(2,4)=-sqrt(d16/d39);A6_2(2,5)= sqrt(d3/d13);  A6_2(2,6)=-sqrt(d3/d13)
        A6_2(3,5)=-sqrt(d1/d52); A6_2(3,6)=-sqrt(d1/d52);A6_2(3,7)=-sqrt(d10/d39);A6_2(3,8)=sqrt(d55/d156);A6_2(3,9)=-sqrt(d55/d156)
        A6_2(4,2)= d1;  A6_2(4,3)=-d1
        A6_2(5,2)= sqrt(d5/d21); A6_2(5,3)=sqrt(d5/d21);  A6_2(5,4)= sqrt(d2/d21); A6_2(5,5)= sqrt(d3/d14);  A6_2(5,6)=-sqrt(d3/d14)
        A6_2(6,5)= sqrt(d5/d28); A6_2(6,6)=sqrt(d5/d28);  A6_2(6,7)= sqrt(d8/d21); A6_2(6,8)= sqrt(d11/d84);A6_2(6,9)=-sqrt(d11/d84)
        A6_2(7,8)= sqrt(d1/d14); A6_2(7,9)= sqrt(d1/d14); A6_2(7,10)= sqrt(d6/d7)
        A6_2(8,1)= sqrt(d7/d13); A6_2(8,2)= sqrt(d6/d13); A6_2(8,3)=-sqrt(d6/d13)
        A6_2(9,2)=-sqrt(d18/d91);A6_2(9,3)=-sqrt(d18/d91);A6_2(9,4)= sqrt(d45/d91);A6_2(9,5)=sqrt(d5/d91); A6_2(9,6)=-sqrt(d5/d91)
        A6_2(10,5)=-sqrt(d55/d182);A6_2(10,6)=-sqrt(d55/d182);A6_2(10,7)=sqrt(d33/d91);A6_2(10,8)=sqrt(d3/d182);
        A6_2(10,9)=-sqrt(d3/d182)
        A6_2(11,8)=-sqrt(d3/d7); A6_2(11,9)=-sqrt(d3/d7); A6_2(11,10)= sqrt(d1/d7)
        do i=1,10; do j=1,11; A6(d61c(i),d61r(j))=A6_1(i,j); enddo; enddo
        do i=1,11; do j=1,10; A6(d62c(i),d62r(j))=A6_2(i,j); enddo; enddo

! divide y by i
        do i=2,8,3; BV2(i)=BV2(i)/(CI); enddo  
        do i=2,14,3; BV4(i)=BV4(i)/(CI); enddo  
        do i=2,20,3; BV6(i)=BV6(i)/(CI); enddo  
        AV2=z; Do i=1,9;  do j=1,9;  
          if ((j.eq.4 .or. j.eq.5).and.(i.eq.1 .or. i.eq.6)) then
            AV2(i)=AV2(i)+A2(i,j)*DReal(BV2(j));
          else if ((j.eq.4 .or. j.eq.5).and.(i.eq.3)) then
            AV2(i)=AV2(i)+A2(i,j)*DImag(BV2(j))*CI;
          else
            AV2(i)=AV2(i)+A2(i,j)*BV2(j)
          endif
        enddo; enddo   
        AV4=z; Do i=1,15; do j=1,15; 
          if ((j.eq.4 .or. j.eq.5).and.(i.eq.1 .or. i.eq.10)) then
            AV4(i)=AV4(i)+A4(i,j)*DReal(BV4(j));
          else if ((j.eq.4 .or. j.eq.5).and.(i.eq.5)) then
            AV4(i)=AV4(i)+A4(i,j)*DImag(BV4(j))*CI;
          else
            AV4(i)=AV4(i)+A4(i,j)*BV4(j); 
          endif
        enddo; enddo   
        AV6=z; Do i=1,21; do j=1,21; 
          if ((j.eq.4 .or. j.eq.5).and.(i.eq.1 .or. i.eq.14)) then
            AV6(i)=AV6(i)+A6(i,j)*DReal(BV6(j));
          else if ((j.eq.4 .or. j.eq.5).and.(i.eq.7)) then
            AV6(i)=AV6(i)+A6(i,j)*DImag(BV6(j))*CI;
          else
            AV6(i)=AV6(i)+A6(i,j)*BV6(j)
          endif
        enddo; enddo   

        do i=1,3
          lamda=2*i
          ii=0
          do j=1,3
            t=lamda-2+j  ! t=lamda-1,lamda,lamda+1 
            do k=1,t+1
              ii=ii+1
              if (i.eq.1) then; AltpR(i,j,k)=DREAL(AV2(ii)); AltpI(i,j,k)=DIMAG(AV2(ii)); endif
              if (i.eq.2) then; AltpR(i,j,k)=DREAL(AV4(ii)); AltpI(i,j,k)=DIMAG(AV4(ii)); endif
              if (i.eq.3) then; AltpR(i,j,k)=DREAL(AV6(ii)); AltpI(i,j,k)=DIMAG(AV6(ii)); endif
            enddo
          enddo
        enddo

        if (check2(4)) then
          write(Idebug,'("A2tp p=   0    +1    +2    +3")')
          write(Idebug,'("t=1    ",3(2F5.1))')(AltpR(1,1,k),AltpI(1,1,k),k=1,2)
          write(Idebug,'("t=2    ",5(2F5.1))')(AltpR(1,2,k),AltpI(1,2,k),k=1,3)
          write(Idebug,'("t=3    ",7(2F5.1))')(AltpR(1,3,k),AltpI(1,3,k),k=1,4)
          write(Idebug,'("A4tp p=   0    +1    +2    +3    +4    +5")')
          write(Idebug,'("t=3    ",7(2F5.1))') (AltpR(2,1,k),AltpI(2,1,k),k=1,4)
          write(Idebug,'("t=4    ",9(2F5.1))') (AltpR(2,2,k),AltpI(2,2,k),k=1,5)
          write(Idebug,'("t=5    ",11(2F5.1))')(AltpR(2,3,k),AltpI(2,3,k),k=1,6)
          write(Idebug,'("A6tp p=   0    +1    +2    +3    +4    +5    +6    +7")')
          write(Idebug,'("t=5    ",11(2F5.1))')(AltpR(3,1,k),AltpI(3,1,k),k=1,6)
          write(Idebug,'("t=6    ",13(2F5.1))')(AltpR(3,2,k),AltpI(3,2,k),k=1,7)
          write(Idebug,'("t=7    ",15(2F5.1))')(AltpR(3,3,k),AltpI(3,3,k),k=1,8)
        endif  ! check2(4)       

      endif ! IntN1.eq.2
!      
      RETURN
      end subroutine  A_B_Mat  
!
!-----------------------------------------------------------------------
!
      SUBROUTINE writeIntPara(iWrite)
!
!  Writes the ALtp, BLli, or Tli intensity parameters.
! 
      IMPLICIT none
      integer i,j,k, ilamda, iWrite
      CHARACTER(LEN=1), PARAMETER :: axes(3)= (/"x","y","z"/)
!      
      if (iWrite.eq.1) then
        WRITE(IO,'("Altp Intensity Parameters")')
        write(IO,'((2(3X,I3,1X,"A21",i1,"R",F10.3)))') (200+2*j-1,j-1,AltpR(1,1,j),j=1,2)
        write(IO,'((2(3X,I3,1X,"A21",i1,"I",F10.3)))') (200+2*j,  j-1,AltpI(1,1,j),j=1,2)
        write(IO,'((3(3X,I3,1X,"A22",i1,"R",F10.3)))') (204+2*j-1,j-1,AltpR(1,2,j),j=1,3)
        write(IO,'((3(3X,I3,1X,"A22",i1,"I",F10.3)))') (204+2*j,  j-1,AltpI(1,2,j),j=1,3)
        write(IO,'((4(3X,I3,1X,"A23",i1,"R",F10.3)))') (210+2*j-1,j-1,AltpR(1,3,j),j=1,4)
        write(IO,'((4(3X,I3,1X,"A23",i1,"I",F10.3)))') (210+2*j,  j-1,AltpI(1,3,j),j=1,4)
        write(IO,'((4(3X,I3,1X,"A43",i1,"R",F10.3)))') (218+2*j-1,j-1,AltpR(2,1,j),j=1,4)
        write(IO,'((4(3X,I3,1X,"A43",i1,"I",F10.3)))') (218+2*j,  j-1,AltpI(2,1,j),j=1,4)
        write(IO,'((5(3X,I3,1X,"A44",i1,"R",F10.3)))') (226+2*j-1,j-1,AltpR(2,2,j),j=1,5)
        write(IO,'((5(3X,I3,1X,"A44",i1,"I",F10.3)))') (226+2*j,  j-1,AltpI(2,2,j),j=1,5)
        write(IO,'((6(3X,I3,1X,"A45",i1,"R",F10.3)))') (236+2*j-1,j-1,AltpR(2,3,j),j=1,6)
        write(IO,'((6(3X,I3,1X,"A45",i1,"I",F10.3)))') (236+2*j,  j-1,AltpI(2,3,j),j=1,6)
        write(IO,'((6(3X,I3,1X,"A65",i1,"R",F10.3)))') (248+2*j-1,j-1,AltpR(3,1,j),j=1,6)
        write(IO,'((6(3X,I3,1X,"A65",i1,"I",F10.3)))') (248+2*j,  j-1,AltpI(3,1,j),j=1,6)
        write(IO,'((7(3X,I3,1X,"A66",i1,"R",F10.3)))') (260+2*j-1,j-1,AltpR(3,2,j),j=1,7)
        write(IO,'((7(3X,I3,1X,"A66",i1,"I",F10.3)))') (260+2*j,  j-1,AltpI(3,2,j),j=1,7)
        write(IO,'((8(3X,I3,1X,"A67",i1,"R",F10.3)))') (274+2*j-1,j-1,AltpR(3,3,j),j=1,8)
        write(IO,'((8(3X,I3,1X,"A67",i1,"I",F10.3)))') (274+2*j  ,j-1,AltpI(3,3,j),j=1,8)
      elseif (iWrite.eq.2) then
        WRITE(IO,'("Blki Intensity Parameters")')
        j=0
        do ilamda=1,3
          do k=0,2*ilamda
            j=j+1
            write(IO,'((3(3X,I3,1X,"B",2i1,A1,"R",F10.3,3X,I3,1X,"B",2i1,A1,"I",F10.3)))')   & 
                   (200+6*(j-1)+(i-1)*2+1,2*ilamda,k,axes(i),BlkiR(j,i),                     &
                    200+6*(j-1)+(i-1)*2+2,2*ilamda,k,axes(i),BlkiI(j,i),i=1,3)
          enddo
        enddo
      elseif (iWrite.eq.3) then
        WRITE(IO,'("Tli Intensity Parameters")')
        if (Lvalue.eq.2) WRITE(IO,'("  PTsig    FTsig    RTsig    PTpiX    FTpiX    RTpiX    PTpiY    FTpiY    RTpiY")')
        if (Lvalue.eq.3) WRITE(IO,'("  DTsig    GTsig    RTsig    DTpiX    GTpiX    RTpiX    DTpiY    GTpiY    RTpiY")')
        write(IO,'(3(3X,(3(I3,1X,F10.3))))') (200+(i-1)*9+j,Tli(i,j),j=1,9)
      elseif (iWrite.eq.4) then
        WRITE(IO,'("Cli Intensity Parameters")')
        if (Lvalue.eq.2) WRITE(IO,'("   C1x      C2x      C1y      C2y      C1z      C2z      C3z")')
        if (Lvalue.eq.3) WRITE(IO,'("   C1x      C2x      C1y      C2y      C1z      C2z      C3z")')
        do i=1,NLigands 
          write(IO,'(3(3X,(3(I3,1X,F10.3))))') (200+(i-1)*7+j,Cli(i,j),j=1,7)
        enddo
      endif
      WRITE(IO,'(106("-"))')
!      
      RETURN
      end subroutine  writeIntPara  
!
!-----------------------------------------------------------------------

END MODULE f_e_LF
!  1777