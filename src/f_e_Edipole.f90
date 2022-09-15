MODULE f_e_Edipole
USE f_e_data
USE f_e_wigner
USE f_e_calculate
!       
!  This module does the electric dipole calculations. 
!

IMPLICIT NONE

CONTAINS
!  Full |SLJMJ> matrix
!  calcEDip() calls buildAMatrix() and calcED(). calcED() calls function transMom(ii,jj,kk).
!
!  Judd-Ofelt (JO)calc in |LSJ> basis
!  calcEDipJO() calls buildJMatrix(), calcU(), diagrs(..) and calcJO()
!
!  Mag Dip calc in |LSJ> basis
!  calcMDipJBasis() calls buildJMatrix(), calcMuJBasis(), diagrs(..) and calcMDJbasis()
!
!-----------------------------------------------------------------------
!
      SUBROUTINE buildAMatrix()
!
!  Builds the ED amplitudes within the original basis:
!  <t,L,S,J,MJ| m(rho) |t',L',S,J',MJ'> for rho = -1,0,+1
!  Notice only states with the same S are connected.
! 
!  The rest must be zero.
!  n may be < n_matrix if symmetry blocking has occurred.
!
!  The ED transition amplitudes are in the units of the A parameters (x10^-12 cm)
!
      IMPLICIT NONE
      integer i,j,ii,jj,k,k1,k2,q,nz,kk,rho,t,pp, allocation_status,NUmax,nele1,nele2,nele3,n1
      INTEGER TwoSi,TwoSj,Li,Lj,TwoJi,TwoJj,TwoMJi,TwoMJj,SENi,SENj
      real*8 D0,D1,D2,X1,decouple,facSen,bit,phase1,phase2,phase3,RH
      real*8 RSi,RSj, RLi,RLj, RJi,RJj, RMJi,RMJj, mat14(14,14)
      REAL*4 Time1,Time2
      COMPLEX*16 CZ,CI,cBkqml,cBkqpl,cBkq0
      PARAMETER(D0=0.0D+00,D1=1.0D+00,D2=2.0D+00,bit=1.0d-6)
      PARAMETER(CZ=(D0,D0), CI=(D0,D1))
!      real*8  Uknq(3,0:6,26)
!      integer Iknq(3,0:6,26),Jknq(3,0:6,26),Nknq(3,0:6)
      logical FirstCall
      save FirstCall
      character label1(20)*7,label2(20)*6
      character(len=99) :: emsg
      data Firstcall/.true./
!
      call CPU_TIME(Time1)
!
!  largest number of non-zero values in the Ukq matrix is 201524, for U61 in f^6 basis.
      NUmax = 201524
      IF (Nelectrons.eq.1) NUmax=26
      IF (Nelectrons.eq.2) NUmax=404
      IF (Nelectrons.eq.3) NUmax=4972
      IF (Nelectrons.eq.4) NUmax=27566
      IF (Nelectrons.eq.5) NUmax=101322
      IF (Nelectrons.eq.6) NUmax=201524
      IF (Nelectrons.eq.7) NUmax=131232
!
      RH=sqrt(D1/D2)
      if (.not.FirstCall) goto 100
!      write(*,'(" Attempting to ALLOCATE Ukq,Ikq,Jkq each (3,0:6,NUmax); NUmax=",I6)') NUmax
      allocate(Ukq(3,0:6,NUmax), stat=allocation_status, errmsg=emsg)
      allocate(Ikq(3,0:6,NUmax), stat=allocation_status, errmsg=emsg)
      allocate(Jkq(3,0:6,NUmax), stat=allocation_status, errmsg=emsg)
      if (allocation_status > 0) then
        write(io,'("***FATAL: error Allocating Ukq(0:2,NUmax):",A)') trim(emsg)
        write(*, '("***FATAL: error Allocating Ukq(0:2,NUmax):",A)') trim(emsg)
        stop
      end if
!      write(*, '(" Attempting to ALLOCATE EDtrans(-1:1,N,N); N=",I5)') n_matrix
      allocate(EDtrans(-1:1,n_matrix,n_matrix), stat=allocation_status, errmsg=emsg)
      if (allocation_status > 0) then
        write(io,'("***FATAL: error Allocating EDtrans(-1:1,N,N):",A)') trim(emsg)
        write(*, '("***FATAL: error Allocating EDtrans(-1:1,N,N):",A)') trim(emsg)
        stop
      end if
!
! zeroise the arrays:
      Ukq=D0; Ikq=0; Jkq=0; Nkq=0
!      Uknq=D0; Iknq=0; Jknq=0; Nknq=0
!
      do i=1,n_matrix
        TwoSi=FullBasis(1,i);  RSi=dble(TwoSi)/D2      ! S
        Li=FullBasis(2,i);     RLi=dble(Li)            ! L
        TwoJi=FullBasis(3,i);  RJi=dble(TwoJi)/D2      ! J
        TwoMJi=FullBasis(4,i); RMJi=dble(TwoMJi)/D2    ! MJ
        ii=FullBasis(5,i)       ! position in N&K order
        SENi=FullBasis(6,i)     ! seniority
        do j=1,n_matrix  ! complete matrix  (not just lower triangle)
          TwoSj=FullBasis(1,j);  RSj=dble(TwoSj)/D2     ! S'
          Lj=FullBasis(2,j);     RLj=dble(Lj)           ! L'
          TwoJj=FullBasis(3,j);  RJj=dble(TwoJj)/D2     ! J'
          TwoMJj=FullBasis(4,j); RMJj=dble(TwoMJj)/D2   ! MJ'
          jj=FullBasis(5,j)
          SENj=FullBasis(6,j)
!          if (Complementary) facSen=DBLE((-1)**(1+(SENi-SENj)/2)) !  (v-v')/2  always integer
!
          if (TwoSi.eq.TwoSj) then  ! delta(S,S')
            do k=1,3
              kk=2*k
              decouple=(-D1)**(RSi+RJi+RLj+dble(kk))*Sqrt((D2*RJi+D1)*(D2*RJj+D1))          &
                        *Wig6J(RJi,  RJj,  DBLE(kk),                      &
                               RLj,  RLi,  RSi)
              do q=0,kk
                X1=(-D1)**(RJi-RMJi)*Wig3J(RJi,   DBLE(kk),  RJj,   &
                                          -RMJi,  DBLE(q),   RMJj)*decouple*UMat(ii,jj,kk)  
                if (abs(X1).gt.1.0d-12) then
                  Nkq(k,q)=Nkq(k,q)+1
                  Ukq(k,q,Nkq(k,q))=X1
                  Ikq(k,q,Nkq(k,q))=i
                  Jkq(k,q,Nkq(k,q))=j
                endif   
!                X1=(-D1)**(RJi-RMJi)*Wig3J(RJi,   DBLE(kk),  RJj,   &
!                                          -RMJi, -DBLE(q),   RMJj)*decouple*UMat(ii,jj,kk)  
!                if (abs(X1).gt.1.0d-12) then
!                  Nknq(k,q)=Nknq(k,q)+1
!                  Uknq(k,q,Nknq(k,q))=X1
!                  Iknq(k,q,Nknq(k,q))=i
!                  Jknq(k,q,Nknq(k,q))=j
!                endif   
              ENDDO ! q
            enddo ! k 
          ENDIF ! S=S'
!
        enddo   !  end j
      enddo  !  end i
!
      if (printMat3(3)) then
        do k=1,3
          do q=0,2*k
            write(imatrix,'(/,"<Ukq> k=",I2," q=",i2,"; Number of non-zero elements:",i5,"; first 100 here")') 2*k,q,Nkq(k,q)
            write(imatrix,'((10(2I4,F8.4)))') (Ikq(k,q,i),Jkq(k,q,i),Ukq(k,q,i),i=1,min(100,Nkq(k,q)))
          enddo  !  q
        enddo  ! k
      endif ! printMat3(3)
!
 100  continue
! zeroise arrays:
      EDtrans=D0
       
      if (IntN1.eq.1) then !  Altp 
! 
        do rho=-1,1,1
          nele2=0
          do k=1,3
            kk=2*k     !  lamda
            do j=1,3
              t=kk-2+j   !  t=kk-1,kk,kk+1
              do pp=-t,t  ! -t,t
                if (abs(AltpR(k,j,abs(pp)+1))+abs(AltpI(k,j,abs(pp)+1)).gt.bit) then

                  nele1=0
                  X1=(-D1)**(rho+pp)*sqrt(D2*t+D1)*Wig3J(     D1,     DBLE(kk),  DBLE(t),    &
                                                        DBLE(rho),-DBLE(rho+pp),  DBLE(pp))
                  phase1=d1; phase3=d1 
                  if (pp.lt.0) then; phase1=(-d1)**(t+pp+1); phase3=-d1; endif   ! phase for Alt-p; must take conjugate (phase3)

                  if (abs(pp+rho).le.kk) then
                    if ((pp+rho).lt.0) then
                      phase2=(-d1)**(pp+rho) ! phase for <|U(l-q)> (& must take transpose)
                      do i=1,Nkq(k,abs(rho+pp))
                        ii=Jkq(k,abs(rho+pp),i); jj=Ikq(k,abs(rho+pp),i)  ! transpose
                        nele1=nele1+1
!  if (i.eq.1) write(io,'("ii,jj=",2I3,"; rho=",I2,"  ltp=",3i2,"; Altp(k,j,|p+1|)*Ukq(k,|rho+pp|),i)=",F8.3,  &
!                         "; Nkq(k,|rho|)=",I5)') ii,jj, rho,kk,t,pp,Altp(k,j,pp+1)*X1*Ukq(k,abs(rho+pp),i),Nkq(k,abs(rho))
                        EDtrans(rho,ii,jj)=EDtrans(rho,ii,jj)+  &
                                phase1*(AltpR(k,j,abs(pp)+1)+phase3*CI*AltpI(k,j,abs(pp)+1))*X1*phase2*Ukq(k,abs(rho+pp),i)
                      enddo  !  i
                    else  
                      do i=1,Nkq(k,abs(rho+pp))
                        ii=Ikq(k,abs(rho+pp),i); jj=Jkq(k,abs(rho+pp),i)
                        nele1=nele1+1
!  if (i.eq.1) write(io,'("ii,jj=",2I3,"; rho=",I2,"  ltp=",3i2,"; Altp(k,j,|p+1|)*Ukq(k,|rho+pp|),i)=",F8.3,  &
!                         "; Nkq(k,|rho|)=",I5)') ii,jj, rho,kk,t,pp,Altp(k,j,pp+1)*X1*Ukq(k,abs(rho+pp),i),Nkq(k,abs(rho))
                        EDtrans(rho,ii,jj)=EDtrans(rho,ii,jj)+  &
                                phase1*(AltpR(k,j,abs(pp)+1)+phase3*CI*AltpI(k,j,abs(pp)+1))*X1*Ukq(k,abs(rho+pp),i)
                      enddo ! i
                    endif ! pp+rho
                  endif  !  |pp+rho| <= lamda
!                write(io,'(" rho=",I2,"  ltp=",3i3,"; number of elements added:",3I6)') rho,kk,t,pp,nele1
                  nele2 = nele2+nele1
                endif  ! Altp # 0
              enddo  !  p
            enddo  !  t
          enddo  !  k
!          write(io,'(" rho=",I2," total number of elements:",I6)') rho,nele2
        enddo  !  rho
      else if (IntN1.eq.2) then   !  Blki
        nele1=0;nele2=0;nele3=0
        do k=1,3
          kk=2*k
          if (k.eq.1) k1=0; if (k.eq.2) k1=3; if (k.eq.3) k1=8;
          do q=0,kk
            k2=k1+q+1
            cBkqml =  RH*(BlkiR(k2,1)+CI*BlkiI(k2,1)) - RH*CI*(BlkiR(k2,2)+CI*BlkiI(k2,2)) 
            cBkqpl = -RH*(BlkiR(k2,1)+CI*BlkiI(k2,1)) - RH*CI*(BlkiR(k2,2)+CI*BlkiI(k2,2)) 
            cBkq0 = BlkiR(k2,3)+CI*BlkiI(k2,3)  ! rho = 0
            if (cdabs(cBkqml).gt.bit) then  ! rho = -1 
              do i=1,Nkq(k,q)
                ii=Ikq(k,q,i); jj=Jkq(k,q,i); nele1=nele1+1
                EDtrans(-1,ii,jj)=EDtrans(-1,ii,jj)+cBkqml*Ukq(k,q,i)
                if (q.ne.0) then
                  phase1=(-D1)**q; phase2=(-D1)**(kk+q+1)
                  EDtrans(+1,jj,ii)=EDtrans(+1,jj,ii)+DCONJG(cBkqml)*phase2*phase1*Ukq(k,q,i)
                endif   
              enddo  !  i
            endif
            if (cdabs(cBkq0).gt.bit) then  ! rho = 0 
              do i=1,Nkq(k,q)
                ii=Ikq(k,q,i); jj=Jkq(k,q,i); nele2=nele2+1 
                EDtrans(0,ii,jj)=EDtrans(0,ii,jj)+cBkq0*Ukq(k,q,i)
                if (q.ne.0) then
                  phase1=(-D1)**q; phase2=(-D1)**(kk+q)
                  EDtrans(0,jj,ii)=EDtrans(0,jj,ii)+DCONJG(cBkq0)*phase2*phase1*Ukq(k,q,i)
                endif   
              enddo  !  i
            endif
            if (cdabs(cBkqpl).gt.bit) then  ! rho = +1
              do i=1,Nkq(k,q)
                ii=Ikq(k,q,i); jj=Jkq(k,q,i); nele3=nele3+1
                EDtrans(+1,ii,jj)=EDtrans(+1,ii,jj)+cBkqpl*Ukq(k,q,i)
                if (q.ne.0) then
                  phase1=(-D1)**q; phase2=(-D1)**(kk+q-1)
                  EDtrans(-1,jj,ii)=EDtrans(-1,jj,ii)+DCONJG(cBkqpl)*phase2*phase1*Ukq(k,q,i)
                endif  
              enddo  !  i
            endif
          enddo  ! q 
        enddo  !  k
        write(io,'(" total number of elements rho=-1:",I6,",  rho=0:",I6,",  rho=+1:",I6)') nele1,nele2,nele3
      endif


      if (printMat3(4).or.printMat3(5)) then
        n1=min(14,n_matrix)
        DO i=1,n1
          TwoJi=FullBasis(3,i);TwoMJi=FullBasis(4,i);ii=FullBasis(5,i)   ! 2J, 2MJ, position in N&K order
          if (mod(TwoJi,2).eq.0) then
            write(label1(I),'(A3,I2)') TS_labels(ii,nelectrons),TwoJi/2  
            write(label2(I),'("(",I2,")")') TwoMJi/2  
          else
            write(label1(I),'(A3,I2,"/2")') TS_labels(ii,nelectrons),TwoJi
            write(label2(I),'("(",I2,"/2)")') TwoMJi  
          endif
        enddo

        if (printMat3(4)) then
          do k=1,3
            do q=0,2*k
              mat14=D0;
              do i=1,min(14*14,Nkq(k,q))  ! 
                ii=Ikq(k,q,i); jj=Jkq(k,q,i)
                if (ii.le.14 .and. jj.le.14) mat14(ii,jj)=Ukq(k,q,i)
              enddo ! i  
              write(imatrix,'(/,"<LSJMJ|U",2i1,"|LSJMJ>")') 2*k,q
              write(imatrix,'(13X,20(2X,A7))') (label1(I),i=1,n1)
              write(imatrix,'(13X,20(3X,A6))') (label2(I),i=1,n1)
              do i=1,n1; write(imatrix,'(A7,A6,(20F9.3))') label1(I),label2(I),(mat14(i,j),j=1,n1); enddo  !  i
            enddo ! q

            phase1=(-D1)**(q) 
            do q=1,2*k
              mat14=D0;
              do i=1,min(14*14,Nkq(k,q))  ! 
                ii=Ikq(k,q,i); jj=Jkq(k,q,i)
!                Li=FullBasis(2,ii);       Lj=FullBasis(2,jj)       ! L
!                TwoJi=FullBasis(3,ii);    TwoJj=FullBasis(3,jj)    ! 2J
!                TwoMJi=FullBasis(4,ii);   TwoMJj=FullBasis(4,jj)   ! 2MJ
!                phase1=(-D1)**((TwoMJi-TwoMJj)/D2) 
                if (ii.le.14 .and. jj.le.14) then
                  mat14(jj,ii)=phase1*Ukq(k,q,i)
                endif  
              enddo ! i  
              write(imatrix,'(/,"<|U",i1,i2,"|>  calc 1")') 2*k,-q
              write(imatrix,'(13X,20(2X,A7))') (label1(I),i=1,n1)
              write(imatrix,'(13X,20(3X,A6))') (label2(I),i=1,n1)
              do i=1,n1; write(imatrix,'(A7,A6,(20F9.3))') label1(I),label2(I),(mat14(i,j),j=1,n1); enddo  !  i
!              mat14=D0;
!              do i=1,min(14*14,Nknq(k,q))  ! 
!                ii=Iknq(k,q,i); jj=Jknq(k,q,i)
!                if (ii.le.14 .and. jj.le.14) mat14(ii,jj)=Uknq(k,q,i)
!              enddo ! i  
!              write(imatrix,'(/,"<|U",i1,i2,"|>  calc 2")') 2*k,-q
!              write(imatrix,'(13X,20(2X,A7))') (label1(I),i=1,n1)
!              write(imatrix,'(13X,20(3X,A6))') (label2(I),i=1,n1)
!              do i=1,n1; write(imatrix,'(A7,A6,(20F9.3))') label1(I),label2(I),(mat14(i,j),j=1,n1); enddo  !  i
            enddo ! q
          enddo !k 
          write(imatrix,'(80("-"))')           
        endif
        
        if (printMat3(5)) then
          write(imatrix,'(/,"EDtrans, X (*PFactor=",F12.4,")")') PFactor
          write(imatrix,'(15X,20(2X,A7,5X))') (label1(I),i=1,n1)
          write(imatrix,'(15X,20(3X,A6,5X))') (label2(I),i=1,n1)
          do i=1,n1
            write(imatrix,'(A7,A6,(14(2F6.1,2X)))') label1(I),label2(I),(PFactor*RH*(EDtrans(-1,i,j)-EDtrans(1,i,j)),j=1,n1)
          enddo  !  i
          write(imatrix,'(/,"EDtrans, Y (*PFactor=",F12.4,")")') PFactor
          write(imatrix,'(15X,20(2X,A7,5X))') (label1(I),i=1,n1)
          write(imatrix,'(15X,20(3X,A6,5X))') (label2(I),i=1,n1)
          do i=1,n1
            write(imatrix,'(A7,A6,(14(2F6.1,2X)))') label1(I),label2(I),(PFactor*RH*(EDtrans(-1,i,j)+EDtrans(1,i,j)),j=1,n1)
          enddo  !  i
          write(imatrix,'(/,"EDtrans, Z (*PFactor=",F12.4,")")') PFactor
          write(imatrix,'(15X,20(2X,A7,5X))') (label1(I),i=1,n1)
          write(imatrix,'(15X,20(3X,A6,5X))') (label2(I),i=1,n1)
          do i=1,n1
            write(imatrix,'(A7,A6,(14(2F6.1,2X)))') label1(I),label2(I),(PFactor*EDtrans(0,i,j),j=1,n1)
          enddo  !  i
          write(imatrix,'(80("-"))')           
        endif
      endif
      
      call CPU_TIME(Time2)
      write( *,'(29x,"Transition matrix in basis:",F8.2," secs")') Time2-Time1
      return  
!
      end subroutine buildAMatrix
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcED()
!
!   This subroutine evaluates electric dipole transition moments
!   using the Altp parameters of Reid & Richardson.
!   The ED transition amplitudes are in the units of the square of the 
!   A parameters (x10^-24 cm^2).
!
      IMPLICIT NONE
      integer*4 i,j,k
      integer*4 IG,IE,IGL,IEL,IIG,IIE
      REAL*8 Rm,Rp,Z,D1
      real*4 time1,time2
      COMPLEX*16 CZ,CI
      parameter(D1=1.0d+00, Z=0.0d+00)
      PARAMETER(CZ=(Z,Z), CI=(Z,D1))
!
      call CPU_TIME(Time1)
      IF (ED(2,1)-ED(1,1).GT.max_GInt .OR. ED(2,2)-ED(1,2).GT.max_EInt) then; WRITE(IO,10); stop; endif
 10   FORMAT(' Array RED overflow in subroutine calcED.')
!
      DO IG=1,ED(2,1)-ED(1,1)+1
        IDEGG(IG)=0
        DO IE=1,ED(2,2)-ED(1,2)+1
          IDEGE(IE)=0
          DO K=1,5
            RED(IG,IE,K)=Z
          Enddo  ! K 
        enddo  ! IE
      enddo  ! IG
!
      EDG=1
      IGL=ED(1,1)
      IIG=0
      DO IG=ED(1,1),ED(2,1)   
        IF (IED.LE.2.AND.ABS(Engs(IG)-Engs(IGL)).GT.Edegen .OR. IED.GE.3.AND.IG.GT.ED(1,1)) EDG=EDG+1
        RENGG(EDG)=Engs(IG)-Engs(1)
        IDEGG(EDG)=IDEGG(EDG)+1
        EDE=1
        IIG=IIG+1
        IEL=ED(1,2)
        IIE=0
        DO IE=ED(1,2),ED(2,2)
          IF (IED.LE.2.AND.ABS(Engs(IE)-Engs(IEL)).GT.Edegen .OR. IED.GE.3.AND.IE.GT.ED(1,2)) EDE=EDE+1
          IF (IG.EQ.ED(1,1)) THEN
            RENGE(EDE)=Engs(IE)-Engs(1)
            IDEGE(EDE)=IDEGE(EDE)+1
          ENDIF
          IIE=IIE+1
          Rp=0.5D+00*(CDABS(transMom(IG,IE,1) + CI*transMom(IG,IE,2)))**2  ! E+^2 = |<Ex> + i.i.<Ey>|^2
          Rm=0.5D+00*(CDABS(transMom(IG,IE,1) - CI*transMom(IG,IE,2)))**2  ! E-^2 = |<Ex> - i.i.<Ey>|^2
          do i=1,3
            CED(IIG,IIE,i)=transMom(IG,IE,i)
            RED(EDG,EDE,i) = RED(EDG,EDE,i) + (CDABS(CED(IIG,IIE,i)))**2 !  Store Ex**2, Ey**2, Ez**2
          enddo ! i
          RED(EDG,EDE,4) = RED(EDG,EDE,4) + Rp  !   E+**2
          RED(EDG,EDE,5) = RED(EDG,EDE,5) + Rm  !   E-**2
          IEL=IE
        enddo   ! IE excited state levels
        IGL=IG
      enddo  ! IG ground state levels
      call CPU_TIME(Time2)
      write( *,'(29x,"Calculating Transitions:   ",F8.2," secs")') Time2-Time1

!
      end subroutine calcED
!
!-----------------------------------------------------------------------
!
      COMPLEX*16 FUNCTION transMom(ii,jj,kk)
!   Evaluates the electric dipole transition moment between eigenvectors ii and jj.
!   kk is 1,2,3 for x,y,z 
!    MAT & CMAT workspace is defined in f_e_data.f90
!
!   The ED transition moments are in the units of the A or B parameters (x10^-12 cm).
!
      IMPLICIT NONE 
      integer ii,jj,kk, i,j, n
      real*8 z,R2,trans,sum
      COMPLEX*16 cTrans,CI
      character*1 lab(3)/"x","y","z"/
      parameter(z=0.0d+00, CI=(z,1.0d+00))
      n=n_matrix
      R2=sqrt(2.0d+00)
!
      if (ii.gt.n_matrix .or. jj.gt.n_matrix) then
        write(io,'("***FATAL: Index out of range in FUNCTION transMom: ii=",i4,", jj=",i4)') ii,jj
        stop
      endif
      cTrans=DCMPLX(z,z)  
      if (MatComplex) then
        if (KK.eq.1) then
          do i=1,n;do j=1,n; cTrans=cTrans+DCONJG(CMAT(i,ii))*CMAT(j,jj)*(EDtrans(-1,i,j)-EDtrans(1,i,j)); enddo;enddo
        else if (kk.eq.2) then
          do i=1,n;do j=1,n; cTrans=cTrans+DCONJG(CMAT(i,ii))*CMAT(j,jj)*(EDtrans(-1,i,j)+EDtrans(1,i,j)); enddo;enddo
        elseif (kk.eq.3) then
          do i=1,n;do j=1,n; cTrans=cTrans+DCONJG(CMAT(i,ii))*CMAT(j,jj)*EDtrans(0,i,j); enddo;enddo
        endif
      else
        ! sum=0;
        if (KK.eq.1) then
          do i=1,n;do j=1,n; cTrans=cTrans+MAT(i,ii)*MAT(j,jj)*(EDtrans(-1,i,j)-EDtrans(1,i,j)); enddo;enddo
        else if (kk.eq.2) then
          do i=1,n;do j=1,n; cTrans=cTrans+MAT(i,ii)*MAT(j,jj)*(EDtrans(-1,i,j)+EDtrans(1,i,j)); enddo;enddo
        elseif (kk.eq.3) then
!        do i=1,n; sum=sum+MAT(i,ii)*MAT(i,jj); enddo;
!        if (Abs(sum).gt.1.0d-06) write(io,'("****Warning: vectors:",I3," and",I3," not orthogonal, V1.V2=",F12.6)') ii,jj,sum
          do i=1,n;do j=1,n; cTrans=cTrans+MAT(i,ii)*MAT(j,jj)*EDtrans(0,i,j); enddo;enddo
        endif
      endif  
!      write(io,'("<",I2,"| E(",A1,") |",i2,"> = ",2E16.5)') ii,lab(kk),jj,cTrans
      if (kk.eq.1) transMom=cTrans/R2; if (kk.eq.2) transMom=CI*cTrans/R2; if (kk.eq.3) transMom=cTrans;
      return
!
      end FUNCTION transMom
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcEDip()
!
!   This subroutine evaluates electric dipole transition moment elements.
!
      IMPLICIT NONE
      integer i,j,rho 
!
      call buildAMatrix()
   !   do rho=-1,1
   !     write(io,'(/,"EDtrans, rho=",I2)') rho
   !     do i=1,min(20,n_matrix)
   !       write(io,'((20F7.4))') (EDtrans(rho,i,j),j=1,min(20,n_matrix))
   !     enddo  !  i
   !   enddo  ! rho
      call calcED()
!
      end subroutine calcEDip
!
!-----------------------------------------------------------------------
!
      subroutine calcU()
!   Evaluates <||U2||>, <||U4||>, <||U6||> in the |SLJ> basis.
!
      IMPLICIT NONE
      character label*9
      logical odd
      integer i,j,k,kk, S,Sp,L,Lp,JJ,Jp, nc,i1,i2,j1,j2
      real*8 z,d1,d2,fac 
      parameter(z=0.0d+00,d1=1.0d+00,d2=2.0d+00)
!
!   Application of Eq(20) in Judd, PR,127,750,(1962).
      do k=1,3
        kk=2*k
        do i=1,n_jmatrix
          S=jbasis(1,i); L=jbasis(2,i); JJ=jbasis(3,i) ! 2S, L, 2J
          do j=1,n_jmatrix
            JO_U(i,j,k)=z
            Sp=jbasis(1,j); Lp=jbasis(2,j); Jp=jbasis(3,j) ! 2S', L', 2J'
            if (S.eq.Sp) then   ! delta(S,S')
              fac=(-d1)**((S+JJ)/2+Lp+kk)*sqrt((JJ+d1)*(Jp+d1))
              fac=fac*WIG6J(dble(L),dble(KK),dble(Lp),dble(Jp)/d2,dble(S)/d2,dble(JJ)/d2)
              if (fac.ne.z) JO_U(i,j,k)=JO_U(i,j,k)+fac*Umat(jbasis(4,i),jbasis(4,j),kk)
            endif  ! delta(S,S')
          enddo ! j
        enddo ! i
        if (printMat2(3)) then
          nc=20; nc=min(nc,n_jmatrix); i1=1; i2=nc; j1=1; j2=nc
          if (printMat2(10)) then; i1=printRng2(1); i2=printRng2(2);  &
                                   j1=printRng2(3); j2=printRng2(4); nc=printRng2(4)-printRng2(3); endif

          odd=.false.; if (mod(Nelectrons,2).eq.1) odd=.true.
          if (odd)      write(iMatrix,'(/,"< ||U^",i1,"|| > ",20(1X,"|",A3,I2,"/2>",1X))') kk,  &
                   (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j),j=j1,j2)
          if (.not.odd) write(iMatrix,'(/,"< ||U^",i1,"|| > ",20(1X,"|",A3,I2,">  ",1X))') kk,  &
                   (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j)/2,j=j1,j2)
          do i=i1,i2
            if (odd)      write(label,'("<",A3,I2,"/2|")') TS_f_labels(jbasis(4,i),Nelectrons),jbasis(3,i)
            if (.not.odd) write(label,'("<",A3,I2,"|  ")') TS_f_labels(jbasis(4,i),Nelectrons),jbasis(3,i)/2
            write(iMatrix,'(1X,A9,1X,20(F9.5,2x))') label,(JO_U(i,j,k),j=j1,j2)
          enddo
        endif  ! printMat2(3)
      enddo ! k
      
      return
!
      end subroutine calcU
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcEDipJO()
!
!   This subroutine evaluates electric dipole transition moment elements.
!
      IMPLICIT NONE
      integer i
      REAL*8 fv1(max_Jstates)
!
      CALL buildJMatrix()
      CALL calcU()
      call diagrs(MAT,max_N,n_jmatrix,engs,fv1,2,io)  ! idiag=2: calc vects
!      write(io,'("n_jmatrix=",I4)') n_jmatrix
!      write(io,'("engs=",/,(10F12.4))') (engs(i),i=1,n_jmatrix)
      call calcJO()

      return
!
      end SUBROUTINE calcEDipJO
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcJO()
!
!   This subroutine evaluates electric dipole transition oscillator strength 
!   f (dimensionless) using the Judd-Ofelt parameters.
!   The values for transitions between states i - j are returned in REDJO(i,j,1).
!   Can be converted to electric dipole transition dipole strength (Debye^2)
!   by D (Debye^2) =  2.127d+6 x f / (Ej - Ei) (cm-1)
!
      IMPLICIT NONE
      integer*4 k,i1,j1,IG,IE
      character*9 label1,label2
      logical odd
      real*8 Z,D1,D2,Fact,UU(3),wf
!  fact = (8*Pi^2*m*c/h) in units of cm-1
      parameter(Z=0.0d+00,D1=1.0d+00,D2=2.0d+00,Fact=1.0847312d+11)
!        write(io,'(" calcJO():",3F12.3)') P(301),P(302),P(303)
!        write(io,'(" calcJO():",3F12.3)') JO_omega(1),JO_omega(2),JO_omega(3)
!
      odd=.false.; if (mod(Nelectrons,2).eq.1) odd=.true.
      if (printMat2(3)) write(iMatrix,'(/,80("-"),/,"Uk matrices in diagonalised |[SL]J> states.",/,  &
                                        34X,"      < ||Uk|| >                    |< ||Uk|| >|^2  ",/, &
                                        34x,"U(2)     U(4)     U(6)         U(2)     U(4)     U(6)")') 
      DO IG=ED(1,1),ED(2,1)  ! start from 1 instead of DO IG=ED(1,1),ED(2,1)
        wfg(ig)=z; wf=z; do i1=1,n_jmatrix; if (mat(i1,ig)**2.gt.wf) then; wfg(ig)=i1; wf=mat(i1,ig)**2; endif; enddo        
        if (printMat2(3)) then 
          if (.not.odd) write(label1,'(A3,"(",I2,")  ")') TS_f_labels(jBasis(4,wfg(ig)),nelectrons),jbasis(3,wfg(ig))/2
          if (     odd) write(label1,'(A3,"(",I2,"/2)")') TS_f_labels(jBasis(4,wfg(ig)),nelectrons),jbasis(3,wfg(ig))
        endif
        DO IE=ED(1,2),ED(2,2)  ! start from 1 instead of DO IE=ED(1,2),ED(2,2)
          wfe(ie)=z; wf=z; do j1=1,n_jmatrix; if (mat(j1,ie)**2.gt.wf) then; wfe(ie)=j1; wf=mat(j1,ie)**2; endif; enddo        
          REDJO(ig,ie,1) = z 
          do k=1,3
            UU(k)=z
            do i1=1,n_jmatrix
              do j1=1,n_jmatrix
                UU(k) = UU(k) + mat(i1,ig)*mat(j1,ie)*JO_U(i1,j1,k)  !  
              enddo   ! j1 basis states
            enddo  ! j2 basis states
            REDJO(ig,ie,1) = REDJO(ig,ie,1) + JO_omega(k)*1.0d-20*UU(k)**2
          enddo  ! k
          if (printMat2(3)) then 
            if (.not.odd)  write(label2,'(A3,"(",I2,")  ")') TS_f_labels(jBasis(4,wfe(ie)),nelectrons),jBasis(3,wfe(ie))/2
            if (    odd)  write(label2,'(A3,"(",I2,"/2)")') TS_f_labels(jBasis(4,wfe(ie)),nelectrons),jBasis(3,wfe(ie))
            write(iMatrix,'("<",A9"||U^k||",A9,"> =",3F9.5,5x,3F9.5)') label1,label2,(UU(k),k=1,3),(UU(k)**2,k=1,3)
          endif
          REDJO(ig,ie,1) = Dielectric*Fact*(Engs(IE)-Engs(IG))*REDJO(ig,ie,1)/(jbasis(3,wfg(ig))+D1)   ! jbasis(3,wfg(i))=2J !??? RED changed to RedJO
          EngJO(ig,ie)=Engs(IE)-Engs(IG)  ! -EAVE
          IntJO(ig,ie)=2.127d+6*108.9d+0*REDJO(ig,ie,1)  
        enddo   ! IE excited state levels
      enddo  ! IG ground state levels
      if (printMat2(3)) write(iMatrix,'(80("-"))') 
!      write(io,'("wfg(ig)=",(20I4))') (wfg(ig),ig=ED(1,1),ED(2,1))
!      write(io,'("wfe(ie)=",(20I4))') (wfe(ie),ie=ED(1,2),ED(2,2))
!
      end subroutine calcJO
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcMDipJBasis()
!
!   This subroutine evaluates magnetic dipole transition moment elements.
!
      IMPLICIT NONE
      integer i
      REAL*8 fv1(max_Jstates)
!
      CALL buildJMatrix()
      CALL calcMuJBasis()
      call diagrs(MAT,max_N,n_jmatrix,engs,fv1,2,io)  ! idiag=2: calc vects
!      write(io,'("n_jmatrix=",I4)') n_jmatrix
!      write(io,'("engs=",/,(10F12.4))') (engs(i),i=1,n_jmatrix)
      call calcMDJbasis()

      return
!
      end SUBROUTINE calcMDipJBasis
!
!-----------------------------------------------------------------------
!
      subroutine calcMuJBasis()
!   Evaluates the <||kL+gS||> matrix elements in the |SLJ> basis.
!
      IMPLICIT NONE
      character label*9, label1*9,label2*9
      logical odd
      integer i,j,k, S,Sp,L,Lp,JJ,Jp, sen,senp, itest,nSLJ, nc,i1,i2,j1,j2
      real*8 z,d1,d2,d3,facL,facS,RS,RL,RJ,g 
      real*8 MD_L(max_Jstates,max_Jstates),MD_S(max_Jstates,max_Jstates),LS(max_Jstates,max_Jstates)  
      parameter(z=0.0d+00,d1=1.0d+00,d2=2.0d+00,d3=3.0d+00)
!
!   Application of Eqs(2)&(4) in Goerller-Walrand,JCP,95,3099(1991). except delta(sen,sen') included.
!              or: Eqs (64) & (65) in  Goerller-Walrand,HPCRE,25,101(1998).
!  In both cases you need an extra (-1)^J factor to convert <LSJ|| mu ||SLJ'> to <LSJ| mu |LSJ'>
!
      nSLJ=0
      do i=1,n_jmatrix
        S=jbasis(1,i); L=jbasis(2,i); JJ=jbasis(3,i); sen=jbasis(5,i) ! 2S, L, 2J, sen
        RS=dble(s)/d2; RL=dble(L);  RJ=dble(JJ)/d2
        do j=1,i
          MD_L(i,j)=z; MD_S(i,j)=z; MD_J(i,j)=z; MD_J(j,i)=z
          Sp=jbasis(1,j); Lp=jbasis(2,j); Jp=jbasis(3,j); senp=jbasis(5,j) ! 2S', L', 2J', sen'
          if (S.eq.Sp .and. L.eq.Lp .and. jbasis(4,i).eq.jbasis(4,j)) then   ! delta(S,S') delta(L,L') delta(sen,sen')
            nSLJ=nSLJ+1
            MD_L(i,j)=(-d1)**((S+JJ)/2+L+1)*sqrt((d2*L+d1)*(JJ+d1)*(Jp+d1)*L*(L+d1)) &
                  *WIG6J(dble(L),dble(JJ)/d2,dble(S)/d2,dble(Jp)/d2,dble(L),d1)
            MD_S(i,j)=(-d1)**((S+Jp)/2+L+1)*sqrt((S+d1)*(JJ+d1)*(Jp+d1)*S/d2*(S/d2+d1)) &
                  *WIG6J(dble(S)/d2,dble(JJ)/d2,dble(L),dble(Jp)/d2,dble(S)/d2,d1)
            MD_J(i,j)=(MD_L(i,j)+d2*MD_S(i,j))
            if (i.ne.j) MD_J(j,i)=-MD_J(i,j) 
!
!            MD_J(i,j)=(MD_L(i,j)+d2*MD_S(i,j))*WIG3J(dble(JJ)/D2,d1,dble(jp)/d2,-dble(JJ)/D2, z,dble(JJ)/d2)
! ! These are equations (66)-(69) of Goerller-Walrand, HPCRE,25,101(1998).
  !          if (Jp.eq.JJ.and.JJ.ne.0) then  
  !            g=D1+(RJ*(RJ+D1)-RL*(RL+D1)+RS*(RS+D1))/(D2*RJ*(RJ+D1))
  !            MD_J(i,j)=g*sqrt(RJ*(RJ+d1)*(d2*RJ+d1))
  !          else if (JP.eq.JJ-2) then  !    (RS+RL-RJ+D1) correct  (RS+RL+RJ-D1) incorrect
  !            MD_J(i,j)=sqrt((RS+RL+RJ+D1)*(RS+RL+RJ-D1)*(RJ+RS-RL)*(RJ+RL-RS)/RJ)/D2
  !          else if (JP.eq.JJ+2) then
  !            MD_J(i,j)=sqrt((RS+RL+RJ+D2)*(RS+RJ+D1-RL)*(RL+RJ+D1-RS)*(RS+RL-RJ)/(RJ+D1))/D2
  !          endif
            
          endif  ! delta(S,S')
        enddo ! j
      enddo ! i
!      
      if (OUTP(7)) then
        write(IO,'("LSJ blocked lower triangle:",I8)') nSLJ 
        write(IO,'("Full matrix lower triangle:",I8,/)') n_jmatrix*(n_jmatrix+1)/2 
      endif   
!
      if (printMat2(4).or.printMat2(5)) then
         write(iMatrix,'(/,"Magnetic matrices in |SLJ> basis")')
         write(iMatrix,'("LSJ blocked lower triangle:",I8)') nSLJ 
         write(iMatrix,'("Full matrix lower triangle:",I8,/)') n_jmatrix*(n_jmatrix+1)/2 
         nc=20; nc=min(nc,TSLJ_f(nelectrons)); i1=1; i2=nc; j1=1; j2=nc
        if (printMat2(10)) then; i1=printRng2(1); i2=printRng2(2);  &
                                 j1=printRng2(3); j2=printRng2(4); nc=printRng2(4)-printRng2(3); endif
      endif ! printMat2(4).or.printMat2(5)                           

     if (printMat2(4)) then
        odd=.false.; if (mod(Nelectrons,2).eq.1) odd=.true.
        if (odd)      write(iMatrix,'(/,"  < ||L|| >   ",20(1X,"|",A3,I2,"/2>"))')   &
                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j),j=j1,j2)
        if (.not.odd) write(iMatrix,'(/,"  < ||L|| >   ",20(1X,"|",A3,I2,">  "))')   &
                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j)/2,j=j1,j2)
        do i=i1,i2
          if (odd)      write(label,'("<",A3,I2,"/2|")') TS_f_labels(jbasis(4,i),Nelectrons),jbasis(3,i)
          if (.not.odd) write(label,'("<",A3,I2,"|  ")') TS_f_labels(jbasis(4,i),Nelectrons),jbasis(3,i)/2
          write(iMatrix,'(3X,A9,1X,20(F8.4,2x))') label,(MD_L(i,j),j=j1,j2)
        enddo
        if (odd)      write(iMatrix,'(/,"  < ||S|| >   ",20(1X,"|",A3,I2,"/2>"))')   &
                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j),j=j1,j2)
        if (.not.odd) write(iMatrix,'(/,"  < ||S|| >   ",20(1X,"|",A3,I2,">  "))')   &
                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j)/2,j=j1,j2)
        do i=i1,i2
          if (odd)      write(label,'("<",A3,I2,"/2|")') TS_f_labels(jbasis(4,i),Nelectrons),jbasis(3,i)
          if (.not.odd) write(label,'("<",A3,I2,"|  ")') TS_f_labels(jbasis(4,i),Nelectrons),jbasis(3,i)/2
          write(iMatrix,'(3X,A9,1X,20(F8.4,2x))') label,(MD_S(i,j),j=j1,j2)
        enddo
      endif ! printMat2(4)
      
      if (printMat2(5)) then
        if (odd)      write(iMatrix,'(/,"<||(L+2S)||> ",20(1X,"|",A3,I2,"/2>"))')   &
                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j),j=j1,j2)
        if (.not.odd) write(iMatrix,'(/,"<||(L+2S)||> ",20(1X,"|",A3,I2,"> "))')   &
                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j)/2,j=j1,j2)
!-------
!  Writes out 7F,5D parts of f6 matrices.
!                 write(iMatrix,'(/,"<||(L+2S)||> ",22(1X,"|",A3,I1,"> "))')   &
!                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j)/2,j=1,7),       &
!                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j)/2,j=12,26)
!-------
        do i=i1,i2
          if (odd)      write(label,'("<",A3,I2,"/2|")') TS_f_labels(jbasis(4,i),Nelectrons),jbasis(3,i)
          if (.not.odd) write(label,'("<",A3,I2,"| ")') TS_f_labels(jbasis(4,i),Nelectrons),jbasis(3,i)/2
          write(iMatrix,'(2X,A9,1X,20(F7.4,2x))') label,(MD_J(i,j),j=j1,j2)
!-------
!          if (i.le.7 .or. (i.ge.12 .and. i.le.26)) write(iMatrix,'(2X,A9,1X,22(F6.3,2x))')   &
!          label,(MD_J(i,j),j=1,7),(MD_J(i,j),j=12,26)
!-------
        enddo
      endif ! printMat2(5)
!      
      if (check2(3)) then
        write(iMatrix,'(/,"SOC matrix elements evaulated by < ||(L)|| >< ||(S)|| > and by using  V11 ")')
        write(iMatrix,'("Disagreeing (if any) matrix elements follow:")')
        itest=0
        do i=1,n_jmatrix
          do j=1,n_jmatrix
            LS(i,j)=z
            do k=1,n_jmatrix
              LS(i,j)=LS(i,j)+MD_L(i,k)*MD_S(k,j)
            enddo ! k
            if (LS(i,j).ne.SOtest(i,j)) then
!            if (LS(i,j).ne.z) then
              itest=itest+1
              if (.not.odd)  write(label1,'(A3,"(",I2,")  ")') TS_f_labels(jBasis(4,i),nelectrons),jBasis(3,i)/2
              if (     odd)  write(label1,'(A3,"(",I2,"/2)")') TS_f_labels(jBasis(4,i),nelectrons),jBasis(3,i)
              if (.not.odd)  write(label2,'(A3,"(",I2,")  ")') TS_f_labels(jBasis(4,j),nelectrons),jBasis(3,j)/2
              if (     odd)  write(label2,'(A3,"(",I2,"/2)")') TS_f_labels(jBasis(4,j),nelectrons),jBasis(3,j)
              write(iMatrix,'("<",A9"||(L.S)||",A9,"> =",2F12.5)') label1,label2,LS(i,j),SOtest(i,j)
            endif 
          enddo ! j
        enddo ! i
        write(iMatrix,'("number of disagreeing matrix elements:",I5)') itest
      endif
!      
      return
!
      end subroutine calcMuJBasis
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcMDJbasis()
!
!   This subroutine evaluates magnetic dipole transition amplitudes 
!   between states i -> j and these are returned in RMD(i,j,1).
!   They are in units of (Bohr Magneton)
!
      IMPLICIT NONE
      integer*4 i1,j1,IG,IE,j
      character*9 label1,label2
      logical odd
      real*8 Z,D1,D2,Fact,wf
!  fact = (8*Pi^2*m*c/h) in units of cm-1
      parameter(Z=0.0d+00,D1=1.0d+00,D2=2.0d+00,Fact=1.0847312d+11)
!
!  Writes out 7F,5D parts of f6 eigenvectors.
!------------------
!                 write(iMatrix,'(/,5X,22(1X,"|",A3,I1,"> "))')   &
!                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j)/2,j=1,7),       &
!                 (TS_f_labels(jbasis(4,j),Nelectrons), jbasis(3,j)/2,j=12,26)
!      write(iMatrix,'("7F0 ",22F8.4)') (mat(j,1),j=1,7), (mat(j,1),j=12,26)           
!      write(iMatrix,'("7F1 ",22F8.4)') (mat(j,2),j=1,7), (mat(j,2),j=12,26)           
!      write(iMatrix,'("5D0 ",22F8.4)') (mat(j,8),j=1,7), (mat(j,8),j=12,26)           
!      write(iMatrix,'("5D1 ",22F8.4)') (mat(j,9),j=1,7), (mat(j,9),j=12,26)           
!      write(iMatrix,'("5D2 ",22F8.4)') (mat(j,10),j=1,7),(mat(j,10),j=12,26)           
!------------------

      odd=.false.; if (mod(Nelectrons,2).eq.1) odd=.true.
      if (printMat2(5)) write(iMatrix,'(/,80("-"),/,"(L+2S) matrices in diagonalised |[SL]J> states.")') 
      MDG=0
      DO IG=MD(1,1),MD(2,1)   
        MDG=MDG+1
        wfg(MDG)=z; wf=z; do i1=1,n_jmatrix; if (mat(i1,IG)**2.gt.wf) then; wfg(MDG)=i1; wf=mat(i1,IG)**2; endif; enddo        
        if (printMat2(5)) then 
          if (.not.odd) write(label1,'(A3,"(",I2,")  ")') TS_f_labels(jBasis(4,wfg(MDG)),nelectrons),jbasis(3,wfg(MDG))/2
          if (     odd) write(label1,'(A3,"(",I2,"/2)")') TS_f_labels(jBasis(4,wfg(MDG)),nelectrons),jbasis(3,wfg(MDG))
        endif
        MDE=0
        DO IE=MD(1,2),MD(2,2)   
          MDE=MDE+1
          wfe(MDE)=z; wf=z; do j1=1,n_jmatrix; if (mat(j1,ie)**2.gt.wf) then; wfe(MDE)=j1; wf=mat(j1,ie)**2; endif; enddo        
          RMD(MDG,MDE,1) = z 
            do i1=1,n_jmatrix
              do j1=1,n_jmatrix
                RMD(MDG,MDE,1) = RMD(MDG,MDE,1) + mat(i1,IG)*mat(j1,ie)*MD_J(i1,j1)  !  
              enddo   ! j1 basis states
            enddo  ! j2 basis states
           if (printMat2(5)) then 
            if (.not.odd)  write(label2,'(A3,"(",I2,")  ")') TS_f_labels(jBasis(4,wfe(MDE)),nelectrons),jBasis(3,wfe(MDE))/2
            if (    odd)  write(label2,'(A3,"(",I2,"/2)")') TS_f_labels(jBasis(4,wfe(MDE)),nelectrons),jBasis(3,wfe(MDE))
            write(iMatrix,'("<",A9"||(L+2S)||",A9,">^2 =",3F9.5)') label1,label2,RMD(MDG,MDE,1)**2
          endif
!          REDJO(MDG,ie,1) = Dielectric*Fact*(Engs(IE)-Engs(ig))*RED(ig,ie,1)/(jbasis(3,wfg(ig))+D1)   ! jbasis(3,wfg(i))=2J
          EngJO(MDG,MDE)=Engs(ie)-Engs(ig)  ! -EAVE
!          IntJO(MDG,ie)=2.127d+6*108.9d+0*REDJO(ig,ie,1)  
          IntJO(MDG,MDE)=RMD(MDG,MDE,1)  
!          WRITE(IO,'("IG,IE,MDG,MDE,EngJO,IntJO=",4I3,F10.2,f10.4)') IG,IE,MDG,MDE,EngJO(MDG,MDE),IntJO(MDG,MDE)
        enddo   ! IE excited state levels
      enddo  ! IG ground state levels
      if (printMat2(5)) write(iMatrix,'(80("-"))') 
!
      end subroutine calcMDJbasis
!
!-----------------------------------------------------------------------

END MODULE f_e_Edipole
!  812 lines 