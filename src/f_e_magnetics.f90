MODULE f_e_magnetics
USE f_e_data
USE f_e_wigner
!       
!  This module does the magnetic calculations.  !  LAST MODIFIED: 1/03/15
!

IMPLICIT NONE

CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcMagMom()
! 
! Calculates the Lx,Ly,Lz and Sx,Sy,Sz matrices into MM(i,j,max_MM).
! Only the non-zero lower triangle is stored. i=1,2 for L,S and j=1,2,3 for x,y,z.
!
! max_MM=36000 (in f_e_data.f90) to cover the maximum of 35,912 for Lx,Ly of f7.
!
! TODO: this needs not to be so big as the delta(NKi,NKj), MEs cannot connect 
! different |SL> term symbols.
!
      IMPLICIT NONE
      integer i,j,NKi,NKj,k,q,nz(6),kk,k1,k2,ik,jk,nSLJ
      INTEGER TwoSi,Li,TwoJi,TwoMJi,SENi, TwoSj,Lj,TwoJj,TwoMJj,SENj
      real*8 Z,D0,D1,D2,RH,bit,facSen,LS(2,3)
      real*8 RSi,RLi,RJi,RMJi, RSj,RLj,RJj,RMJj
      real*8 fac_m1,fac_0,fac_p1,Lm,Lm1,Lp1,Lx,Ly,Lz,Sm,Sm1,Sp1,Sx,Sy,Sz
      Real*8 MAT2(max_N,max_N),MAT3(max_N,max_N),MAT4(max_N,max_N)
      character*1 lab(3)/"x","y","z"/
      PARAMETER(Z=0.0d+00,D0=0.0D+00,D1=1.0D+00,D2=2.0D+00, bit=1.0D-10)
!
      RH=sqrt(D1/D2); nSLJ=0 
      do k1=1,2; do k2=1,3; nMu(k1,k2)=0; enddo; enddo
!      facSen=D1
!
      do i=1,n_matrix
        TwoSi=FullBasis(1,i);  RSi=DBLE(TwoSi)/D2
        Li=FullBasis(2,i);     RLi=DBLE(Li) 
        TwoJi=FullBasis(3,i);  RJi=DBLE(TwoJi)/D2
        TwoMJi=FullBasis(4,i);RMJi=DBLE(TwoMJi)/D2
        NKi=FullBasis(5,i) ! N&K number
!        SENi=FullBasis(6,i)
        do j=1,i
          TwoSj=FullBasis(1,j);  RSj=DBLE(TwoSj)/D2
          Lj=FullBasis(2,j);     RLj=DBLE(Lj)
          TwoJj=FullBasis(3,j);  RJj=DBLE(TwoJj)/D2
          TwoMJj=FullBasis(4,j);RMJj=DBLE(TwoMJj)/D2
          NKj=FullBasis(5,j)
!          SENj=FullBasis(6,j)
!          if (Complementary) facSen=DBLE((-1)**(1+(SENi-SENj)/2)) !  (v-v')/2  always integer

          if (Li.eq.Lj .and. TwoSi.eq.TwoSj .and. NKi.eq.NKj) then
            nSLJ=nSLJ+1           
            fac_m1=(-D1)**(RJi-RMJi)*Wig3J(RJi,D1,RJj, -RMJi,-D1,RMJj)
            fac_0 =(-D1)**(RJi-RMJi)*Wig3J(RJi,D1,RJj, -RMJi, D0,RMJj)
            fac_p1=(-D1)**(RJi-RMJi)*Wig3J(RJi,D1,RJj, -RMJi,+D1,RMJj)
            Lm=sqrt((D2*RJi+D1)*(D2*RJj+D1)*RLi*(RLi+D1)*(D2*RLi+D1))*(-D1)**(RLi+RSi+RJi+D1)* &
                 Wig6J(RLi,RJi,RSi,RJj,RLi,D1)
            Sm=sqrt((D2*RJi+D1)*(D2*RJj+D1)*RSi*(RSi+D1)*(D2*RSi+D1))*(-D1)**(RLi+RSi+RJj+D1)* &
                 Wig6J(RSi,RJi,RLi,RJj,RSi,D1)
            Lm1=fac_m1*Lm; Lz=fac_0*Lm; Lp1=fac_p1*Lm
            Sm1=fac_m1*Sm; Sz=fac_0*Sm; Sp1=fac_p1*Sm
            LS(1,1)=RH*(Lm1-Lp1); LS(1,2)=RH*(Lm1+Lp1); LS(1,3)=Lz   ! Lx, Ly/i, Lz
            LS(2,1)=RH*(Sm1-Sp1); LS(2,2)=RH*(Sm1+Sp1); LS(2,3)=Sz   ! Sx, Sy/i, Sz
!            if(abs(LS(1,1)).gt.bit) then
!              write(iMatrix,'("Ji,MJi=",2F4.1,", Jj,MJj=",2F4.1,", fac_m1,fac_p1,Lm,Lm1,Lp1,L="2F4.1,4F5.2)') &
!                 RJi,RMJi,RJj,RMJj,fac_m1,fac_p1,Lm,Lm1,Lp1,LS(1,1)
!            endif 
            do k1=1,2
              do k2=1,3  
                if (abs(LS(k1,k2)).gt.bit) then
                  nMu(k1,k2)=nMu(k1,k2)+1
                  ImuIndex(k1,k2,nMu(k1,k2))=i
                  JmuIndex(k1,k2,nMu(k1,k2))=j
                  MM(k1,k2,nMu(k1,k2))=LS(k1,k2)
                endif
              enddo ! k2=x,y,z
            enddo ! k1=L,S
            
          endif ! Li.eq.Lj .and. TwoSi.eq.TwoSj .and. SENi.eq.SENj
!
        enddo   !  end j
      enddo  !  end i
!
      if (OUTP(7)) then
        write(IO,'(/," Magnetic matrices, non-zero elements:")')
        write(IO,'("Lx:",I8,",   Ly/i:",I8,",   Lz:",I8)') (nMu(1,k),k=1,3)
        write(IO,'("Sx:",I8,",   Sy/i:",I8,",   Sz:",I8)') (nMu(2,k),k=1,3)
        write(IO,'("All tLS blocked lower triangle:",I8)') nSLJ 
        write(IO,'("All Full matrix lower triangle:",I8,/)') n_matrix*(n_matrix+1)/2 
      endif
!      
      if (printMat3(6)) then
        write(iMatrix,'(/," Magnetic matrices, non-zero elements:")')
        write(iMatrix,'("Lx:",I8,",   Ly/i:",I8,",   Lz:",I8)') (nMu(1,k),k=1,3)
        write(iMatrix,'("Sx:",I8,",   Sy/i:",I8,",   Sz:",I8)') (nMu(2,k),k=1,3)
        write(iMatrix,'("All tLS blocked lower triangle:",I8)') nSLJ 
        write(iMatrix,'("All Full matrix lower triangle:",I8,/)') n_matrix*(n_matrix+1)/2 
        do i=1,n_matrix; do j=1,n_matrix; MAT4(i,j)=z; enddo; enddo
        do kk=1,3
          write(iMatrix,'(/,"L",A1,":")') lab(kk)
          do i=1,n_matrix; do j=1,n_matrix; MAT1(i,j)=z; enddo; enddo
          do k=1,nMu(1,kk)  ! L  
            ik=ImuIndex(1,kk,k); jk=JmuIndex(1,kk,k)         
            MAT1(ik,jk)=MM(1,kk,k)
            if (kk.eq.2 .and. ik.ne.jk) MAT1(jk,ik)=-MAT1(ik,jk)  
            if (kk.ne.2 .and. ik.ne.jk) MAT1(jk,ik)= MAT1(ik,jk)  
          enddo
          do i=1,n_matrix; do j=1,n_matrix; MAT2(i,j)=MAT1(i,j); enddo; enddo
          do i=1,min(n_matrix,20)
            write(iMatrix,'((20F5.2))') (Mat2(i,j),j=1,min(n_matrix,20))
          enddo
          write(iMatrix,'(/,"S",A1,":")') lab(kk)
          do i=1,n_matrix; do j=1,n_matrix; MAT1(i,j)=z; enddo; enddo
          do k=1,nMu(2,kk)  ! S  
            ik=ImuIndex(2,kk,k); jk=JmuIndex(2,kk,k)         
            MAT1(ik,jk)=MM(2,kk,k)
            if (kk.eq.2 .and. ik.ne.jk) MAT1(jk,ik)=-MM(2,kk,k) 
            if (kk.ne.2 .and. ik.ne.jk) MAT1(jk,ik)= MM(2,kk,k) 
          enddo
          do i=1,min(n_matrix,20)
            write(iMatrix,'((20F5.2))') (Mat1(i,j),j=1,min(n_matrix,20))
          enddo
          do i=1,n_matrix; do j=1,n_matrix; MAT3(i,j)=z; enddo; enddo
          if (kk.eq.2) then  
            do i=1,n_matrix; do j=1,n_matrix; do k=1,n_matrix; MAT3(i,j)=MAT3(i,j)-MAT1(i,k)*MAT2(k,j); enddo; enddo; enddo
          else
            do i=1,n_matrix; do j=1,n_matrix; do k=1,n_matrix; MAT3(i,j)=MAT3(i,j)+MAT1(i,k)*MAT2(k,j); enddo; enddo; enddo
          endif
          write(iMatrix,'(/,"L",A1,"S",A1,":")')  lab(kk),lab(kk)
          do i=1,min(n_matrix,20)
            write(iMatrix,'((20F5.2))') (Mat3(i,j),j=1,min(n_matrix,20))
          enddo
          do i=1,n_matrix; do j=1,n_matrix; MAT4(i,j)=MAT4(i,j)+MAT3(i,j); enddo; enddo
        enddo ! kk=x,y,z
        write(iMatrix,'(/,"LxSx-LySy+LzSz:")')   
        do i=1,min(n_matrix,20)
          write(iMatrix,'((20F5.2))') (Mat4(i,j),j=1,min(n_matrix,20))
        enddo
      ENDIF ! printMat3(6)

      return
!
      end subroutine calcMagMom
!
!-----------------------------------------------------------------------
!
      COMPLEX*16 FUNCTION calcMu(ii,jj,kk)
!
!   Evaluates <ii|mu(kk)|jj> = K(kk) L(kk) + ge S(kk) between eigenvectors ii and jj.
!   kk is 1,2,3 for x,y,z 
!    MAT1 is workspace defined in f_e_data.f90
!   Note that calcMu can be complex if the WFs are complex.
!   Both Ly and Sy are /i, so calcMu(ii,jj,2) is /i (divided by i).
!
      IMPLICIT NONE 
      integer ii,jj,kk, i,j,k,ik,jk
      real*8 z,mu
      COMPLEX*16 cmu
      character*1 lab(3)/"x","y","z"/
      parameter(z=0.0d+00)
!
      if (ii.gt.n_matrix .or. jj.gt.n_matrix) then
        write(io,'("***FATAL: Index out of range in FUNCTION calcMu: ii=",i4,", jj=",i4)') ii,jj
        stop
      endif
      do i=1,n_matrix
        do j=1,n_matrix
          MAT1(i,j)=z
        enddo
      enddo
      do k=1,nMu(1,kk)  ! L  
        ik=ImuIndex(1,kk,k); jk=JmuIndex(1,kk,k)         
        MAT1(ik,jk)=RK(kk)*MM(1,kk,k)
        if (kk.eq.2 .and. ik.ne.jk) MAT1(jk,ik)=-MAT1(ik,jk)  
        if (kk.ne.2 .and. ik.ne.jk) MAT1(jk,ik)= MAT1(ik,jk)  
      enddo
      do k=1,nMu(2,kk)  ! S  
        ik=ImuIndex(2,kk,k); jk=JmuIndex(2,kk,k)         
        MAT1(ik,jk)=MAT1(ik,jk)+Ge*MM(2,kk,k)
        if (kk.eq.2 .and. ik.ne.jk) MAT1(jk,ik)=MAT1(jk,ik)-Ge*MM(2,kk,k) 
        if (kk.ne.2 .and. ik.ne.jk) MAT1(jk,ik)=MAT1(jk,ik)+Ge*MM(2,kk,k) 
      enddo
!      
      if (printMat3(7) .and. ii.eq.1 .and. jj.eq.1) then
        write(iMatrix,'(/,"Mat1 = k",A1," L(",A1,") + ge S(",A1,")")') lab(kk),lab(kk),lab(kk)
        do i=1,min(20,n_matrix)
          write(iMatrix,'(20F8.4)') (MAT1(i,j),j=1,min(20,n_matrix))
        enddo  
      endif

      mu=z; cmu=DCMPLX(z,z) 
      if (MatComplex) then
        do i=1,n_matrix
          do j=1,n_matrix
            cmu = cmu + DCONJG(CMAT(i,ii))*MAT1(i,j)*CMAT(j,jj)
          enddo
        enddo
      else
        do i=1,n_matrix
          do j=1,n_matrix
            mu = mu + MAT(i,ii)*MAT1(i,j)*MAT(j,jj)
          enddo
        enddo
        cmu=DCMPLX(mu,0.0d+00)
      endif  
      if (printMat3(7) .and. ii.le.min(20,n_matrix) .and. jj.le.min(20,n_matrix)) then
        write(iMatrix,'("<",I2,"| mu(",A1,") |",i2,"> = ",2F10.5)') ii,lab(kk),jj,cmu
      endif
      calcMu=cmu
      return
!
      end FUNCTION calcMu
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcGval()
!
!   Evaluates the g-values for a simple doublet (Kramers or not) given by states i & j.
!   Principla axes are always calculated (regardless of g_axes value)
!
      IMPLICIT NONE
      integer i,j,k,i1,j1
      real*8 Z,d1,d2,test 
      complex*16 gx11,gx12,gx21,gx22, gy11,gy12,gy21,gy22, gz11,gz12,gz21,gz22, t,ci
      real*8 G2(3,3),eng(3),fv1(3)
      parameter(z=0.0d+00, D1=1.0d+00, D2=2.0d+00)
!
      ci=(z,d1)
      do i=1,ngexs  ! (gexs(2)-gexs(1)+1)/2
        j=gexs(2*i-1)   ! gexs(1)+2*(i-1)
        gx11=   calcMu(j,j,1); gx12=   calcMu(j,j+1,1); gx21=   calcMu(j+1,j,1); gx22=   calcMu(j+1,j+1,1);
        gy11=ci*calcMu(j,j,2); gy12=ci*calcMu(j,j+1,2); gy21=ci*calcMu(j+1,j,2); gy22=ci*calcMu(j+1,j+1,2);
        gz11=   calcMu(j,j,3); gz12=   calcMu(j,j+1,3); gz21=   calcMu(j+1,j,3); gz22=   calcMu(j+1,j+1,3);
!          write(io,'("j=",I3)') j
!          write(io,'("gx11,gx12,gx21,gx22=",4(2X,2F8.2))') gx11,gx12,gx21,gx22
!          write(io,'("gy11,gy12,gy21,gy22=",4(2X,2F8.2))') gy11,gy12,gy21,gy22
!          write(io,'("gz11,gz21,gz12,gz22=",4(2X,2F8.2))') gz11,gz12,gz21,gz22

        t=gx11*gx11 + gx12*gx21 + gx21*gx12 + gx22*gx22; G2(1,1)=dble(t)
        if (imag(t) .gt. 0.00001) write(io,'("j=",I3," G2(",i1,",",i1,") complex=",2F8.2)') j,1,1,t
        t=gx11*gy11 + gx12*gy21 + gx21*gy12 + gx22*gy22; G2(2,1)=dble(t)
        if (imag(t) .gt. 0.00001) write(io,'("j=",I3," G2(",i1,",",i1,") complex")') j,2,1
        t=gx11*gz11 + gx12*gz21 + gx21*gz12 + gx22*gz22; G2(3,1)=dble(t)
        if (imag(t) .gt. 0.00001) write(io,'("j=",I3," G2(",i1,",",i1,") complex")') j,3,1
        t=gy11*gy11 + gy12*gy21 + gy21*gy12 + gy22*gy22; G2(2,2)=dble(t) 
        if (imag(t) .gt. 0.00001) write(io,'("j=",I3," G2(",i1,",",i1,") complex")') j,2,2
        t=gy11*gz11 + gy12*gz21 + gy21*gz12 + gy22*gz22; G2(3,2)=dble(t) 
        if (imag(t) .gt. 0.00001) write(io,'("j=",I3," G2(",i1,",",i1,") complex")') j,3,2
        t=gz11*gz11 + gz12*gz21 + gz21*gz12 + gz22*gz22; G2(3,3)=dble(t) 
        if (imag(t) .gt. 0.00001) write(io,'("j=",I3," G2(",i1,",",i1,") complex")') j,3,3
 !       G2(1,1) = gx11*gx11 + gx12*gx21 + gx21*gx21 + gx22*gx22
 !       G2(2,1) = gx11*gy11 + gx12*gy21 + gx21*gy21 + gx22*gy22
 !       G2(3,1) = gx11*gz11 + gx12*gz21 + gx21*gz21 + gx22*gz22
 !       G2(2,2) = gy11*gy11 + gy12*gy21 + gy21*gy21 + gy22*gy22
 !       G2(3,2) = gy11*gz11 + gy12*gz21 + gy21*gz21 + gy22*gz22
 !       G2(3,3) = gz11*gz11 + gz12*gz21 + gz21*gz21 + gz22*gz22
        call diagrs(G2,3,3,eng,fv1,2,io)
! If   Z.(X x Y) is 1, then a right handed system, otherwise reverse Z.
!  X x Y = (X2 Y3 - X3 Y2)i + (X3 Y1 - X1 Y3)j + (X1 Y2 - X2 Y1)k
        test = G2(1,3)*(G2(2,1)*G2(3,2)-G2(3,1)*G2(2,2))  &
             + G2(2,3)*(G2(3,1)*G2(1,2)-G2(1,1)*G2(3,2))  &
             + G2(3,3)*(G2(1,1)*G2(2,2)-G2(2,1)*G2(1,2))  
!        write(io,'("i=",I4,",  Test=",F12.8)') i,test
        if (abs(test+1.d+00).lt.0.00001) then
          do i1=1,3
            G2(3,i1) = -G2(3,i1)
          enddo
        elseif (abs(test-d1).lt.0.00001) then
        ELSE
          write(io,'("Incorrect coordinate system. Z.(X x Y)=",F8.3)') test
        endif  
        do i1=1,3
          if (eng(i1).ge.z) then
            gval(i1,i)=sqrt(d2*eng(i1))
          else
            if (abs(eng(i1)).gt.1.0d-6) write(io,'("**** WARNING: Diagonalized g2<0; eng=",F12.8)') eng(i1)
            gval(i1,i)=z
          endif          
          do j1=1,3
            g_prin(j1,i1,i)=G2(j1,i1)
          enddo
        enddo    
!        write(io,'(2i3,2F9.2,", |gx| =",F8.4,", |gy| =",F8.4,", |gz| =",F8.4)') &
!                                   j,j+1,engs(j),engs(j+1),gval(1,i),gval(2,i),gval(3,i)
      enddo ! ! i
!
      end subroutine calcGval
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcMDip()
!
!   This subroutine evaluates magnetic dipole transition moment elements between states I and J:
!    
!   <I|Ux|J>, TMD(2)=<I|Uy|J>, TMD(3)=<I|Uz|J> where Ux = RK*Lx + ge*Sx, etc.  
!   RKx is an orbital reduction parameter. 
!
      IMPLICIT NONE
      integer*4 i,j,k
      integer*4 IG,IE,IGL,IEL,IIG,IIE
      REAL*8 Rm,Rp
      real*8 Z,D1
      COMPLEX*16 CZ,CI
      parameter(D1=1.0d+00, Z=0.0d+00)
      PARAMETER(CZ=(Z,Z), CI=(Z,D1))
!
      IF (MD(2,1)-MD(1,1).GT.50 .OR. MD(2,2)-MD(1,2).GT.100) WRITE(IO,10)
 10   FORMAT(' Array RMD overflow in subroutine calcMDip.')
!
      DO IG=1,MD(2,1)-MD(1,1)+1
        IDEGG(IG)=0
        DO IE=1,MD(2,2)-MD(1,2)+1
          IDEGE(IE)=0
          DO K=1,5
            RMD(IG,IE,K)=Z
          Enddo  ! K 
        enddo  ! IE
      enddo  ! IG
!
      MDG=1
      IGL=MD(1,1)
      IIG=0
      DO IG=MD(1,1),MD(2,1)   
        IF (IMD.LE.2.AND.ABS(Engs(IG)-Engs(IGL)).GT.Edegen .OR. IMD.GE.3.AND.IG.GT.MD(1,1)) MDG=MDG+1
        RENGG(MDG)=Engs(IG)-Engs(1)
        IDEGG(MDG)=IDEGG(MDG)+1
        MDE=1
        IIG=IIG+1
        IEL=MD(1,2)
        IIE=0
        DO IE=MD(1,2),MD(2,2)
          IF (IMD.LE.2.AND.ABS(Engs(IE)-Engs(IEL)).GT.Edegen .OR. IMD.GE.3.AND.IE.GT.MD(1,2)) MDE=MDE+1
          IF (IG.EQ.MD(1,1)) THEN
            RENGE(MDE)=Engs(IE)-Engs(1)
            IDEGE(MDE)=IDEGE(MDE)+1
          ENDIF
          IIE=IIE+1
          Rp=0.5D+00*(CDABS(calcMu(IG,IE,1) - calcMu(IG,IE,2)))**2  ! R+ = |<Mx> + i.i.<My>|^2
          Rm=0.5D+00*(CDABS(calcMu(IG,IE,1) + calcMu(IG,IE,2)))**2  ! R- = |<Mx> - i.i.<My>|^2
          do i=1,3
            CMD(IIG,IIE,i)=calcMu(IG,IE,i)
            RMD(MDG,MDE,i) = RMD(MDG,MDE,i) + (CDABS(CMD(IIG,IIE,i)))**2 !  Store MX**2, MY**2, MZ**2
          enddo ! i
          RMD(MDG,MDE,4) = RMD(MDG,MDE,4) + Rp  !   M+**2
          RMD(MDG,MDE,5) = RMD(MDG,MDE,5) + Rm  !   M-**2
          IEL=IE
        enddo   ! IE excited state levels
        IGL=IG
      enddo  ! IG ground state levels
!
      end subroutine calcMDip
!
!-----------------------------------------------------------------------

END MODULE f_e_magnetics
!  359 lines