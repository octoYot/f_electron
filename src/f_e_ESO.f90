MODULE f_e_ESO
USE f_e_data
USE f_e_parameters
USE f_e_calculate ! calls diagonalise

!
!  This module calculates the energy levels of a multiplet where the crystal field 
!  has been quantified in terms of the "extended Steven Operators" (ESOs) equivalent 
!  operators.

IMPLICIT NONE

CONTAINS
!
!-----------------------------------------------------------------------
!
      Subroutine calc_ESO()
      integer*4 nbas,doVects,i,j
      real*8 D0,Eng(max_N)
      Logical mC,cV 
      character*100 cAllMJ(41)
      parameter(D0=0.0D+00)
      
      nbas=ESO2J+1
      call build_ESO_Mat()
      do i=1,nbas; do j=1,nbas; AR(i,j)=DREAL(ESOcmat(i,j)); AI(i,j)=DIMAG(ESOcmat(i,j)); VR(i,j)=D0; VI(i,j)=D0; enddo; enddo
      
      doVects=2   ! doVects=1 gives an error. (TODO: what does this mean?)
      mC=MatComplex; cV=calcVecs; MatComplex=.True.; calcVecs=.True.
      call diagonalise(nbas,Eng)
      MatComplex=mC; calcVecs=cV
!      
      call getWFs(cAllMJ)
      write(io,'(/,100("-"),/," Energies calculated using Extended Stevens Operators:")')
      do i=1,nbas
        write(io,'(F12.4,2x,"|",2x,A100)') Eng(i)-Eng(1),cAllMJ(i)
      enddo
      WRITE(IO,'(100("-"))'); 
!      
      return
      End Subroutine calc_ESO
!
!-----------------------------------------------------------------------
!
      Subroutine build_ESO_Mat()
!
!  This subroutine writes the matrix of a J multiplet where the ZFS/LF is  expressed in terms of 
!  Extended Steven Operators (ESO) operators.
!      
      INTEGER*4 m,n,nn,k,q,i,j,im,nbas,nhalf
      REAL*8 RS, D0,D1,D2,  RNkq(40,0:40) 
      real*8 Akqm(40,0:40,0:40),ESOmat(18,18)
      complex*16 C0,C1,CI
      Logical neven
      Integer(KIND=8) Fkq(14,0:14)
!      Integer(KIND=8) Fkq(12,0:12)
      data Fkq/                                                                    &
        2,4,24,48,480,  2880,40320,80640,1451520,14515200,  319334400,1916006400,49816166400,697426329600,  &
        1,2,6, 24,240,  1440,5040, 40320,725760, 7257600,   79833600, 958003200, 24908083200,348713164800,  &
        0,1,6, 8, 240,  360, 1680, 40320,725760, 1209600,   79833600, 958003200, 8302694400, 43589145600,   &
        0,0,1, 4, 10,   60,  168,  6720, 60480,  604800,    13305600, 31933440,  518918400,  7264857600,    &
        0,0,0, 1, 10,   12,  168,  672,  60480,  86400,     2661120,  3991680,   103783680,  7264857600,    &
        0,0,0, 0, 1,    6,   14,   336,  864,    2880,      23760,    1995840,   5765760,    726485760,     &
        0,0,0, 0, 0,    1,   14,   16,   288,    360,       7920,     31680,     5765760,    40360320,      &
        0,0,0, 0, 0,    0,   1,    8,    18,     180,       1320,     15840,     205920,     2882880,       &
        0,0,0, 0, 0,    0,   0,    1,    18,     20,        1320,     1584,      68640,      262080,        &         
        0,0,0, 0, 0,    0,   0,    0,    1,      10,        22,       264,       624,        43680,         &    
        0,0,0, 0, 0,    0,   0,    0,    0,      1,         22,       24,        624,        2184,          &       
        0,0,0, 0, 0,    0,   0,    0,    0,      0,         1,        12,        26,         1092,          &       
        0,0,0, 0, 0,    0,   0,    0,    0,      0,         0,        1,         26,         28,            &       
        0,0,0, 0, 0,    0,   0,    0,    0,      0,         0,        0,         1,          14,            &                         
        0,0,0, 0, 0,    0,   0,    0,    0,      0,         0,        0,         0,          1  /
!        0,1,6, 8, 240,  360, 1680, 40320,725760, 1209600,   79833600, 958003200,   &
!        0,0,1, 4, 10,   60,  168,  6720, 60480,  604800,    13305600, 31933440,    &
!        0,0,0, 1, 10,   12,  168,  672,  60480,  86400,     2661120,  3991680,     &
!        0,0,0, 0, 1,    6,   14,   336,  864,    2880,      23760,    1995840,     &
!        0,0,0, 0, 0,    1,   14,   16,   288,    360,       7920,     31680,       &
!        0,0,0, 0, 0,    0,   1,    8,    18,     180,       1320,     15840,       &
!        0,0,0, 0, 0,    0,   0,    1,    18,     20,        1320,     1584,        &         
!        0,0,0, 0, 0,    0,   0,    0,    1,      10,        22,       264,         &    
!        0,0,0, 0, 0,    0,   0,    0,    0,      1,         22,       24,          &       
!        0,0,0, 0, 0,    0,   0,    0,    0,      0,         1,        12,          &       
!        0,0,0, 0, 0,    0,   0,    0,    0,      0,         0,        1,           &       
!        0,0,0, 0, 0,    0,   0,    0,    0,      0,         0,        0,           &                         
!        0,0,0, 0, 0,    0,   0,    0,    0,      0,         0,        0           /          
                                                                                       
      parameter(D0=0.0D+00, D1=1.0d+00, D2=2.0d+00)
      parameter(C0=(D0,D0), C1=(D1,D0), CI=(D0,D1))
!
      nbas=ESO2J+1; RS=ESO2J/D2
!      
      do k=1,ESO2J
        do q=0,k
          do m=0,k
             Akqm(k,q,m)=D0
          enddo
        enddo
      enddo
!
      do k=1,ESOkmax   ! k=1,2,3,4,....
        do q=0,k ! 
          do m=0,k
            Akqm(k,q,m)=D0  ! 
          enddo 
        enddo
        RNkq(k,k)=(-1)**k*sqrt(dble(fact(2*k)))/dble(2**k*fact(k))
        Akqm(k,k,0)=D1  ! a(k,k,0)
        if (printMat4(1)) then
          write(imatrix,'(/,"  k  q  m     Akqm/Fkq      Fkq        Nkq")')
          write(imatrix,'(3i3,F16.4,I10,3x,E16.8)') k,k,0,Akqm(k,k,0)/dble(Fkq(k,k)),Fkq(k,k),RNkq(k,k)
        endif
        do q=k-1,0,-1 !  k-1, ... 1,0
          im=0
          RNkq(k,q)=(-d1)**q/(D2**k*dble(fact(k)))*sqrt(dble(fact(k+q)/dble(fact(k-q))))
!          if (printMat4(1)) write(imatrix,'(/)')
          do m=0,k-q 
            
            if (q+m.le.k) then
              if (m-1.ge.0) Akqm(k,q,m)=Akqm(k,q,m)+dble(2*q+m+1)*Akqm(k,q+1,m-1)
              Akqm(k,q,m)=Akqm(k,q,m)+dble((q+1)*q-m*(m+1)/2)*Akqm(k,q+1,m)
              nn=int(k-q-1-m)
              if (nn.gt.0) then
                do n=1,nn                    ! now use real RS*(RS+1) 
                  Akqm(k,q,m)=Akqm(k,q,m)+(-D1)**n*(comb(m+n,m)*RS*(RS+D1)-comb(m+n,m-1)-comb(m+n,m-2))*Akqm(k,q+1,m+n)
                enddo
              endif
            endif
            if (printMat4(1)) then
              if (m.eq.0) then
                write(imatrix,'(3i3,F16.4,I10,3x,E16.8)') k,q,m,Akqm(k,q,m)/DBLE(Fkq(k,q)),Fkq(k,q),RNkq(k,q)
              else
                if (Abs(Akqm(k,q,m)).gt.1.0d-10) write(imatrix,'(3i3,F16.4,I10,3x,E16.8,i4)') k,q,m,Akqm(k,q,m)/DBLE(Fkq(k,q))
              endif                                    
            endif ! printMat4(1)                                   
            
          enddo  ! m
        enddo  ! q
      enddo  ! k
       
!      call GetKfactors(RKmn2)

      do k=2,ESOkmax,2
        do q=-k,k
          if (abs(ESOBkq(k/2,k+q+1)).gt.1.0d-12) then
!            write(io,'("k,q,ESOBkq(k,q)=",2I3,4X,F14.12)') k,q,ESOBkq(k,2*k+q)
            call buildMat(Akqm,Fkq,k,q,ESOmat)
            if (printMat4(1)) then
              write(imatrix,'(" Matrix for ESO Bkq*Okq k=",I2,", q=",I2)') k,q
              write(imatrix,'(18("|",i3,"/2> "))') (-ESO2J+2*(j-1),j=1,nbas)
              do i=1,nbas;  write(imatrix,'((18F8.3))') (ESOBkq(k/2,k+q+1)*ESOmat(i,j),j=1,nbas); enddo  
            endif
            if (q.ge.0) then
              do i=1,nbas; DO j=1,nbas; ESOcmat(i,j)=ESOcmat(i,j)+ESOBkq(k/2,k+q+1)*ESOmat(i,j)*C1; enddo; enddo
            else
              do i=1,nbas; DO j=1,nbas; ESOcmat(i,j)=ESOcmat(i,j)+ESOBkq(k/2,k+q+1)*ESOmat(i,j)*CI; enddo; enddo
            endif
          endif ! abs(ESOBkq(k,q)).gt.1.0d-12
        enddo ! q
      enddo ! k 
      
      if (printMat4(2)) then  !     ! 2J+1 even; 2J odd              2J+1 odd; 2J even
        if (mod(nbas,2).eq.0) then; neven=.true.;  nhalf=nbas/2;  else ;  neven=.false.; nhalf=ESO2J/2; endif 
        write(imatrix,'(/,"(Half) Matrix in ESO basis",/," Real part of Matrix ")') 
        if (neven) then
          write(imatrix,'(8X,18(" |",i3,"/2> "))') (-ESO2J+2*(j-1),j=1,nhalf)
          do i=1,nbas;  write(imatrix,'("<",i3,"/2|",(18F9.3))') -ESO2J+2*(i-1),(Dreal(ESOcmat(i,j)),j=1,nhalf); enddo  
        else
          write(imatrix,'(8X,18("  |",i3,">  "))') (-ESO2J/2+(j-1),j=1,nhalf)
          do i=1,nbas;  write(imatrix,'("<",i3,"|",(18F9.3))') -ESO2J/2+(i-1),(Dreal(ESOcmat(i,j)),j=1,nhalf); enddo  
        endif  
        write(imatrix,'(/," Imaginary part of Matrix ")') 
        if (neven) then
          write(imatrix,'(8X,18(" |",i3,"/2> "))') (-ESO2J+2*(j-1),j=1,nhalf)
          do i=1,nbas;  write(imatrix,'("<",i3,"/2|",(18F9.3))') -ESO2J+2*(i-1),(Dimag(ESOcmat(i,j)),j=1,nhalf); enddo  
        else
          write(imatrix,'(8X,18("  |",i3,">  "))') (-ESO2J/2+(j-1),j=1,nhalf)
          do i=1,nbas;  write(imatrix,'("<",i3,"|",(18F9.3))') -ESO2J/2+(i-1),(Dimag(ESOcmat(i,j)),j=1,nhalf); enddo  
        endif  
      enDIF  ! printMat4(2)
     
!             
      end subroutine build_ESO_Mat
!
!------------------------------------------------------------------------------
      integer function comb(n,m)
      integer m,n
      comb=0
      if (m.lt.0 .or. (n-m).lt.0) return
      comb=int(fact(n)/(fact(m)*fact(n-m)))
      return
      end function comb
!------------------------------------------------------------------------------
      real*8 function fact(n)
      integer done_once/0/, i,n
      real*8 rfact(40)
!      real*8 Ri
      save done_once,rfact
! write(*,'("n=",i3,";  done_once=",I2)') n,done_once
!
      if (n.gt.40) then
        write(io,'("FUNCTION fact(n) failed for n=",I4)') n
        stop
      endif  
      if (n.lt.0) then; fact=0.0d+00; return; endif
      if (n.eq.1 .or. n.eq.0) then; fact=1.0d+00; return; endif
      if (done_once.eq.1) goto 10
      done_once=1
      rfact(1)=1.0d+00
! write(*,'("n=",i3)') n
      do i=2,40
        rfact(i)=rfact(i-1)*dble(i)
!        write(*,'(i2,"! too big for integer. integer!-real!=",E16.8)') i,dble(ifact(i))-Ri
      enddo  
 10   fact=rfact(n)
!      write(*,'("done_once=",I2,"   ",I2,"!=",F16.2)') done_once,n,fact
!
      end function fact
!-----------------------------------------------------------------------

      SUBROUTINE buildMat(Akqm,Fkq,k,q,ESOmat)
      integer*4 k,q,nbas,m,mm,TwoM,TwoMp,i,j
      integer(Kind=8) Fkq(14,0:14)
      real*8 Akqm(40,0:40,0:40),ESOmat(18,18),s1,alpha,D0,D1,D2
      parameter(D0=0.0d+00, D1=1.0D+00, D2=2.0D+00)
!                                    k even              q odd  (doesn't work for q negative!)     
      nbas=ESO2J+1; alpha=D1; if (mod(k,2).eq.0 .and. mod(abs(q),2).eq.1) alpha=0.5d+00
!      write(io,'("k,q,alpha=",2I3,4X,F6.2)') k,q,alpha
      do i=1,nbas; do j=1,nbas; ESOmat(i,j)=D0; enddo; enddo
      do i=1,nbas
        TwoMp=-ESO2J+2*(i-1)
        do j=1,nbas
          TwoM=-ESO2J+2*(j-1)
          mm=k-abs(q)
          do m=0,mm
            S1=0.0d+00
            if (q.lt.0) then
             if (TwoMp.eq.TwoM+2*abs(q)) s1=s1+sqrt(fact((ESO2J-TwoM)/2)*fact((ESO2J+TwoM)/2+abs(q))/  &
                                                   (fact((ESO2J+TwoM)/2)*fact((ESO2J-TwoM)/2-abs(q))))
             if (TwoMp.eq.TwoM-2*abs(q)) s1=s1-sqrt(fact((ESO2J+TwoM)/2)*fact((ESO2J-TwoM)/2+abs(q))/  &
                                                   (fact((ESO2J-TwoM)/2)*fact((ESO2J+TwoM)/2-abs(q))))*(-D1)**(k-abs(q)-m)
            else IF (q.eq.0) then
              if (TwoMp.eq.TwoM) s1=(d1+(-D1)**(k-abs(q)-m))*sqrt(fact((ESO2J-TwoM)/2)*fact((ESO2J+TwoM)/2+abs(q))/  &
                                                                 (fact((ESO2J+TwoM)/2)*fact((ESO2J-TwoM)/2-abs(q))))
            else IF (q.gt.0) then         
            if (TwoMp.eq.TwoM+2*abs(q)) s1=s1+sqrt(fact((ESO2J-TwoM)/2)*fact((ESO2J+TwoM)/2+abs(q))/  &
                                                  (fact((ESO2J+TwoM)/2)*fact((ESO2J-TwoM)/2-abs(q))))
            if (TwoMp.eq.TwoM-2*abs(q)) s1=s1+sqrt(fact((ESO2J+TwoM)/2)*fact((ESO2J-TwoM)/2+abs(q))/  &
                                                  (fact((ESO2J-TwoM)/2)*fact((ESO2J+TwoM)/2-abs(q))))*(-D1)**(k-abs(q)-m)
            endif            
            ESOmat(i,j)=ESOmat(i,j)+s1*(dble(TwoM)/d2)**m*Akqm(k,abs(q),m)
          enddo  ! m
          ESOmat(i,j)=ESOmat(i,j)*alpha/(D2*Fkq(k,abs(q)))
        enddo  ! j
      enddo  ! i
            
      end subroutine buildMat
!-----------------------------------------------------------------------

      SUBROUTINE GetKfactors(RKmn2)
      integer k,q
      real*8 RKmn2(20,-20:20)
      
      do k=1,ESOkmax
        do q=-k,k
          RKmn2(k,q)=1.0d+00
        enddo
      enddo

      RKmn2(2,-2)=1.5d+00
      RKmn2(2,-1)=6.0d+00
      RKmn2(2, 0)=1.0d+00
      RKmn2(2, 1)=6.0d+00
      RKmn2(2, 2)=1.5d+00
!-------------
      RKmn2(4,-4)=17.5d+00
      RKmn2(4,-3)=140.0d+00
      RKmn2(4,-2)=10.0d+00
      RKmn2(4,-1)=20.0d+00
      RKmn2(4, 0)=1.0d+00
      RKmn2(4, 1)=20.0d+00
      RKmn2(4, 2)=10.0d+00
      RKmn2(4, 3)=140.0d+00
      RKmn2(4, 4)=17.5d+00
!---------2----
      RKmn2(6,-6)=57.75d+00
      RKmn2(6,-5)=693.0d+00
      RKmn2(6,-4)=31.5d+00
      RKmn2(6,-3)=105.0d+00
      RKmn2(6,-2)=26.25d+00
      RKmn2(6,-1)=42.0d+00
      RKmn2(6, 0)=1.0d+00
      RKmn2(6, 1)=42.0d+00
      RKmn2(6, 2)=26.25d+00
      RKmn2(6, 3)=105.0d+00
      RKmn2(6, 4)=31.5d+00
      RKmn2(6, 5)=693.0d+00
      RKmn2(6, 6)=57.75d+00
!---------2----
      RKmn2(8,-8)=3217.5d+00
      RKmn2(8,-7)=51480.0d+00
      RKmn2(8,-6)=1716.0d+00
      RKmn2(8,-5)=72072.0d+00
      RKmn2(8,-4)=1386.0d+00
      RKmn2(8,-3)=9240.0d+00
      RKmn2(8,-2)=1260.0d+00
      RKmn2(8,-1)=72.0d+00
      RKmn2(8, 0)=1.0d+00
      RKmn2(8, 1)=72.0d+00
      RKmn2(8, 2)=1260.0d+00
      RKmn2(8, 3)=9240.0d+00
      RKmn2(8, 4)=1386.0d+00
      RKmn2(8, 5)=72072.0d+00
      RKmn2(8, 6)=1716.0d+00
      RKmn2(8, 7)=51480.0d+00
      RKmn2(8, 8)=3217.5d+00
!---------2--
      RKmn2(10,-10)=11547.25d+00
      RKmn2(10, -9)=230945.0d+00
      RKmn2(10, -8)=6077.5d+00
      RKmn2(10, -7)=36465.0d+00
      RKmn2(10, -6)=536.25d+00
      RKmn2(10, -5)=1716.0d+00
      RKmn2(10, -4)=4290.0d+00
      RKmn2(10, -3)=8580.0d+00
      RKmn2(10, -2)=82.5d+00
      RKmn2(10, -1)=110.0d+00
      RKmn2(10,  0)=1.0d+00
      RKmn2(10,  1)=110.0d+00
      RKmn2(10,  2)=82.5d+00
      RKmn2(10,  3)=8580.0d+00
      RKmn2(10,  4)=4290.0d+00
      RKmn2(10,  5)=1716.0d+00
      RKmn2(10,  6)=536.25d+00
      RKmn2(10,  7)=36465.0d+00
      RKmn2(10,  8)=6077.5d+00
      RKmn2(10,  9)=230945.0d+00
      RKmn2(10, 10)=11547.25d+00
!---------2--   
      RKmn2(12,-12)=169009.75d+00
      RKmn2(12,-11)=4056234.0d+00
      RKmn2(12,-10)=88179.0d+00
      RKmn2(12, -9)=646646.0d+00
      RKmn2(12, -8)=69283.5d+00
      RKmn2(12, -7)=277134.0d+00
      RKmn2(12, -6)=2431.0d+00
      RKmn2(12, -5)=306306.0d+00
      RKmn2(12, -4)=2252.25d+00
      RKmn2(12, -3)=4004.0d+00
      RKmn2(12, -2)=6006.0d+00
      RKmn2(12, -1)=156.0d+00
      RKmn2(12,  0)=1.0d+00
      RKmn2(12,  1)=156.0d+00
      RKmn2(12,  2)=6006.0d+00
      RKmn2(12,  3)=4004.0d+00
      RKmn2(12,  4)=2252.25d+00
      RKmn2(12,  5)=306306.0d+00
      RKmn2(12,  6)=2431.0d+00
      RKmn2(12,  7)=277134.0d+00
      RKmn2(12,  8)=69283.5d+00
      RKmn2(12,  9)=646646.0d+00
      RKmn2(12, 10)=88179.0d+00
      RKmn2(12, 11)=4056234.0d+00
      RKmn2(12, 12)=169009.75d+0
      
      end subroutine GetKfactors
!-----------------------------------------------------------------------
!
      character*100 function makeESOMJLabel(pMJ,cMJ,nMJ)
!
!  Returns the MJ projections as a single label.
!  Order so the components are listed in descending order.
!
      IMPLICIT none
      integer i,j,nMJ,inMJ
      integer pMJ(41),p,sumMJ 
      character cMJ(41)*7, c*7, lab1*4, label*100
! Check% add up
      sumMJ=0
      do i=1,ESO2J+1
        sumMJ=sumMJ+pMJ(i)
      enddo
      if (abs(sumMJ-100).gt.25) write(io,'("***WARNING: Sum in makeMJLabel far from 100%; sumMJ=",I3)') sumMJ
!  Make them descending order:
 10   continue
      do i=1,ESO2J+1
        do j=i+1,ESO2J+1
          if (pMJ(j).gt.pMJ(i)) then
            p=pMJ(j)
            pMJ(j)=pMJ(i)
            pMJ(i)=p
            c=cMJ(j)
            cMJ(j)=cMJ(i)
            cMJ(i)=c
            goto 10
          endif
        enddo ! j
      enddo ! i
                  
      write(lab1,'(I3)') pMJ(1)
      label=trim(lab1)//"*"
      do i=1,len(trim(cMJ(1)))
        if (cMJ(1)(i:i).ne." ") label=trim(label)//cMJ(1)(i:i)
      enddo
      if (nMJ.gt.1) then
        do inMJ=2,nMJ
          write(lab1,'(" +",I2)') pMJ(inMJ)
          label=trim(label)//trim(lab1)//"*"
          do i=1,len(trim(cMJ(inMJ)))
            if (cMJ(inMJ)(i:i).ne." ") label=trim(label)//cMJ(inMJ)(i:i)
          enddo
        enddo
      endif
!      if (len(label).eq.len(TRIM(label))) write(IO,'("***Warning label in makeMJLabel is full")') 
!      write(IO,'("***Label in makeMJLabel:",A80)') label 
      makeESOMJLabel = label
      RETURN 
      END function makeESOMJLabel
!      
!-----------------------------------------------------------------------
!
      SUBROUTINE getWFs(cAllMJ)
!
!  Calculates: The MJ projections of each wavefunction in cMJ. 
!
      IMPLICIT none
      integer i,ii,is,nMJ,pMJ(41),TwoMJ
      REAL*8 wf2,limit
      character*7 label,allLabels(41),cMJ(41)
      character*100 cAllMJ(41)
      PARAMETER(limit=0.02d+00)
!
      DO is=1,ESO2J+1
        TwoMJ = -ESO2J + (is-1)*2
        if (mod(TwoMJ,2).eq.0) then
          write(label,'("|",I3,">")') TwoMJ/2   ! TwoMJ even 
        else
          write(label,'("|",I3,"/2>")') TwoMJ   ! TwoMJ odd
        endif
        allLabels(is)=label
      enddo
!      write(IO,'("allLabels(k)=",40A9)') (allLabels(i),i=1,ESO2J+1)

!  Loop for every eigenvector:
      DO II=1,ESO2J+1
        DO i=1,ESO2J+1       !
          pMJ(i)=0; cMJ(i)="      "
        enddo
        nMJ=0
        do i=1,ESO2J+1
          wf2=VR(i,ii)**2+VI(i,ii)**2
          if (wf2.gt.limit) then
            nMJ=nMJ+1
            pMJ(nMJ)=INT(100*wf2)
            cMJ(nMJ)=allLabels(i)
          endif  
        enddo  ! i
        cAllMJ(ii)=makeESOMJLabel(pMJ,cMJ,nMJ)
!
      enddo  !  II  loop over the neng eigenvectors
!
      RETURN
      END Subroutine getWFs
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcESOfit(m,fvec)
!
!  Calculate the square of the differences between Calc MEs and ESO MEs in the vector fvec.
!
      IMPLICIT none
      integer nbas,nhalf,jj,i,j,k,m
      logical neven
      real*8 fvec(m)
!
      nbas=ESO2J+1
      if (mod(nbas,2).eq.0) then   ! 2J+1 even; 2J odd
        neven=.true.
        nhalf=nbas/2
      else                         ! 2J+1 odd; 2J even
        neven=.false.
        nhalf=ESO2J/2
      endif       
!
      k=0
      if (mod(nbas,2).eq.0) then  ! 2J+1 even; 2J odd
        do i=1,nbas/2
          do j=1,i
            k=k+1
            fvec(k)=(AR(i,j)-Dreal(ESOcmat(i,j)))**2
            if (j.ne.i) then
              k=k+1
              fvec(k)=(AI(i,j)-Dimag(ESOcmat(i,j)))**2
            endif
          enddo
        enddo
        do i=nbas/2+1,nbas
          do j=1,nbas/2-(i-nbas/2-1)
            k=k+1
            fvec(k)=(AR(i,j)-Dreal(ESOcmat(i,j)))**2
            k=k+1
            fvec(k)=(AI(i,j)-Dimag(ESOcmat(i,j)))**2
          enddo
        enddo
      else  ! 2J+1 odd; 2J even
        jj=ESO2J/2
        do i=1,jj+1
          do j=1,i
            k=k+1
            fvec(k)=(AR(i,j)-Dreal(ESOcmat(i,j)))**2
            if (j.ne.i) then
              k=k+1
              fvec(k)=(AI(i,j)-Dimag(ESOcmat(i,j)))**2
            endif
          enddo
        enddo
        do i=jj+2,nbas
          do j=1,j+1-(i-(j+1))
            k=k+1
            fvec(k)=(AR(i,j)-Dreal(ESOcmat(i,j)))**2
            k=k+1
            fvec(k)=(AI(i,j)-Dimag(ESOcmat(i,j)))**2
          enddo
        enddo
      endif
      write(imatrix,'("calcESOfit: P(21)=",f18.12,";  P(41)=",f18.12)') P(21),P(41)
      write(imatrix,'("calcESOfit: AR(1,1)=",f18.12,";  Dreal(ESOcmat(1,1))=",f18.12)') AR(1,1),Dreal(ESOcmat(1,1))    
      write(imatrix,'("calcESOfit: AR(2,2)=",f18.12,";  Dreal(ESOcmat(2,2))=",f18.12)') AR(2,2),Dreal(ESOcmat(2,2))        
      write(imatrix,'("calcESOfit: fvec:  k=",i3)')  k      
      write(imatrix,'((20E10.2))') (fvec(i),i=1,k)
      if (k.ne.m .or. k.ne. nexp) then
        write(io,'("****FATAL: problem loading ESO MEs in calcESOfit; k=",i4,"; m=",i4,"; nexp=",i4)') k,m,nexp
        stop
      endif
!
      RETURN
      END Subroutine calcESOfit
!
!-----------------------------------------------------------------------
END MODULE f_e_ESO
! 540