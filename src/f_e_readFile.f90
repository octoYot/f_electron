MODULE f_e_readFile

USE f_e_data
USE f_e_parameters
!
!  This module reads the input file.
!

IMPLICIT NONE

! PRIVATE

PUBLIC 

CONTAINS
!
!
!-----------------------------------------------------------------------
!
      INTEGER FUNCTION TS_LookupC(N,Lab)
 
!  Returns the index of the term symbol for the d^N or f^N system in the standard order
!  used by Nielson & Koster, given the three character label "Lab" for the term symbol.
 
      IMPLICIT none
      integer n,i 
      character*3 Lab 
!
      do 50 i=1,TS_n(N)
        if (Lab.eq.TS_labels(i,N)) then
          TS_LookupC=i
          return
        endif
 50   continue
      
      write(*,'("***FATAL: Index for ",A1,"^",I1," not found for label= ",A3)') e_lab,N,Lab 
      stop
      end function TS_LookupC
!
!-----------------------------------------------------------------------
!
      real*8 FUNCTION decode1(Lab)
 
!  Returns the real number coded as a product of powers of prime numbers 
!  in the format of  Neilson & Koster. 
!      A = r* Sqrt( Prod(i=1,12) p(i)^n(i) )
!  ie the square root of the product of the first 12 primes (2,3,5...37) raised
!  to some integer power. The remainder r can contain higher primes.
!  The smallest and biggest powers that are encountered for f^n system 
!  is -11(-B) and 29( T) respectively.
!  f^5  <4K1| S4 |2M2>  ( P)
!  f^6  <4F3| S4 |2P4>  ( T)
!  
      IMPLICIT none
      integer i,rem,n(12)
      real*8 p(12),x
      character  Lab*40, ln(12)*2
      data p/ 2.0d+00, 3.0d+00, 5.0d+00, 7.0d+00,11.0d+00,13.0d+00, &
             17.0d+00,19.0d+00,23.0d+00,29.0d+00,31.0d+00,37.0d+00/
      read(lab,'(I10,3(2x,4A2))') rem,(ln(i),i=1,12)
      x=1.0d+00
      do 10 i=1,12
        if (ln(i).eq." A") then 
          n(i)=10
        else if (ln(i).eq."-A") then
          n(i)=-10
        else if (ln(i).eq." B") then
          n(i)=11
        else if (ln(i).eq."-B") then
          n(i)=-11
        else if (ln(i).eq." C") then
          n(i)=12
        else if (ln(i).eq." D") then
          n(i)=13
        else if (ln(i).eq." E") then
          n(i)=14
        else if (ln(i).eq." F") then
          n(i)=15
        else if (ln(i).eq." G") then
          n(i)=16
        else if (ln(i).eq." H") then
          n(i)=17
        else if (ln(i).eq." I") then
          n(i)=18
        else if (ln(i).eq." J") then
          n(i)=19
        else if (ln(i).eq." K") then
          n(i)=20
        else if (ln(i).eq." L") then
          n(i)=21
        else if (ln(i).eq." M") then
          n(i)=22
        else if (ln(i).eq." N") then
          n(i)=23
        else if (ln(i).eq." O") then
          n(i)=24
        else if (ln(i).eq." P") then
          n(i)=25
        else if (ln(i).eq." T") then
          n(i)=29
        else
          read(ln(i),'(i2)') n(i)
          if (n(i).lt.-9 .or. n(i).gt.9) write(io,'("*** Invalid value found in FUNCTION decode1: ln, n=",A2,I3)') ln(i),n(i) 
        endif
        if (n(i).ne.0) x=x*(p(i)**n(i))
 10   continue
      decode1=rem*sqrt(x)
!      write(idebug,'("LAB=",A40,"decode1=",F10.8)') LAB,decode1 
!      write(idebug,'("    ",I6,3(2x,4I2))') REM,(n(i),i=1,12)
      return
      end function decode1
!
!-----------------------------------------------------------------------
!
      SUBROUTINE decode2(ii,jj,Lab,diag,EEpar)
 
!  Returns the real number values for F2, F4, F6 electron repulsion parameters 
!  coded in a character string. The value of F0 is ignored.
!
!  Note the string is scanned for +/- which is assumed to be infront of either 
!  a number (which also contains the symbols F2,F4,F6), the symbol F2,F4,F6 itself,
!  or a "(" which is the start of a sqrt sign. 
!  
      IMPLICIT none
      logical diag
      integer i,t1,t2,t3,n1,n2,n3,n4,n,  ii,jj
      real*8 EEpar(3),F2,F4,F6
      character  Lab*80,LAB1*25,LAB2*25,LAB3*25,num*15,dem*15
      do  i=1,3
        EEpar(i)=0.0d+00
      enddo
      If (Nelectrons.eq.1) return
      n=len(lab)
! n1 start of F0      
      t1=scan(lab,"+"); t2=scan(lab,"-")
      n1=0
      if (t1.lt.t2 .and. t1.ne.0) n1=t1 
      if (t2.lt.t1 .and. t2.ne.0) n1=t2 
      if (t1.eq.0) n1=t2
      if (t2.eq.0) n1=t1
! n2 start of F2     
      t1=n1+scan(lab(n1+1:n),"+"); t2=n1+scan(lab(n1+1:n),"-")
      n2=0      
      if (t1.lt.t2 .and. t1.ne.n1) n2=t1 
      if (t2.lt.t1 .and. t2.ne.n1) n2=t2
      if (t1.eq.n1) n2=t2
      if (t2.eq.n1) n2=t1
!  n3 start of F4      
      t1=n2+scan(lab(n2+1:n),"+"); t2=n2+scan(lab(n2+1:n),"-") 
      n3=0      
      if (t1.lt.t2 .and. t1.ne.n2) n3=t1 
      if (t2.lt.t1 .and. t2.ne.n2) n3=t2 
      if (t1.eq.n2) n3=t2
      if (t2.eq.n2) n3=t1
!  n4 start of F6      
      t1=n3+scan(lab(n3+1:n),"+"); t2=n3+scan(lab(n3+1:n),"-")
      n4=0
      if (t1.lt.t2 .and. t1.ne.n3) n4=t1 
      if (t2.lt.t1 .and. t2.ne.n3) n4=t2 
      if (t1.eq.n3) n4=t2
      if (t2.eq.n3) n4=t1
!      write(idebug,'("n1,n2,n3,n4=",4i4,";  diag=",L)') n1,n2,n3,n4,diag 
      if (diag) then  ! ignoring F0 on diagonal
        n1=n2
        n2=n3
        n3=n4
        n4=0
      endif
      LAB1=LAB(n1:n2-1);  
      if (Lvalue.eq.2) then
        LAB2=LAB(n2:n)       ! only F2, F4
        LAB3=""
      else
        LAB2=LAB(n2:n3-1)    ! F2, F4 
        LAB3=LAB(n3:n)       ! F6
      endif  
!      write(idebug,'("LAB=",A80,"L1,L2,L3",3(2X,A25),"n1,n2,n3,n4=",4i4)') LAB,LAB1,LAB2,LAB3,n1,n2,n3,n4 
      n1=scan(lab1,"F"); num=lab1(1:n1-1); dem=lab1(n1+3:25)   
      call decode3(lab1,num,dem,n1,t1,t2,t3)
      if (t2.gt.0) then
        EEpar(1)=dble(t1)*sqrt(dble(t2))/dble(t3)
       else
        write(IO,'("***FATAL: sqrt of a negative number in decode2")'); stop
      endif
!      write(idebug,'("LAB1=",A25,4x,i15,"(",i15,")F2/",i15,"  F2=",F12.6)') LAB1,t1,t2,t3,EEpar(1) 
      n1=scan(lab2,"F"); num=lab2(1:n1-1); dem=lab2(n1+3:25)
      call decode3(lab2,num,dem,n1,t1,t2,t3)
      if (t2.gt.0) then
        EEpar(2)=dble(t1)*sqrt(dble(t2))/dble(t3)
       else
        write(IO,'("***FATAL: sqrt of a negative number in decode2")'); stop
      endif
!      write(idebug,'("LAB2=",A25,4x,i15,"(",i15,")F4/",i15,"  F4=",F12.6)') LAB2,t1,t2,t3,EEpar(2) 
      n1=scan(lab3,"F")  
      if (n1.ne.0) then  ! F6 missing for d electrons
        num=lab3(1:n1-1); dem=lab3(n1+3:25)
        call decode3(lab3,num,dem,n1,t1,t2,t3)
        if (t2.gt.0) then
          EEpar(3)=dble(t1)*sqrt(dble(t2))/dble(t3)
         else
          write(IO,'("***FATAL: sqrt of a negative number in decode2")'); stop
        endif
!        write(idebug,'("LAB3=",A25,4x,i15,"(",i15,")F6/",i15,"  F6=",F12.6)') LAB3,t1,t2,t3,EEpar(3) 
      endif

      return
      end SUBROUTINE decode2
!
!-----------------------------------------------------------------------
!
      SUBROUTINE decode3(lab1,num,dem,n1,t1,t2,t3)
! Has to deal with:
!   num being strings like: "+", "-5", "-54(233)", "-(56)"
!   dem being strings like: "34"
      IMPLICIT NONE
      integer i,n1,n2,t1,t2,t3
      character LAB1*25,num*15,dem*15,sr*15 
!      
      if (n1.eq.2 .and. num(1:1).eq."-") num="-1"
      if (n1.eq.2 .and. num(1:1).eq."+") num="+1"
      n1=scan(lab1,"("); n2=scan(lab1,")")
!      write(idebug,'("num=",A15," dem=",A15,"n1,n2=",2i4)') num,dem,n1,n2 
      if (n1.ne.0 .and. n2.ne.0) then   !  sqrt ( )
        num=lab1(1:n1-1); sr=lab1(n1+1:n2-1)
        if (n1.eq.2 .and. num(1:1).eq."-") THEN
          num="-1"
        else if (n1.eq.2 .and. num(1:1).eq."+") then
          num="+1"
        endif
        read(num,'(I15)') t1
        READ(sr, '(I15)') t2
      else
        read(num,'(I15)') t1
        t2=1
      endif
      read(dem,'(I15)') t3      
      return
      end SUBROUTINE decode3
!
!-----------------------------------------------------------------------
!
      SUBROUTINE readUMat(nline)
!
!  Reads the Uk, k=1,6 (for f electrons) reduced matrix elements. 
! 
      IMPLICIT none
      integer i,j,k, nelem,nline, Lij
      character lab1*3, lab2*3, Lab*40
!----------
!  Positions the data file containing the matrix elements to the start of the Nelectrons position.
!        write(idebug,'("Inside readUMat")')  
      nelem=0; k=0
 5    READ(Idata,'(1X,A3,1X,A3)') LAB1,LAB2
!      write(idebug,'(A3,"|",A3)') LAB1,LAB2 
      nline=nline+1
      if (LAB1.eq."***") goto 5
!      write(idebug,'("Nelectrons=",I4)') Nelectrons 
      if (Nelectrons.eq.1) goto 35
      DO 30 I=1,Nelectrons-1
 10     nline=nline+1
        READ(Idata,'(1X,A3,1X,A3)') LAB1,LAB2
!        write(idebug,'(A3,"|",A3,"|",A40)') LAB1,LAB2,LAB
        IF (LAB1.eq."   ".and.LAB2.eq."   ") GOTO 30
        GOTO 10
 30   CONTINUE
!----------
!
 35   CONTINUE
      if (check1(2)) then
        write(idebug,'(80("-"),/,"Uk(i,j,k) reduced matrix elements, positioned at line:",I5)') nline
      endif
 50   READ(Idata,'(1X,A3,1X,A3,1X,A40)',end=100) LAB1,LAB2,LAB
!      write(idebug,'(2A3,A40)') LAB1,LAB2,LAB
      nline=nline+1
      if (LAB1.eq."V11") then  ! f^n
        if (check1(2)) write(idebug,'(I5," Uk reduced matrix elements read.")') nelem
        return
      else if (LAB1(1:1).eq."U") then
        k=ichar(LAB1(2:2))-48
        if (k.lt.1 .or. k.gt.6) Then
          write(io,'("***FATAL: Invalid UMat index=",I4,". Program stopped in readUMat at line:",  &
          i4,/," line:",1X,A3,1X,A3,3X,A40)') k,nline,LAB1,LAB2,LAB
          STOP
        endif
      else
        nelem=nelem+1
        if (LAB1.ne."   ") i=TS_LookupC(Nelectrons,LAB1)
        j=TS_LookupC(Nelectrons,LAB2)
        Lij=(TS_Bases(2,i,Nelectrons)-TS_Bases(2,j,Nelectrons))    !  (L - L') 
        UMAT(i,j,k)=decode1(Lab)
        UMAT(j,i,k)=DBLE((-1)**Lij)*UMAT(i,j,k)
        if (check1(2)) write(idebug,'(1X,A3,1X,A3,3I4,3x,F12.8)')LAB1,LAB2,i,j,k,UMAT(i,j,k)
      endif
      goto 50
 100  write(*,'("**** End of the data file reached in readUMat")')   
      return
      end subroutine readUMat
!
!-----------------------------------------------------------------------
!
      SUBROUTINE readVMat(nline)
!
!  Reads the V11 reduced matrix elements used to calculate the spin-orbit coupling.
!  Lvalue=2 for d-electrons and 3 for f-electrons.
! 
      IMPLICIT none
      integer nline
      integer i,j, nelem, Sij,Lij
      character lab1*3, lab2*3, Lab*40
!----------
!  Should already be positioned at the "ZV11..." line in the data file.
      nelem=0
      if (check1(3)) then
        write(idebug,'(80("-"),/,"V11 reduced matrix elements, positioned at line:",I5)') nline
      endif
 50   READ(idata,'(1X,A3,1X,A3,1X,A40)',end=100) LAB1,LAB2,LAB
      nline=nline+1
!      write(idebug,'(2A3,A40)') LAB1,LAB2,LAB
      if (LAB1.eq."ELE".or.LAB1.eq."M0 ") then  !  f^1 or f^(n>1)
        if (check1(3)) write(idebug,'(I5," V11 reduced matrix elements read.")') nelem
        return
      else
        nelem=nelem+1
        if (LAB1.ne."   ") i=TS_LookupC(Nelectrons,LAB1)
        j=TS_LookupC(Nelectrons,LAB2)
        Sij=(TS_Bases(1,i,Nelectrons)-TS_Bases(1,j,Nelectrons))/2  !  (S - S')  always integer 
        Lij=(TS_Bases(2,i,Nelectrons)-TS_Bases(2,j,Nelectrons))    !  (L - L')  always integer 
        VMAT(i,j)=decode1(Lab)
        VMAT(j,i)=DBLE((-1)**(Lij-Sij))*VMAT(i,j)
        if (check1(3)) write(idebug,'(1X,A3,1X,A3,2I4,3X,F12.8)')LAB1,LAB2,i,j,VMAT(i,j)
      endif
      goto 50
 100  write(*,'("**** End of the data file reached in readVMat")')   
      return
      end subroutine readVMat
!
!-----------------------------------------------------------------------
!
      SUBROUTINE readEEMat(nline)
!
!  Reads the matrices of the Slater F^k electron-repulsion reduced matrix elements.
! 
      IMPLICIT none
      logical diag
      integer nline,nelem,i,j,k,totstates
      REAL*8 EEpar(3), sumX(3)
      character lab1*3, lab2*3, Lab*80
!----------
!  Should already be positioned at the "ZELECTROSTATIC MATRIX FOR..." line in the data file.
      if (Nelectrons.lt.2) return
      nelem=0
      do k=1,3
        sumX(k)=0.0d+00
      enddo
      if (check1(4)) then
        write(idebug,'(80("-"),/,"E-E repulsion matrix elements, positioned at line:",I5)') nline
        write(idebug,'(25X,"F2          F4          F6          e3")')
      endif
 50   READ(idata,'(1X,A3,1X,A3,2X,A80)',end=100) LAB1,LAB2,LAB
      nline=nline+1
!      write(idebug,'(2A3,A80)') LAB1,LAB2,LAB
      if (LAB1.eq."Z  ") then   !  f^n
        if (check1(4)) write(idebug,'(I5," E-E matrix elements read.")') nelem
        totstates=0
        DO i=1,TS_n(Nelectrons)
          totstates=totstates+(TS_Bases(1,i,Nelectrons)+1.0d+00)*(2.0d+00*TS_Bases(2,i,Nelectrons)+1.0d+00)
        enddo
        DO i=1,TS_n(Nelectrons)
          do k=1,3
            EEMAT(i,i,k)=EEMAT(i,i,k)-sumX(k)/dble(totstates)
          enddo
        enddo
        if (check1(4)) then
          write(idebug,'("Subtracting mean values from diagonal:",3F12.6," (/",I4,")")') (sumX(k),k=1,3),totstates
          DO i=1,TS_n(Nelectrons)
            LAB1=TS_labels(i,Nelectrons)   
            write(idebug,'(1X,A3,1X,A3,2I4,3X,3F12.6)')LAB1,LAB1,i,i,(EEMAT(i,i,k),k=1,3)
          enddo
        endif
        return
      else
        nelem=nelem+1
        if (LAB1.ne."   ") i=TS_LookupC(Nelectrons,LAB1)
        j=TS_LookupC(Nelectrons,LAB2)
        diag = i.eq.j
        call decode2(i,j,LAB,diag,EEpar)
        do k=1,3
          EEMAT(i,j,k)=EEpar(k)
          EEMAT(j,i,k)=EEpar(k)
        enddo
        EEMAT(i,j,4)=825.0d+00/14.0d+00*EEpar(1)+396.0d+00/7.0d+00*EEpar(2)-5577.0d+00/50.0d+00*EEpar(3)
        EEMAT(j,i,4)=825.0d+00/14.0d+00*EEpar(1)+396.0d+00/7.0d+00*EEpar(2)-5577.0d+00/50.0d+00*EEpar(3)
        if (i.eq.j) then
          do k=1,3
            sumX(k)=sumX(k)+EEMAT(i,j,k)*(TS_Bases(1,i,Nelectrons)+1.0d+00)*(2.0d+00*TS_Bases(2,i,Nelectrons)+1.0d+00)
          enddo
        endif 
        if (check1(4)) write(idebug,'(1X,A3,1X,A3,2I4,3X,4F12.6)')LAB1,LAB2,i,j,(EEMAT(i,j,k),k=1,4)
      endif
      goto 50
 100  write(*,'("**** End of the data file reached in readEEMat")')   
      return
      end subroutine readEEMat
!
!-----------------------------------------------------------------------
!
      SUBROUTINE readABGMat()
!
!  "Reads" the diagonal matrices for the orbit-orbit CI 2e- matrix elements
!  parameterised by ALPHA, BETA, GAMMA.  Also called Electrostatic CI parameters. 
!  Not read from file, generated from the QNs given by L, G2 and R7 respectively.
! 
      IMPLICIT none
      integer i,j,k, Li,totstates 
      real*8 sumX(3)
!
      if (Nelectrons.lt.2) return  
      totstates=0
      do j=1,3
        sumX(j)=0.0d+00
        DO i=1,TS_n(Nelectrons)
          ABGMat(i,j)=0.0d+00
        Enddo  
      enddo
      
      if (Lvalue.eq.2) then  ! trees correction only

        DO i=1,TS_n(Nelectrons)
          Li=TS_Bases(2,i,Nelectrons)
          ABGMat(i,1)=DBLE(Li)*(DBLE(Li+1))
        enddo
        goto 10
      
      else if (Lvalue.eq.3) then  ! alpha, beta, gamma only for f electrons


        if (Nelectrons.eq.1.or.Nelectrons.eq.13) goto 10
        DO i=1,TS_n(Nelectrons)
          Li=TS_Bases(2,i,Nelectrons)
          ABGMat(i,1)=DBLE(Li)*(DBLE(Li+1))
          ABGMat(i,2)=G2R7lookup(TS_f_QNs(2,i,Nelectrons),2)
          ABGMat(i,3)=G2R7lookup(TS_f_QNs(1,i,Nelectrons),3)
          totstates=totstates+(TS_Bases(1,i,Nelectrons)+1.0d+00)*(2.0d+00*TS_Bases(2,i,Nelectrons)+1.0d+00)
          do j=1,3
            sumX(j)=sumX(j)+ABGMat(i,j)*(TS_Bases(1,i,Nelectrons)+1.0d+00)*(2.0d+00*TS_Bases(2,i,Nelectrons)+1.0d+00)
          enddo
        enddo
        DO i=1,TS_n(Nelectrons)
          do j=1,3
            ABGMat(i,j)=ABGMat(i,j)-sumX(j)/dble(totstates)
          enddo
        enddo
        
      endif  

 10   if (check1(5)) then
        write(idebug,'(80("-"),/,"alpha RMEs: (Centre of gravity:",F10.4," subtracted)")') sumX(1)/dble(totstates)
        write(idebug,'(10F10.6)') (ABGMat(i,1),i=1,TS_n(Nelectrons))
        write(idebug,'("beta  RMEs: (Centre of gravity:",F10.4," subtracted)")') sumX(2)/dble(totstates)
        write(idebug,'(10F10.6)') (ABGMat(i,2),i=1,TS_n(Nelectrons))
        write(idebug,'("gamma RMEs: (Centre of gravity:",F10.4," subtracted)")') sumX(3)/dble(totstates)
        write(idebug,'(10F10.6)') (ABGMat(i,3),i=1,TS_n(Nelectrons))
      endif
      return
      end subroutine readABGMat
!
!-----------------------------------------------------------------------
!
      SUBROUTINE readMnPnSnMat(nline)
!
!  Reads the matrices for the spin-other orbit Mn matrix elements
!  and the electrostatically correlated spin-orbit Pn matrix elements.
!  and the spin-spin matrix elements (uses Mn parameters too).
      IMPLICIT none
      integer i,j,k,ii,jj,kk, nelem,nline, Sij,Lij,mp
      real*8 sum1
      character lab1*3, lab2*3, Lab*40
      logical OKtoStart
!----------
      nelem=0; OKtoStart=.true.; mp=1; k=1
      if (Nelectrons.lt.2) return   
      if (Lvalue.ne.3) return   !  Only for f electrons
      if (check1(6)) then
        write(idebug,'(80("-"),/,"Mn,Pn,Sn reduced matrix elements, positioned at line:",I5)') nline
      endif
 50   READ(Idata,'(1X,A3,1X,A3,1X,A40)',end=100) LAB1,LAB2,LAB
!      write(idebug,'(A3,"|",A3,"|",A40)') LAB1,LAB2,LAB
      nline=nline+1
      if (.not.OKtoStart .and. LAB1.eq."   ") then
        write(io,'(" No Mk,Pk,Sk reduced matrix elements in data file")')
        return
      endif
      if (LAB1.eq."ELE" .or. LAB1.eq."T2 ") then   !  f^2 or f^(n>2)
        if (check1(6)) write(idebug,'(I5," Mk,Pk,Sk reduced matrix elements read.")') nelem
        if (check1(6)) then
          do ii=1,TSL_f(Nelectrons)
            do jj=ii,TSL_f(Nelectrons)
              sum1=0.0d+00
              do kk=1,3
                sum1=sum1+MnMAT(ii,jj,kk)+PnMAT(ii,jj,kk)+SnMAT(ii,jj,kk)
              enddo 
              if (sum1.ne.0.0d+00) then
                write(idebug,'(1X,A3,1X,A3,3(3x,3F12.7))') TS_labels(ii,Nelectrons),TS_labels(jj,Nelectrons),   &
                 (MnMAT(ii,jj,kk),kk=1,3),(PnMAT(ii,jj,kk),kk=1,3),(SnMAT(ii,jj,kk),kk=1,3)
              endif
            enddo
          enddo
        endif      
        return
      else if (LAB1(1:1).eq."M".or.LAB1(1:1).eq."P".or.LAB1(1:1).eq."S") then
        OKtoStart=.true.
        if (LAB1(1:1).eq."M") then
          mp=1
          k=(ichar(LAB1(2:2))-48)/2+1  ! M0,M2,M4: k=1,2,3
 !         write(idebug,'("M, LAB1=",A3,", k=",I2)') LAB1,k
        endif  
        if (LAB1(1:1).eq."P") then
          mp=2
          k=(ichar(LAB1(2:2))-48)/2    ! P2,P4,P6: k=1,2,3
!          write(idebug,'("P, LAB1=",A3,", k=",I2)') LAB1,k
        endif
        if (LAB1(1:1).eq."S") then
          mp=3
          k=(ichar(LAB1(2:2))-48)/2+1   ! S0,S2,S4: k=1,2,3
!          write(idebug,'("S, LAB1=",A3,", k=",I2,", lab=",A40)') LAB1,k,lab
        endif
        if (k.lt.1 .or. k.gt.3) Then
          write(io,'("***FATAL: Invalid MnPnSnMat index=",I4,". Program stopped in readMnPnSnMat at line:",  &
          i4,/," line:",1X,A3,1X,A3,3X,A40)') k,nline,LAB1,LAB2,LAB
          STOP
        endif
      else
        if (OKtoStart) then
          nelem=nelem+1
          if (LAB1.ne."   ") i=TS_LookupC(Nelectrons,LAB1)
          j=TS_LookupC(Nelectrons,LAB2)
          Sij=(TS_f_Bases(1,i,Nelectrons)-TS_f_Bases(1,j,Nelectrons))/2  !  (S - S')  always integer
          Lij=(TS_f_Bases(2,i,Nelectrons)-TS_f_Bases(2,j,Nelectrons))    !  (L - L')  always integer
          if (mp.eq.1) then
            MnMAT(i,j,k)=decode1(Lab)
            MnMAT(j,i,k)=DBLE((-1)**(Lij-Sij))*MnMAT(i,j,k)  !  Eq. after 4 Judd, Phys.Rev.169,130(1968).
!            if (check1(6)) write(idebug,'(1X,A3,1X,A3,3I4,3x,F12.8)')LAB1,LAB2,i,j,k,MnMAT(i,j,k)
          else if (mp.eq.2) then
            PnMAT(i,j,k)=decode1(Lab)
            PnMAT(j,i,k)=DBLE((-1)**(Lij-Sij))*PnMAT(i,j,k)
          else if (mp.eq.3 .and. .not.(option(5))) then
            SnMAT(i,j,k)=decode1(Lab)
            SnMAT(j,i,k)=DBLE((-1)**(Lij-Sij))*SnMAT(i,j,k)
          endif  ! MP
        endif  ! OKtoStart
      endif  ! LAB1
      goto 50
 100  write(*, '("***FATAL: End of the data file reached in readMnPnSnMat")')   
      write(io,'("***FATAL: End of the data file reached in readMnPnSnMat")') 
      stop      
!
      return
      end subroutine readMnPnSnMat
!
!-----------------------------------------------------------------------
!
      SUBROUTINE readTnMat(nline)
!
!  Reads the matrices for the three particle CI Tn matrix elements.
!  T2,T3,T4,T6,T7,T8
! 
      IMPLICIT none
      integer i,j,k, nelem,nline,Tindex(8),Lij,Sij,DelSeniority2
      character lab1*3, lab2*3, Lab*40
      data Tindex/-1,1,2,3,-1,4,5,6/
      logical OKtoStart
!----------
      nelem=0; OKtoStart=.true.; k=1
      if (Nelectrons.lt.3) return
      if (Lvalue.ne.3) return   !  Only for f electrons
      if (check1(7)) then
        write(idebug,'(80("-"),/,"Tn reduced matrix elements, positioned at line:",I5)') nline
      endif
 50   READ(Idata,'(1X,A3,1X,A3,1X,A40)',end=100) LAB1,LAB2,LAB
!      write(idebug,'(2A3,A40)') LAB1,LAB2,LAB
      nline=nline+1
      if (.not.OKtoStart .and. LAB1.eq."U2 ") then
        write(io,'(" No Tk reduced matrix elements in data file")')
        return
      endif
      if (LAB1.eq."ELE") then  ! f^(n>3)
        if (check1(7)) write(idebug,'(I5," Tn reduced matrix elements read.")') nelem
        return
      else if (LAB1(1:1).eq."T") then
        OKtoStart=.true.
        i=ichar(LAB1(2:2))-48
        if (i.lt.2 .or. i.gt.8 .or. i.eq.5) Then
          write(io,'("***FATAL: Invalid TnMat index=",I4,". Program stopped in readTnMat at line:",  &
          i4,/," line:",1X,A3,1X,A3,3X,A40)') k,nline,LAB1,LAB2,LAB
          STOP
        endif
        k=Tindex(i)
        if (k.lt.1 .or. k.gt.6) Then
          write(io,'("***FATAL: Invalid TnMat index=",I4,". Program stopped in readTnMat at line:",  &
          i4,/," line:",1X,A3,1X,A3,3X,A40)') k,nline,LAB1,LAB2,LAB
          STOP
        endif
      else
        if (OKtoStart) then
          nelem=nelem+1
          if (LAB1.ne."   ") i=TS_LookupC(Nelectrons,LAB1)
          j=TS_LookupC(Nelectrons,LAB2)
          TnMAT(i,j,k)=decode1(Lab)
!          Sij=(TS_f_Bases(1,i,Nelectrons)-TS_f_Bases(1,j,Nelectrons))/2  !  (S - S')  always integer
!          Lij=(TS_f_Bases(2,i,Nelectrons)-TS_f_Bases(2,j,Nelectrons))    !  (L - L')  always integer
!          TnMAT(j,i,k)=DBLE((-1)**(Lij-Sij))*TnMAT(i,j,k)
!          DelSeniority2=(TS_f_Bases(3,i,Nelectrons)-TS_f_Bases(3,j,Nelectrons))/2    !  (v - v')/2  always integer
!          if (DelSeniority2*2.ne.(TS_f_Bases(3,i,Nelectrons)-TS_f_Bases(3,j,Nelectrons))) then
!            write(io,'("i,j,k=",3I2,";  Sen1=",I2,";  Sen2=",I2,";  (Sen1-Sen2)/2=",I2)')     &
!                            i,j,k,TS_f_Bases(3,i,Nelectrons),TS_f_Bases(3,j,Nelectrons),DelSeniority2
!          endif
!          TnMAT(j,i,k)=DBLE((-1)**(1+DelSeniority2))*TnMAT(i,j,k)
          if (complementary .and. k.gt.1) TnMAT(i,j,k)=-TnMAT(i,j,k)  ! take care of T2 TnMAT(i,j,1) elements in t2MatComp below.
           TnMAT(j,i,k)=TnMAT(i,j,k)
          if (check1(7)) write(idebug,'(1X,A3,1X,A3,3I4,3x,F12.8)')LAB1,LAB2,i,j,k,TnMAT(i,j,k)
        endif ! OKtoStart
      endif
      goto 50
 100  write(*,'("***FATAL: End of the data file reached in readTnMat")')
      stop 
!
      return
      end subroutine readTnMat
!
!-----------------------------------------------------------------------
!
      SUBROUTINE readCCFMat()
!
!  Reads the matrices for CCF matrix elements from file CCF.dat
!  CCFmat(i1,j1,i,k)
!  i1,j1 is the N&K ordered SL states
!  i=1,12 corresponding to 1,2,3,4,5,6,7,8,9,10A,10B,11
!  k=1,6 corresponding to k=2,4,6,8,10,12
!
!  for delta function model, 
! 
      IMPLICIT none
      real*8 RME
      integer i,j,k, i1,j1, nelem,nc
      character lab1*4, lab2*4, lab3*5, Lab*40, ltest*40, FMTR*20
      logical OKtoStart,writeMat
!----------
      nelem=0; k=1
      if (Nelectrons.lt.2) return
      if (Lvalue.ne.3) return   !  Only for f electrons
      if (check2(5)) write(idebug,'(80("-"),/,"CCF reduced matrix elements")') 

      write(ltest,'("FOR f",i1)') Nelectrons
40   READ(iccf,'(A4,A4,A5,A40)',end=100) LAB1,LAB2,LAB3,LAB
!      write(idebug,'("LAB1,LAB2,LAB3,LAB:|",A4,"|",A4,"|",A5,"|",A40,"|")') LAB1,LAB2,LAB3,LAB
      if (LAB.ne.ltest)  goto 40 ! position at the appropriate part of the file
      
!      
 50   READ(iccf,'(A4,A4,A5,A40)',end=100) LAB1,LAB2,LAB3,LAB
!      write(idebug,'("LAB1,LAB2,LAB3,LAB:|",A4,"|",A4,"|",A5,"|",A40,"|")') LAB1,LAB2,LAB3,LAB
      if (LAB1.eq."****") goto 50
      if (LAB1.eq."CCF ") then    ! finish
        if (check2(5)) write(idebug,'(I5," CCF reduced matrix elements read.")') nelem
          nc=20; nc=min(nc,TS_n(nelectrons))
          WRITE(FMTR,'("(",I2,"(F7.2))")') nc
!

        if (CCFtype.eq.1) then
          if (check2(5)) write(idebug,'("  i  j",16x,"k=2        k=4")') 
          do i1=1,TS_n(Nelectrons)
            do j1=i1,TS_n(Nelectrons)
              if (complementary) then
                CCFmat(i1,j1,1,1)=-CCFmat(i1,j1,1,1)  
                CCFmat(i1,j1,1,2)=-CCFmat(i1,j1,1,2)  
              endif  
              CCFmat(j1,i1,1,1)= CCFmat(i1,j1,1,1)
              CCFmat(j1,i1,1,2)= CCFmat(i1,j1,1,2)
              if (check2(5)) then
                 do k=1,3; CCFtest(j1,i1,k)= CCFtest(i1,j1,k); enddo
                 write(idebug,'(2I3,1X,A3,1X,A3,3x,2F12.8)') i1,j1,TS_labels(i1,Nelectrons), &
                                                TS_labels(j1,Nelectrons),CCFmat(i1,j1,1,1),CCFmat(i1,j1,1,2)
              endif                         
            enddo ! i1
          enddo ! j1
          if (check2(5)) then
            do k=1,3
            if (k.eq.1) write(idebug,'(80("-"),/,"CCF delta model g20 matrix elements:",/,  &
            " 7/18*sr(105)*g12 + 35/6*sr(14)*g22 -35/22*sr(154)*g32 -28/143*sr(15015)*g102", &
            /,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=1,nc)
            if (k.eq.2) write(idebug,'(80("-"),/,"CCF delta model g40 matrix elements:",/,  &
            " -21/44*sr(154)*g14 -21/22*sr(1155)*g24 -63/22*sr(105)*g34 -84/715*sr(30030)*g10A4 -8232/12155*sr(3315)*g10B4", &
            /,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=1,nc)
            if (k.eq.3) write(idebug,'(80("-"),/,"CCF delta model g60 matrix elements:",/,  &
            " -35/198*sr(3003)*g16 -35/66*sr(10010)*g26 -35/22*sr(910)*g36 +56/6479*sr(8843835)*g10A6 +588/5797*sr(34255)*g10B6", &
            /,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=1,nc)
            do i=1,nc; write(idebug,FMTR)  (CCFtest(i,j,k)*PFactor,j=1,nc); enddo
            enddo ! k
          endif ! check2(5)
!
        else if (CCFtype.eq.2) then
          write(idebug,'(80("-"),/,"SCCF model matrix elements: VK=",i2,/,20(" |",A3,"> "))') 2*k,(TS_labels(j,Nelectrons),j=1,nc)
          do i=1,nc; write(idebug,FMTR)  (CCFmat(i,j,1,k)*PFactor,j=1,nc); enddo  
        else if (CCFtype.eq.3) then
          do k=1,6 
            do i=1,12
              writeMat=.false.
              do j=1,N_CCF
                if (i.eq.CCFindex(j,1) .and. 2*k.eq.CCFindex(j,2)) then
		   writeMat=.true.
                   write(idebug,'("j,CCFindex(j,1),CCFindex(j,2)=",3i3)') j,CCFindex(j,1),CCFindex(j,2)
		   exit
		endif   
              enddo 
              if (writeMat) then
                  if (check2(5)) write(idebug,'(80("-"),/,"CCF general matrix elements: g_",A3,"_",i2,/,20(" |",A3,"> "))') &
                                                                       CCF_label(i),2*k,(TS_labels(j,Nelectrons),j=1,nc)
                  do i1=1,nc; do j1=i1,nc; CCFmat(j1,i1,i,k)=CCFmat(i1,j1,i,k); enddo; enddo
                  if (complementary) then
                    do i1=1,nc; do j1=1,nc; CCFmat(i1,j1,i,k)=-CCFmat(i1,j1,i,k); CCFmat(j1,i1,i,k)=CCFmat(i1,j1,i,k);enddo; enddo
                  endif  
                  if (check2(5)) then; do i1=1,nc; write(idebug,FMTR)  (CCFmat(i1,j1,i,k)*PFactor,j1=1,nc); enddo; endif
              endif ! writeMat
            enddo ! i
          enddo ! k            
        endif  !  CCFtype
        return
      else if (LAB1(1:3).eq."Zg_") then   ! Zg_10B_12 MATRIX FOR f3  longest
        OKtoStart=.true.
        ltest=LAB1(4:4)//LAB2//LAB3(1:1)
!        write(idebug,'("g_=",A40)') ltest
        call getGik(ltest,i,k)  ! returns: i,k where k=k/2
      else
        if (OKtoStart) then
          nelem=nelem+1
          if (CCFtype.eq.1.and.k.gt.2) then  ! k#2,4
            OKtoStart=.false.; goto 50
          endif
          if (LAB1.ne."    ") i1=TS_LookupC(Nelectrons,LAB1(2:4))
          j1=TS_LookupC(Nelectrons,LAB2(2:4))
          RME=decode1(Lab)
 !         if (check2(5)) write(idebug,'("I1,J1=",2I4,3x,F12.8)') i1,j1,RME

          if (CCFtype.eq.1) then
            if (k.eq.1) then  ! k=2
              if (i.eq.2) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)+(35.0d+00/3.0d+00)*SQRT(3.5d+00)*RME
              if (i.eq.3) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)-(35.0d+00)*SQRT(7.0d+00/22.0d+00)*RME
              if (i.eq.10) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)-(28.0d+00)*SQRT(105.0d+00/143.0d+00)*RME 
            else if (k.eq.2) then  ! k=4
              if (i.eq.2) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)-(10.5d+00)*SQRT(105.0d+00/11.d+00)*RME
              if (i.eq.3) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)+(63.0d+00/22.00d+00)*SQRT(105.0d+00)*RME
              if (i.eq.10) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)+(84.0d+00)*SQRT(42.0d+00/715.0d+00)*RME    
              if (i.eq.11) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)+(8232.0d+00/11.0d+00)*SQRT(3.0d+00/1105.0d+00)*RME   
            endif
            if (check2(5)) then
              if (k.eq.1) then  ! k=2
                if (i.eq.1)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)+(7.0d+00/18.0d+00)*SQRT(105.0d+00)*RME 
                if (i.eq.2)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)+(35.0d+00/3.0d+00)*SQRT(3.5d+00)*RME
                if (i.eq.3)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)-(35.0d+00)*SQRT(7.0d+00/22.0d+00)*RME
                if (i.eq.10) CCFtest(i1,j1,k)=CCFtest(i1,j1,k)-(28.0d+00)*SQRT(105.0d+00/143.0d+00)*RME 
              else if (k.eq.2) then  ! k=4
                if (i.eq.1)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)-(21.0d+00/44.0d+00)*SQRT(154.0d+00)*RME  
                if (i.eq.2)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)-(10.5d+00)*SQRT(105.0d+00/11.d+00)*RME
                if (i.eq.3)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)+(63.0d+00/22.00d+00)*SQRT(105.0d+00)*RME
                if (i.eq.10) CCFtest(i1,j1,k)=CCFtest(i1,j1,k)+(84.0d+00)*SQRT(42.0d+00/715.0d+00)*RME  
                if (i.eq.11) CCFtest(i1,j1,k)=CCFtest(i1,j1,k)+(8232.0d+00/11.0d+00)*SQRT(3.0d+00/1105.0d+00)*RME   
              else if (k.eq.3) then  ! k=6
                if (i.eq.1)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)-(35.0d+00/198.0d+00)*SQRT(3003.0d+00)*RME  
                if (i.eq.2)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)-(35d+00/66.0d+00)*SQRT(10010.0d+00)*RME
                if (i.eq.3)  CCFtest(i1,j1,k)=CCFtest(i1,j1,k)-(35.0d+00/22.00d+00)*SQRT(910.0d+00)*RME
                if (i.eq.10) CCFtest(i1,j1,k)=CCFtest(i1,j1,k)+(56.0d+00/6479.0d+00)*SQRT(8843835.0d+00)*RME  
                if (i.eq.11) CCFtest(i1,j1,k)=CCFtest(i1,j1,k)+(588.0d+00/5797.0d+00)*SQRT(34255.0d+00)*RME   
              endif ! k
            endif ! check2(5)
          else if (CCFtype.eq.2) then
              if (i.eq.1) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)+((7.0d+00-Nelectrons)/8.0d+00)*RME
              if (i.eq.2) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)-SQRT(30.0d+00)/8.d+00*RME
              if (i.eq.3) CCFmat(i1,j1,1,k)=CCFmat(i1,j1,1,k)-SQRT(330.0d+00)/8.d+00*RME
          else if (CCFtype.eq.3) then
            CCFmat(i1,j1,i,k)=RME
          endif  ! CCFtype

        endif ! OKtoStart
        goto 50
      endif ! LAB1
!
      goto 50
 100  write(*,'("***FATAL: End of the data file ""ccf.dat"" reached in readCCFMat")')
      stop 
!
      return
      end subroutine readCCFMat
!
!-----------------------------------------------------------------------
!
      SUBROUTINE getGik(ltest,i,k)
!
!  Reads the label for the CCF matrix elements.
!  examples: 1_2, 9_12,10B_12
!  i will return as 1-12 for i=1,2,3,4,5,6,7,8,9,10A,10B,11
!  k will return as 1-6  for k=2,4,6,8,10,12
! 
      IMPLICIT none
      integer i,i1,k
      character ltest*40,lab1*3,lab2*2, labs1(12)*3,labs2(6)*2
      data labs1/"1","2","3","4","5","6","7","8","9","10A","10B","11"/
      data labs2/"2","4","6","8","10","12"/
! 
      
      i=scan(ltest,"_")
      if (i.eq.0) then
        write( *,'("***FATAL: getGik; Error reading file ""ccf.dat""; a ""_"" not found, ltest=",A44)') ltest
        write(io,'("***FATAL: getGik; Error reading file ""ccf.dat""; a ""_"" not found, ltest=",A44)') ltest
        stop
      endif
      i1=i-3; if (i1.lt.1) i1=1
      lab1=ltest(i1:i-1); lab2=ltest(i+1:i+3)
      i=0; k=0
      do i1=1,12; if (lab1.eq.labs1(i1)) i=i1; enddo
      do i1=1,6;  if (lab2.eq.labs2(i1)) k=i1; enddo
!      write(idebug,'("ltest=",A44,", lab1,lab2=",A3,1x,A2,",  i,k=",2i3)') ltest,lab1,lab2,i,k
      if (i.eq.0 .or. k.eq.0) then
        write( *,'("***FATAL: getGik; Error reading file ""ccf.dat""; i,k=",2i3)')  i,k
        write(io,'("***FATAL: getGik; Error reading file ""ccf.dat""; i,k=",2i3)')  i,k
        stop
      endif
!      
      return
      end subroutine getGik
!
!-----------------------------------------------------------------------
!
      SUBROUTINE t2MatComp(Wrongt2)
!
!  Adjusts the t2 matrix elements for complementary configurations.
!  Taking into account the 2-electron component. See: 
!  Hansen, Judd, Crosswhite, Atomic Data & Nucl. Data Tables, 62, 1, (1996).
!  Yeung, Tanner, J. Alloys & Comps, 575, 54, (2013).
! 
      IMPLICIT none
      logical Wrongt2
      integer i,j,ncomp,N
      real*8 C,T2,T2p,T2c,bit
      parameter(bit=1.0d-10)
      if (Lvalue.ne.3) return   !  Only for f electrons
      n=nelectrons
      ncomp = 14 - nelectrons
      C=70.0d+00*sqrt(2.0d+00)
      if (Wrongt2) write(io,110)
 110  FORMAT('**** WARNING: Errors have been introduced into the T2 MEs for complimentary configurations.')
      if (check1(7)) then
        if (Wrongt2) write(idebug,110)
        if (.not. Wrongt2) write(idebug,'(80("-"),/,"SUBROUTINE t2MatComp",/,    &
                                   17X,"f^",I1," t2",7X,"f^",I1," t2''",7X,"f^",I2," t2")') n,n,ncomp
      endif                             
      do i=1,TS_n(n)
        do j=i,TS_n(n)
          if (Wrongt2) then           ! FNCR(7) = T
            TnMAT(i,j,1) =-TnMAT(i,j,1)
            TnMAT(j,i,1) = TnMAT(i,j,1)
          else
            T2=TnMAT(i,j,1)  ! non-complementary t2
            TnMAT(i,j,1) = T2 - (n-2)*EEMat(i,j,4)/C ! non-complementary t2'
            t2p=TnMAT(i,j,1)
            TnMAT(i,j,1) =-T2p + (ncomp-2)*EEMat(i,j,4)/C   ! complementary t2
            t2c=TnMAT(i,j,1)
            TnMAT(j,i,1) = TnMAT(i,j,1)
            if (check1(7).and. abs(t2).gt.bit) write(idebug,'("<",A3,"|T2|",A3,">:",3F13.8)') &
                   TS_labels(i,n),TS_labels(j,n),T2,T2p,T2c
          endif
        enddo ! j 
      enddo ! i  
      if (check1(7)) write(idebug,'(80("-"))')
!
      return
      end subroutine t2MatComp
!
!-----------------------------------------------------------------------
!
      SUBROUTINE testTnMat()
!
!  Tests the orthogonality conditions for the Tn matrix elements.
!  Equation (3) in Hansen et al, Atomic Data Nucl Data,62,1,(1996).
!  Note: it tests t2' rather than t2. 
!  Uses eq(1) in above paper, together with e3 stored in EEmat(i,j,4).
! 
      IMPLICIT none
      integer i,j,k1,k2,nerr,n, TSp1i,TLp1i,TSp1j,TLp1j
      Real*8 F(6),SumT,c1,rFactorial(14)
!
      rFactorial(1)=1.0d+00
      do i=2,14; rFactorial(i)=rFactorial(i-1)*dble(i); enddo   

      if (Lvalue.ne.3) return   !  Only for f electrons
      F(1)=63.0d+00  ! (This value for t2') 540.0d+00/7.0d+00   !  These values from table 1
      F(2)=90.0d+00
      F(3)=90.0d+00
      F(4)=30.0d+00
      F(5)=30.0d+00
      F(6)=30.0d+00
      if (nelectrons.lt.3) return
      if (nelectrons.gt.3) then
        do i=1,6    ! Equation (4)
          F(i)=F(i)*rFactorial(8)/(rFactorial(nelectrons-3)*rFactorial(11-nelectrons))
        enddo
      endif
      if (check1(7)) write(idebug,'("SUBROUTINE testTnMat")')
!      if (idebug.gt.0) write(idebug,'(17X,"t2          e3         t2'' ")')
      n=TSL_f(Nelectrons); nerr=0
      c1=dble(nelectrons-2)/(70.0d+00*sqrt(2.0d+00))      
      if (complementary) c1=dble((14-nelectrons)-2)/(70.0d+00*sqrt(2.0d+00))      
      do k1=1,6      
        do k2=1,6      
          SumT=0.0d+00
          do i=1,n
            TSp1i=  TS_f_Bases(1,i,nelectrons)+1
            TLp1i=2*TS_f_Bases(2,i,nelectrons)+1
            do j=1,n   
              TSp1j=  TS_f_Bases(1,j,nelectrons)+1
              TLp1j=2*TS_f_Bases(2,j,nelectrons)+1
              if (TSp1i.eq.TSp1j .and. TLp1i.eq.TLp1j) then
      if (k1.eq.1 .and. k2.eq.1) then
!        if (TSp1i.eq.TSp1j.and.TLp1i.eq.TLp1j .and. j.ge.i) write(idebug,'(2(1X,A3,1X),3F12.6)') TS_labels(i,Nelectrons),  &
!           TS_labels(j,Nelectrons),TnMAT(i,j,k1),EEmat(i,j,4),(TnMAT(i,j,k1)-c1*EEmat(i,j,4))
        SumT = SumT + (TnMAT(i,j,k1)-c1*EEmat(i,j,4))*(TnMAT(i,j,k2)-c1*EEmat(i,j,4))*TSp1i*TLp1i
      endif
      if (k1.eq.1 .and. k2.ne.1) SumT = SumT + (TnMAT(i,j,k1)-c1*EEmat(i,j,4))*TnMAT(i,j,k2)*TSp1i*TLp1i
      if (k1.ne.1 .and. k2.eq.1) SumT = SumT +  TnMAT(i,j,k1)*(TnMAT(i,j,k2)-c1*EEmat(i,j,4))*TSp1i*TLp1i
      if (k1.ne.1 .and. k2.ne.1) SumT = SumT +  TnMAT(i,j,k1)*TnMAT(i,j,k2)*TSp1i*TLp1i
              ENDIF
            enddo ! end j
          enddo ! end i
          if (k1.eq.1 .and. k2.eq.1 .and. check1(7)) write(idebug,'("t2'' Sum=",F12.6)') SumT
          if (k1.eq.k2) then
            if (abs(F(k1)-SumT).gt.1.0d-10) then
              write(io,'("**** F(",I1,")=",F12.2," does not equal SumT=",F12.2)') k1,F(k1),SumT
              if (check1(7)) write(idebug,'("**** F(",I1,")=",F12.2," does not equal SumT=",F12.2)') k1,F(k1),SumT
              nerr=nerr+1
            endif  
          else
            if (abs(SumT).gt.1.0d-10) then
              write(io,'("**** T(",I1,") not orthogonal to T(",i1,");  SumT=",F12.2)') k1,k2,SumT
              if (check1(7)) write(idebug,'("**** T(",I1,") not orthogonal to T(",i1,");  SumT=",F12.2)') k1,k2,SumT
              nerr=nerr+1
            endif  
          endif  ! k1.eq.k2
        endDO  ! end k2
      Enddo  ! end k1
      if (check1(7)) then
        if (nerr.eq.0 .and. check1(7)) write(idebug,'("Passed Tn test")')
        if (nerr.gt.0 .and. check1(7)) write(idebug,'("**** Did not pass Tn test; nerr=",I4)') nerr
        write(idebug,'(80("-"))')
      endif
!
      return
      end subroutine testTnMat
!
!-----------------------------------------------------------------------
!
      SUBROUTINE getLSBasis()
! Fill the arrays for d or f orbital states:
! TS_n(i)        The number of states for a d^i or f^i system.
! TS_labels(j,i) The labels j=1,TS_n(i) for each d(f)^i configuration.
! TS_bases(k,j,i) The q.n. k=L,S,sen Quant. No. for each j state of the d(f)^i configuration. 
      IMPLICIT none
      INTEGER*4 i,j,k,ii
      INTEGER*8 factorial(14)
!
      factorial(1)=1
      do i=2,14; factorial(i)=factorial(i-1)*i; enddo   
!
      if (Lvalue.eq.2) then 
        do i=1,2*Lvalue+1   ! i=1,..5
          TS_n(i)=TSL_d(i)
          TJp1(i)=TJp1_d(i)
          TS_full(i)=factorial(10)/(factorial(i)*factorial(10-i))
          do j=1,n_states
            TS_labels(j,i)=TS_d_labels(j,i)
            do k=1,3
              TS_Bases(k,j,i)=TS_d_Bases(k,j,i)
            enddo
          enddo  
        enddo 
      else
        do i=1,2*Lvalue+1 ! i=1,..7
          TS_n(i)=TSL_f(i)
          TJp1(i)=TJp1_f(i)
          TS_full(i)=factorial(14)/(factorial(i)*factorial(14-i))
          do j=1,n_states
            TS_labels(j,i)=TS_f_labels(j,i)
            do k=1,3
              TS_Bases(k,j,i)=TS_f_Bases(k,j,i)
            enddo
          enddo  
        enddo 
      Endif
!
      return
      end subroutine getLSBasis
!
!-----------------------------------------------------------------------
!
      SUBROUTINE writeLSBasis()
! Write the basis quantum numbers into "debug.dat"
      IMPLICIT none
      INTEGER*4 irepeat(max_LSstates),i,j,k,ii
!
      k=Nelectrons
      write(idebug,'(80("-"),/,A1,"^",I2," configuration, SL basis in terms of:")') e_lab,Nelectrons 
      write(idebug,'("  2S L sen N&K     R7  G2 other")')
      DO i=1,n_states
        if (TS_Bases(1,i,k).ne.-1) then
          irepeat(i)=0
          do ii=1,i
            if (TS_Bases(1,ii,k).eq.TS_Bases(1,i,k).and.TS_Bases(2,ii,k).eq.TS_Bases(2,i,k)) irepeat(i)=irepeat(i)+1
          enddo ! ii
        endif
      enddo ! i
      DO i=1,n_states-1
        if (irepeat(i).eq.1.and.irepeat(i+1).ne.2) irepeat(i)=0
      enddo ! i 
      DO i=1,n_states
        if (TS_Bases(1,i,k).ne.-1) then
          if (irepeat(i).eq.0) write(idebug,'(3I3,3X,I1,A1,1X,2X,3I4)')    &
            (TS_Bases(j,i,k),j=1,3),TS_Bases(1,i,k)+1,orbital(TS_Bases(2,i,k)),(TS_f_QNs(j,i,k),j=1,3)
          if (irepeat(i).ne.0) write(idebug,'(3I3,3X,I1,A1,I1,2X,3I4)')    &
            (TS_Bases(j,i,k),j=1,3),TS_Bases(1,i,k)+1,orbital(TS_Bases(2,i,k)),irepeat(i),(TS_f_QNs(j,i,k),j=1,3)
        endif
      enddo ! i
!      
      return
      end subroutine writeLSBasis
!
!-----------------------------------------------------------------------
!
      SUBROUTINE getFullBasis(nExplicitbasis,nExBasisNums)
!
!  Make the full J,MJ basis from the appropriate LS basis.
!  The QNs in the full basis are 2S, L, 2J, 2MJ, K, seniority.
!  Where K is the number in the "standard" order used by Nielson & Koster.
!
!  If explicitMin=explicitMax then all the basis functions are used. 
!  Otherwise only those that fall in the range [explicitMin,explicitMax],
!  using the standard order, are included. 
!
!  The number of basis functions:
!    n  !   SL  |  SLJ  |  SLJMJ
!  -------------------------------
!    1  |    1  |    2  |    14 
!    2  |   17  |   13  |    91 
!    3  |   17  |   41  |   364 
!    4  |   47  |  107  |  1001 
!    5  |   73  |  295  |  2002 
!    6  |  119  |  295  |  3003 
!    7  |  119  |  327  |  3432 
!
!  max_Jstates = 327
!
      IMPLICIT none
      INTEGER*4 i,ii,it,j,k, idiag
      integer*4 iLS,iTwoS,iL,iTwoJ,iTwoMJ,iJ1,iJ2,iLSJ, jLS,jTwoS,jL,jTwoJ,jTwoMJ,jJ1,jJ2,jLSJ
      INTEGER*4 nExplicitbasis,nExBasisNums(max_Jstates)
      REAL*8 D1,D2,SOfac,W6J,CV,SOC,ER
      logical found,putIn
      Parameter(D1=1.0d+00, D2=2.0d+00)

      found=.true. 
      n_states=TS_n(Nelectrons)
      i=0; ii=0; it=0; Nspin=0

      if (nExplicitbasis.eq.0) then  ! then full basis
        if (check1(1)) write(idebug,'(80("-"),/," Full SLJMJ basis",/,"Number  2S   L   2J   2MJ index Sen  SL term ")')
        do iLS=1,n_states
          iTwoS=TS_Bases(1,iLS,Nelectrons)
          iL   =TS_Bases(2,iLS,Nelectrons)
          iJ1=abs(2*iL-iTwoS)
          iJ2=abs(2*iL+iTwoS)
 !       NJ=(J2-J1)/2
          do iTwoJ=iJ1,iJ2,2
            ii=ii+1; JBasis(1,ii)=iTwoS;JBasis(2,ii)=iL;
            JBasis(3,ii)=iTwoJ;JBasis(4,ii)=iLS;JBasis(5,ii)=TS_Bases(3,iLS,Nelectrons)
            do iTwoMJ=-iTwoJ,iTwoJ,2
              i=i+1
              FullBasis(1,i) = iTwoS                       ! 2S
              FullBasis(2,i) = iL                          ! L
              FullBasis(3,i) = iTwoJ                       ! 2J
              FullBasis(4,i) = iTwoMJ                      ! 2MJ
              FullBasis(5,i) = iLS                         ! Place in the LS reduce matrix elements of Nielson&Koster
              FullBasis(6,i) = TS_Bases(3,iLS,Nelectrons)  ! seniority
              if (check1(1)) write(idebug,'(I4,6I5,6X,A3)') i,(FullBasis(k,i),k=1,6),TS_labels(iLS,Nelectrons)
            enddo ! MJ
          enddo ! J
        enddo !iSL 
        if (i.ne.TS_full(Nelectrons)) then
          write(io,'("***FATAL: Invalid basis size=",I4,"; should be=",i4)') i,TS_full(Nelectrons)
          stop
        endif 
      endif  ! .not.explicitBasis
!
      if (nExplicitbasis.gt.0) then  ! then use only basis functions that fall in the range [explicitMin:explicitMax]
        if (check1(1)) write(idebug,'(80("-"),/," Explicit SLJMJ basis",/,"Number  2S   L   2J   2MJ index Sen  SL term ")')
        do iLS=1,n_states
          iTwoS=TS_Bases(1,iLS,Nelectrons)
          iL   =TS_Bases(2,iLS,Nelectrons)
          iJ1=abs(2*iL-iTwoS)
          iJ2=abs(2*iL+iTwoS)
 !       NJ=(J2-J1)/2
          do iTwoJ=iJ1,iJ2,2
            ii=ii+1
            JBasis(1,ii)=iTwoS;JBasis(2,ii)=iL;
            JBasis(3,ii)=iTwoJ;JBasis(4,ii)=iLS;JBasis(5,ii)=TS_Bases(3,iLS,Nelectrons)
            putIn=.false.
            do k=1,nExplicitbasis; if (ii.eq.nExBasisNums(k)) putIn=.true.; enddo
            do iTwoMJ=-iTwoJ,iTwoJ,2
              it=it+1     ! counts total
              if (putIn) then
                i=i+1     ! counts included
                FullBasis(1,i) = iTwoS                       ! 2S
                FullBasis(2,i) = iL                          ! L
                FullBasis(3,i) = iTwoJ                       ! 2J
                FullBasis(4,i) = iTwoMJ                      ! 2MJ
                FullBasis(5,i) = iLS                         ! Place in the LS reduce matrix elements of Nielson&Koster
                FullBasis(6,i) = TS_Bases(3,iLS,Nelectrons)  ! seniority
                if (check1(1)) write(idebug,'(I4,6I5,6X,A3)') i,(FullBasis(k,i),k=1,6),TS_labels(iLS,Nelectrons)
              endif  
            enddo ! MJ
          enddo ! J
        enddo !iSL 
!        if (i.ne.nExplicitbasis) then
!          write(io,'("***FATAL: Inconsistent explicit basis specified; nExplicitbasis=",I4,"; found=",i4)') nExplicitbasis,i
!          stop
!        endif 
        if (it.ne.TS_full(Nelectrons)) then
          write(io,'("***FATAL: Invalid basis size=",I4,"; should be=",i4)') i,TS_full(Nelectrons)
          stop
        endif 
      endif  ! nExplicitbasis.gt.0
!
      n_matrix=i; n_jmatrix=ii
      if (check1(1)) write(idebug,'("A total of ",I4," |J,MJ> states in basis.")') i
      if (n_matrix.gt.max_N) then
        write(io,'("***FATAL: Problem too big, must redimension program from",I5," to",I5)') max_N,n_matrix
        STOP
      endif
      
      Nspin=1
      S_mult(Nspin)=FullBasis(1,1)+1   ! 2S+1
      do iLS=2,n_matrix
        found=.false.
        do i=1,Nspin
          if ((FullBasis(1,iLS)+1).eq.S_mult(i)) found=.true.
!          write(IO,'("ILS,I,FullBasis(1,iLS)+1, S_mult(i),found=",2I3,2X,2I4,L4)') iLS,i,FullBasis(1,iLS)+1,S_mult(i),found
        enddo
        if (.not.found) then  
          Nspin=Nspin+1
          S_mult(Nspin)=FullBasis(1,iLS)+1
        endif
      enddo ! iLS
      
      if (check1(1)) then
        write(idebug,'("Spin multiplicities:",4I4)') (S_mult(i),i=1,Nspin)
        write(idebug,'(80("-"),/," SLJ basis",/,"Number  2S   L   2J   index Sen  SL term ")')
        do i=1,n_jmatrix
          write(idebug,'(I4,5I5,6X,A3)') i,(JBasis(k,i),k=1,5),TS_labels(JBasis(4,i),Nelectrons)
        enddo
      endif
      
      if (nExplicitbasis.gt.0) then
        write(io,'("Explicit basis; instead of the full ",I4," basis, only ",I4," |SLJMJ> basis functions are used,")') it,n_matrix
        write(io,'("that are in the explicitly specified |SLJ> basis: ")')  
        if (N_odd) then
          write(io,'(10(I4,2X,A3,"(",I2,"/2)"  ))') (nExBasisNums(k),TS_labels(JBasis(4,nExBasisNums(k)),Nelectrons), &
                                                                       JBasis(3,nExBasisNums(k)),k=1,nExplicitbasis)
        else
          write(io,'(10(I4,2X,A3,"(",I2,")"  ))') (nExBasisNums(k),TS_labels(JBasis(4,nExBasisNums(k)),Nelectrons),  &
                                                                   JBasis(3,nExBasisNums(k))/2,k=1,nExplicitbasis)
        endif
        write(io,'("Look at file ""debug.dat"" to see what the complete SLJ basis is.")')
        write(io,'(100("-"))')
        write(idebug,'(80("-"),/," SLJ basis")')
        if (N_odd) then
          write(idebug,'(10(I4,2X,A3,"(",I2,"/2)"  ))') (i,TS_labels(JBasis(4,i),Nelectrons),JBasis(3,i),i=1,n_jmatrix)
        else
          write(idebug,'(10(I4,2X,A3,"(",I2,")"  ))') (i,TS_labels(JBasis(4,i),Nelectrons),JBasis(3,i)/2,i=1,n_jmatrix)
        endif
      endif
      
      return
      end subroutine getFullBasis
!
!-----------------------------------------------------------------------
!
      real*8 FUNCTION G2R7lookup(QN,Ident)
 
!  Returns the eigenvalues of the G2 (Ident=2) and R7 (Ident=3) operators.
!
      IMPLICIT none
      integer QN,Ident
      G2R7lookup=0.0d+00
      if (Ident.eq.2) then
        select case (QN)
        case(0); G2R7lookup=0.0d+00
        case(10); G2R7lookup=0.5d+00
        case(11); G2R7lookup=1.0d+00
        case(20); G2R7lookup=14.0d+00/12.0d+00
        case(21); G2R7lookup=21.0d+00/12.0d+00
        case(22); G2R7lookup=30.0d+00/12.0d+00
        case(30); G2R7lookup=24.0d+00/12.0d+00
        case(31); G2R7lookup=32.0d+00/12.0d+00
        case(40); G2R7lookup=36.0d+00/12.0d+00
        case default; Write(io,'("Unknown eigenvalue of G2")') 
        END Select
      ELSE IF (Ident.eq.3) then
        select case (QN)
        case(0);   G2R7lookup=0.0d+00
        case(100); G2R7lookup=0.6d+00
        case(110); G2R7lookup=1.0d+00
        case(111); G2R7lookup=1.2d+00
        case(200); G2R7lookup=1.4d+00
        case(210); G2R7lookup=1.8d+00
        case(211); G2R7lookup=2.0d+00
        case(220); G2R7lookup=2.4d+00
        case(221); G2R7lookup=2.6d+00
        case(222); G2R7lookup=3.0d+00
        case default; Write(io,'("Unknown eigenvalue of R7")') 
        END Select
      ENDIF
      Return
      end function G2R7lookup
!      
!-----------------------------------------------------------------------
!
END MODULE f_e_readFile
!  1233