MODULE f_e_Block
USE f_e_data
USE f_e_parameters
!       
!  This module does the symmetry blocking to reduce the size of the matrix to be diagonalised.
!

IMPLICIT NONE

CONTAINS
!
!-----------------------------------------------------------------------
!
      subroutine BlockBasis(qM)
!
!  Determines what sort of blocking can be made based on the non-zero Bkq values.
!  Separates the basis functions into separate blocks.
!  See: Goeller-Walrand, "Handbook" PCREs,23, 121, (1996). pp 156-7.
!  
!            even             number(size)          odd                number(size)    
!   qmin   M(mod(qmin))       of matrices       2M(mod(qmin))         of matrices  
!     2      0,1                 2(~N/2)    (-1/2,+1/2)                 1(~N/2)
!     3      0,(-1,+1)           2(~N/3)    (-1/2,+1/2),3/2             2(~N/4)
!     4      0,(-1,+1),2         3(~N/4)    (-1/2,+1/2),(-3/2,+3/2),
!     5      0,(-1,+1),(-2,+2)              (-1/2,+1/2),(-3/2,+3/2),5/2
!     6      0,(-1,+1),(-2,+2),3            (-1/2,+1/2),(-3/2,+3/2),(-5/2,+5/2)
!
!     7      MJ blocking
!
!  If qMin is non-zero on input, then the symmetry is forced (neglecting the actual values of the crystal field).
!  This happens, when PRED is used and diagonalisation of the atom terms only is made, before the crystal field is applied.
!  It is important to get the correct WFs if the calculation is subsequently going to use the matrix elements to compare
!  with the Extended Stevens Operator calculation. The degenerate levels will have the MJ values mixed.  
!
      IMPLICIT NONE
      character*5 lab,CryQN(2,6,6) ! crystal quantum numbers: (odd,even),qMin,mu
      DATA CryQN/72*"XXXXX"/ 
      logical odd,q1,q2,q3,q4,q5,q6, b3,b4,b6,b7,b8,b9,b11,b12,b13,b14,forceSym/.false./,TR/.false./
      integer i,j,TwoMJi,TMJ, IB,qMin,qM,MM,iodd, maxn
      real*8 z,s 
      parameter(z=0.0d+00)
      s = CFbit
!
      CryQN(1,2,1)=" -1/2"; CryQN(1,2,2)=" +1/2" 
      CryQN(1,3,1)=" -1/2"; CryQN(1,3,2)=" +1/2"; CryQN(1,3,3)="  3/2"
      CryQN(1,4,1)=" -1/2"; CryQN(1,4,2)=" +1/2"; CryQN(1,4,3)=" +3/2"; CryQN(1,4,4)=" -3/2"
      CryQN(1,5,1)=" -1/2"; CryQN(1,5,2)=" +1/2"; CryQN(1,5,3)=" +3/2"; CryQN(1,5,4)=" +5/2"; CryQN(1,5,5)=" -3/2"
      CryQN(1,6,1)=" -1/2"; CryQN(1,6,2)=" +1/2"; CryQN(1,6,3)=" +3/2"; CryQN(1,6,4)=" +5/2"; CryQN(1,6,5)=" -5/2"
      CryQN(1,6,6)=" -3/2"
      CryQN(2,2,1)="  0  "; CryQN(2,2,2)="  1  "; 
      CryQN(2,3,1)="  0  "; CryQN(2,3,2)=" +1  "; CryQN(2,3,3)=" -1  "
      CryQN(2,4,1)="  0  "; CryQN(2,4,2)=" +1  "; CryQN(2,4,3)="  2  "; CryQN(2,4,4)=" -1  "
      CryQN(2,5,1)="  0  "; CryQN(2,5,2)=" +1  "; CryQN(2,5,3)=" +2  "; CryQN(2,5,4)=" -2  "; CryQN(2,5,5)=" -1  "
      CryQN(2,6,1)="  0  "; CryQN(2,6,2)=" +1  "; CryQN(2,6,3)=" +2  "; CryQN(2,6,4)="  3  "; CryQN(2,6,5)=" -2  "
      CryQN(2,6,6)=" -1  "
      odd=.false.; if (mod(Nelectrons,2).eq.1) odd=.true.

      qMin=qM
      if (qMin.ne.0) goto 100  ! force symmetry
      
! If there is no symmtery blocking, important for some things to treat as if there is 1 block.
      if (.not.block) then
        nblocks=1; Nblock(1)=n_matrix
        do i=1,n_matrix; Do j=1,6; BBasis(1,j,i)=FullBasis(j,i); enddo; BBasis(1,7,i)=i; enddo
        return
      endif
!      
!                   1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16 
! DATA BkqR_label/"B00 ","B20 ","B21 ","B22 ","B40 ","B41 ","B42 ","B43 ","B44 ","B60 ","B61 ","B62 ","B63 ","B64 ","B65 ","B66 "/
      b3=(abs(BkqR(3)).lt.s.and.abs(BkqI(3)).lt.s)
      b4=(abs(BkqR(4)).lt.s.and.abs(BkqI(4)).lt.s)
      b6=(abs(BkqR(6)).lt.s.and.abs(BkqI(6)).lt.s)
      b7=(abs(BkqR(7)).lt.s.and.abs(BkqI(7)).lt.s)
      b8=(abs(BkqR(8)).lt.s.and.abs(BkqI(8)).lt.s)
      b9=(abs(BkqR(9)).lt.s.and.abs(BkqI(9)).lt.s)
      b11=(abs(BkqR(11)).lt.s.and.abs(BkqI(11)).lt.s)
      b12=(abs(BkqR(12)).lt.s.and.abs(BkqI(12)).lt.s)
      b13=(abs(BkqR(13)).lt.s.and.abs(BkqI(13)).lt.s)
      b14=(abs(BkqR(14)).lt.s.and.abs(BkqI(14)).lt.s)
!      write(io,'("b3,b4,b6,b7,b8,b9,b11,b12,b13,b14=",10(2X,L1))') b3,b4,b6,b7,b8,b9,b11,b12,b13,b14
      q1=(b3 .and. b6 .and. b11)
      q2=(b4 .and. b7 .and. b12)
      q3=(b8 .and. b13)
      q4=(b9 .and. b14)
      q5=(abs(BkqR(15)).lt.s.and.abs(BkqI(15)).lt.s)
      q6=(abs(BkqR(16)).lt.s.and.abs(BkqI(16)).lt.s)
!      write(io,'("q1,q2,q3,q4,q5,q6=",6(2X,L1))') q1,q2,q3,q4,q5,q6
      qMin=1
      if (.not.q1) qMin=1 
      if (q1.and. .not.q2) qMin=2 
      if (q1.and.q2.and. .not.q3) qMin=3 
      if (q1.and.q2.and.q3 .and. .not.q4) qMin=4 
      if (q1.and.q2.and.q3.and.q4 .and. .not.q5) qMin=5 
      if (q1.and.q2.and.q3.and.q4.and.q5 .and. .not.q6) qMin=6 
      if (q1.and.q2.and.q3.and.q4.and.q5 .and. q6) qMin=7      ! All B2q,B4q,B6q zero. 
! Forces an artificial (high symmetry) ligand field (only for a fit) to limit the array sizes generate during a fit.
      if (qMin.eq.7 .and. nFit.ne.0) then; forceSym=.true.; qMin=6; endif
! Check Time reversal symmetry exists again      
      TR=odd.and.(q1.and.q3.and.q5) 
      if (TR_exists .and. .not.TR) write(io,'("****WARNING: TR_exists=",L,"; TR=",L)') TR_exists,TR
!      
 100  continue

      if (wordy.or.PrintCQN) write(io,'("Symmetry blocking:",/,115("-"))') 
!      write(io,'("Nelectrons=",I4)') Nelectrons
      if (qMin.eq.1) then
        nblocks=1; Nblock(1)=n_matrix
        do i=1,n_matrix; Do j=1,6; BBasis(1,j,i)=FullBasis(j,i); enddo; BBasis(1,7,i)=i; enddo
        if (wordy.or.PrintCQN) write(io,'("No symmetry blocking can be done; nBlocks=",i2)') nBlocks
        return
      else if (qMin.eq.7) then  ! Arrays TJp1_d and TJp1_f are the 2Jmax+1, where Jmax is the maximum J value.
        if (Lvalue.eq.2) then; nBlocks=TJp1_d(Nelectrons); else; nBlocks=TJp1_f(Nelectrons); endif
        if (TR_exists) nBlocks=nBlocks/2
        if (wordy.or.PrintCQN) write(io,'("MJ blocking of atomic calculation; nblocks=",i2)') nBlocks
      ELSE      
        nBlocks=qMin; if (TR_exists) nBlocks=nBlocks/2
        if (wordy.or.PrintCQN) write(io,'("The minimum q amongst non-zero Bkq is qMin=",i2)') qMin 
        if (forceSym .and. (wordy.or.PrintCQN))   &
        write(io,'("The MJ blocking is limited to an artificial symmetry to reduce array sizes during fit")') 
      endif
      do i=1,nBlocks; Nblock(i)=0; enddo
!
      do i=1,n_matrix
        TwoMJi=FullBasis(4,i)
        IB=1 
        if (odd) then
          MM=(TwoMJi+1)/2
          if (qmin.eq.2) then
            if (MODULO(MM,2).eq.0) IB=1; if (MODULO(MM,2).eq.1) IB=2 
          else if (qMin.eq.3) then
            if (MODULO(MM,3).eq.0) IB=1; if (MODULO(MM,3).eq.1) IB=2; if (MODULO(MM,3).eq.2) IB=3 
          else if (qMin.eq.4) then
            if (MODULO(MM,4).eq.0) IB=1; if (MODULO(MM,4).eq.1) IB=2; if (MODULO(MM,4).eq.2) IB=3 
            if (MODULO(MM,4).eq.3) IB=4 
          else if (qMin.eq.5) then
            if (MODULO(MM,5).eq.0) IB=1; if (MODULO(MM,5).eq.1) IB=2; if (MODULO(MM,5).eq.2) IB=3 
            if (MODULO(MM,5).eq.3) IB=4; if (MODULO(MM,5).eq.4) IB=5 
          else if (qMin.eq.6) then
            if (MODULO(MM,6).eq.0) IB=1; if (MODULO(MM,6).eq.1) IB=2; if (MODULO(MM,6).eq.2) IB=3 
            if (MODULO(MM,6).eq.3) IB=4; if (MODULO(MM,6).eq.4) IB=5; if (MODULO(MM,6).eq.5) IB=6 
          else if (qMin.eq.7) then  ! Arrays TJp1_d and TJp1_f are the 2Jmax+1, where Jmax is the maximum J value. 
                    ! IB is then the position of the MJ value.
            if (Lvalue.eq.2) then; 
              IB=(TJp1_d(Nelectrons)-1+TwoMJi)/2+1
              if (TR_exists) then
                if (MODULO(MM,2).eq.1) IB=(TJp1_d(Nelectrons)-1+TwoMJi)/4+1  ! Only store 1 Kramers doublet componant
                if (MODULO(MM,2).ne.1) IB=0                                  ! Don't store the other
              endif  
            else 
              IB=(TJp1_f(Nelectrons)-1+TwoMJi)/2+1
              if (TR_exists) then
                if (MODULO(MM,2).eq.1) IB=(TJp1_f(Nelectrons)-1+TwoMJi)/4+1  ! Only store 1 Kramers doublet componant
                if (MODULO(MM,2).ne.1) IB=0                                  ! Don't store the other
              endif  
            endif
          endif
        else if (.not.odd) then 
          MM=TwoMJi/2        
          if (qMin.eq.2) then
            if (MODULO(MM,2).eq.0) IB=1; if (MODULO(MM,2).eq.1) IB=2
          else if (qMin.eq.3) then
            if (MODULO(MM,3).eq.0) IB=1; if (MODULO(MM,3).eq.1) IB=2; if (MODULO(MM,3).eq.2) IB=3
          else if (qMin.eq.4) then
            if (MODULO(MM,4).eq.0) IB=1; if (MODULO(MM,4).eq.1) IB=2; if (MODULO(MM,4).eq.2) IB=3
            if (MODULO(MM,4).eq.3) IB=4
          else if (qMin.eq.5) then
            if (MODULO(MM,5).eq.0) IB=1; if (MODULO(MM,5).eq.1) IB=2; if (MODULO(MM,5).eq.2) IB=3
            if (MODULO(MM,5).eq.3) IB=4; if (MODULO(MM,5).eq.4) IB=5
          else if (qMin.eq.6) then
            if (MODULO(MM,6).eq.0) IB=1; if (MODULO(MM,6).eq.1) IB=2; if (MODULO(MM,6).eq.2) IB=3
            if (MODULO(MM,6).eq.3) IB=4; if (MODULO(MM,6).eq.4) IB=5; if (MODULO(MM,6).eq.5) IB=6
          else if (qMin.eq.7) then  ! Arrays TJp1_d and TJp1_f are the 2Jmax+1, where Jmax is the maximum J value.
            if (Lvalue.eq.2) then; IB=(TJp1_d(Nelectrons)-1+TwoMJi)/2+1; else; IB=(TJp1_f(Nelectrons)-1+TwoMJi)/2+1; endif
          endif  
        endif
        if (IB.gt.0) then
          Nblock(IB)=NBlock(IB)+1;
!        write(io,'("IB=",I4,";   TwoMJi=",i4,"; Nblock(IB)=",I4)') IB,TwoMJi,Nblock(IB)
          do j=1,6; Bbasis(IB,j,Nblock(IB))=FullBasis(j,i); enddo
          Bbasis(IB,7,Nblock(IB))=i ! position in full basis.
        endif
      enddo ! i
      
      iodd=2; if (odd) iodd=1
      if (qMin.ne.7) then  !  (qMin=7 MJ blocking)
        if (wordy.or.PrintCQN) then
          write(io,'("The (",I4,"x",I4,") matrix was blocked (Crystal quantum numbers below):")') n_matrix,n_matrix
          write(io,'(6("(",i4,"x",i4,") "))') (Nblock(i),Nblock(i),i=1,nBlocks)
          write(io,'(6(2x,A5,5X))') (CryQN(iodd,qMin,i),i=1,nBlocks)
        endif  
        do i=1,nBlocks; CQN(i)= CryQN(iodd,qMin,i); enddo
      else
        do ib=1,nBlocks
!          if (Lvalue.eq.2) then; TMJ=(IB-1)*2+1-TJp1_d(Nelectrons); else; TMJ=(IB-1)*2+1-TJp1_f(Nelectrons); endif
          TMJ=Bbasis(IB,4,1) ! Just take the first one, they should all be the same
          if (odd.and.TMJ.le.0) write(lab,'(I3,"/2")') TMJ
          if (odd.and.TMJ.gt.0) write(lab,'("+",I2,"/2")') TMJ
          if (.not.odd.and.TMJ.le.0) write(lab,'(" ",I3," ")') TMJ/2
          if (.not.odd.and.TMJ.gt.0) write(lab,'(" +",I2," ")') TMJ/2
          CQN(ib) = lab
        enddo
        if (assignMethod.eq.1) then
          write(io,'("Assignment of experimental energy levels to Crystal Quantum numbers:",/, &
                     "Block:",26(3x,i2,2X))') (i,i=1,nBlocks)
          write(io,'("  CQN:",26(1x,A5,1X))') (CQN(i),i=1,nBlocks)
        endif 
        if (wordy.or.PrintCQN) then
          write(io,'("The (",I4,"x",I4") matrix was blocked into:",26("(",i3,"x",i3,")"))') &
                                     n_matrix,n_matrix,(Nblock(i),Nblock(i),i=1,Min(nBlocks,13))
          write(io,'("Crystal quantum numbers:",16X,26(2x,A5,2X))') (CQN(i),i=1,Min(nBlocks,13))
          if (nBlocks.ge.14) then
            write(io,'(48X,26("(",i3,"x",i3,")"))') (Nblock(i),Nblock(i),i=14,nBlocks)
            write(io,'(48X,26(2x,A5,2X))') (CQN(i),i=14,nBlocks)
          endif  
        endif
      endif      
!
      if (check1(1)) then
        maxn=nblock(1); do i=2,nBlocks; if (Nblock(i).gt.maxn) maxn=Nblock(i); enddo
        write(idebug,'("Blocked symmetry numbers",/," 2S L 2J 2MJ N&K sen Original")')
        do i=1,maxn
          write(idebug,'(6(6I3,I4,3X))') ((Bbasis(IB,j,i),j=1,7),IB=1,nBlocks)
        enddo
      endif
!      
      return
!
      end subroutine BlockBasis
!
!-----------------------------------------------------------------------
END MODULE f_e_Block
!  232 lines