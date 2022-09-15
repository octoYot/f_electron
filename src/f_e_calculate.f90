MODULE f_e_calculate
USE f_e_data
USE f_e_parameters
USE f_e_wigner
USE f_e_lf             ! calls function Bkq_index(k,q)
!
!  This module does the main calculations.
!
IMPLICIT NONE

CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE buildAndDiagonalise(iblock,iPRED)
!
!  Build & diagonalise matrix.      
!  The symmetry blocking happens in here. This is the iblock symmetry block.
!  iPRED=0 if normal calculation; energies for all iblock appended to same ENGS(i)
!  iPRED=1 if the preliminary part of a PRED calculation. LF set to zero. 
!          Must be real, Eigenvectors must be calculated.
!          Energies are saved in PReng(iblock,i); Eigenvectors saved in PRvec(iblock,i,j)
! not used iPRED=2 if the preliminary part of a PRED calculation. LF matrix only in original basis. 
!                  No diagonalisation needed.
! not used iPRED=3 if PRED calculation in the prediagonalised basis.
!
      IMPLICIT NONE
      integer iblock,iPRED, i,j,k, n, Nsum, basis(7,Max_N) 
      real*8 E(max_N),d0
      real t1,t2
      LOGICAL MatComplex1
      parameter(d0=0.0d+00)
      save Nsum     
      call CPU_TIME(T1)
!   
!     if (iblock.eq.1) Nsum=0  ! First time called. 

!  iPRED=1, LF=0, calculate vects, matrix real.
      if (iPRED.eq.1) then; MatComplex1=MatComplex; MatComplex=.false.; endif
      if (iblock.eq.1) then  ! First time called.  ! (Unnecessary??, AR,AI zeroed now in BuildMatrix)
        Nsum=0
        if (calcVecs) then
          if (MatComplex) then
            do i=1,n_matrix; do j=1,n_matrix; CMAT(i,j)=DCMPLX(d0,d0); enddo; enddo
          else 
            do i=1,n_matrix; do j=1,n_matrix; MAT(i,j)=d0; enddo; enddo
          endif
        endif 
      endif
      call getBBasis(basis,iblock,n)
!      if (n.eq.0) write(io,'("iblock,iPRED,n,nBlock(iblock)=",4i3)') iblock,iPRED,n,nBlock(iblock)
!      if (n.eq.0) goto 100  ! Should never be zero, as this is the full, not truncated basis (even if PRED, full atomic ini this call) 
      CALL buildMatrix(basis,n,iPRED)
!
      if (iPRED.eq.2)  goto 100
      if (option(9)) goto 100  ! Write matrix, skip diagonalisation
!
      CALL diagonalise(n,e)
!
      if (iPRED.eq.0) then
        if (idebug.gt.0) write(idebug,'("In buildAndDiagonalise; iblock,n,Nsum=",3i4)') iblock,n,Nsum
        if (idebug.gt.0) write(idebug,'("Energies:",/,(10F10.1))') (e(i),i=1,n)
        do i=1,n 
          engs(i+Nsum)=e(i)
          allCQN(i+Nsum)=CQN(iblock)
          iCQN(i+Nsum)=iblock
        enddo
        if (calcVecs) then
          if (nblocks.eq.1 .and. .not.fastMat(2)) then
            if (MatComplex) then
              do i=1,N; do j=1,N; CMAT(i,j)=DCMPLX(AR(i,j),AI(i,j)); enddo; enddo
            else
              do i=1,N; do j=1,N; MAT(i,j)=AR(i,j); enddo; enddo
            endif      
          else  ! nblocks >1 of time-reversal symmetry blocking
            if (MatComplex) then
              do i=1,N; do j=1,N; CMAT(Bbasis(iblock,7,i),j+Nsum)=DCMPLX(AR(i,j),AI(i,j)); enddo; enddo
            else ! .not.MatComplex  
              do i=1,N; do j=1,N; MAT(Bbasis(iblock,7,i),j+Nsum)=AR(i,j); enddo; enddo
            endif      
          endif ! nblocks 
        endif ! calcVecs  
      elseif (iPRED.eq.1) then
        nPRmat(iblock)=n
        do i=1,n; PReng(iblock,i)=e(i); enddo      
        do i=1,N; do j=1,N; PRvec(iblock,i,j)=AR(i,j); enddo; enddo
      endif
!      write(io,'(I5," energies and vectors added to arrays, now there are:",I5)') Nsum,Nsum+n
 100  if (iPRED.eq.1) then; MatComplex=MatComplex1; endif
      Nsum=Nsum+n
      if (iTime) then
        call CPU_TIME(T2)
        if (iPRED.eq.0) write(io,'("**TIME**: buildAndDiagonalise: Normal calculation CPU=",F6.1,"s")') T2-T1
        if (iPRED.eq.1) write(io,'("**TIME**: buildAndDiagonalise: PRED calc; LF=0    CPU=",F6.1,"s")') T2-T1
        if (iPRED.eq.2) write(io,'("**TIME**: buildAndDiagonalise: PRED calc; LF only CPU=",F6.1,"s")') T2-T1
        if (iPRED.eq.3) write(io,'("**TIME**: buildAndDiagonalise: PRED calc; iPRED.eq.3  CPU=",F6.1,"s")') T2-T1
      endif
!  
      return
      end subroutine buildAndDiagonalise
!
!-----------------------------------------------------------------------
!
      SUBROUTINE buildAndDiagonalisePRED()
!  Build & diagonalise matrix using the Prediagonalised basis.
!  nPRED is the dimension of the restricted basis.
!  PRvec(iblock,i,k) contains the nPRED prediagonalised basis-vectors in terms of the original basis.
!  AR/AI conatins the real or complex eigenvectors in terms of the prediagonalised basis (for each iblock).
!  MAT/CMAT contains the real or complex eigenvectors in terms of the original basis (for each iblock).
!
      IMPLICIT NONE
      integer i,j,k,n,n1,iblock,Nsum
      real*8 E(max_N),d0
      COMPLEX*16 c0
      real t1,t2
      parameter(d0=0.0d+00, c0=(0.0d+00,0.0d+00))  
!      
      if (MatComplex) then
        do i=1,n_matrix; do j=1,nPRED; CMAT(i,j)=c0; enddo; enddo
      else  
        do i=1,n_matrix; do j=1,nPRED;  MAT(i,j)=d0; enddo; enddo
      endif 
!      
      call CPU_TIME(T1)
!       
      Nsum=0
      do iblock=1,nblocks
        n=nPRmat(iblock)
!        if (n.le.0) write(io,'("iblock,n,nPRmat(iblock)=",4i3)') iblock,n,nPRmat(iblock)
        if (n.le.0) goto 100  !     can be 0 if using PRED & BLOC, atomic calculation, some high MJ blocks 0 elements
        call transformPRED(iblock)
        if (printMat3(1)) then 
          write(Imatrix,'(/,"Full matrix written in buildAndDiagonalisePRED; iblock=",i2)') iblock
          call printMatrix3(1,n)
        endif  
        if (printMat3(2)) then
          write(Imatrix,'(/,"Full zero/nonzero (./*) matrix written in buildAndDiagonalisePRED; iblock=",i2)') iblock
          If (printMat3(2)) call printMatrix3(2,n)
        endif
!
        CALL diagonalise(n,E)
        do i=1,n; engs(i+Nsum)=e(i); allCQN(i+Nsum)=CQN(iblock); iCQN(i+Nsum)=iblock; enddo
!
        if (.not.calcVecs) goto 100
!      
! transfrom eigenvectors back into original basis
!       write(imatrix,'(/,"buildAndDiagonalisePRED iblock:",i2," Energy & vects:")') iblock
!       write(imatrix,'(20F9.2)') (E(j),j=1,n) 
!       do i=1,n; write(imatrix,'(20F9.4)') (AR(i,j),j=1,n); enddo 
!       write(imatrix,'(/,"buildAndDiagonalisePRED iblock:",i2," PRvec:")') iblock
!       do i=1,nblock(iblock); write(imatrix,'(20F9.4)') (PRvec(iblock,i,j),j=1,n); enddo 

        if (MatComplex) then
          do i=1,nblock(iblock)
            n1=BBasis(iblock,7,i)
            do j=1,n
              do k=1,n; CMAT(n1,j+Nsum) = CMAT(n1,j+Nsum) + PRvec(iblock,i,k)*DCMPLX(AR(k,j),AI(k,j)); enddo 
            enddo
          enddo
        else
          do i=1,nblock(iblock)
            n1=BBasis(iblock,7,i)
!            write(imatrix,'("n,i,n1="3i3)') n,i,n1
            do j=1,n
              do k=1,n; MAT(n1,j+Nsum) = MAT(n1,j+Nsum) + PRvec(iblock,i,k)*AR(k,j); enddo  
            enddo
          enddo
        endif      
 !       write(imatrix,'(/,"After backtransform:")') 
 !       do i=1,n_matrix; write(imatrix,'(20F9.4)') (MAT(i,j+Nsum),j=1,n); enddo 
!  
 100    CONTINUE
!
        Nsum=Nsum+nPRmat(iblock)
      enddo ! iblock
      
      call orderE_V(calcVecs)
!      
      if (printMat3(8)) then 
        write(imatrix,'(/,"buildAndDiagonalisePRED full nPRED:",I4,"(maximum of 20 Eng/vectors shown")') nPRED
        write(imatrix,'(20F9.2)') (Engs(j),j=1,min(20,nPRED)) 
        do i=1,n_matrix; write(imatrix,'(20F9.4)') (MAT(i,j),j=1,min(20,nPRED)); enddo 
      endif
!
      if (iTime) then
        call CPU_TIME(T2)
        write(io,'("**TIME**: buildAndDiagonalisePRED: Calc in prediagonalised basis CPU=",F6.1,"s")') T2-T1
      endif
!      
      return
      end subroutine buildAndDiagonalisePRED
!
!-----------------------------------------------------------------------
!
      subroutine doPRED1() 
!  Do a prediagonalisation in the full basis using atomic only parameters.
!
!  If the matrix elements in a prediagonalised matrix are to be compared to the MEs
!  of an EXSO calculation, then you have to make sure that the states are ordered 
!  the same way: MJ = -J, ...., +J.
!
      IMPLICIT NONE
      integer i,j,k,TwoMJval(Max_N),TwoMJ,mj,n1,allocation_status,iBlock,n2,n3,i3
      real*8 maxMJ,binMJ(-30:30),c2
      logical MJfound(-30:30),found,calcVecs1
      character(len=99) :: emsg
!
      calcVecs1=calcVecs; calcVecs=.true. ! have to calculate eigenvectors, override main calculation
!  this will have already been called.
!     call BlockBasis(0)   ! No ligand field; but do the symmetry blocking that is appropriate for a non-zero LF.
!
      if (.not.ALLOCATED(PRvec)) then
        n1 = (max_N+24)/nblocks  ! nblocks=1,2,3,4,6 divides evenly into max_N=3432; for nblocks>6, MJ blocking below.
                                 ! +24 for f^7
        if (Lvalue.eq.3 .and. nblocks.gt.6) n1=TSLJ_f(Nelectrons)  ! needs to be bigger for the case of f electrons MJ blocking
        if (wordy) write(io,'("*wordy: About to Allocate PRvec(nblocks,n1,n1); nblocks=",I4,", n1=",i4)') nblocks,n1
        allocate (PRvec(nblocks,n1,n1), stat=allocation_status, errmsg=emsg)
        allocate (LFMAT(nblocks,n1,n1), stat=allocation_status, errmsg=emsg)
        allocate (CLFMAT(nblocks,n1,n1), stat=allocation_status, errmsg=emsg)
        if (allocation_status > 0) then
          write(io,'("***FATAL: error Allocating PRvec/LFMAT/CLFMAT in subroutine doPRED1:",A)') trim(emsg)
          write(*, '("***FATAL: error Allocating PRvec/LFMAT/CLFMAT in subroutine doPRED1:",A)') trim(emsg)
          stop
        end if
      endif
!
      if (doAtomic) call unloadP1() ! different atomic parameters than main calulation. 
      do i=1,nBlocks
        CALL buildAndDiagonalise(i,1)   ! no LF
      enddo
      if (doAtomic) call unloadP() ! unload parameters for main calculation. 
      call orderPRED()  ! find nPRvec(iblock) & nPRmult(iblock,i) for a particular nPRED
!
! The energies & eigenvectors are in separate blocks 
! The eigenvectors are in PRvec. They must be real; as the atomic matrix elements are all real.

!  Offset
      if (offsetType.eq.1) then
        DO iBlock=1,nBlocks
          DO i=1,nPRmat(iblock)
!            PReng(iBlock,i)=Engs(i)
            do j=1,offsetN
              if (ioffset(j,1).eq.nPRmult(iblock,i)) then
                PReng(iblock,i)=PReng(iBlock,i)+Roffset(j)
              endif
            enddo ! offsetN 
          enddo ! nPRmat(iblock)
        enddo ! nBlocks
        call orderE_V(calcVecs)
      endif ! offsetType.eq.1
!
      if (ESO2J.gt.0) then    
!  Loop for every eigenvector:
        DO i=1,nPRED
          do j=-30,30; binMJ(j)=z; enddo
          DO j=1,n_matrix
            TwoMJ=FullBasis(4,j)
            c2 = MAT(j,i)**2
            if (TwoMJ.lt.-30 .or. TwoMJ.gt.30) then
              write(io,'("TwoMJ=",I4," out of range in FreeIonPR1")') TwoMJ
              stop
            else
              binMJ(TwoMJ) = binMJ(TwoMJ)+c2
            endif
          enddo  !  j
          maxMJ=0.0d+00 
          do j=-30,30; if (binMJ(j).gt.maxMJ) then; maxMJ=binMJ(j); TwoMJval(I)=j; endif; enddo
          write(IO,'("i=",I4,"; Eng=",F9.1,"; 2MJ=",i3,"; binMJ(j)=",31F5.2)') i,PReng(1,i),TwoMJval(i),(binMJ(j),j=-15,15)
        enddo  !  i  loop over the neng eigenvectors
        
        k=0
        do mj=-ESO2J,ESO2J,2; MJfound(mj)=.false.; enddo
        do mj=-ESO2J,ESO2J,2
          found=.false.
          do i=1,nPRED
            if (TwoMJval(i).eq.mj) then
              k=k+1
              DO j=1,n_matrix; MAT1(J,k)=MAT(J,I); ENDDO; 
              found=.true.
            endif
          enddo ! i
          if (found .and. MJfound(mj)) then
            write(io,'("There is more than one MJ=",i3," in the lowest nPRED (=",I4,") states")') mj,nPRED
            stop
          endif
          MJfound(mj)=found          
        enddo ! mj
        do mj=-ESO2J,ESO2J,2 
          if (.not.MJfound(mj)) then
            write(io,'("There is no MJ=",i3," in the lowest nPRED (=",I4,") states")') mj,nPRED
            stop
          endif  
        enddo  
      endif ! ESO2J.gt.0     
!  
      if (printMat3(1)) then
      write(imatrix,'(/,"Atomic parameters only eigenvalues/eigenvectors in doPRED1")')
        write(imatrix,'(80("-"))')
        DO iBlock=1,nBlocks
          n3=(nPRmat(iblock)-1)/20        
          n2=0
          if (n3.gt.0) then
            do i3=1,n3
              n1=n2+1; n2=n1-1+20
              if (block) write(imatrix,'(20(2X,A5,3X))') (CQN(iblock),j=n1,n2) 
              write(imatrix,'(20(3X,I4,3X))') (j,j=n1,n2)
              write(imatrix,'((20F10.3))')    (PReng(iblock,j),j=n1,n2)
              do i=1,nBlock(iblock); write(imatrix,'(20F10.4)') (PRvec(iblock,i,j),j=n1,n2); enddo 
            enddo 
          endif
          n1=n2+1; n2=nPRmat(iblock)
          if (block) write(imatrix,'(20(2X,A5,3X))') (CQN(iblock),j=n1,n2) 
          write(imatrix,'(20(3X,I4,3X))') (j,j=n1,n2)
          write(imatrix,'((20F10.3))')    (PReng(iblock,j),j=n1,n2)
          do i=1,nBlock(iblock); write(imatrix,'(20F10.4)') (PRvec(iblock,i,j),j=n1,n2); enddo 
        enddo
        write(imatrix,'(80("-"))')
      endif
!
      calcVecs=calcVecs1
      return
      end subroutine doPRED1
!     
!-----------------------------------------------------------------------
!
      subroutine doPRED2()
!  Calcuate the LF matrix in the full original blocked basis      
      IMPLICIT NONE
      integer i,j,iblock, basis(7,Max_N),n 
!
      do iblock=1,nBlocks
        n=nblock(iblock)  ! note, must use nblock(iblock); nPRmat(iblock) < nblock(iblock)
        do i=1,n; do j=1,7; basis(j,i)=BBasis(iblock,j,i); enddo; enddo  
        CALL buildMatrix(basis,n,2)  ! only LF
        if (CFcomplex) then  ! Store LF matrix in original basis
          DO i=1,n; DO j=1,n; CLFMAT(iblock,J,I)=COMPLEX(AR(J,I),AI(J,I)); ENDDO; ENDDO   ! Store LF matrix in original basis
        ELSE
          DO i=1,n; DO j=1,n; LFMAT(iblock,J,I)=AR(J,I); ENDDO; ENDDO   ! Store LF matrix in original basis
        ENDIF
!
        if (printMat3(1)) then 
          write(imatrix,'(/,"Ligand field matrix in original blocked basis (max of 20 columns shown)")')
          write(imatrix,'(20(2X,I2,3X))') (j,j=1,min(20,n))
          do i=1,n; write(imatrix,'(20F7.1)') (LFMAT(iblock,i,j),j=1,min(20,n)); enddo 
        endif
!
      enddo ! iblock
  
      return
      end subroutine doPRED2
!     
!-----------------------------------------------------------------------
!
      SUBROUTINE transformPRED(iblock)
! transforms the LF in the original basis, to the PRED basis.  TODO symmetry blocking
!
! Form (H') = (V)^-1.(H).(V) + (Ediag)
!             where V(N_matrix,nPRED) contains columns of eigenvectors.
!                   H(N_matrix,N_matrix) is the LF matrix in original basis.
!                   Ediag(nPRED,nPRED) is the diagonal energies of diagonalised atomic only matrix.
!                   H'(nPRED,nPRED) is the reduced matrix in the prediagonalised basis.
!
      IMPLICIT NONE
      integer i,j,k,iblock,n
      real*8 d0,h,ave
      parameter(d0=0.0d+00, h=0.5d+00)
!      
      n=nblock(iblock)
      if (CFcomplex) then
        do i=1,n
          Do j=1,nPRmat(iblock)
            VR(i,j)=d0;  VI(i,j)=d0
            do k=1,n; VR(i,j) = VR(i,j) + DREAL(CLFMAT(iblock,i,k))*PRvec(iblock,k,j); enddo  
            do k=1,n; VI(i,j) = VI(i,j) + DIMAG(CLFMAT(iblock,i,k))*PRvec(iblock,k,j); enddo  
          enddo  
        enddo  
!      
        do i=1,nPRmat(iblock)
           Do j=1,nPRmat(iblock)
             AR(i,j)=d0;  AI(i,j)=d0
             do k=1,n; AR(i,j) = AR(i,j) + PRvec(iblock,k,i)*VR(k,j); enddo  
             do k=1,n; AI(i,j) = AI(i,j) + PRvec(iblock,k,i)*VI(k,j); enddo  
          enddo  
        enddo
!  Average in forward diagonal
!        do i=1,nPRED
!           Do j=1,nPRED+1-i
!             ave = h*(AR(i,j)+AR(nPRED+1-j,nPRED+1-i))
!             AR(i,j) = ave; AR(nPRED+1-j,nPRED+1-i) = ave  
!             AI(i,j) = ave; AI(nPRED+1-j,nPRED+1-i) = ave  
!          enddo  
!        enddo

      ELSE   ! not CFcomplex
        do i=1,n   ! 6/12/17  n_matrix
          Do j=1,nPRmat(iblock)
            VR(i,j)=d0
            do k=1,n; VR(i,j) = VR(i,j) + LFMAT(iblock,i,k)*PRvec(iblock,k,j); enddo  
          enddo  
        enddo  
 !       write(imatrix,'(/,"PRvec (n x nPRmat(iblock)):")')
 !       do i=1,n; write(imatrix,'(20F8.4)') (PRvec(iblock,i,j),j=1,nPRmat(iblock)); enddo 
 !       write(imatrix,'(/,"LFMAT (n x n):")')
 !       do i=1,n; write(imatrix,'(20F9.2)') (LFMAT(iblock,i,j),j=1,n); enddo 
 !       write(imatrix,'(/,"VR (n x nPRmat(iblock)):")')
 !       do i=1,n; write(imatrix,'(20F9.2)') (VR(i,j),j=1,nPRmat(iblock)); enddo 

        do i=1,nPRmat(iblock)
          Do j=1,nPRmat(iblock)
            AR(i,j)=d0
            do k=1,n; AR(i,j) = AR(i,j) + PRvec(iblock,k,i)*VR(k,j); enddo  
          enddo  
        enddo
 !       write(imatrix,'(/,"AR (nPRmat(iblock) x nPRmat(iblock)):")')
 !       do i=1,nPRmat(iblock); write(imatrix,'(20F9.2)') (AR(i,j),j=1,nPRmat(iblock)); enddo 

!  Average in forward diagonal
!        do i=1,nPRED
!           Do j=1,nPRED+1-i
!             ave = h*(AR(i,j)+AR(nPRED+1-j,nPRED+1-i))
!             AR(i,j) = ave; AR(nPRED+1-j,nPRED+1-i) = ave  
!          enddo  
!        enddo
      endif  ! CFcomplex
! 
      do i=1,nPRmat(iblock); AR(i,i) = AR(i,i) +  PReng(iblock,i); enddo 
!
      if (printMat3(8)) then  
        write(imatrix,'(/,"Ligand field matrix in original basis")')
        if (CFcomplex) then
          write(imatrix,'(20(3X,I3,2X))') (j,j=1,min(20,n))
          do i=1,min(20,n); write(imatrix,'(20F8.1)') (DREAL(CLFMAT(iblock,i,j)),j=1,min(20,n)); enddo 
          write(imatrix,'(20(3X,I3,2X))') (j,j=1,min(20,n))
          do i=1,min(20,n); write(imatrix,'(20F8.1)') (DIMAG(CLFMAT(iblock,i,j)),j=1,min(20,n)); enddo 
        else
          write(imatrix,'(20(3X,I3,2X))') (j,j=1,min(20,n))
          do i=1,min(20,n); write(imatrix,'(20F8.1)') (LFMAT(iblock,i,j),j=1,min(20,n)); enddo 
        endif !  (CFcomplex) 
        write(imatrix,'(/,"Ligand field matrix in PRED basis; iblock=",i2,", nPRmat(iblock)=",i4,   &
                          ", nblock(iblock)=",i4)') iblock,nPRmat(iblock),n
        write(imatrix,'(20(3X,I3,2X))') (j,j=1,min(20,nPRmat(iblock)))
        do i=1,min(20,nPRmat(iblock)); write(imatrix,'(20F8.1)') (AR(i,j),j=1,min(20,nPRmat(iblock))); enddo
        if (CFcomplex) then
          write(imatrix,'(20(3X,I3,2X))') (j,j=1,min(20,nPRmat(iblock)))
          do i=1,min(20,nPRmat(iblock)); write(imatrix,'(20F8.1)') (AI(i,j),j=1,min(20,nPRmat(iblock))); enddo
        endif
      endif  ! printMat3(8)
!      
      return
      end subroutine transformPRED
!
!-----------------------------------------------------------------------
!
      SUBROUTINE buildMatrix(basis,n,iPRED)
!
!  Writes the LOWER triangle.
!  The rest must be zero.
!  n may be < n_matrix if symmetry blocking has occurred.
!
!  IPRED=0 normal
!  IPRED=1 no LF
!  IPRED=2 only LF
!  IPRED=3 only a particular parameter will be non-zero (Call from calcExpect())
!
      IMPLICIT NONE
      integer n,iPRED, i,j,ii,jj,k,q,nz,noff,kk, i1
      INTEGER TwoSi,TwoSj,Li,Lj,TwoJi,TwoJj,TwoMJi,TwoMJj,SENi,SENj
      INTEGER basis(7,Max_N) 
      real*8 D0,D1,D2,RL,XX,X1,X2,ThreeJ,SixJ,CKK,decouple
      real*8 facC,facSen,SOfac,bit,RCCF
      COMPLEX*16 CI,CX,C0,C1,CC
      logical magFld
      PARAMETER(D0=0.0D+00,D1=1.0D+00,D2=2.0D+00,bit=1.0d-6)
      PARAMETER(C0=(D0,D0),C1=(D1,D0),CI=(D0,D1))
!
      magFld=.false.; if (abs(MAGF(1))+abs(MAGF(2))+abs(MAGF(3)).gt.bit) magFld=.true.
!      write(io,'("P(306),P(401),MAGF(3)=",2F14.8)')P(306),P(401),MAGF(3)
!
      AR=D0; AI=D0;   ! zeroise matrices
!
!      write(io,'("CFcomplex, MatComplex=",2L2)') CFcomplex,MatComplex
!      write(io,'("P(121),P(22),P(42)=",3F14.8," ",A4,"=",F14.8," ",A4,"=",F14.8)') &
!                     P(121),P(22),P(42),BkqR_label(3),BkqR(3),BkqI_label(3),BkqI(3)

      facC=D1
      if (Complementary) facC=-D1
      facSen=D1
      SOfac=SQRT(Dble(Lvalue)*(Lvalue+D1)*(D2*Lvalue+D1))
!
      if (Racah) then
        F2 = 49.0d+00*Bracah + 7.0d+00*Cracah
        F4 = 63.0d+00/5.0d+00*Cracah
        F6 = z
      endif
!
      do i=1,n
        TwoSi=basis(1,i)    ! 2S
        Li=basis(2,i)       ! L
        TwoJi=basis(3,i)    ! 2J
        TwoMJi=basis(4,i)   ! 2MJ
        ii=basis(5,i)       ! position in N&K order
        SENi=basis(6,i)     ! seniority
        do j=1,i
          TwoSj=basis(1,j)
          Lj=basis(2,j)
          TwoJj=basis(3,j)
          TwoMJj=basis(4,j)
          jj=basis(5,j)
          SENj=basis(6,j)
          if (Complementary) facSen=DBLE((-1)**(1+(SENi-SENj)/2)) !  (v-v')/2  always integer
!
!  E-E repulsion elements
!                      
          if (iPRED.eq.2) goto 100  !  skip atomic parameters
          
!          if (i.eq.j) AR(i,j)=EAVE  ! 11/12/17  added EAVE in SUBROUTINE zeroE(e,n)
          
          if (TwoJi.eq.TwoJj .and. TwoMJi.eq.TwoMJj) then
            if (Li.eq.Lj .and. TwoSi.eq.TwoSj) then           
!              if (Complementary) then
!                AR(i,j)=AR(i,j)-facSen*(F2*EEmat(ii,jj,1)+F4*EEmat(ii,jj,2)+F6*EEmat(ii,jj,3))
!              else
                AR(i,j)=AR(i,j)+F2*EEmat(ii,jj,1)+F4*EEmat(ii,jj,2)+F6*EEmat(ii,jj,3)
!              endif
            endif
          endif
!
!  Spin-orbit coupling
!
          if (TwoJi.eq.TwoJj .and. TwoMJi.eq.TwoMJj) then
            AR(i,j)=AR(i,j)+facC*zeta*(-1)**((TwoJi+TwoSj)/2+Li)*SOfac                         &
                    *Wig6J(DBLE(Li),DBLE(Lj),D1,DBLE(TwoSj)/D2,DBLE(TwoSi)/D2,DBLE(TwoJi)/D2)  &
                    *VMat(ii,jj)
          endif
          if (Lvalue.eq.2) goto 100
!
!  M0,M2,M4 spin-other-orbit elements
!  P2,P4,P6 electrostatically correlated spin-orbit elements
!   Application of Eq.3 Judd, et al, Phys.Rev., 169, 130, (1968).
          if (TwoJi.eq.TwoJj .and. TwoMJi.eq.TwoMJj) then
            AR(i,j)=AR(i,j)+(-D1)**((TwoJi+TwoSj)/2+Li)                                        &
                    *Wig6J(DBLE(TwoSj)/D2,DBLE(Lj),DBLE(TwoJi)/D2,DBLE(Li),DBLE(TwoSi)/D2,D1)  &
                    *(M0*MnMAT(ii,jj,1)+M2*MnMAT(ii,jj,2)+M4*MnMAT(ii,jj,3)                    &
                     +P2*PnMAT(ii,jj,1)+P4*PnMAT(ii,jj,2)+P6*PnMAT(ii,jj,3))
          Endif
!
!  M0,M2,M4 spin-spin elements
!   Application of Eq.5 Yeung&Tanner, J.Alloys&Comps., 575, 54, (2013).
          if (.not.(option(5))) then
            if (TwoJi.eq.TwoJj .and. TwoMJi.eq.TwoMJj) then
              AR(i,j)=AR(i,j)+(-D1)**((TwoJi+TwoSj)/2+Li)                                      &
                    *Wig6J(DBLE(TwoSj)/D2,DBLE(Lj),DBLE(TwoJi)/D2,DBLE(Li),DBLE(TwoSi)/D2,D2)  &
                    *(M0*SnMAT(ii,jj,1)+M2*SnMAT(ii,jj,2)+M4*SnMAT(ii,jj,3))
            Endif
          endif
!
!  T2,T3,T4,T6,T7,T8 three body electrostatic elements
!
          if (TwoJi.eq.TwoJj .and. TwoMJi.eq.TwoMJj) then
            if (Li.eq.Lj .and. TwoSi.eq.TwoSj) then           
                AR(i,j)=AR(i,j)+(T2*TnMAT(ii,jj,1)+T3*TnMAT(ii,jj,2)+T4*TnMAT(ii,jj,3)   &
                                +T6*TnMAT(ii,jj,4)+T7*TnMAT(ii,jj,5)+T8*TnMAT(ii,jj,6))
            endif
          Endif
!
 100      continue
          if (iPRED.eq.1) goto 200 ! no LF
!  Ligand field
!
!       if (i.eq.1.and.j.eq.1) write(io,'("3J((3,2,3),(0,0,0))=",F8.6)') Wig3J(3.0d+00,2.0d+00,3.0d+00,D0,D0,D0)
          if (CFcomplex) then
!            CMAT(i,j)=DCMPLX(Mat(I,J))
            if (TwoSi.eq.TwoSj) then
              do k=2,2*Lvalue,2
                RL=DBLE(Lvalue)
                CKK=(-D1)**RL*(D2*RL+D1)*Wig3J(RL,DBLE(k),RL,D0,D0,D0)
                do q=0,k
                  CX=C0
                  if (ABS(BkqR(Bkq_index(k,q)))+ABS(BkqI(Bkq_index(k,q))).gt.CFbit) then
                    decouple=(-D1)**((TwoJi+TwoSi)/2+Lj+k)*Sqrt((TwoJi+D1)*(TwoJj+D1))           &
                              *Wig6J(DBLE(TwoJi)/D2,DBLE(TwoJj)/D2,DBLE(k),DBLE(Lj),DBLE(Li),DBLE(TwoSi)/D2)
                    X1=(-D1)**q*CKK*(-D1)**((TwoJi-TwoMJi)/2)*decouple*UMat(ii,jj,k)           &
                     *Wig3J(DBLE(TwoJi)/D2,DBLE(k),DBLE(TwoJj)/D2,-DBLE(TwoMJi)/D2, DBLE(q),DBLE(TwoMJj)/D2) 
                    if (q.eq.0) then
                      CX=BkqR(Bkq_index(k,q))*X1 
                    else
                      X2=CKK*(-D1)**((TwoJi-TwoMJi)/2)*decouple*UMat(ii,jj,k)                  &
                      *Wig3J(DBLE(TwoJi)/D2,DBLE(k),DBLE(TwoJj)/D2,-DBLE(TwoMJi)/D2,-DBLE(q),DBLE(TwoMJj)/D2) 
                      CX=BkqR(Bkq_index(k,q))*(X1+X2)+CI*BkqI(Bkq_index(k,q))*(-X1+X2)
                    endif
                    AR(i,j)=AR(i,j)+DREAL(facC*CX)
                    AI(i,j)=AI(i,j)+DIMAG(facC*CX)
                  endif ! BkqR/BkqI nonzero
                enddo ! q
              enddo ! k
            endif ! TwoSi.eq.TwoSj
          ELSE  ! Real Ligand field
  !    if (i.eq.1.and.j.eq.1) write(io,'(" Umat(ii,jj,1)=",F8.4,";  Umat(ii,jj,2)=",F8.4,";  Umat(ii,jj,3)=",F8.4)') &
  !                         UMat(ii,jj,1),UMat(ii,jj,2),UMat(ii,jj,3)
            if (TwoSi.eq.TwoSj) then
              do k=2,2*Lvalue,2
                RL=DBLE(Lvalue)
                CKK=(-D1)**RL*(D2*RL+D1)*Wig3J(RL,DBLE(k),RL,D0,D0,D0)
!                if (i.eq.1.and.j.eq.1) write(io,'("RL=",F4.1,"; k=",i2,"; 3J=",F8.6)') RL,K,Wig3J(RL,DBLE(k),RL,D0,D0,D0)
                do q=0,k
                  XX=D0
                  if (ABS(BkqR(Bkq_index(k,q))).gt.CFbit) then
                    decouple=(-D1)**((TwoJi+TwoSi)/2+Lj+k)*Sqrt((TwoJi+D1)*(TwoJj+D1))           &
                              *Wig6J(DBLE(TwoJi)/D2,DBLE(TwoJj)/D2,DBLE(k),DBLE(Lj),DBLE(Li),DBLE(TwoSi)/D2)
                    X1=CKK*(-D1)**((TwoJi-TwoMJi)/2)*decouple*UMat(ii,jj,k)     
                    if (q.eq.0) then
                      XX=BkqR(Bkq_index(k,q))*X1*Wig3J(DBLE(TwoJi)/D2, DBLE(k),DBLE(TwoJj)/D2,  &
                                                      -DBLE(TwoMJi)/D2,DBLE(q),DBLE(TwoMJj)/D2) 
                    else
                      XX=BkqR(Bkq_index(k,q))*X1*((-D1)**q*Wig3J(DBLE(TwoJi)/D2,  DBLE(k),DBLE(TwoJj)/D2,  &
                                                                -DBLE(TwoMJi)/D2,-DBLE(q),DBLE(TwoMJj)/D2) &
                                                          +Wig3J(DBLE(TwoJi)/D2,  DBLE(k),DBLE(TwoJj)/D2,  &
                                                                -DBLE(TwoMJi)/D2, DBLE(q),DBLE(TwoMJj)/D2))
                    endif ! CF
  !!    if (i.eq.8.and.j.eq.1) then
  !!      X2=((-D1)**q*Wig3J(DBLE(TwoJi)/D2,  DBLE(k),DBLE(TwoJj)/D2,-DBLE(TwoMJi)/D2,-DBLE(q),DBLE(TwoMJj)/D2) &
  !!                  +Wig3J(DBLE(TwoJi)/D2,  DBLE(k),DBLE(TwoJj)/D2,-DBLE(TwoMJi)/D2, DBLE(q),DBLE(TwoMJj)/D2))
  !!      ThreeJ=Wig3J(DBLE(TwoJi)/D2,DBLE(k),DBLE(TwoJj)/D2,-DBLE(TwoMJi)/D2, DBLE(q),DBLE(TwoMJj)/D2)
  !!      SixJ=Wig6J(DBLE(TwoJi)/D2,DBLE(TwoJj)/D2,DBLE(k),DBLE(Lj),DBLE(Li),DBLE(TwoSi)/D2)
  !!      write(io,'("i,j=",2I4,";  k,q=",2I3,";  XX=",F10.4,";  3J=",F10.4,";  6J=",F10.4,";  x1=",F10.4, & 
  !!                 ";  x2=",F10.4,";  Uk=",F10.4,";  BkqR=",F10.4)') i,j,k,q,XX,ThreeJ,SixJ,x1,x2,  &
  !!                  UMat(ii,jj,k),BkqR(Bkq_index(k,q))
  !!      write(io,'(" Uk(ii,jj,k)=",F10.4,"; ii,jj=",2I4"; k=",I4)') & 
  !!                   UMat(ii,jj,k),ii,jj,k
  !!    endif
!                    if (Complementary) then
!                      mat(i,j)=mat(i,j)-facSen*facC*XX
!                    else
                      AR(i,j)=AR(i,j)+facC*XX
!                    endif 
                  endif  ! |BkqR|>0
                enddo   ! q
              enddo  ! k
            ENDIF  ! S=S'
          endif  ! CF Complex/Real
!
 200      CONTINUE
!
        enddo   !  end j
!
!  CI two body electrostatic parameters (orbit-orbit interactions)
        if (iPRED.NE.2) AR(i,i)=AR(i,i)+ALPHA*ABGMat(ii,1)+BETA*ABGMat(ii,2)+GAMMA*ABGMat(ii,3)

      enddo  !  end i
!
!  CCF
!
      if (N_CCF.gt.0) then
!
        do i1=1,12
          do k=2,2*Lvalue,2
            do q=0,k
              RCCF=CCFvalue(i1,k,q) ! *BkqR(Bkq_index(k,q))/BkqR(Bkq_index(k,0))
              if (ABS(RCCF).gt.CFbit) then
                 write(io,'("i,k,q=",3I3,"; CCF=",F12.4)') i1,k,q,RCCF
        
        nz=0
        do i=1,n
          TwoSi=basis(1,i)    ! 2S
          Li=basis(2,i)       ! L
          TwoJi=basis(3,i)    ! 2J
          TwoMJi=basis(4,i)   ! 2MJ
          ii=basis(5,i)       ! position in N&K order
          do j=1,i
            TwoSj=basis(1,j)
            Lj=basis(2,j)
            TwoJj=basis(3,j)
            TwoMJj=basis(4,j)
            jj=basis(5,j)
            
            if (TwoSi.eq.TwoSj) then
!              RL=DBLE(Lvalue)
!              CKK=(-D1)**RL*(D2*RL+D1)*Wig3J(RL,DBLE(k),RL,D0,D0,D0)
              XX=D0; 
              decouple=(-D1)**((TwoJi+TwoSi)/2+Lj+k)*Sqrt((TwoJi+D1)*(TwoJj+D1))           &
                        *Wig6J(DBLE(TwoJi)/D2,DBLE(TwoJj)/D2,DBLE(k),DBLE(Lj),DBLE(Li),DBLE(TwoSi)/D2)
!! ??              X1=CKK*(-D1)**((TwoJi-TwoMJi)/2)*decouple*CCFmat(ii,jj,i1,k/2)   
              X1=(-D1)**((TwoJi-TwoMJi)/2)*decouple*CCFmat(ii,jj,i1,k/2)   
              if (q.eq.0) then
                XX=RCCF*X1*Wig3J(DBLE(TwoJi)/D2, DBLE(k),DBLE(TwoJj)/D2,  &
                                -DBLE(TwoMJi)/D2,DBLE(q),DBLE(TwoMJj)/D2) 
              else
                XX=RCCF*X1*((-D1)**q*Wig3J(DBLE(TwoJi)/D2,  DBLE(k),DBLE(TwoJj)/D2,  &
                                          -DBLE(TwoMJi)/D2,-DBLE(q),DBLE(TwoMJj)/D2) &
                                    +Wig3J(DBLE(TwoJi)/D2,  DBLE(k),DBLE(TwoJj)/D2,  &
                                          -DBLE(TwoMJi)/D2, DBLE(q),DBLE(TwoMJj)/D2))
              endif
              if (abs(XX).gt.1.0d-8) nz=nz+1
              if (abs(XX).gt.1.0d-8) write(idebug,'("i,k,q=",3i3,"<",A3,2i2,"| |",A3,2i2,">=",F12.4)') &
                   i1,k,q, TS_labels(ii,Nelectrons),TwoJi,TwoMJi, TS_labels(jj,Nelectrons),TwoJj,TwoMJj,XX
              AR(i,j) = AR(i,j)+ XX
            ENDIF  ! S=S'
          enddo   !  end j
        enddo  !  end i
!        write(io,'("Number of non-zero lower triangle CCF(i,k,q)=",3i2," elements:",i6)') i1,k,q,nz
        
              endif  ! |CCFvalue|>0
            enddo   ! q
          enddo  ! k
        enddo ! i1
      endif  ! CCF
!
!  Magnetic Field
!
      if (magFld) then
        if (MatComplex) then
          do k=1,3
            CC=C1; if (k.eq.2) CC=CI
!            WRITE(imatrix,'("k=",I2,", MagF(k)=",F8.4)') k,MagF(k)
!            WRITE(imatrix,'("L, nmu=",i4)') nmu(1,k)
!            WRITE(imatrix,'("S, nmu=",i4)') nmu(2,k)
            if (MagF(k).gt.bit) then
! expand each into full (unBLOCked) matrix, (not n)
              do i=1,n_matrix; do j=1,n_matrix; VR(i,j)=d0; VI(i,j)=d0; enddo; enddo
              do i=1,nmu(1,k)
                ii=ImuIndex(1,k,i)
                jj=JmuIndex(1,k,i)
                VR(ii,jj)=DREAL(CC*Bohr*MAGF(k)*RK(k)*MM(1,k,i))
                VI(ii,jj)=DIMAG(CC*Bohr*MAGF(k)*RK(k)*MM(1,k,i))
              enddo  
              do i=1,nmu(2,k)
                ii=ImuIndex(2,k,i)
                jj=JmuIndex(2,k,i)           
                VR(ii,jj)=VR(ii,jj)+DREAL(CC*Bohr*MAGF(k)*Ge*MM(2,k,i))
                VI(ii,jj)=VI(ii,jj)+DIMAG(CC*Bohr*MAGF(k)*Ge*MM(2,k,i))
              enddo
              do i=1,n
                do j=1,i
                  AR(i,j)=AR(i,j)+VR(basis(7,i),basis(7,j))
                  AI(i,j)=AI(i,j)+VI(basis(7,i),basis(7,j))
!                WRITE(imatrix,'(" II,JJ=",2I3,", AR,AI=",2F10.2)') ii,jj,AR(ii,jj),AI(ii,jj)
                enddo
              enddo
            endif !  MagF(k).gt.bit
          enddo ! k
        endif ! MatComplex
        if (.not.MatComplex) then
          do k=1,3,2 ! skip k=2   Ly & Sy not needed
            if (MagF(k).gt.bit) then
              do i=1,n_matrix; do j=1,n_matrix; VR(i,j)=d0; enddo; enddo
              do i=1,nmu(1,k)
                ii=ImuIndex(1,k,i)
                jj=JmuIndex(1,k,i)           
                VR(ii,jj)=Bohr*MAGF(k)*RK(k)*MM(1,k,i)
              enddo
              do i=1,nmu(2,k)
                ii=ImuIndex(2,k,i)
                jj=JmuIndex(2,k,i)           
                VR(ii,jj)=VR(ii,jj)+Bohr*MAGF(k)*Ge*MM(2,k,i)
              enddo
              do i=1,n
                do j=1,i
                  AR(i,j)=AR(i,j)+VR(basis(7,i),basis(7,j))
!                WRITE(imatrix,'(" II,JJ=",2I3,", AR,AI=",2F10.2)') ii,jj,AR(ii,jj),AI(ii,jj)
                enddo
              enddo
            endif ! MagF(k).gt.bit
          enddo ! k
        endif ! .not.MatComplex
      Endif ! magFld
!
      if (OUTP(7)) then
        nz=0; noff=0
        do i=1,n
          do j=1,i
            if (Abs(AR(i,j))+Abs(AI(i,j)).gt.1.0D-10) then
              nz=nz+1
              if (i.ne.j) noff=noff+1
            endif  
          enddo   !  end j
        enddo  !  end i
        write(io,850) n,n,nz,noff,n*(n+1)/2
 850    format(/," Matrix size:(",I4,"x",I4,").   Number of non-zero matrix elements:",I8, /,  &
               " Number of non-zero lower triangle MEs=",I8," (out of a possible lower triangle:",I8,")",/)
      endif   !  OUTP(7)
      
      if (iPRED.eq.0) then
        if (printMat3(1)) write(Imatrix,'("Full matrix written in BuildMatrix")')
        if (printMat3(2)) write(Imatrix,'("Full zero/nonzero (./*) matrix written in BuildMatrix")')
      elseif (iPRED.eq.1) then
        if (printMat3(1)) write(Imatrix,'("Atomic Terms only matrix written in BuildMatrix")')
        if (printMat3(2)) write(Imatrix,'("Atomic Terms only matrix zero/nonzero (./*) written in BuildMatrix")')
      elseif (iPRED.eq.2.or.iPRED.eq.3) then  ! Need full matrix.
        do i=1,n-1
          do j=i+1,n
            AR(i,j) = AR(j,i)
            AI(i,j) =-AI(j,i)
          enddo
        enddo
        if (printMat3(1)) write(Imatrix,'("LF matrix only written in BuildMatrix")')
        if (printMat3(2)) write(Imatrix,'("LF matrix only zero/nonzero (./*) written in BuildMatrix")')
      endif
      if (printMat3(1)) call printMatrix3(1,n)
      if (printMat3(2)) call printMatrix3(2,n)

      call printMatrix1()
      call printMatrix4()

      return
!
  900 write(io,'("***FATAL: Error in Subroutine buildMatrix()")')   
      stop
!
      end subroutine buildMatrix
!
!-----------------------------------------------------------------------
!
      SUBROUTINE buildJMatrix()
! Writes the LOWER triangle matrix in the |SLJ> basis.
! No ligand field matrix elelemnts, so matrix is real.
      IMPLICIT NONE
      integer i,j,ii,jj,k,q,nz,noff,kk
      INTEGER TwoSi,TwoSj,Li,Lj,TwoJi,TwoJj,SENi,SENj
      real*8 D0,D1,D2,RL,XX,X1,X2,ThreeJ,SixJ,CKK,decouple
      real*8 facC,facSen,SOfac,bit
      COMPLEX*16 CI,CX,C0,C1,CC
      logical magFld
      PARAMETER(D0=0.0D+00,D1=1.0D+00,D2=2.0D+00,bit=1.0d-6)
      PARAMETER(C0=(D0,D0),C1=(D1,D0),CI=(D0,D1))
!
      facC=D1
      if (Complementary) facC=-D1
      facSen=D1
      SOfac=SQRT(Dble(Lvalue)*(Lvalue+D1)*(D2*Lvalue+D1))
!
      do i=1,n_jmatrix
        TwoSi=jbasis(1,i)
        Li=jbasis(2,i)
        TwoJi=jbasis(3,i)
        ii=jbasis(4,i)
        SENi=jbasis(5,i)
        do j=1,i
          TwoSj=jbasis(1,j)
          Lj=jbasis(2,j)
          TwoJj=jbasis(3,j)
          jj=jbasis(4,j)
          SENj=jbasis(5,j)
          if (Complementary) facSen=DBLE((-1)**(1+(SENi-SENj)/2)) !  (v-v')/2  always integer

          mat(i,j)=D0; SOtest(i,j)=D0
!
!  E-E repulsion elements
!                      
          if (TwoJi.eq.TwoJj) then
            if (Li.eq.Lj .and. TwoSi.eq.TwoSj) then           
!              if (Complementary) then
!                mat(i,j)=mat(i,j)-facSen*(F2*EEmat(ii,jj,1)+F4*EEmat(ii,jj,2)+F6*EEmat(ii,jj,3))
!              else
                mat(i,j)=mat(i,j)+F2*EEmat(ii,jj,1)+F4*EEmat(ii,jj,2)+F6*EEmat(ii,jj,3)
!              endif
            endif
          endif
!
!  Spin-orbit coupling
!
          if (TwoJi.eq.TwoJj) then
            mat(i,j)=mat(i,j)+facC*zeta*(-1)**((TwoJi+TwoSj)/2+Li)*SOfac                         &
                      *Wig6J(DBLE(Li),DBLE(Lj),D1,DBLE(TwoSj)/D2,DBLE(TwoSi)/D2,DBLE(TwoJi)/D2)  &
                      *VMat(ii,jj)
            if (check2(3)) then
              SOtest(i,j)=facC*(-1)**((TwoJi+TwoSj)/2+Li)*SOfac                         &
                        *Wig6J(DBLE(Li),DBLE(Lj),D1,DBLE(TwoSj)/D2,DBLE(TwoSi)/D2,DBLE(TwoJi)/D2)  &
                        *VMat(ii,jj)
              SOtest(j,i)=SOtest(i,j)
            endif 
          endif
!
!  M0,M2,M4 spin-other-orbit elements
!  P2,P4,P6 electrostatically correlated spin-orbit elements
!   Application of Eq.3 Judd, et al, Phys.Rev., 169, 130, (1968).
          if (TwoJi.eq.TwoJj) then
            mat(i,j)=mat(i,j)+(-D1)**((TwoJi+TwoSj)/2+Li)                                        &
                      *Wig6J(DBLE(TwoSj)/D2,DBLE(Lj),DBLE(TwoJi)/D2,DBLE(Li),DBLE(TwoSi)/D2,D1)  &
                      *(M0*MnMAT(ii,jj,1)+M2*MnMAT(ii,jj,2)+M4*MnMAT(ii,jj,3)                    &
                       +P2*PnMAT(ii,jj,1)+P4*PnMAT(ii,jj,2)+P6*PnMAT(ii,jj,3))
          Endif
!
!  M0,M2,M4 spin-spin elements
!   Application of Eq.5 Yeung&Tanner, J.Alloys&Comps., 575, 54, (2013).
          if (.not.(option(5))) then
            if (TwoJi.eq.TwoJj) then
              mat(i,j)=mat(i,j)+(-D1)**((TwoJi+TwoSj)/2+Li)                                      &
                      *Wig6J(DBLE(TwoSj)/D2,DBLE(Lj),DBLE(TwoJi)/D2,DBLE(Li),DBLE(TwoSi)/D2,D2)  &
                      *(M0*SnMAT(ii,jj,1)+M2*SnMAT(ii,jj,2)+M4*SnMAT(ii,jj,3))
            Endif
          endif
!
!  T2,T3,T4,T6,T7,T8 three body electrostatic elements
!
          if (TwoJi.eq.TwoJj) then
            if (Li.eq.Lj .and. TwoSi.eq.TwoSj) then           
                mat(i,j)=mat(i,j)+(T2*TnMAT(ii,jj,1)+T3*TnMAT(ii,jj,2)+T4*TnMAT(ii,jj,3)   &
                                  +T6*TnMAT(ii,jj,4)+T7*TnMAT(ii,jj,5)+T8*TnMAT(ii,jj,6))
            endif
          Endif
        enddo  ! j   
!
!  CI two body electrostatic parameters (orbit-orbit interactions)
!
        if (abs(ALPHA)+ABS(BETA)+abs(GAMMA).gt.CFbit) then
          mat(i,i)=mat(i,i)+ALPHA*ABGMat(ii,1)+BETA*ABGMat(ii,2)+GAMMA*ABGMat(ii,3)
        endif
      enddo  !  end i
!
      nz=0; noff=0
       do i=1,n_jmatrix
        do j=1,i
          IF (Abs(mat(i,j)).gt.1.0D-10) then
            nz=nz+1       
            if (i.ne.j) noff=noff+1
          endif  
        enddo   !  end j
      enddo  !  end i
      if (OUTP(7)) write(io,'(/," Matrix size:(",I4,"x",I4,").   Number of non-zero matrix elements:",I8,  &
                 /," Number of non-zero lower triangle MEs=",I8," (out of a possible lower triangle:",  &
                   I8,")",/)') n_jmatrix,n_jmatrix,nz,noff,n_jmatrix*(n_jmatrix+1)/2 
      do i=1,9
        if (printMat2(i)) call printMatrix2(i)
      enddo
      return
!
  900 write(io,'("***FATAL: Error in Subroutine buildJMatrix()")')   
      stop
!
      end subroutine buildJMatrix
!
!-----------------------------------------------------------------------
!
      SUBROUTINE diagonalise(N,E)
!      
!  N: The dimension of the matrix to be diagonalised.
!  E: The N Energies calculated.
!  The real (AR) or complex (AR,AI) matrix is diagonalised.
!  Eigenvectors (if required) are returned in the original matrix (AR) or (AR,AI). 
!  Lanczos options:
!  IFL  : 0 STARTING VECTOR PROVIDED.
!       : 1 RANDOM VECTOR USED.
!  LP1  : 0 LANCZOS VECTORS NOT STORED.
!              (MUST BE REGENERATED IF THE EIGENVECTORS ARE WANTED.)
!       :10 LANCZOS VECTORS ARE STORED ON DISK FILE 10.
!  IMESS:   DIAGNOSTIC MESSAGES TO BE PRINTED? (Y/N:1/0)
!
      IMPLICIT NONE
      integer n,i,j,idiag,ifl,lp1,imess
      real*8 fv1(max_N)
      real*8 wk1(max_N),wk2(max_N),wk3(max_N),E(max_N)
      idiag=1 ! just energies
      if (calcVecs) idiag=2
      ifl=1; lp1=0; imess=0
!      write(idebug,'("** idiag=",I4"; calcVecs=",L2,"; N=",I4,"; MaxN=",I4)') idiag, calcVecs, N,max_N
      if (MatComplex) then
        call DIAGCH(AR,max_N,AI,max_N,N,E,VR,max_N,VI,max_N,wk1,wk2,wk3,IDIAG,IO)
      else
        if (lancIt.ne.0) then 
!          call LANCZOS(N,lancIt,E,idiag,ifl,lp1,msome,rlb,rub,imess,io)
        ELSE  
          call diagrs(AR,max_N,N,E,fv1,idiag,io)  ! AR overwritten by eigenvectors here.
        ENDIF  
      endif
! 
      if (check1(9)) then
        write(idebug,'(80("-"),/,"Eigenvalues:")')
        do i=1,N
          write(idebug,'(F12.4)') E(i)
        enddo
      endif

      return
!
      end subroutine diagonalise
!
!-----------------------------------------------------------------------
!
      SUBROUTINE doubleKramers()
!      
!  If time-reversal symmetry is exploited, the other componant of the Kramers doublet 
!  if regenerated here. fastMat(2)=T, fastMat(3)=T
!
!  neng = n_matrix or nPRED (nPRED.ne.0) 
!
!  The eigenvectors have already been expanded into their correct positions in the
!  standard basis order in either subroutines: buildAndDiagonalise() or buildAndDiagonalisePRED()
!
      IMPLICIT NONE
      integer i,j, base,twoMJ, old2S,oldL,old2J,oldLS, minusMJ(max_N),phase(max_N)
      integer base1(max_N)
!
      do i=n_matrix/2,1,-1; engs(2*i)=engs(i); engs(2*i-1)=engs(i); enddo  
!      
      if (.not.calcVecs) return
      minusMJ=0; base=0; old2S=Fullbasis(1,1); oldL=Fullbasis(2,1); old2J=Fullbasis(3,1); oldLS=Fullbasis(5,1)
      do i=1,n_matrix
        if (Fullbasis(1,i).ne.old2S.or.Fullbasis(2,i).ne.oldL.or.Fullbasis(3,i).ne.old2J.or.Fullbasis(5,i).ne.oldLS) then
          base = base + old2J+1  ! 2J+1
          old2S=Fullbasis(1,i); oldL=Fullbasis(2,i); old2J=Fullbasis(3,i); oldLS=Fullbasis(5,i)
        endif 
        twoMJ=Fullbasis(4,i)
        if (twoMJ.lt.0) then
          minusMJ(i)=base+(Fullbasis(3,i)+1)/2-(twoMJ+1)/2+1
        else  
          minusMJ(i)=base+(Fullbasis(3,i)+1)/2-(twoMJ+1)/2+1
        endif  
        base1(i)=base
        phase(i)=(Fullbasis(3,i)-Fullbasis(4,i))/2  ! J-MJ
        if (mod(phase(i),2).eq.0) then; phase(i)=1; else; phase(i)=-1; endif
      enddo
!      write(io,'(//,"    i  base  J    MJ   -MJ    phase  Eng")')
!      do i=1,n_matrix
!        write(io,'(6I5,f12.2)') i,base1(i),Fullbasis(3,i),Fullbasis(4,i),minusMJ(i),phase(i),engs(i)-engs(1)      
!      enddo
!       
      do i=n_matrix/2,1,-1
        if (MatComplex) then
          do j=1,n_matrix; cmat(minusMJ(j),2*i)= phase(j)*DCONJG(cmat(j,i)); enddo          
          do j=1,n_matrix; cmat(j,2*i-1)       = cmat(j,i); enddo
        else
          do j=1,n_matrix;  mat(minusMJ(j),2*i)= phase(j)*mat(j,i); enddo          
          do j=1,n_matrix;  mat(j,2*i-1)       = mat(j,i); enddo
        endif
      enddo  
!
      return
!
      end subroutine doubleKramers
!
!-----------------------------------------------------------------------
!
      SUBROUTINE FreeIonPR1(neng)
!
!  Calculates: 1) The free ion term projections of each wavefunction as a Label in cFreeIon.
!              2) The MJ projections of each wavefunction in cMJ. 
!  The internal constants FImax and MJmax can be changed in the input file.
!
      IMPLICIT none
      integer neng,i,ii,is,inFI,inMJ,k,ktot, TwoS,L,TwoJ,TwoMJ,LS,old2S,oldL,old2J,oldLS
      integer, parameter :: maxFI=8
      integer pFI(maxFI),nFI, pMJ(maxFI),nMJ
      real*8 z,bin(1000),binmax,c2,bin1MJ(-30:30)
      character*9 label,allLabels(1000),cFI(maxFI),cMJ(maxFI),MJlab
      PARAMETER(Z=0.0D+00)
!
      k=0; old2S=-1; oldL=-1; old2J=-1; oldLS=-1
      DO is=1,n_matrix
        TwoS =FullBasis(1,is)   ! 2S
        L    =FullBasis(2,is)   ! L
        TwoJ =FullBasis(3,is)   ! 2J
        LS   =FullBasis(5,is)   !  Place in the LS reduce matrix elements of Nielson&Koster
        if (TwoS.ne.old2S .or. L.ne.oldL .or.  TwoJ.ne.old2J .or. LS.ne.oldLS) then
          k=k+1
          if (mod(TwoJ,2).eq.0) then
            write(label,'(A3,"(",I2,")")') TS_labels(LS,nelectrons),TwoJ/2 
          else
            write(label,'(A3,"(",I2,"/2)")') TS_labels(LS,nelectrons),TwoJ 
          endif
          allLabels(k)=label
          old2S=TwoS; oldL=L; old2J=TwoJ; oldLS=LS
        endif 
      enddo
!      write(IO,'("allLabels(k)=",40A9)') (allLabels(i),i=1,k)

!  Loop for every eigenvector:
      DO II=1,neng
        DO inFI=1,maxFI       ! Each state can have up to maxFI different free-ion components
          pFI(inFI)=0;           pMJ(inFI)=0 
          cFI(inFI)="         "; cMJ(inFI)="         "
        enddo
        do k=-30,30
          bin1MJ(k)=z
        enddo
!
        k=0; old2S=-1; oldL=-1; old2J=-1; oldLS=-1
        DO is=1,n_matrix
          TwoS=FullBasis(1,is); L=FullBasis(2,is); TwoJ=FullBasis(3,is); TwoMJ=FullBasis(4,is); LS=FullBasis(5,is)
          if (TwoS.ne.old2S .or. L.ne.oldL .or.  TwoJ.ne.old2J .or. LS.ne.oldLS) then
            k=k+1
            bin(k) = 0.0d+00
          endif
          IF (.not.MatComplex) c2 = MAT(is,II)**2
          IF (     MatComplex) c2=CDABS(CMAT(is,II))**2
          if (TwoMJ.lt.-30 .or. TwoMJ.gt.30) then
            write(io,'("TwoMJ=",I4," out of range in FreeIonPR1")') TwoMJ
            stop
          else
            bin1MJ(TwoMJ) = bin1MJ(TwoMJ)+c2
          endif
          bin(k) = bin(k)+ c2
          old2S=TwoS; oldL=L; old2J=TwoJ; oldLS=LS
        enddo  !  is
        ktot=k
        do i=-30,30; binMJ(II,i) = bin1MJ(i); enddo
!        write(IO,'("II=",I4,"; binMJ(II,i)=",9F6.3)') II,(binMJ(II,i),i=-4,4)
!
        nFI=0
        do inFI=1,min(maxFI,ktot)
          binmax=FImax
          do k=1,ktot
            if (bin(k).gt.binmax) then
              binmax=bin(k); bin(k)=z; nFI=nFI+1
              pFI(inFI)=INT(100.0*binmax); cFI(inFI)=allLabels(k)
              goto 20
            endif
          enddo  ! ktot
 20       continue
        enddo  ! inFI
        do k=1,ktot
          bin(k)=z
        enddo
!        
        nMJ=0
        do inMJ=1,min(maxFI,20)
          binmax=MJmax
          do k=-30,30
            if (bin1MJ(k).gt.binmax) then
!              binmax=bin1MJ(k); 
              pMJ(inMJ)=INT(100.0*bin1MJ(k)); 
              bin1MJ(k)=z; nMJ=nMJ+1
              if (     N_odd) write(MJlab,'("(",I3,"/2)")') k
              if (.not.N_odd) write(MJlab,'("(",I3,")")')   k/2
              cMJ(inMJ)=MJlab
              goto 30
            endif
          enddo  ! k
 30       continue
        enddo  ! inMJ

!        write(IO,'("SUBROUTINE FreeIonPR: ii,nMJ=",2I4,";  pMJ,cMJ:",4(I3,A9))') ii,nMJ,(pMJ(inMJ),cMJ(inMJ),inMJ=1,nMJ)
!        write(IO,'("SUBROUTINE FreeIonPR: ",4(I3,A9))') (pFI(inFI),cFI(inFI),inFI=1,nFI)
        cFreeIon(ii)=makeFullLabel(pFI,cFI,nFI)
        cAllMJ(ii)=makeMJLabel(pMJ,cMJ,nMJ)
!        write(IO,'("SUBROUTINE FreeIonPR: cAllMJ(ii)=",A60)') cAllMJ(ii)
!        write(IO,'("SUBROUTINE FreeIonPR: cFreeIon(ii)=",A60)') cFreeIon(ii)
!
      enddo  !  II  loop over the neng eigenvectors
!
      RETURN
      END Subroutine FreeIonPR1
!
!-----------------------------------------------------------------------
!
      SUBROUTINE getBBasis(basis,ib,n)
!      
!  Returns basis,n:
!  basis:  the basis functions of block ib.
!  n:      the number of basis functions in basis.
!
      implicit none
      INTEGER i,j,n,ib, basis(7,Max_N)
            
      if (nblocks.eq.1 .and. .not.fastMat(2)) then  ! Time reversal symmetry, basis in Bbasis even if nblocks=1
        n=n_matrix
        do i=1,n_matrix
          do j=1,6
            basis(j,i)=FullBasis(j,i)
          enddo
          basis(7,i)=i
        enddo
      else
        n=Nblock(ib) 
        do i=1,n
          do j=1,7
            basis(j,i)=Bbasis(ib,j,i)
          enddo
        enddo        
      endif

      return
      end subroutine getBBasis
!       
!-----------------------------------------------------------------------
!
      character*80 function makeFullLabel(pFI,cFI,nFI)
!
!  Returns the free ion term projections as a single label.
!  Order so the components are listed in descending order.
!
      IMPLICIT none
      integer i,j,nFI,inFI
      integer, parameter :: maxFI=8
      integer pFI(maxFI),p 
      character cFI(maxFI)*9, c*9, lab1*5, label*80
!
      do i=1,maxFI
        do j=i+1,maxFI
          if (pFI(j).gt.pFI(i)) then
            p=pFI(j)
            pFI(j)=pFI(i)
            pFI(i)=p
            c=cFI(j)
            cFI(j)=cFI(i)
            cFI(i)=c
            goto 10
          endif
        enddo
 10     continue
      enddo
                  
      write(lab1,'(I3)') pFI(1)
      label=trim(lab1)//"*"
      do i=1,len(trim(cFI(1)))
        if (cFI(1)(i:i).ne." ") label=trim(label)//cFI(1)(i:i)
      enddo
      if (nFI.gt.1) then
        do inFI=2,nFI
          write(lab1,'(" +",I3)') pFI(inFI)
          label=trim(label)//trim(lab1)//"*"
          do i=1,len(trim(cFI(inFI)))
            if (cFI(inFI)(i:i).ne." ") label=trim(label)//cFI(inFI)(i:i)
          enddo
        enddo
      endif
!      if (len(label).eq.len(TRIM(label))) write(IO,'("***Warning label in makeFullLabel is full")') 
      if (option(3)) makeFullLabel=label(1:40)
      if (option(4)) makeFullLabel=label 
      RETURN 
      END function makeFullLabel
!
!-----------------------------------------------------------------------
!
      character*80 function makeMJLabel(pMJ,cMJ,nMJ)
!
!  Returns the MJ projections as a single label.
!  Order so the components are listed in descending order.
!
      IMPLICIT none
      integer i,j,nMJ,inMJ
      integer, parameter :: maxFI=8
      integer pMJ(maxFI),p,sumMJ 
      character cMJ(maxFI)*9, c*9, lab1*5, label*80
! Check% add up
      sumMJ=0
      do i=1,maxFI
        sumMJ=sumMJ+pMJ(i)
      enddo
!      if (abs(sumMJ-100).gt.25) write(io,'("***WARNING: Sum in makeMJLabel far from 100%; sumMJ=",I3)') sumMJ
!  Make them descending order:
 10   continue
      do i=1,maxFI
        do j=i+1,maxFI
          if (pMJ(j).gt.pMJ(i)) then
            p=pMJ(j)
            pMJ(j)=pMJ(i)
            pMJ(i)=p
            c=cMJ(j)
            cMJ(j)=cMJ(i)
            cMJ(i)=c
            goto 10
          endif
        enddo
      enddo
                  
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
      if (option(6)) makeMJLabel=label(1:40)
      if (option(7)) makeMJLabel=label 
      RETURN 
      END function makeMJLabel
!      
!-----------------------------------------------------------------------
!
      SUBROUTINE orderE_V(orderVecs)
! Order the eigenvalues, crystal quantum numbers and (IF (orderVecs) ) eigenvectors
      IMPLICIT NONE
      logical orderVecs
      integer i,ii,j,k,n,iq
      real*8 p
      complex*16 CP
      character*5 qn
      
      n=n_matrix; if (fastMat(2)) n=n/2
      if (nPRED.gt.0) n= nPRED
      DO II=2,n
        I=II-1
        K=I
        P=Engs(I)
        qn=allCQN(i)
        iq=iCQN(i)
!
        DO J=II,n
          IF (Engs(J).GE.P) GOTO 260
          K=J
          P=Engs(J)
          qn=allCQN(j) 
          iq=iCQN(j) 
 260      CONTINUE
        enddo ! J
!
        IF (K.EQ.I) GOTO 300
        Engs(K)=Engs(I)
        Engs(I)=P
        allCQN(K)=allCQN(I)
        allCQN(I)=qn
        iCQN(K)=iCQN(I)
        iCQN(I)=iq
!
        if (orderVecs) then
          DO J=1,n_matrix   ! full length.
            if (MatComplex) then
              CP=CMAT(J,I)
              CMAT(J,I)=CMAT(J,K)
              CMAT(J,K)=CP
            else
              P=MAT(J,I)
              MAT(J,I)=MAT(J,K)
              MAT(J,K)=P
            endif
          enddo ! J
        endif ! orderVecs
!
 300    CONTINUE
! 
      enddo  ! II
!
      RETURN 
      END SUBROUTINE orderE_V
!      
!-----------------------------------------------------------------------
!
      SUBROUTINE orderPRED()
!      
! In a PRED caclualtion, the lowest nPRED energy levels are cacluated in a PREDiagonalised (atomic) basis. 
! This subroutines determines the number of levels, and therefore matrix dimensions in each symmetry 
! block that this nPRED corresponds to. These numbers are set in nPRmat(6) and must sum to nPRED. 
! 
! Order the eigenvalues, crystal quantum numbers and (IF (orderVecs) ) eigenvectors
!
      IMPLICIT NONE
      integer i,i1,i2,ii,j,k,iblock,sum,mult,ib(max_N),ibi(max_N)
      real*8 p, e1(max_N)
!     
      sum=0
      do iblock=1,nblocks
        do i=1,nPRmat(iblock); e1(i+sum)=PReng(iblock,i); ib(i+sum)=iblock; ibi(i+sum)=i; enddo
        sum=sum+nPRmat(iblock)
      enddo
!      
      DO II=2,n_matrix
        I=II-1
        K=I
        P=e1(i)
        i1=ib(i)
        i2=ibi(i)
!
        DO J=II,n_matrix
          IF (e1(J).lt.P) then
            K=J
            P=e1(J)
            i1=ib(j)
            i2=ibi(j)
          endif  
        enddo ! J
!
        IF (K.ne.I) then
          e1(K)=e1(I); e1(I)=P
          ib(K)=ib(I); ib(I)=i1
          ibi(K)=ibi(I); ibi(I)=i2
        endif
      enddo  ! II
!
      if (wordy) write(io,'("*wordy: Before truncation nPRED=",i4,",  block sizes:",6i4)') nPRED,(nPRmat(i),i=1,nblocks)
      mult=1
      nPRmult(ib(1),ibi(1))=mult
      do iblock=1,nblocks; nPRmat(iblock)=0; enddo
      nPRmat(ib(1))=1
      do i=2,nPred
        if (e1(i)-e1(i-1).gt.Edegen) mult=mult+1
        nPRmat(ib(i))= nPRmat(ib(i))+1     
        nPRmult(ib(i),ibi(i))=mult
      enddo ! i
      sum=0
      do iblock=1,nblocks; sum=sum+nPRmat(iblock); enddo
!
      if (wordy) then
        write(io,'("*wordy: After truncation nPRED=",i4,",  block sizes:",6i4)') nPRED,(nPRmat(i),i=1,nblocks)
        do iblock=1,nblocks
          write(io,'("*wordy: Block:",i2," truncated to:",I4)') iblock,nPRmat(iblock)
          if (nPRmat(iblock).gt.0) then
            write(io,'("*wordy: Engs:",10F12.2,/,(13x,10F12.2))') (PReng(iblock,i),i=1,nPRmat(iblock))
            write(io,'("*wordy: Mult:",10(4X,I4,4X),/,(13x,10(4X,I4,4X)))') (nPRmult(iblock,i),i=1,nPRmat(iblock))
          endif
        enddo 
      endif      
      if (sum.ne.nPRED) then
        write(io,'("****FATAL: The symmetry block numbers in orderPRED do not sum to nPRED; sum,nPRED=",2I4)') sum,nPRED
        stop
      endif  
!
      RETURN 
      END SUBROUTINE orderPRED
!      
!-----------------------------------------------------------------------
!
      subroutine prepCalc(ipar) 
! Sets the parameters before a calculation
! 
      IMPLICIT NONE
      integer i,ipar
      real DTR
!
!      write(io,'(" prepCalc 1: P(21-26),P(121)=",7F7.1)') (P(i),i=21,26),P(121) 
      call setpEq(ipar)  ! can set AOMchanged to .true.
      call makeLinks()
!      write( *,'(" LFtype:",A4,"; AOMchanged:",L4)') LFtype,AOMchanged 
      CALL unloadP()   ! Everything out of P
!      write(io,'(" prepCalc 1: BKQR(2)=",F12.3,",  BKQR(16)=",F12.3)')BKQR(2),BKQR(16)
!      write(io,'("P(21),P(22):")'); write(io,'((10F12.3))') P(21),P(22)
      if (LFtype.eq."AOM " .or. LFtype.eq."LF1e") then
        if (AOMchanged) then
          if (Lvalue.eq.2) call AOMmatrixD()
          if (Lvalue.eq.3) call AOMmatrixF()
        endif  
      endif 
      if ((abs(RotLF(1))+abs(RotLF(2))+abs(RotLF(3))).gt.1.0d-12) then
        if (fitType.eq.3) then
          do i=1,16; BKQR(i)=savedBKQR(i); BKQI(i)=savedBKQI(i); enddo
        endif
!        write(io,'(" prepCalc: about to call RotateLF; a,b,g=",3F12.4,",  BKQI(3)=",F16.9,",  BKQI(4)=",F16.9)') &
!           RotLF(1),RotLF(2),RotLF(3),BKQI(3),BKQI(4)
        call RotateLF()
      endif  
!      write(io,'(" prepCalc 2: BKQI(15)=",F12.3,",  BKQI(16)=",F12.3,",  MAGF(3)=",F6.2)')BKQI(15),BKQI(16),MAGF(3)
!        write(io,'("Contents of P:")'); write(io,'((10F12.3))') (P(i),i=1,160)
!
      return
      end subroutine prepCalc
!
!-----------------------------------------------------------------------
!
      SUBROUTINE printMatrix1()
! Called at the end of BuildMatrix.
!
! Prints matrices in the |SL> basis.
!  ip = 1  UMat
!       2  Vmat
!       3  F2,F4,F6
!       4  Mn 
!       5  Pn 
!       6  Sn 
!       7  Tn 
!       8  gik
!
      IMPLICIT NONE
      integer i,j,k,n,nc, i1,i2,j1,j2 
      character FMTR*20, FMTC*20
      logical doOut/.false./
      
      do i=1,8; if (printMat1(i)) doOut=.true.; enddo
      if (.not.doOut) return
      
      nc=20; nc=min(nc,TS_n(nelectrons))
      i1=1; i2=nc; j1=1; j2=nc
      if (printMat1(10)) then; i1=printRng1(1); i2=printRng1(2);  &
                               j1=printRng1(3); j2=printRng1(4); nc=printRng1(4)-printRng1(3)+1; endif
!      n=TSL_f(Nelectrons)      
!
      write(iMatrix,'(80("-"),/,"Matrix elements in |SL> basis multiplied by: ",F12.8)') PFactor
      if (printMat1(10)) write(iMatrix,'(" Subset of Matrices printed: i1,i2, j1,j2:",4i3)')  i1,i2, j1,j2
            
      if (printMat1(1)) then
        WRITE(FMTR,'("(",I2,"(F7.4))")') nc
        do k=1,3
          if (k.eq.1) write(iMatrix,'(80("-"),/,"U2 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.2) write(iMatrix,'(/,"U4 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.3) write(iMatrix,'(/,"U6 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          do i=i1,i2
            write(iMatrix,FMTR)  (UMAT(j,i,k)*PFactor,j=j1,j2)           
          enddo
        enddo 
      endif  
      if (printMat1(2)) then
        WRITE(FMTR,'("(",I2,"(F7.3))")') nc
        write(iMatrix,'(80("-"),/,"V11 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
        do i=i1,i2
          write(iMatrix,FMTR)  (Vmat(i,j)*PFactor,j=j1,j2)           
        enddo      
      endif  
      if (printMat1(3)) then
        WRITE(FMTR,'("(",I2,"(F7.4))")') nc
        do k=1,3
          if (k.eq.1) write(iMatrix,'(80("-"),/,"F2 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.2) write(iMatrix,'(/,"F4 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.3) write(iMatrix,'(/,"F6 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          do i=i1,i2
            write(iMatrix,FMTR)  (EEmat(i,j,k)*PFactor,j=j1,j2)           
          enddo
        enddo      
      endif  
      if (printMat1(4)) then
        WRITE(FMTR,'("(",I2,"(F7.2))")') nc
        do k=1,3
          if (k.eq.1) write(iMatrix,'(80("-"),/,"M0 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.2) write(iMatrix,'(/,"M2 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.3) write(iMatrix,'(/,"M4 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          do i=i1,i2
            write(iMatrix,FMTR)  (MnMat(i,j,k)*PFactor,j=j1,j2)           
          enddo
        enddo      
      endif  
      if (printMat1(5)) then
        WRITE(FMTR,'("(",I2,"(F7.4))")') nc
        do k=1,3
          if (k.eq.1) write(iMatrix,'(80("-"),/,"P2 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.2) write(iMatrix,'(/,"P4 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.3) write(iMatrix,'(/,"P6 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          do i=i1,i2
            write(iMatrix,FMTR)  (PnMat(i,j,k)*PFactor,j=j1,j2)           
          enddo
        enddo      
      endif  
      if (printMat1(6)) then
        WRITE(FMTR,'("(",I2,"(F7.4))")') nc
        do k=1,3
          if (k.eq.1) write(iMatrix,'(80("-"),/,"S0 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.2) write(iMatrix,'(/,"S2 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.3) write(iMatrix,'(/,"S4 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          do i=i1,i2
            write(iMatrix,FMTR)  (SnMat(i,j,k)*PFactor,j=j1,j2)           
          enddo
        enddo      
      endif  
      if (printMat1(7)) then
        WRITE(FMTR,'("(",I2,"(F7.4))")') nc
        do k=1,6
          if (k.eq.1) write(iMatrix,'(80("-"),/,"T2 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.2) write(iMatrix,'(/,"T3 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.3) write(iMatrix,'(/,"T4 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.4) write(iMatrix,'(/,"T6 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.5) write(iMatrix,'(/,"T7 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          if (k.eq.6) write(iMatrix,'(/,"T8 matrix elements",/,20(" |",A3,"> "))') (TS_labels(j,Nelectrons),j=j1,j2)
          do i=i1,i2
            write(iMatrix,FMTR)  (TnMat(i,j,k)*PFactor,j=j1,j2)           
          enddo
        enddo      
      endif  
      if (printMat1(8)) then
        WRITE(FMTR,'("(",I2,"(F7.2))")') nc
        if (CCFtype.eq.1) then
          do k=1,N_CCF
            write(iMatrix,'(80("-"),/,"CCF delta model g",i1,"0 matrix elements",/,20(" |",A3,"> "))') 2*k,  &
                                                                               (TS_labels(j,Nelectrons),j=j1,j2)
            do i=i1,i2 
              write(iMatrix,FMTR)  (CCFmat(i,j,1,k)*PFactor,j=j1,j2)           
            enddo
          enddo
        else if (CCFtype.eq.2) then
        else if (CCFtype.eq.3) then
          do k=1,N_CCF
            write(iMatrix,'(80("-"),/,"CCF general model g",A3,"_",i2," matrix elements",/,20(" |",A3,"> "))') &
                                       CCF_label(CCFindex(k,1)),CCFindex(k,2), (TS_labels(j,Nelectrons),j=j1,j2)
            do i=i1,i2 
              write(iMatrix,FMTR)  (CCFmat(i,j,CCFindex(k,1),CCFindex(k,2)/2)*PFactor,j=j1,j2)           
            enddo
          enddo

        endif
      endif
!
      return
      end subroutine printMatrix1
!
!-----------------------------------------------------------------------
!
      SUBROUTINE printMatrix2(ipr)
!
! Prints matrices in |SLJ> basis.
!
!  ip = 1  Full matrix (max nc x nc) 
!       2  Full matrix */. (max 100x100)
!       3  <||Uk||> (printed elsewhere)
!       4  L/S matrices (printed elsewhere)
!       5  <||kL+gS||> (printed elsewhere)
!       6  Info on ESOs in the (2J+1)basis of ground state. (printed elsewhere)
!       7  Full ESO and multiplet matrices (printed elsewhere)
!
      IMPLICIT NONE
      integer i,j,k,n,nc,ipr, i1,i2,j1,j2 
      character map(100,100)*1, FMTR*20, FMTC*20
      
      nc=20; nc=min(nc,TSLJ_f(nelectrons))
      i1=1; i2=nc; j1=1; j2=nc
      if (printMat2(10)) then; i1=printRng2(1); i2=printRng2(2);  &
                               j1=printRng2(3); j2=printRng2(4); nc=printRng2(4)-printRng2(3)+1; endif

      write(iMatrix,'(80("-"),/,"Matrix elements in |SLJ> basis ",/, &
      " Range printed: (",I2,"-",I3," x ",I2,"-",I3,")")')  i1,i2, j1,j2
    
      WRITE(FMTR,'("(",I2,"(    F9.3))")') nc
! 
      select case (ipr)
      CASE(1)
        write(iMatrix,'(80("-"),/," Subset of Full Matrix; multiplied by: ",F12.8)') PFactor
        do i=i1,i2
          write(iMatrix,FMTR) (mat(i,j)*PFactor,j=j1,j2)           
        enddo
      CASE(2)
        n=min(TSLJ_f(nelectrons),100)      
        do i=1,n
          do j=1,n
            map(i,j)='.'
            if (abs(mat(i,j)).gt.1.0D-10) Map(i,j)="*"
          enddo
        enddo
        write(iMatrix,'(80("-"),/,i3,"x",i3," picture of the non-zero elements")') n,n
        do i=1,n
          write(iMatrix,'(100A2)')  (Map(i,j),j=1,n)           
        enddo 
      CASE(3)  ! Uk matrices printed elsewhere
      CASE(4)  ! L,S matrices printed elsewhere
      CASE(5)  ! <L+2S> matrices printed elsewhere
      case default;  write(iMatrix,'("Invalid value for ipr=",I2," in Subroutine PrintMatrix2")') ipr
      END SELECT 
!
      return
      end subroutine printMatrix2
!
!-----------------------------------------------------------------------
!
      SUBROUTINE printMatrix3(ipr,nn)
!
! Prints matrices in full |SLJMJ> basis.
!  ipr = 1  Full matrix (max nc x nc) 
!        2  Full matrix */. (max 100x100)
!        3  non-zero <Ukq> (max 100)
!        4  <Ukq> matrix (max 14x14)
!        5  <ED>  matrix
!        6  <L/S>
!        7  <kL+gS>
!        8  LF in original & PRED basis
!        9  NU
!        10 selection 
!
      IMPLICIT NONE
      integer i,j,k,n,nc,ipr,nn, i1,i2,j1,j2 
      character map(100,100)*1, FMTR*20, FMTC*20

      nc=20; nc=min(nc,nn)
      i1=1; i2=nc; j1=1; j2=nc
      if (printMat3(10)) then; i1=printRng3(1); i2=printRng3(2);  &
                               j1=printRng3(3); j2=printRng3(4); nc=printRng3(4)-printRng3(3)+1; endif
    
      WRITE(FMTR,'("(",I2,"(    F7.1))")') nc
      WRITE(FMTC,'("(",I2,"(2x,2F8.3))")') nc

      select case (ipr)
      CASE(1)
        write(iMatrix,'(80("-"),/," Subset of Full Matrix; multiplied by: ",F12.8)') PFactor
        do i=i1,i2
          if (.not.MatComplex) write(iMatrix,FMTR) (AR(i,j)*PFactor,j=j1,j2)           
          if (     MatComplex) write(iMatrix,FMTC) (AR(i,j)*PFactor,AI(i,j)*PFactor,j=j1,j2)           
        enddo
      CASE(2)
        n=min(nn,100)      
        do i=1,n
          do j=1,n
            map(i,j)='.'
            if (.not.MatComplex.and.abs(AR(i,j)).gt.1.0D-10) Map(i,j)="*"
            if (     MatComplex.and.abs(AR(i,j))+abs(AI(i,j)).gt.1.0D-10) Map(i,j)="*"
          enddo
        enddo
        write(iMatrix,'(80("-"),/,"100x100 picture of the lower-diagonal non-zero elements")')
        do i=1,n
          write(iMatrix,'(100A2)')  (Map(i,j),j=1,n)           
        enddo 
      CASE(3)  ! non-zero Uk printed elsewhere
      CASE(4)  ! Uk matrices printed elsewhere
      CASE(5)  ! ED matrices printed elsewhere
      CASE(6)  ! L,S matrices printed elsewhere
      CASE(7)  ! <L+2S> matrices printed elsewhere
      CASE(8)  ! LF matrices printed elsewhere
      case default;  write(iMatrix,'("Invalid value for ipr=",I2," in Subroutine PrintMatrix3")') ipr
      END SELECT 
!
      return
      end subroutine printMatrix3
!
!-----------------------------------------------------------------------
!
      SUBROUTINE printMatrix4()
!
! This is called at the end of BuildMatrix
! Prints various matrices.
!  ip = 1 Info on ESOs in (2J+1) |MJ> basis. (printed elsewhere) 
!       2 Half ESO matrices in (2J+1) |MJ> basis. (ESO part printed elsewhere) 
!       3 Matrix in uncoupled |SLMlMs> basis.
!       4 Matrix in uncoupled |SLMlMs=S> basis.
!       5 Matrix in uncoupled |Sfi(i=1,7)> basis.
!       6 Lx,Ly,Lz matrices in uncoupled |SLMlMs=S> basis. (printed elsewhere) 
!
!   ???TODO, This subroutine not finished/checked
!
      IMPLICIT NONE
      integer i,j,k,n,nc,nn, i1,i2,j1,j2, Ns 
      integer IT,ITD,II,ISTART,L,basis(2,max_n),nhalf
      real*8 PHASE,XX, Z,D1,D2,BIT, v(max_n)
      character FMT*30
      PARAMETER(Z=0.0D+00, D1=1.0D+00, D2=2.0D+00, BIT=0.001D+00)
      logical doOut/.false./, neven
      
      do i=1,8; if (printMat4(i)) doOut=.true.; enddo
      if (.not.doOut) return

      nc=20; nc=min(nc,n_matrix)
      i1=1; i2=nc; j1=1; j2=nc
      if (printMat4(10)) then; i1=printRng4(1); i2=printRng4(2);  &
                               j1=printRng4(3); j2=printRng4(4); nc=printRng4(4)-printRng4(3)+1; endif    
                             
!          if (n_odd)      write(imatrix,'(20(" |",2i3,"/2> "))') (basis(2,j),basis(1,j),j=j1,j2)
!          if (.not.n_odd) write(imatrix,'(20(" |",2i3,">   "))') (basis(2,j),basis(1,j)/2,j=j1,j2)
      if (n_odd)      WRITE(FMT,'("""<"",2i3,""/2| "",(",I2,"(    F7.1))")') nc
      if (.not.n_odd) WRITE(FMT,'("""  <"",2i3,""| "",(",I2,"(    F7.1))")') nc
!      WRITE(FMTC,'("(",I2,"(2x,2F8.3))")') nc

       ! DO I=1,n_matrix
       !   DO J=1,n_matrix
       !     VR(I,J)=Z
       !   enddo
       ! enddo  
      VR=Z
      
      if (printMat4(2)) then  !     ! 2J+1 even; 2J odd              2J+1 odd; 2J even
        if (mod(ESO2J+1,2).eq.0) then; neven=.true.;  nhalf=(ESO2J+1)/2;  else ;  neven=.false.; nhalf=ESO2J/2; endif 
        write(imatrix,'(/,"(Half) Matrix in full basis",/," Real part of Matrix ")') 
        if (neven) then
          write(imatrix,'(8X,18(" |",i3,"/2> "))') (-ESO2J+2*(j-1),j=1,nhalf)
          do i=1,ESO2J+1;  write(imatrix,'("<",i3,"/2|",(18F9.3))') -ESO2J+2*(i-1),(AR(i,j)*PFactor,j=1,nhalf); enddo  
        else
          write(imatrix,'(8X,18("  |",i3,">  "))') (-ESO2J/2+(j-1),j=1,nhalf)
          do i=1,ESO2J+1;  write(imatrix,'("<",i3,"|",(18F9.3))') -ESO2J/2+(i-1),(AR(i,j)*PFactor,j=1,nhalf); enddo  
        endif  
        write(imatrix,'(/," Imaginary part of Matrix ")') 
        if (neven) then
          write(imatrix,'(8X,18(" |",i3,"/2> "))') (-ESO2J+2*(j-1),j=1,nhalf)
          do i=1,ESO2J+1;  write(imatrix,'("<",i3,"/2|",(18F9.3))') -ESO2J+2*(i-1),(AI(i,j)*PFactor,j=1,nhalf); enddo  
        else
          write(imatrix,'(8X,18("  |",i3,">  "))') (-ESO2J/2+(j-1),j=1,nhalf)
          do i=1,ESO2J+1;  write(imatrix,'("<",i3,"|",(18F9.3))') -ESO2J/2+(i-1),(AI(i,j)*PFactor,j=1,nhalf); enddo  
        endif  
      enDIF  ! printMat4(2)
        
      if (printMat4(3)) then
        ISTART=0
        DO IT=1,n_states
          ITD=(2*TS_Bases(2,IT,nelectrons) + 1)*(TS_Bases(1,IT,nelectrons) + 1)  ! ( 2L + 1)*(2S + 1 )
          IF (IT.NE.1) ISTART=ISTART+(2*TS_Bases(2,IT-1,nelectrons)+1)*(TS_Bases(1,IT-1,nelectrons)+1)
          DO I=1,ITD
           II=I+ISTART
           PHASE=SQRT(DBLE(FullBasis(3,II))+D1)  ! (2J + 1)
            !        (-1)^(S + J - L)
            IF (MOD((FullBasis(1,II)+FullBasis(3,II))/2-FullBasis(2,II),2).GT.0.02D+00) PHASE=-PHASE
            L=ISTART
            DO K=1,TS_Bases(1,IT,nelectrons) + 1  ! (2S + 1)
              DO J=1,2*TS_Bases(2,IT,nelectrons) + 1  ! (2L + 1)
                L=L+1
                basis(1,L)=FullBasis(1,II) ! 2S
                basis(2,L)=FullBasis(2,II) ! L
                VR(II,L)=PHASE*WIG3J(DBLE(FullBasis(2,II)),     DBLE(FullBasis(1,II))/D2,        DBLE(FullBasis(3,II))/D2,  &      !  (  L  S  J  )
                                     DBLE(FullBasis(2,II)-J)+D1,DBLE(FullBasis(1,II)-2*K)/D2+D1,-DBLE(FullBasis(3,II)-2*J)/D2+D1)  !  ( ML MS -MJ )
              enddo ! j                                                                                              
            enddo ! k  
          enddo  ! i
        enddo  !  it
        write(iMatrix,'(80("-"),/,"Full Matrix given in |ML,MS> basis; multiplied by: ",F12.8)') PFactor
        if (printMat4(10)) write(iMatrix,'(" Full matrix size:",I4,"; the (",I3,"-",I3,",",I3"-",I3") subset is printed")') &
                                                                                                                  i1,i2,j1,j2
        if (n_odd)      write(imatrix,'(20(" |",2i3,"/2> "))') (TS_Bases(2,j,nelectrons),TS_Bases(1,j,nelectrons),j=j1,j2)
        if (.not.n_odd) write(imatrix,'(20(" |",2i3,">   "))') (TS_Bases(2,j,nelectrons),TS_Bases(1,j,nelectrons)/2,j=j1,j2)
        Ns=n_matrix
      endif
        
      if (printMat4(4)) then
        ISTART=0
        DO IT=1,n_states
          ITD=(2*TS_Bases(2,IT,nelectrons) + 1)*(TS_Bases(1,IT,nelectrons) + 1)  ! ( 2L + 1)*(2S + 1 )
          IF (IT.NE.1) ISTART=ISTART+(2*TS_Bases(2,IT-1,nelectrons)+1)*(TS_Bases(1,IT-1,nelectrons)+1)
          DO I=1,ITD
            II=I+ISTART
            PHASE=SQRT(DBLE(FullBasis(3,II))+D1)  ! (2J + 1)
            !        (-1)^(S + J - L)
            IF (MOD((FullBasis(1,II)+FullBasis(3,II))/2-FullBasis(2,II),2).GT.0.02D+00) PHASE=-PHASE
            L=ISTART
!            DO K=1,TS_Bases(1,IT,nelectrons) + 1  ! (2S + 1)
              DO J=1,2*TS_Bases(2,IT,nelectrons) + 1  ! (2L + 1)
                L=L+1
                basis(1,L)=FullBasis(1,II) ! 2S
                basis(2,L)=FullBasis(2,II) ! L
                VR(II,L)=PHASE*WIG3J(DBLE(FullBasis(2,II)),     DBLE(FullBasis(1,II))/D2,  DBLE(FullBasis(3,II))/D2,  &       !  (  L   S   J  ) 
                                     DBLE(FullBasis(2,II)-J)+D1,DBLE(FullBasis(1,II))/D2, -DBLE(FullBasis(3,II)-2*J)/D2+D1)   !  ( ML MS=S -MJ )
              enddo ! j
!            enddo ! k  
          enddo  ! i
        enddo  !  it
        write(iMatrix,'(80("-"),/,"Full Matrix given in |ML,MS=S> basis; multiplied by: ",F12.8)') PFactor
        if (printMat4(10)) write(iMatrix,'(" Full matrix size:",I4,"; the (",I3,"-",I3,",",I3"-",I3") subset is printed")') &
                                                                                                                  i1,i2,j1,j2
        if (n_odd)      write(imatrix,'(20(" |",2i3,"/2> "))') (basis(2,j),basis(1,j),j=j1,j2)
        if (.not.n_odd) write(imatrix,'(20(" |",2i3,">   "))') (basis(2,j),basis(1,j)/2,j=j1,j2)
        Ns=TMsS_f(nelectrons)
      endif
        
      if (printMat4(3) .or. printMat4(4)) then
        if (MatComplex) write(imatrix,'("Real part of complex matrix")')  
        ! Mat = AR*VR
        Mat=z
        do i=1,n_matrix
          do j=1,Ns
            do k=1,n_matrix
              Mat(i,j)=Mat(i,j)+AR(i,k)*VR(k,j)
            enddo
          enddo 
        enddo
        do i=1,Ns
          V=Z
          do j=1,Ns
            do k=1,n_matrix
              V(j)=V(j)+VR(k,i)*Mat(k,j)
            enddo
          enddo 
          if (i.ge.i1 .and. i.le.i2) then
            if (n_odd)      write(imatrix,FMT) (basis(2,i),basis(1,i),j=j1,j2),(V(j)*PFactor,j=j1,j2)
            if (.not.n_odd) write(imatrix,FMT) (basis(2,i),basis(1,i)/2,j=j1,j2),(V(j)*PFactor,j=j1,j2)
          endif
        enddo

        if (MatComplex) then
          write(imatrix,'("Imaginary part of complex matrix")')  
          Mat=z
          do i=1,n_matrix
            do j=1,Ns
              do k=1,n_matrix
                Mat(i,j)=Mat(i,j)+AI(i,k)*VR(k,j)
              enddo
            enddo 
          enddo
          do i=1,Ns
            V=Z
            do j=1,Ns
              do k=1,n_matrix
                V(j)=V(j)+VR(k,i)*Mat(k,j)
              enddo
            enddo 
            if (i.ge.i1 .and. i.le.i2) then
              if (n_odd)      write(imatrix,FMT) (basis(2,i),basis(1,i),j=j1,j2),(V(j)*PFactor,j=j1,j2)
              if (.not.n_odd) write(imatrix,FMT) (basis(2,i),basis(1,i)/2,j=j1,j2),(V(j)*PFactor,j=j1,j2)
            endif
          enddo
        endif  ! MatComplex
          
      endif ! printMat4(3) .or. printMat4(4)
!
      return
      end subroutine printMatrix4
!
!-----------------------------------------------------------------------
!
      subroutine calcExpect() 
!  Calculates the expectation values of certain parameters.
!  The eigenvectors are in matrices Mat(Cmat) and the matrix elements 
!  of the parameter operator are load into AR,(AR,AI) 
!
      IMPLICIT NONE
      integer in_e,ie,i,j,k,n, basis(7,Max_N) 
      real*8 d0,d1, v1(Max_N),SaveP(maxTotPar)
      PARAMETER(d0=0.0D+00, d1=1.0D+00)
!
      saveP=P ! save the parameters, before resetting them 
!
      call getBBasis(basis,1,n)
      do in_e=1,n_expect
        p=d0;  p(P_expect(in_e))=d1
        write(io,'(i2,2X,i2,2x,A9,"P:",20F8.1)') in_e,P_expect(in_e),Plab(P_expect(in_e)), &
                                                              (p(P_expect(i)),i=1,n_expect)
        call unloadP() ! only P_expect(in_e) parameter non-zero. 
        CALL buildMatrix(basis,n,3) ! Ipred=3
        do ie = 1, n
          if (.not.MatComplex) then
            do k=1,n
              v1(k)=d0
              do j=1,n
                v1(k)=v1(k)+ar(j,k)*mat(j,ie)
              enddo
            enddo
            Expect(ie,in_e)=d0
            do k=1,n
              Expect(ie,in_e) = Expect(ie,in_e) + mat(ie,k)*v1(k)
            enddo
          else
            write(io,'("Complex matrices not yet implemented in calcExpect")')          
          endif  
        enddo  
      enddo
      P=saveP; call unloadP()   ! restore parameters
!
      return
      end subroutine calcExpect
!     
!-----------------------------------------------------------------------
!
      SUBROUTINE calcSpinAllowed(neng)
!
!  Returns the "spin-allowed"ness of the transitions from [SpAllow(1),SpAllow(2)] 
!  to the higher excited states [SpAllow(2)+1,n_matrix].
!  These values will be summed over both inital and final degeneracies, and then 
!  divided by the total degeneracy. Values are only given once for degenerate states
!  (other values set to -1). The values for the inital levels are set to 1.0
!  When there is no spin-orbit coupling, the spin-allowedness should be exactly 0 or 1.
 
      IMPLICIT none
      COMPLEX*16 VECg ,VECe      
      integer neng, i,j, ig,ie, ISG,ISE, dege,degg
      real*8 Z, Spin(8), sumSA
      PARAMETER(Z=0.0D+00)
!
      DO I=1,n_matrix; SpinAllowed(I)=Z; enddo
      do I=SpAllow(1),SpAllow(2); SpinAllowed(I)=1.0d+00; enddo
! 
      write(IO,'("neng=",I5,";  SpAllow(1),SpAllow(2)=",2I4)') neng,SpAllow(1),SpAllow(2)
      write(io,'(/,"SpinAllowed:")')
      write(io,'((20F8.3))') (SpinAllowed(j), j=1,n_matrix)

!  Loop for every transition:
      DO IE=SpAllow(2)+1,neng
        DO IG=SpAllow(1),SpAllow(2)
          do I=1,8; Spin(I)=z; enddo
          DO I=1,n_matrix
            ISG=  FullBasis(1,I) + 1  !  2S + 1 
            IF (.not.MatComplex) VECg = DCMPLX(MAT(I,IG),Z)
            IF (     MatComplex) VECg = CMAT(I,IG)
            DO J=1,n_matrix
              ISE=  FullBasis(1,J) + 1  !  2S + 1 
              if (ISG.eq.ISE) then
                IF (.not.MatComplex) VECe = DCMPLX(MAT(J,IE),Z)
                IF (     MatComplex) VECe = CMAT(J,IE)
                 Spin(ISG)=Spin(ISG)+DBLE(DCONJG(VECg)*VECe)
              endif
            enddo  ! J
          enddo  ! I
        enddo  ! IG
        do I=1,8; SpinAllowed(IE) = SpinAllowed(IE) + SPIN(I)**2; enddo
      enddo  ! IE
      write(io,'(/,"SpinAllowed:")')
      write(io,'((20F8.3))') (SpinAllowed(j), j=1,n_matrix)
!  Sum degeneracies:
      degg=SpAllow(2)-SpAllow(1)
      i=SpAllow(1); dege=1; sumSA=SpinAllowed(IG)
      DO IG=SpAllow(1),n_matrix-1
        
        if (abs(engs(ig+1)-engs(ig)).lt.Edegen) then
          sumSA=sumSA+SpinAllowed(ig); dege=dege+1
        else
          if (dege.eq.1) then 
            SpinAllowed(ig)=sumSA/dble(degg)
          else
            do j=i,ig
              if (IG.le.SpAllow(2)) SpinAllowed(j)=sumSA/dble(degg) 
              if (IG.gt.SpAllow(2)) SpinAllowed(j)=sumSA/dble(degg*dege) 
            enddo
            dege=1; sumSA=z; i=ig
          endif
        endif
        write(IO,'("IG,dege,engs(ig),sumSA=",2I4,F10.1,F8.3)') IG,dege,engs(ig),sumSA

      enddo   ! IG
!
      write(io,'(/,"SpinAllowed:")')
      write(io,'((20F8.3))') (SpinAllowed(j), j=1,n_matrix)

!
      RETURN
      END Subroutine calcSpinAllowed
!
!-----------------------------------------------------------------------
!
      SUBROUTINE calcSpinPR(neng)
!
!  Returns the spin projection of each wavefunction in SPIN.
!
      IMPLICIT none
      COMPLEX*16 VEC      
      integer neng,i,j,ii,ISD,ISP
      real*8 Z, SMALL, PROJ
      PARAMETER(Z=0.0D+00, SMALL=1.0D-12)
!
      DO I=1,n_matrix
        DO J=1,4       ! Each state can have up to 4 different spin components
          SPIN(I,J)=Z
        enddo
      enddo
! 
!  Loop for every eigenvector:
      DO II=1,neng
        DO I=1,n_matrix
          ISD=  FullBasis(1,I) + 1  !  2S + 1 
          IF (.not.MatComplex) VEC = DCMPLX(MAT(I,II),Z)
          IF (     MatComplex) VEC = CMAT(I,II)
          IF (ABS(VEC).gt.SMALL) then
            PROJ = DBLE(VEC*DCONJG(VEC))
!  Find the spin projections:
            DO ISP=1,NSPIN
              IF (ISD.EQ.S_mult(ISP)) SPIN(II,ISP)=SPIN(II,ISP)+PROJ
            enddo
          endif ! VEC>small
        enddo ! i
!        write(Idebug,'("SUBROUTINE SPINPR: SPIN=",4F16.12)') (SPIN(II,ISP),ISP=1,NSPIN)
!
      enddo ! ii
!
      RETURN
      END Subroutine calcSpinPR
!
!-----------------------------------------------------------------------

END MODULE f_e_calculate
!  2031