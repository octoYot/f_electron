MODULE f_e_fit !  

USE f_e_data
USE f_e_eso
USE f_e_parameters
USE f_e_calculate
USE f_e_edipole
USE f_e_magnetics

!
!  This module contains the fitting routines.
!

IMPLICIT NONE

PUBLIC 

CONTAINS
!
!------------------------------------------------------------------------------
!
      SUBROUTINE prepareFit()
!  This prepares the matrices for the fitting subroutine.
!  fitType: 1 Eng/Int;   2 J-O Eng/Int;  3 Bkq values; 4 EXSO MEs
!           5 g-values;  6 1e MEs
!  If a prediagonalised (PRED) basis is used, most of this is bypassed.
!
      implicit none
      INTEGER i,j,ifit,n,ib, basis(7,Max_N),ii,n1, allocation_status,nr,ni 
      real*8 saveP(maxTotPar)
      character(len=99) :: emsg
      logical op3,op9,c18
!
      op3=outp(3); op9=outp(9); c18=check1(8); outp(3)=.false.; outp(9)=.false.; check1(8)=.false.
! stop unnecessary printing in AOMmatrixD/F via call to prepCalc

      call loadP()
      if (fitType.eq.3 .or. fitType.eq.4 .or. fitType.eq.6 .or. nPRED.ne.0) goto 100 ! Fits that don't need the matrix diagonalization
      
      do i=1,maxTotPar; saveP(i) = P(i); enddo  ! Save original parameters
!
      if (linearFit) then ! All parameters are linear.

        n1 = (max_N+24)/nblocks  ! nblocks=1,2,3,4,6 divides evenly into max_N=(3432+24)
! The +24 is because for f7, nblocks=6, the blocks are 2x570+2x572+2x574; so 3432/6=572 is too small.
! Note for MJ blocking, nblocks set to 6.
        if (FitOpt(1)) then 
          write(io,'("FitOpt(1)=T; Fitting info will be printed")')
!          write(io,'(" Attempting to ALLOCATE nblocks=",I2,", n1=",I4,", n1=",I4)') nblocks,n1,n1
!          write(*, '(" Attempting to ALLOCATE nblocks=",I2,", n1=",I4,", n1=",I4)') nblocks,n1,n1
        endif        
        allocate (ARfit(nblocks,n1,n1), stat=allocation_status, errmsg=emsg)
        if (allocation_status > 0) then
          write(io,'("***FATAL: error Allocating ARfit:",A)') trim(emsg)
          write(*, '("***FATAL: error Allocating ARfit:",A)') trim(emsg)
          stop
        end if
        allocate (AIfit(nblocks,n1,n1), stat=allocation_status, errmsg=emsg)
        if (allocation_status > 0) then
          write(io,'("***FATAL: error Allocating AIfit:",A)') trim(emsg)
          write(*, '("***FATAL: error Allocating AIfit:",A)') trim(emsg)
          stop
        end if
        if (FitOpt(1)) then
          write(io,'(2X,"Matrices ARfit,AIfit memory ALLOCATE: (",i2,",",i4,",",i4,") total:",i8)') nblocks,n1,n1,nblocks*n1**2
!          write(*, '(2X,"Matrices ARfit,AIfit memory ALLOCATE: (",i2,",",i4,",",i4,") total:",i8)') nblocks,n1,n1,nblocks*n1**2
!          write(io,'(2X,"ARfit,AIfit SHAPE =",(10I6))') shape(ARfit)
!          write(*, '(2X,"ARfit,AIfit SHAPE =",(10I6))') shape(ARfit)
        endif
        n1=maxNF/nblocks ! nblocks=1,2,3,4,6 divides evenly into maxNF=630000  ! Note for MJ blocking, nblocks set to 6.
        allocate (FitMatR(NFit,nblocks,n1), stat=allocation_status, errmsg=emsg)
        if (allocation_status > 0) then
          write(io,'("***FATAL: error Allocating FitMatR:",A)') trim(emsg)
          write(*, '("***FATAL: error Allocating FitMatR:",A)') trim(emsg)
          stop
        end if
        allocate (iIndexR(NFit,nblocks,n1), stat=allocation_status, errmsg=emsg)
        if (allocation_status > 0) then
          write(io,'("***FATAL: error Allocating iIndexR:",A)') trim(emsg)
          write(*, '("***FATAL: error Allocating iIndexR:",A)') trim(emsg)
          stop
        end if
        allocate (jIndexR(NFit,nblocks,n1), stat=allocation_status, errmsg=emsg)
        if (allocation_status > 0) then
          write(io,'("***FATAL: error Allocating jIndexR:",A)') trim(emsg)
          write(*, '("***FATAL: error Allocating jIndexR:",A)') trim(emsg)
          stop
        end if
        if (FitOpt(1)) then
          i=NFit*nblocks*n1
          write(io,'(2X,"Matrices FitMatR,iIndexR,jIndexR memory ALLOCATE: (",i2,",",i2,",",i6") total:",i8)') NFit,nblocks,n1,i
!          write(io,'(2X,"FitMatR,iIndexR,jIndexR SHAPE =",(10I6))') shape(FitMatR)
        endif  
        if (MatComplex) then
          allocate (FitMatI(NFit,nblocks,n1), stat=allocation_status, errmsg=emsg)
          if (allocation_status > 0) then
            write(io,'("***FATAL: error Allocating FitMatI:",A)') trim(emsg)
            write(*, '("***FATAL: error Allocating FitMatI:",A)') trim(emsg)
            stop
          end if
          allocate (iIndexI(NFit,nblocks,n1), stat=allocation_status, errmsg=emsg)
          if (allocation_status > 0) then
            write(io,'("***FATAL: error Allocating iIndexI:",A)') trim(emsg)
            write(*, '("***FATAL: error Allocating iIndexI:",A)') trim(emsg)
            stop
          end if
          allocate (jIndexI(NFit,nblocks,n1), stat=allocation_status, errmsg=emsg)
          if (allocation_status > 0) then
            write(io,'("***FATAL: error Allocating jIndexI:",A)') trim(emsg)
            write(*, '("***FATAL: error Allocating jIndexI:",A)') trim(emsg)
            stop
          end if
          if (FitOpt(1)) then
            i=NFit*nblocks*n1
            write(io,'(2X,"Matrices FitMatI,iIndexI,jIndexI memory ALLOCATE: (",i2,",",i2,",",i6") total:",i8)') NFit,nblocks,n1,i
!            write(io,'(2X,"FitMatI,iIndexI,jIndexI SHAPE =",(10I6))') shape(FitMatI)
          endif  
        endif
!       
        do ib = 1,nblocks
          call getBBasis(basis,ib,n)  ! n and basis are returned        
!  Build matrix with fitted parameters all zero.
          do i=1,maxTotPar; P(i)=saveP(i); enddo
          do ifit=1,nFit
            P(FitN(ifit)) = 0.0d+00
            if (NLinks.gt.0) then
              do j=1,NLinks
                if (nlink(2,j).eq.FitN(ifit)) then
                  P(nlink(1,j)) = 0.0d+00
                endif
              enddo
            endif  
          enddo
          
          call prepCalc(0)        
          call makeFitMat(ifit,ib,n,basis,0,1)
!
!  Make a compressed matrix for each fit parameter.
          do ifit=1,nFit   
            do i=1,80; P(i)=0.0d+00; enddo  ! Want theta,phi,chi to take their correct values.
            do i=141,450; P(i)=0.0d+00; enddo
            P(FitN(ifit)) = 1.0d+00
            call prepCalc(FitN(ifit))        
!          call PrintParameters(1)
            call makeFitMat(ifit,ib,n,basis,0,2)
            if (NLinks.gt.0) then
              do j=1,NLinks
                if (nlink(2,j).eq.FitN(ifit)) then
                  do i=1,80; P(i) = 0.0d+00; enddo
                  P(nlink(1,j)) = 1.0d+00
                  call unloadP()
                  if ((LFtype.eq."AOM ".or. LFtype.eq."LF1e") .and. FitN(ifit).gt.20 .and. FitN(ifit).le.140) then
                    if (Lvalue.eq.2) call AOMmatrixD()
                    if (Lvalue.eq.3) call AOMmatrixF()
                  endif 
!          call PrintParameters(1)
                  call makeFitMat(ifit,ib,n,basis,j,3)
                endif
              enddo
            endif
          enddo ! ifit
!
          call flush(IO)
          if (FitOpt(1)) then
            n1=min(n,14)
            write(io,'(100("-"),/,"Matrix of nonvarying parameters")')
            do i=1,n1
              write(io,'((14F10.1))') (ARfit(ib,i,j),j=1,n1)           
            enddo
            do ifit=1,nFit   
              do i=1,n1; do j=1,n1; mat1(i,j)=0.0d+00; enddo; enddo          
              do i=1,nMatFR(ifit,ib)
                mat1(iIndexR(ifit,ib,i),jIndexR(ifit,ib,i)) =  FitMatR(ifit,ib,i)    ! unpack
              enddo
              write(io,'(/,100("-"),/,"Real Matrix of parameter=",i3,2X,I3,2X,A)') ifit,FitN(ifit),Plab(FitN(ifit))
              do i=1,n1
                Write(io,'((14F10.4))') (mat1(i,j),j=1,n1)           
              enddo
              
              if (nMatC(ifit)) then
                do i=1,n1; do j=1,n1; mat1(i,j)=0.0d+00; enddo; enddo          
                do i=1,nMatFI(ifit,ib)
                   mat1(iIndexI(ifit,ib,i),jIndexI(ifit,ib,i)) =  FitMatI(ifit,ib,i)    ! unpack
                enddo
                write(io,'(/,100("-"),/,"Imag Matrix of parameter=",i3,2X,I3,2X,A)') ifit,FitN(ifit),Plab(FitN(ifit))
                do i=1,n1
                  Write(io,'((14F10.4))') (mat1(i,j),j=1,n1)           
                enddo
              endif  
            enddo  
            write(io,'(/,"Finished writing FIT matrices.")')

          endif ! FitOpt(1)
!        
        enddo ! ib
      endif ! linearFit

      if (.not.linearFit) write(io,'("The fit is using non-linear parameter/variables")') 
!
!  reset starting parameters & load initial parameters to be fitted.
      do i=1,maxTotPar; P(i)=saveP(i); enddo
!      call unloadP()   ! puts P back into FI, Bkq, AOM parameter

 100  if (fitType.eq.3 .and. (abs(RotLF(1))+abs(RotLF(2))+abs(RotLF(3))).gt.1.0d-12) then
        do i=1,16; savedBKQR(i)=BKQR(i); savedBKQI(i)=BKQI(i); enddo
      endif
            
!      if (FitType.eq.6) then
!        WRITE(IO,'(115("-"))'); 
!        write(io,'(" The 7x7 real 1-electron f-orbital matrix elements to be fitted (x indicates zero weighting)")') 
!        write(io,'(14X,"|sig>       |piS>       |piC>       |delS>      |delC>      |psiS>      |psiC>  ")') 
!        write(io,'(2X," <sig|",7(F11.2,A1))')  EXP_E( 1),  WGTchar( 1)
!        write(io,'(2X," <piS|",7(F11.2,A1))') (EXP_E( 1+i),WGTchar( 1+i),i=1,2)
!        write(io,'(2X," <piC|",7(F11.2,A1))') (EXP_E( 3+i),WGTchar( 3+i),i=1,3)
!        write(io,'(2X,"<delS|",7(F11.2,A1))') (EXP_E( 6+i),WGTchar( 6+i),i=1,4)
!        write(io,'(2X,"<delC|",7(F11.2,A1))') (EXP_E(10+i),WGTchar(10+i),i=1,5)
!        write(io,'(2X,"<psiS|",7(F11.2,A1))') (EXP_E(15+i),WGTchar(15+i),i=1,6)
!        write(io,'(2X,"<psiC|",7(F11.2,A1))') (EXP_E(21+i),WGTchar(21+i),i=1,7)
!      endif

      call prepCalc(0)        
      do ifit=1,nFit; fitP(ifit)=P(FitN(ifit)); enddo  
!      
      write(io,'(/,"Parameters to be fitted:  (maxIts=",I4,", tol=",E12.4,")")') maxIts,FitTol
      write(io,'("    Number Parameter   Value       Low        High  ")') 
      do i=1,nFit
        write(io,'(I2,2X,I3,2X,A9,2A1,3F11.2)') i,FitN(i),Plab(FitN(i)),V1(FitN(i)),V2(FitN(i)),  &
                                                  fitP(i),FitL(i),FitH(i)
      enddo  
      outp(3)=op3; outp(9)=op9; check1(8)=c18
      call flush(IO)
!
      return
!
      END SUBROUTINE prepareFit
!
!------------------------------------------------------------------------------
!
      SUBROUTINE prepareFitJO()
!   This prepares for fitting the JO parameters.
!
      implicit none
      INTEGER i,ifit 
!
      call loadP()
      call prepCalc(0)        
      do ifit=1,nFit; fitP(ifit)=P(FitN(ifit)); enddo  
!      
      write(io,'(/,"Parameters to be fitted:  (maxIts=",I4,", tol=",E12.4,")  FitOpt=",8L2)') &
                                                                  maxIts,FitTol,(FitOpt(i),i=1,8)
      write(io,'("    Number Parameter   Value       Low        High  ")') 
      do i=1,nFit
        write(io,'(I2,2X,I3,2X,A9,2A1,3F11.2)') i,FitN(i),Plab(FitN(i)),V1(FitN(i)),V2(FitN(i)),  &
                                                  fitP(i),FitL(i),FitH(i)
      enddo  
!
      return
!
      END SUBROUTINE prepareFitJO
!
!------------------------------------------------------------------------------
!
      SUBROUTINE makeFitMat(ifit,ib,n,basis,ilink,ii)
!
!  The subroutine makes the symmetry adapted matrices
!  ii=1 Builds the matrix with the parameters to be fitted set to zero. 
!       This full matrix contains all of the things not being changed in the fit. 
!  ii=2 Builds the separate matrices for each of the fitted parameters
!       which are "linear" functions.
!       During the fit these matrices can just be multiplied by the parameter being varied.
!  ii=3 A parameter P2 (which is being fitted) is linked to another parameter P1.
!         P1 = a*P2 + b  
!       The parameter P2 has a matrix Mat2. It will now have the matrix (Mat2 + a*Mat1), 
!       where Mat1 is the matrix for P1. The matrix of non-varying parameters (ii=1) will 
!       have b*Mat1 added to it. Note that P1 cannot be varied.
!  ii=4 A variable is being fitted (values stored in the i=401-450 of the P(i) matrix), and 
!       this variable could be used in the definition of other parameters P(i<401) using formula.
!       Note that these parameters (which depend on a fitted variable) themselves cannot be varied.
!  ***ii=4 not used.
!
!  Each separate matrix for the parameters being fitted are stored in a sparse storage format.
!
!  Before each call to this subroutine, the appropriate values for the parameters 
!  and the basis must be set.
!
!  Note that each parameter cannot have both real and imaginary elements.  
!
      implicit none
      INTEGER ifit,ib,n,ilink,ii, i,j,kr,ki, basis(7,Max_N)
      real*8 B
      
      CALL buildMatrix(basis,n,0)
!      write(io,'("II=",i2"; sigma(1-6),piX(1-6),sigma=",2(6F7.1,2X),F7.1)') ii,(eSigma(i),i=1,6), &
!          (ePiX(i),i=1,6),variablesvalues(1)
      if (ii.eq.1) then
        do i=1,n; do j=1,i; ARfit(ib,i,j)=AR(i,j); enddo; enddo
        if (MatComplex) then
          do i=1,n; do j=1,i; AIfit(ib,i,j)=AI(i,j); enddo; enddo
        endif
        return
      endif  !   ii=1
      
      if (ii.eq.2) then
        if (FitN(ifit).eq.1) return  ! Parameter is EAVE
        kr=0; ki=0
        do i=1,n
          do j=1,i
            if (AR(i,j).ne.z) then 
              kr=kr+1
              FitMatR(ifit,ib,kr) = AR(i,j)
              iIndexR(ifit,ib,kr) = i
              jIndexR(ifit,ib,kr) = j 
            endif
           enddo ! j
        enddo ! i
        
        if (ii.eq.2 .and. MatComplex) then
          do i=1,n
            do j=1,i
            if (AI(i,j).ne.z) then 
                ki=ki+1
                FitMatI(ifit,ib,ki) = AI(i,j)
                iIndexI(ifit,ib,ki) = i
                jIndexI(ifit,ib,ki) = j
              endif
            enddo ! j
          enddo ! i
        endif  
 !       write(io,'("**** ifit,kr,ki=",3I4)') ifit,kr,ki
        if (kr.eq.0 .and. ki.eq.0) then
          write(io,'("***FATAL: Parameter to be fitted:",I4," has no elements in matrix, kr,ki=0,0")') ifit
          STOP
        elseif (kr.ne.0 .and. ki.eq.0) then 
          nMatFR(ifit,ib)=kr; nMatC(ifit)=.false.  ! just real
        elseif (ki.ne.0) then
          nMatFR(ifit,ib)=kr; nMatFI(ifit,ib)=ki; nMatC(ifit)=.true.  ! Complex
        endif
!
        if (FitOpt(1)) then
          if (nMatC(ifit)) then
            write(io,20) ib,ifit,FitN(ifit),PLab(FitN(ifit)),nMatFR(ifit,ib),nMatFI(ifit,ib)
          else  
            write(io,21) ib,ifit,FitN(ifit),PLab(FitN(ifit)),nMatFR(ifit,ib)
          endif  
 20       format(2X,"Block:",I2," Parameter to be fitted:",I2,I4,2X,A9,"has",I6," real and",I6," imaginary non-zero elements")
 21       format(2X,"Block:",I2," Parameter to be fitted:",I2,I4,2X,A9,"has",I6," real non-zero elements")
        endif ! (FitOpt(1))
      endif  !  ii=2
      
      if (ii.eq.3) then
        do i=1,n; do j=1,i; mat(i,j)=0.0d+00; enddo; enddo          
        do i=1,nMatFR(ifit,ib)
          mat(iIndexR(ifit,ib,i),jIndexR(ifit,ib,i)) =  FitMatR(ifit,ib,i)    ! unpack
        enddo
        if (nMatC(ifit)) then
          do i=1,n; do j=1,i; mat1(i,j)=0.0d+00; enddo; enddo          
          do i=1,nMatFI(ifit,ib)
            mat1(iIndexI(ifit,ib,i),jIndexI(ifit,ib,i)) =  FitMatI(ifit,ib,i)    ! unpack
          enddo
        Endif
        
        B=rlink(2,ilink)
! Add the B from the link to the original (non-varying) matrix. 
        do i=1,n
          do j=1,i
            if (AR(i,j).ne.z) then 
              ARfit(ib,i,j) = ARfit(ib,i,j) + B*AR(i,j)       ! the b in P1 = a*p2 + b
              mat(i,j) = mat(i,j) + rlink(1,ilink)*AR(i,j)    ! the a 
            endif
            if (AI(i,j).ne.z) then 
              AIfit(ib,i,j) = AIfit(ib,i,j) + B*AI(i,j)
              mat1(i,j) = mat1(i,j) + rlink(1,ilink)*AI(i,j) 
            endif
          enddo ! j
        enddo ! i      
        kr=0
        do i=1,n
          do j=1,i
            if (mat(i,j).ne.z) then 
              kr=kr+1
              if (kr.gt.maxNF) then
                write(io,'("***FATAL: kr=",I6,", greater than maxNF=",I6," in Subroutine makeFitMat")') kr,maxNF
                write(io,'("You need either a use a higher symmetry, or to recompile the program with a larger maxNF")')
                stop
              endif
              FitMatR(ifit,ib,kr) = mat(i,j)    ! 
              iIndexR(ifit,ib,kr) = i
              jIndexR(ifit,ib,kr) = j 
            endif
            if (mat1(i,j).ne.z) then 
              ki=ki+1
              if (ki.gt.maxNF) then
                write(io,'("***FATAL: ki=",I6,", greater than maxNF=",I6," in Subroutine makeFitMat")') ki,maxNF
                write(io,'("You need either a use a higher symmetry, or to recompile the program with a larger maxNF")')
                stop
              endif
              FitMatI(ifit,ib,kr) = mat1(i,j)    ! 
              iIndexI(ifit,ib,kr) = i
              jIndexI(ifit,ib,kr) = j 
            endif
          enddo ! j
        enddo ! i
        nMatFR(ifit,ib)=kr
        nMatFI(ifit,ib)=ki
!
        if (FitOpt(1)) write(io,30) nlink(1,ilink),FitN(ifit),nMatFR(ifit,ib),nMatFI(ifit,ib)
 30     format(" Adding linked parameter:",I3,"; the parameter",I3," to be fitted now has",I6," real, and",  &
                                                                                           I6," imaginary non-zero elements")
      endif  ! ii=3
!      
      return
!
      END SUBROUTINE makeFitMat
!
!------------------------------------------------------------------------------
!
      subroutine fcn(m,n,x,fvec,iflag)
! ----------
! Calculate the function at x and return this vector in fvec.
! ----------
! Input:
! m       The number of experimental values to be used.
! n       The number of parameters to be varied. The parameters are in X(N) 
! x(n) =fitP   The n parameters. 
! Output:
! fvec(M) The m calculated values corresponding to the experimental values.
! iflag   Don't change, but set to negative if you wish to terminate. 
!
!  If JuddOfelt .true. assume a fit for JO intensities.
!
! fitType=1   Fit Eng/Int 
! fitType=2   Fit Judd-Ofelt Eng/Int 
! fitType=3   Fit AOM parameters to Crystal Field Coefficients
! fitType=4   Fit calculated MEs to ESO MEs
! fitType=5   Fit g-values
! fitType=6   Fit calculated MEs to one electron LF MEs
! 
      implicit none
      integer m,mm,n,iflag,m1
      real*8 x(n),fvec(m),fnorm,addInt
! 
      INTEGER i,j,j1,k,ifit,ib,Nsum,nn, basis(7,Max_N),nb,neng
      real*8 E(max_N),penalty 
      character line1*120, line2*120, bit1*12,bit2*12, CL(maxFit)*1,CH(maxFit)*1
      logical AtomChanged, LFChanged, OffsetOnly, cV
!
      Its=Its+1;  AtomChanged=.False.; LFChanged=.False.
      penalty=0.0d+00  ! To constrain fit parameters to be within certain limits.
      OffsetOnly=.true.; cV=calcVecs
      calcVecs=.false.; if (fitType.eq.5) calcVecs=.true.
!      
      do ifit = 1,n  ! = Nfit
        if (offsetType.ne.2 .and. FitN(ifit).ne.1)  OffsetOnly=.false.
        if (offsetType.eq.2  .and. FitN(ifit).ne.1 .and. .not.(FitN(ifit).ge.351.and.FitN(ifit).le.400)) OffsetOnly=.false.
        CL(ifit)=" "; CH(ifit)=" "
        if (x(ifit).lt.FitL(ifit)) then
          CL(ifit)="*"
          if (FitOpt(6)) x(ifit)=FitL(ifit)  !  FitOpt(6) means hard limits
          if (.not.FitOpt(6)) penalty=penalty+100.0d+00*(x(ifit)-FitL(ifit))**2
        endif  
        if (x(ifit).gt.FitH(ifit)) then
          CH(ifit)="*"
          if (FitOpt(6)) x(ifit)=FitH(ifit)
          if (.not.FitOpt(6)) penalty=penalty+100.0d+00*(x(ifit)-FitH(ifit))**2 
        endif  
        if ((FitN(ifit).ge.2 .and. FitN(ifit).le. 20) .or.                                         &
           (offsetType.eq.1.and.FitN(ifit).ge.351.and.FitN(ifit).le.400)) AtomChanged=.True.
        if (FitN(ifit).ge.21.and. FitN(ifit).le.140) LFChanged=.True.
      enddo
      if (penalty.ne.0.0d+00) write(IO,'("****penalty=",E12.6,(10(2X,F9.3,2(F9.1,A1))))') penalty, &
                                               (x(ifit),FitL(ifit),CL(ifit),FitH(ifit),CH(ifit),ifit=1,n) 
      if (penalty.ne.0.0d+00) write(* ,'("****penalty=",E12.6,(10(2X,F9.3,2(F9.1,A1))))') penalty, &
                                               (x(ifit),FitL(ifit),CL(ifit),FitH(ifit),CH(ifit),ifit=1,n) 
!
!      write(IO,'("fitType,linearFit,nPRED=",I2,L2,I4,";  its,AOMchanged=",I3,L2)') fitType,linearFit,nPRED,  &
!                                                                                                 its,AOMchanged
!     FitType = 1 Normal Energy Fit or FitType = 5  g-value fit (both need matrix diagonalisation)
!     ***********                      ***********
      if (fitType.eq.1 .or. fitType.eq.5) then 

        if (.not.linearFit .or. nPRED.ne.0) then
          do ifit=1,n  ! =nFit   
            P(FitN(ifit)) = x(ifit)
          enddo
          call prepCalc(0)        
!          write(io,'("P(62),P(121),theta(2)=",3F12.8," ",A4,"=",F12.4)')P(62),P(121),theta(2),BkqI_label(3),BkqI(3)
        endif
!        
        if (nPRED.ne.0) then; neng=nPRED; else; neng=n_matrix; endif
        if (its.gt.1 .and. OffsetOnly) then   ! only EAVE/OFFSETs are being varied. Do quick a quick fit by missing diagonalisation; use previous eigenvalues.
          do i=1,neng; engs(i)=save_engs(i); enddo
          goto 100
        endif
!
        if (nPRED.ne.0) then
          if (its.eq.1 .or. AtomChanged) call doPRED1()  ! Atomic parameters changed
          if (its.eq.1 .or. LFChanged)   call doPRED2()  ! Ligand Field parameters only using the full |LSJMJ> basis; save LF 
          call buildAndDiagonalisePRED()
        else
          Nsum=0
          do ib = 1,Nblocks
            nn=Nblock(ib)
            if (nn.le.0) goto 90            
            if (.not.linearFit) then
              call getBBasis(basis,ib,nb)
              if (nb.ne.nn) then; WRITE(io,'("Sanity fail nn,nb="2i4)') nn,nb; STOP; endif;
              CALL buildMatrix(basis,nn,0)
            else  ! Linear: The fit can be decomposed into a sum of matrices for each fit parameter.
!            write(io,'("**nn=",I4,",  ib=",i4,",  nMatFR(1,ib)=",i4)') nn,ib,nMatFR(1,ib)
              do i=1,nn; do j=1,nn; AR(i,j)=0.0D+00; enddo; enddo
              do i=1,nn; do j=1,i; AR(i,j)=ARfit(ib,i,j); enddo; enddo
              if (MatComplex) then
                do i=1,nn; do j=1,nn; AI(i,j)=0.0D+00; enddo; enddo
                do i=1,nn; do j=1,i; AI(i,j)=AIfit(ib,i,j); enddo; enddo
              endif
              do ifit = 1,n  ! = Nfit
!            write(IO,'("ifit,ib,fitN(ifit),nMatFR(ifit,ib)=",3I3,I6)') ifit,ib,fitN(ifit),nMatFR(ifit,ib)
                IF (fitN(ifit).eq.1) EAVE=x(ifit)   ! So CALL zeroE() will work. 
                do i = 1,nMatFR(ifit,ib)
                  AR(iIndexR(ifit,ib,i),jIndexR(ifit,ib,i)) = AR(iIndexR(ifit,ib,i),jIndexR(ifit,ib,i))+x(ifit)*FitMatR(ifit,ib,i)
                enddo ! i
                if (nMatC(ifit)) then
                  do i = 1,nMatFI(ifit,ib)
                    AI(iIndexI(ifit,ib,i),jIndexI(ifit,ib,i)) = AI(iIndexI(ifit,ib,i),jIndexI(ifit,ib,i))+x(ifit)*FitMatI(ifit,ib,i)
                  enddo ! i
                endif
              enddo ! ifit
            endif ! .not.linearFit 
!          write(io,'("**fcn: x=",10F9.2)') (x(ifit),ifit=1,n)
            if (FitOpt(5)) then
              write(iMatrix,'("Full Fit Matrix to be diagonalised")')
              call printMatrix3(1,nn)
            endif
!          write(io,'("**MatComplex,CFcomplex=",2L2)') MatComplex,CFcomplex
            CALL diagonalise(nn,e)
!        write(io,'("**nsum=",I4)') nsum
            do i=1,nn; engs(i+Nsum)=e(i); enddo
!            
            if (calcVecs) then
              if (nblocks.eq.1) then
                if (MatComplex) then
                  do i=1,nn; do j=1,nn; CMAT(i,j)=DCMPLX(AR(i,j),AI(i,j)); enddo; enddo
                else
                  do i=1,nn; do j=1,nn; MAT(i,j)=AR(i,j); enddo; enddo
                endif      
              else  ! nblocks >1
                if (MatComplex) then
                  do i=1,nn; do j=1,nn; CMAT(Bbasis(ib,7,i),j+Nsum)=DCMPLX(AR(i,j),AI(i,j)); enddo; enddo
                else ! .not.MatComplex  
                  do i=1,nn; do j=1,nn; MAT(Bbasis(ib,7,i),j+Nsum)=AR(i,j); enddo; enddo
                endif      
              endif ! nblocks 
            endif ! calcVecs  
!    
 90         nsum=nsum+nn
          enddo ! ib
          if (nsum.ne.n_matrix) write(io,'("***Warning: nsum=",I4," not equal n_matrix=",I4)') nsum,n_matrix
          if (nBlocks.ne.1) call orderE_V(calcVecs) 
        endif  ! nPRED  
!  Offset
 100    continue
        if (offsetType.eq.2) then
          DO i=1,neng  !  n_matrix
            do j=1,offsetN
             if (i.ge.ioffset(j,1) .and. i.le.ioffset(j,2)) Engs(i)=Engs(i)+Roffset(j)
            enddo
          enddo  
          if (nBlocks.ne.1) call orderE_V(calcVecs)
        endif 
        if (its.eq.1) then
          do i=1,neng; save_engs(i)=engs(i); enddo
          if (EaveOnly) save_EAVE=x(1) 
        endif
        call zeroE(engs,neng)   ! make relative to E0 if option(1)
        if (fitType.eq.5) goto 110 
!
!        write(io,'("CalcE:    ",10F9.1,/,(10X,10F9.1))') (engs(i),i=1,min(m,100))  
        if (EaveOnly) then
!          write(io,'("EaveOnly:",L1,";  Eave=",F14.8,"; save_EAVE=",F14.8)') EaveOnly,x(1),save_EAVE
          DO i=1,neng; Engs(i)=Engs(i)+x(1)-save_EAVE; enddo  
        endif 
!
        do i=1,m  ! m=Nexp
          if (assignMethod.eq.0) then
            fvec(i) = Wgt_E(i)*(EXP_E(i)-engs(NAssign(i,1)))
            if (Nass(i).gt.1) then
              do j=1,Nass(i)
!            write(io,'("i,j=",2I4,",  NAssign(i,j)=",I4,", ExpE=",F10.2,", calcE(Nassign)=",F10.2 )') &
!                                                            i,j,NAssign(i,j),EXP_E(i),engs(NAssign(i,j))
                fvec(i) = fvec(i) + Wgt_E(i)*(EXP_E(i)-engs(NAssign(i,j)))
              enddo
            endif
          else if (assignMethod.eq.1) then  ! Assign to levels in a CQN block.
            k=0
            do j=1,neng
              if (iCQN(j).eq.NAssign(i,2)) k=k+1
              if (k.eq.NAssign(i,1)) goto 105
            enddo
            write(io,'("Level:",i3," of block",i2," not found.")') NAssign(i,1),NAssign(i,2)
            stop
 105        fvec(i)=Wgt_E(i)*(EXP_E(i)-engs(k))           
          else if (assignMethod.eq.2) then  ! Assign to levels to spin.
            k=0
            do j=1,neng
              do j1=1,NSPIN  ! SPIN(J,K),K=1,NSPIN
                if (SPIN(j,j1).eq.NAssign(i,2)) k=k+1
                if (k.eq.NAssign(i,1)) goto 106
              enddo  ! j1
            enddo  ! j
            write(io,'("Level:",i3," of spin",i2," not found.")') NAssign(i,1),NAssign(i,2)
            stop
 106        fvec(i)=Wgt_E(i)*(EXP_E(i)-engs(k))          
          endif  ! assignMethod
        enddo  ! m
        fnorm = enorm(m,fvec)/sqrt(dble(m))
!        write(*,'(5X,"Calc:",10F12.4,/,(10X,10F9.1))') (engs(i),i=1,nsum)        
        if (Its.eq.1) then
          write(*,'("     Exp: ",10F9.1,/,(10X,10F9.1))') (EXP_E(i),i=1,m)        
          write(*,'("Exp-Calc: ",10F9.1,/,(10X,10F9.1))') (fvec(i),i=1,m)
        endif 
        if (PLOTout1(7)) write(ip1,'(I4,E16.8,10F20.8)') Its,fnorm,(x(i),i=1,n)
 110    continue       
      endif  ! fitType.eq.1 .or. fitType.eq.5
!
!     FitType = 2  JuddOfelt
!     **********************
      if (fitType.eq.2) then  
        mm=m/2; if (2*mm.ne.m) then; write(io,'("m=",i2," must be even for JuddOfelt")') m; endif
        do ifit=1,n  ! =nFit   
          P(FitN(ifit)) = x(ifit)
        enddo
        call prepCalc(0)        
        CALL calcEDipJO()
        do i=1,m/2  ! m/2=Nexp
!          write(io,'("EXP_E(i),(EngJO(1,NAssign(i,j)),j=1,Nass(i)),EXP_I(i),(IntJO(1,NAssign(i,j)),j=1,Nass(i))=",20F12.2)')  &
!                      EXP_E(i),(EngJO(1,NAssign(i,j)),j=1,Nass(i)),EXP_I(i),(IntJO(1,NAssign(i,j)),j=1,Nass(i))
          addInt=0.0d+00 
          do j=1,Nass(i)  
            fvec(i)    =Wgt_E(i)*(EXP_E(i)-EngJO(1,NAssign(i,j)))
            if (IED.EQ.1) addInt=addInt+REDJO(1,NAssign(i,j),1)*1d+6
            if (IED.EQ.2) addInt=addInt+2.127d+12*REDJO(1,NAssign(i,j),1)/EngJO(1,NAssign(i,j))
            if (IED.EQ.3) addInt=addInt+IntJO(1,NAssign(i,j))
          enddo ! k1  
          fvec(m/2+i)=Wgt_I(i)*(EXP_I(i)-addInt)
        enddo ! i
        fnorm = enorm(m,fvec)/sqrt(dble(m))+penalty
        if (Its.eq.1) then
          write(*,'(5X,"ExpE:   ",10F9.1,/,(10X,10F9.1))') (EXP_E(i),i=1,m/2) 
          write(*,'(5X,"W*(E-C):",10F9.1,/,(10X,10F9.1))') (fvec(i),i=1,m/2)
          write(*,'(5X,"ExpI:   ",10F9.1,/,(10X,10F9.1))') (EXP_I(i),i=1,m/2)        
          write(*,'(5X,"W*(I-C):",10F9.1,/,(10X,10F9.1))') (fvec(m/2+i),i=1,m/2)
        endif 
        if (PLOTout1(7)) write(ip1,'(I4,E16.8,10F20.8)') Its,fnorm,(x(i),i=1,n)
      endif  ! fitType.eq.2  JuddOfelt
!                  
!     FitType = 3  Fitting AOM parameters to Bkq Crystal Field Coefficients
!     *********************************************************************
      if (fitType.eq.3) then 
        do ifit=1,n  ! =nFit   
          P(FitN(ifit)) = x(ifit)
        enddo
        call prepCalc(0)
        do i=1,m  ! m=Nexp
          if (NAssign(i,1).lt.0) fvec(i)=Wgt_E(i)*(EXP_E(i)-BKQI(abs(NAssign(i,1))))
          if (NAssign(i,1).gt.0) fvec(i)=Wgt_E(i)*(EXP_E(i)-BKQR(NAssign(i,1)))
        enddo
        fnorm = enorm(m,fvec)/sqrt(dble(m))+penalty
        if (Its.eq.1) then
          line1=""; line2=""
          do i=1,m
            if (NAssign(i,1).gt.0) then
              write(bit1,'(8X,A4)') BkqR_label(NAssign(i,1))
              write(bit2,'(F12.4)') BKQR(NAssign(i,1))
            else
              write(bit1,'(8X,A4)') BkqI_label(abs(NAssign(i,1)))
              write(bit2,'(F12.4)') BKQI(abs(NAssign(i,1)))
            endif
            line1=trim(line1)//" "//bit1
            line2=trim(line2)//bit2
          enddo  
          if (Its.eq.1) then
            write(*,'(7X,A120)') line1     
            write(*,'(2X,"Exp: ",10F12.4)') (EXP_E(i),i=1,m) 
          endif
          write(*,  '(2X,"Calc:",A120)') line2     
        endif 
!      write(*,'(5X,"diff:",10F9.1,/,(10X,10F9.1))') (fvec(i),i=1,m)
        if (PLOTout1(7)) write(ip1,'(I4,E16.8,10F20.8)') Its,fnorm,(x(i),i=1,n)
      endif   ! fitType=3
      
!                  
!     FitType = 4  Fitting calculated MEs to ESO MEs
!     **********************************************
      if (fitType.eq.4) then 
        do ifit=1,n; P(FitN(ifit)) = x(ifit); enddo  ! n=nFit   
        call prepCalc(0)
        call doPRED2()
!        call getBBasis(basis,1,nb)
!        CALL buildMatrix(basis,nb,2)   ! 2 for  PRED
        call transformPRED(Its)
        call calcESOFit(m,fvec)
        fnorm = enorm(m,fvec)/sqrt(dble(m))
        write(*,'("in fcn; fnorm=",e16.4)') fnorm
      endif ! fitType.eq.4

      if (fitType.eq.5) then  !  fitting calculated g-values to experimental
!        call calcMagMom()
        call calcGval()
        do i=1,m   !
          j=(NAssignG(i)+1)/2
          if (j.lt.1 .or. j.gt.m) then
            write(io,'("invalid j=",I4,",  NAssignG(",i3,")=",I4)') j,i,NAssignG(i)
            stop
          endif
          fvec(i) = Wgt_G(1,i)*(EXP_G(1,i)-gval(1,j)) + Wgt_G(2,i)*(EXP_G(2,i)-gval(2,j)) + Wgt_G(3,i)*(EXP_G(3,i)-gval(3,j))
!          write(io,'("i,j=",2I4,",  NAssignG(j)=",I4,", ExpG=",3F8.4,", calcG(Nassign)=",3F8.4,", fvec=",f8.4)') &
!               i,j,NAssignG(j),EXP_G(1,i),EXP_G(2,i),EXP_G(3,i),gval(1,j),gval(2,j),gval(3,j),fvec(i)
        enddo
        fnorm = enorm(m,fvec)/sqrt(dble(m))
!  ToDo, this only works for odd electron systems        
      endif ! fitType.eq.5

!                  
!     FitType = 6  Fitting calculated 1-e LF MEs to a given set
!     *********************************************************
      if (fitType.eq.6) then
        do ifit=1,n; P(FitN(ifit)) = x(ifit); enddo  ! n=nFit   
        call prepCalc(0)  ! Calculates LF1eMat here
        fnorm=0.0d+00; m=0 
        do i=1,2*Lvalue+1
          do j=1,i
            m=m+1
            fvec(m) = Wgt_E(m)*(EXP_E(m)-LF1eMat(i,j))
          enddo
        enddo  
        fnorm = enorm(m,fvec)/sqrt(dble(m))  ! m=15(=1+2+3+4+5) or 28(=1+2+3+4+5+6+7)
!        write(*,'("in fcn; fnorm=",e16.8,";  m=",i3)') fnorm,m
      endif ! fitType.eq.6
!      
      calcVecs=cV
!
      if (FitOpt(3)) then
        write(io,'("it=",I4,"; sq(del^2/N)+pen=",E14.8,"; pen=",E14.8,"; P="(10F16.8))') Its,fnorm,penalty,(x(i),i=1,n)   
        write(io,'("fvec=",10F9.2,/,(5X,10F9.2))') (fvec(i),i=1,m)
      endif
      write(*,'("it=",I4,"; sq(del^2/N)+pen=",E14.8,"; pen=",E14.8,"; P="(10F16.8))') Its,fnorm,penalty,(x(i),i=1,n)      
!
      return
!
      END SUBROUTINE fcn
!
!------------------------------------------------------------------------------
!
      subroutine doFit()
! ----------
! Calls the Levenberg-Marquardt SUBROUTINE
! ----------
! Input:
! m     (=Nexp)  The number of experimental values to be used. 
! n     (=Nfit)  The number of parameters to be varied. The parameters are in X(N) 
! x(n)  (=FitP)  The n parameters. 
! Output:
! fvec(m) The m calculated values corresponding to the experimental values.
! info    
! lwa   is a positive integer input variable not less than m*n+5*n+m.
! nfev  the total number of function evaluations.
      implicit none
      integer i,info,lwa, iwa(maxFit),nfev,m
      REAL*8 x(maxFit),fvec(max_exp),wa(6400)  ! max_exp*maxFit+5*maxFit+max_exp = 300*20+5*20+300
! 
      if ((Nexp.le.0 .and. Nexpg.le.0) .or. NFit.le.0) then
        write(io,'("Cannot have Nexp&NexpG=0 or Nfit=0 in doFit() call")') 
        return
      endif
      m=Nexp
      if (JuddOfelt) m=2*Nexp
      if (fitType.eq.5) m=NexpG
      info = 0
      lwa = m*Nfit + 5*Nfit + m
      if (JuddOfelt) lwa = 2*m*Nfit + 5*Nfit + 2*m      
!      
!      write(io,'("fitType,m,NFit=",3i4)') fitType,m,NFit
      call lmdif1(m,NFit,fitP,fvec,FitTol,info,iwa,wa,lwa,nfev)
      if (info.lt.1 .or. info.gt.3) then
        write(io,'("Back from fitting algorithm after",i3," function evaluations; info=",i2)') nfev,info
        if (info.eq.0) write(io,'("Incorrect input parameters to lmdif1")') 
        if (info.eq.4) write(io,'("fvec is orthogonal to the columns of the jacobian to machine precision")') 
        if (info.eq.5) write(io,'("Number of calls to fcn has reached or exceeded 200*(n+1)")') 
        if (info.eq.6) write(io,'("tol is too small, no further reduction in the sum of squares is possible")') 
        if (info.eq.7) write(io,'("tol is too small, no further improvement in the approximate solution x is possible")') 
      endif
      call loadP()
      do i=1,Nfit
        P(FitN(i))=fitP(i)
      enddo
      call prepCalc(0)
!
      return
      END SUBROUTINE doFit
!
!------------------------------------------------------------------------------
!
!      subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
      subroutine lmdif1(m,n,x,fvec,tol,info,iwa,wa,lwa,nfev)
      implicit none
      integer m,n,info,lwa, iwa(n)
      REAL*8 tol, x(n),fvec(m),wa(lwa)
!      external fcn
!
!   The purpose of lmdif1 is to minimize the sum of the squares of m nonlinear functions 
!   in n variables by a modification of the Levenberg-Marquardt algorithm. This is done 
!   by using the more general least-squares solver lmdif. The user must provide a
!   subroutine which calculates the functions. The jacobian is then calculated by a 
!   forward-difference approximation.
!
!   The subroutine statement is
!
!       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         real*8 x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif1.
!         in this case set iflag to a negative integer.
!
!  m is a positive integer input variable set to the number of functions.
!
!  n is a positive integer input variable set to the number of variables. n must not exceed m.
!
!  x is an array of length n. on input x must contain an initial estimate of the solution vector. 
!       On output x contains the final estimate of the solution vector.
!
!  fvec is an output array of length m which contains the functions evaluated at the output x.
!
!  tol is a nonnegative input variable. Termination occurs when the algorithm estimates either that 
!      the relative error in the sum of squares is at most tol or that the relative error between 
!      x and the solution is at most tol.
!
!  info is an integer output variable. If the user has terminated execution, info is set to the 
!     (negative) value of iflag. See description of fcn. Otherwise, info is set as follows.
!
!     info = 0  improper input parameters.
!     info = 1  algorithm estimates that the relative error in the sum of squares is at most tol.
!     info = 2  algorithm estimates that the relative error between x and the solution is at most tol.
!     info = 3  conditions for info = 1 and info = 2 both hold.
!     info = 4  fvec is orthogonal to the columns of the jacobian to machine precision.
!     info = 5  number of calls to fcn has reached or exceeded 200*(n+1).
!     info = 6  tol is too small. no further reduction in the sum of squares is possible.
!     info = 7  tol is too small. no further improvement in the approximate solution x is possible.
!
!  iwa is an integer work array of length n.
!
!  wa is a work array of length lwa.
!
!  lwa is a positive integer input variable not less than m*n+5*n+m.
!
!     subprograms called
!       user-supplied ...... fcn
!       minpack-supplied ... lmdif
!
!   argonne national laboratory. minpack project. march 1980.
!   burton s. garbow, kenneth e. hillstrom, jorge j. more
!
      integer maxfev,mode,mp5n,nfev,nprint
      real*8 epsfcn,ftol,gtol,xtol,zero,factor
      Parameter(zero=0.0d0)
      info = 0
!
!  Check the input parameters for errors.
      if (n.le.0 .or. m.lt.n .or. tol.lt.zero .or. lwa.lt.m*n+5*n+m) return
!
!  Call lmdif.
!      maxfev = 200*(n + 1)  ! maximum number of function evaluations
      maxfev = maxIts
      ftol = tol  ! Termination occurs when the relative error between x and the solution is at most tol. 
      xtol = tol  ! Termination occurs when the relative error in the sum of squares is at most tol.
      gtol = zero
      factor=100d+00  ! factor is a positive input variable used in determining the initial L-M step bound. 
                      ! In most cases factor should lie in the range (.1,100.), 100. is a generally recommended value.
      epsfcn = 1.0d-14
      mode = 1  !  If mode=1, the variables will be scaled internally, if mode=2, the scaling is specified by the input diag. 
      nprint = 0
      mp5n = m + 5*n
!      call lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa(1) 
      call lmdif(m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa(1),  &
                 mode,factor,nprint,info,nfev,wa(mp5n+1),m,iwa,  &
                 wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info .eq. 8) info = 4
!      write(io,'("Back from lmdif; info,nfev=",2I4)') info,nfev
      return
!
      end subroutine lmdif1
!
!------------------------------------------------------------------------------
!
      subroutine covar1(n,r,ldr,ipvt,tol,wa) 
      implicit none
      integer n,ldr 
      integer ipvt(n) 
      real*8 tol 
      real*8 r(ldr,n),wa(n) 
! 
!     subroutine covar 
! 
!     given an m by n matrix a, the problem is to determine 
!     the covariance matrix corresponding to a, defined as 
! 
!                    t 
!           inverse(a *a) . 
! 
!     this subroutine completes the solution of the problem 
!     if it is provided with the necessary information from the 
!     qr factorization, with column pivoting, of a. that is, if 
!     a*p = q*r, where p is a permutation matrix, q has orthogonal 
!     columns, and r is an upper triangular matrix with diagonal 
!     elements of nonincreasing magnitude, then covar expects 
!     the full upper triangle of r and the permutation matrix p. 
!     the covariance matrix is then computed as 
! 
!                      t     t 
!           p*inverse(r *r)*p  . 
! 
!     if a is nearly rank deficient, it may be desirable to compute 
!     the covariance matrix corresponding to the linearly independent 
!     columns of a. to define the numerical rank of a, covar uses 
!     the tolerance tol. if l is the largest integer such that 
! 
!           abs(r(l,l)) .gt. tol*abs(r(1,1)) , 
! 
!     then covar computes the covariance matrix corresponding to 
!     the first l columns of r. for k greater than l, column 
!     and row ipvt(k) of the covariance matrix are set to zero. 
! 
!     the subroutine statement is 
! 
!       subroutine covar(n,r,ldr,ipvt,tol,wa) 
! 
!     where 
! 
!       n is a positive integer input variable set to the order of r. 
! 
!       r is an n by n array. on input the full upper triangle must 
!         contain the full upper triangle of the matrix r. on output 
!         r contains the square symmetric covariance matrix. 
! 
!       ldr is a positive integer input variable not less than n 
!         which specifies the leading dimension of the array r. 
! 
!       ipvt is an integer input array of length n which defines the 
!         permutation matrix p such that a*p = q*r. column j of p 
!         is column ipvt(j) of the identity matrix. 
! 
!       tol is a nonnegative input variable used to define the 
!         numerical rank of a in the manner described above. 
! 
!       wa is a work array of length n. 
! 
!     subprograms called 
! 
!       fortran-supplied ... dabs 
! 
!     argonne national laboratory. minpack project. august 1980. 
!     burton s. garbow, kenneth e. hillstrom, jorge j. more 
! 
! 
!     form the inverse of r in the full upper triangle of r. 
! 
      integer i,ii,j,jj,k,km1,l
      real*8 temp,tolr,zero,one
      logical sing
      
      zero=0.0d+00; one=1.0d+00
!      
      tolr = tol*dabs(r(1,1)) 
      l = 0 
      do 40 k = 1, n 
         if (dabs(r(k,k)) .le. tolr) go to 50 
         r(k,k) = one/r(k,k) 
         km1 = k - 1 
         if (km1 .lt. 1) go to 30 
         do 20 j = 1, km1 
            temp = r(k,k)*r(j,k) 
            r(j,k) = zero 
            do 10 i = 1, j 
               r(i,k) = r(i,k) - temp*r(i,j) 
   10          continue 
   20       continue 
   30    continue 
         l = k 
   40    continue 
   50 continue 
! 
!  form the full upper triangle of the inverse of (r transpose)*r 
!  in the full upper triangle of r. 
! 
      if (l .lt. 1) go to 110 
      do 100 k = 1, l 
         km1 = k - 1 
         if (km1 .lt. 1) go to 80 
         do 70 j = 1, km1 
            temp = r(j,k) 
            do 60 i = 1, j 
               r(i,j) = r(i,j) + temp*r(i,k) 
   60          continue 
   70       continue 
   80    continue 
         temp = r(k,k) 
         do 90 i = 1, k 
            r(i,k) = temp*r(i,k) 
   90       continue 
  100    continue 
  110 continue 
! 
!  form the full lower triangle of the covariance matrix 
!  in the strict lower triangle of r and in wa. 
! 
      do 130 j = 1, n 
         jj = ipvt(j) 
         sing = j .gt. l 
         do 120 i = 1, j 
            if (sing) r(i,j) = zero 
            ii = ipvt(i) 
            if (ii .gt. jj) r(ii,jj) = r(i,j) 
            if (ii .lt. jj) r(jj,ii) = r(i,j) 
  120       continue 
         wa(jj) = r(j,j) 
  130    continue 
! 
! symmetrize the covariance matrix in r. 
! 
      do 150 j = 1, n 
         do 140 i = 1, j 
            r(i,j) = r(j,i) 
  140       continue 
         r(j,j) = wa(j) 
  150    continue 
      return 
! 
! last card of subroutine covar. 
! 
      end subroutine covar1
!
!------------------------------------------------------------------------------
!
!      subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,   
      subroutine lmdif(m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn, &
                       diag,mode,factor,nprint,info,nfev,fjac,  &
                       ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
      implicit none
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac
      integer ipvt(n)
      real*8 ftol,xtol,gtol,epsfcn,factor
      real*8 x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),wa1(n),wa2(n),wa3(n),wa4(m)
!      external fcn
!
!  The purpose of lmdif is to minimize the sum of the squares of m nonlinear functions 
!  in n variables by a modification of the Levenberg-Marquardt algorithm. The user must 
!  provide a subroutine which calculates the functions. The Jacobian is then calculated 
!  by a forward-difference approximation.
!
!  The subroutine statement is
!
!    subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
!                     diag,mode,factor,nprint,info,nfev,fjac,
!                     ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!    where
!
!  fcn is the name of the user-supplied subroutine which calculates the functions. 
!      fcn must be declared  in an external statement in the user calling program, 
!      and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!       The value of iflag should not be changed by fcn unless the user wants to 
!         terminate execution of lmdif. In this case set iflag to a negative integer.
!
!   m is a positive integer input variable set to the number of functions.
!
!   n is a positive integer input variable set to the number of variables. 
!     n must not exceed m.
!
!   x is an array of length n. on input x must contain an initial estimate of the 
!     solution vector. on output x contains the final estimate of the solution vector.
!
!   fvec is an output array of length m which contains the functions evaluated at the output x.
!
!   ftol is a nonnegative input variable. Termination occurs when both the actual and 
!     predicted relative reductions in the sum of squares are at most ftol.
!     Therefore, ftol measures the relative error desired in the sum of squares.
!
!   xtol is a nonnegative input variable. Termination occurs when the relative 
!     error between two consecutive iterates is at most xtol. Therefore, xtol measures the
!     relative error desired in the approximate solution.
!
!   gtol is a nonnegative input variable. Termination occurs when the cosine of the 
!     angle between fvec and any column of the jacobian is at most gtol in absolute
!     value. Therefore, gtol measures the orthogonality desired between the function 
!     vector and the columns of the jacobian.
!
!   maxfev is a positive integer input variable. Termination occurs when the number of 
!     calls to fcn is at least maxfev by the end of an iteration.
!
!   epsfcn is an input variable used in determining a suitable step length for the 
!     forward-difference approximation. This approximation assumes that the relative errors 
!     in the functions are of the order of epsfcn. If epsfcn is less than the machine precision, 
!     it is assumed that the relative errors in the functions are of the order of the machine precision.
!
!   diag is an array of length n. if mode = 1 (see below), diag is internally set. If mode = 2, diag
!     must contain positive entries that serve as multiplicative scale factors for the variables.
!
!   mode is an integer input variable. 
!        If mode = 1, the variables will be scaled internally. 
!        If mode = 2, the scaling is specified by the input diag. 
!        Other values of mode are equivalent to mode = 1.
!
!   factor is a positive input variable used in determining the initial step bound. 
!       This bound is set to the product of factor and the euclidean norm of diag*x 
!       if nonzero, or else to factor itself. In most cases factor should lie in the
!       interval (.1,100.). 100. is a generally recommended value.
!
!   nprint is an integer input variable that enables controlled printing of iterates 
!       if it is positive. in this case, fcn is called with iflag = 0 at the beginning
!       of the first iteration and every nprint iterations thereafter and immediately 
!       prior to return, with x and fvec available for printing. if nprint is not positive, 
!       no special calls of fcn with iflag = 0 are made.
!
!  info is an integer output variable. if the user has terminated execution, info is set 
!       to the (negative) value of iflag. See description of fcn. Otherwise,
!         info is set as follows:
!
!   info = 0  improper input parameters.
!   info = 1  both actual and predicted relative reductions in the sum of squares are at most ftol.
!   info = 2  relative error between two consecutive iterates is at most xtol.
!   info = 3  conditions for info = 1 and info = 2 both hold.
!   info = 4  the cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.
!   info = 5  number of calls to fcn has reached or exceeded maxfev.
!   info = 6  ftol is too small. no further reduction in the sum of squares is possible.
!   info = 7  xtol is too small. no further improvement in the approximate solution x is possible.
!   info = 8  gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision.
!
!  nfev is an integer output variable set to the number of calls to fcn.
!
!  fjac is an output m by n array. the upper n by n submatrix of fjac contains an 
!       upper triangular matrix r with diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!        where p is a permutation matrix and jac is the final calculated jacobian.
!        Column j of p is column ipvt(j) (see below) of the identity matrix. The lower trapezoidal
!        part of fjac contains information generated during the computation of r.
!
!  ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.
!
!  ipvt is an integer output array of length n. ipvt defines a permutation matrix p 
!      such that jac*p = q*r, where jac is the final calculated jacobian, q is orthogonal 
!      (not stored), and r is upper triangular with diagonal elements of nonincreasing 
!      magnitude. Column j of p is column ipvt(j) of the identity matrix.
!
!  qtf is an output array of length n which contains
!      the first n elements of the vector (q transpose)*fvec.
!
!  wa1, wa2, and wa3 are work arrays of length n.
!
!  wa4 is a work array of length m.
!
!     subprograms called
!       user-supplied ...... fcn
!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
      integer i,iflag,iter,j,l
      real*8 actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,one,par,pnorm
      real*8 prered,p1,p5,p25,p75,p0001,ratio,sum,temp,temp1,temp2,xnorm,zero
!      real*8 dpmpar,enorm
      data one,p1,p5,p25,p75,p0001,zero/1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
!      external enorm
!
!  epsmch is the machine precision.
!      epsmch = dpmpar(1)
      epsmch = 2.22044604926d-16
!
      info = 0;   iflag = 0;     nfev = 0
!
!  Check the input parameters for errors.
      if (n.le.0 .or. m.lt.n .or. ldfjac.lt.m .or. ftol.lt.zero .or. xtol.lt.zero .or.   &
          gtol.lt.zero .or. maxfev.le.0 .or. factor.le.zero) return
      if (mode .eq. 2) then
        do j = 1, n
          if (diag(j) .le. zero) go to 300
        enddo 
      endif
!
!  Evaluate the function at the starting point and calculate its norm.
      iflag = 1
      call fcn(m,n,x,fvec,iflag)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(m,fvec)
      fnorm_init=fnorm/sqrt(dble(m)) 
!
!  Initialize levenberg-marquardt parameter and iteration counter.
      par = zero
      iter = 1
!
!  Beginning of the outer loop.
   30 continue
!
!  Calculate the jacobian matrix.
        iflag = 2
!        call fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
        call fdjac2(m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
        nfev = nfev + n
        if (iflag .lt. 0) go to 300
!
!  If requested, call fcn to enable printing of iterates.
        if (nprint .gt. 0) then
          iflag = 0
          if (mod(iter-1,nprint) .eq. 0) call fcn(m,n,x,fvec,iflag)
          if (iflag .lt. 0) go to 300
        endif
!
!  Compute the qr factorization of the jacobian.
        call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
!
!  On the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
        if (iter .eq. 1) then
          if (mode .ne. 2) then
            do j = 1, n
              diag(j) = wa2(j)
              if (wa2(j) .eq. zero) diag(j) = one
            enddo
!            write(io,'("iter,diag= ",I3,10F10.4)') iter,(diag(i),i=1,n)
          endif
!
!  On the first iteration, calculate the norm of the scaled x
!  and initialize the step bound delta.
!            write(io,'(" x= ",10F10.4)') (x(i),i=1,n)
          do j = 1, n
            wa3(j) = diag(j)*x(j)
          enddo
          xnorm = enorm(n,wa3)
          delta = factor*xnorm
          if (delta .eq. zero) delta = factor
!          write(io,'("xnorm,factor,delta=",10F10.4)') xnorm,factor,delta
        endif
!
!  Form (q transpose)*fvec and store the first n components in qtf.
        do i = 1, m
          wa4(i) = fvec(i)
        enddo
        do j = 1, n
          if (fjac(j,j) .ne. zero) then
            sum = zero
            do i = j, m
              sum = sum + fjac(i,j)*wa4(i)
            enddo
            temp = -sum/fjac(j,j)
            do i = j, m
              wa4(i) = wa4(i) + fjac(i,j)*temp
            enddo
          endif
          fjac(j,j) = wa1(j)
          qtf(j) = wa4(j)
        enddo
!
!  Compute the norm of the scaled gradient.
        gnorm = zero
        if (fnorm .ne. zero) then
          do j = 1, n
            l = ipvt(j)
            if (wa2(l) .ne. zero) then
              sum = zero
              do i = 1, j
                sum = sum + fjac(i,j)*(qtf(i)/fnorm)
              enddo 
              gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
            endif
          enddo
        endif
!
!  Test for convergence of the gradient norm.
        if (gnorm .le. gtol) info = 4
        if (info .ne. 0) go to 300
!
!  Rescale if necessary.
        if (mode .ne. 2) then
          do j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
          enddo
!            write(io,'("iter,diag= ",I3,10F10.4)') iter,(diag(i),i=1,n)
        endif
!
!  Beginning of the inner loop.
 200    continue
!
!  Determine the levenberg-marquardt parameter.
          call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,wa3,wa4)
!
!  Store the direction p and x + p. calculate the norm of p.
            do 210 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  210          continue
            pnorm = enorm(n,wa3)
!
!  On the first iteration, adjust the initial step bound.
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
!            write(io,'("iter,delta,pnorm= ",I3,10F10.4)') iter,delta,pnorm
!
!  Evaluate the function at x + p and calculate its norm.
            iflag = 1
            call fcn(m,n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(m,wa4)
!
!  Compute the scaled actual reduction.
            actred = -one
            if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
!
!  Compute the scaled predicted reduction and the scaled directional derivative.
            do j = 1, n
              wa3(j) = zero
              l = ipvt(j)
              temp = wa1(l)
              do i = 1, j
                wa3(i) = wa3(i) + fjac(i,j)*temp
              enddo
            enddo
            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
!
!  Compute the ratio of the actual to the predicted reduction.
            ratio = zero
            if (prered .ne. zero) ratio = actred/prered
!
!  Update the step bound.
            if (ratio .gt. p25) go to 240
               if (actred .ge. zero) temp = p5
               if (actred .lt. zero) temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
               go to 260
  240       continue
               if (par .ne. zero .and. ratio .lt. p75) go to 250
               delta = pnorm/p5
               par = p5*par
  250          continue
  260       continue
!
!  Test for successful iteration.
            if (ratio .lt. p0001) go to 290
!
!  Successful iteration. update x, fvec, and their norms.
            do j = 1, n
              x(j) = wa2(j)
              wa2(j) = diag(j)*x(j)
            enddo
            do i = 1, m
              fvec(i) = wa4(i)
            enddo
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
 290      continue
!
!  Tests for convergence.
          if (dabs(actred).le.ftol .and. prered.le.ftol .and. p5*ratio.le.one) info = 1
          if (delta .le. xtol*xnorm) info = 2
          if (dabs(actred).le.ftol .and. prered.le.ftol .and. p5*ratio.le.one .and. info.eq.2) info = 3
          if (info .ne. 0) go to 300
!
!  Tests for termination and stringent tolerances.
          if (nfev .ge. maxfev) info = 5
          if (dabs(actred).le.epsmch .and. prered.le.epsmch .and. p5*ratio.le.one) info = 6
          if (delta .le. epsmch*xnorm) info = 7
          if (gnorm .le. epsmch) info = 8
          if (info .ne. 0) go to 300
!
!  End of the inner loop. repeat if iteration unsuccessful.
          if (ratio .lt. p0001) go to 200
!
!  End of the outer loop.
        go to 30
 300  continue
!
!  Termination, either normal or user imposed.
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(m,n,x,fvec,iflag)
!     
!      call covar(n,r,ldr,ipvt,tol,wa) 
      call covar1(n,fjac,ldfjac,ipvt,xtol,wa1) 
      do i=1,n
        do j=1,n
          cov(i,j)=fjac(i,j)
        enddo
      enddo        
!       
      return
!
      end subroutine lmdif
!      
!------------------------------------------------------------------------------
!
!      subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
      subroutine fdjac2(m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
      implicit none
      integer m,n,ldfjac,iflag
      real*8 epsfcn
      real*8 x(n),fvec(m),fjac(ldfjac,n),wa(m)
!
!  this subroutine computes a forward-difference approximation
!  to the m by n jacobian matrix associated with a specified
!  problem of m functions in n variables.
!
!  the subroutine statement is
!
!    subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
!
!  where
!
!    fcn is the name of the user-supplied subroutine which
!      calculates the functions. fcn must be declared
!      in an external statement in the user calling
!      program, and should be written as follows.
!
!      subroutine fcn(m,n,x,fvec,iflag)
!      integer m,n,iflag
!      double precision x(n),fvec(m)
!      ----------
!      calculate the functions at x and
!      return this vector in fvec.
!      ----------
!      return
!      end
!
!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of fdjac2.
!      in this case set iflag to a negative integer.
!
!   m is a positive integer input variable set to the number
!     of functions.
!
!   n is a positive integer input variable set to the number
!     of variables. n must not exceed m.
!
!   x is an input array of length n.
!
!   fvec is an input array of length m which must contain the
!     functions evaluated at x.
!
!   fjac is an output m by n array which contains the
!     approximation to the jacobian matrix evaluated at x.
!
!   ldfjac is a positive integer input variable not less than m
!     which specifies the leading dimension of the array fjac.
!
!   iflag is an integer variable which can be used to terminate
!     the execution of fdjac2. see description of fcn.
!
!   epsfcn is an input variable used in determining a suitable
!     step length for the forward-difference approximation. this
!     approximation assumes that the relative errors in the
!     functions are of the order of epsfcn. if epsfcn is less
!     than the machine precision, it is assumed that the relative
!     errors in the functions are of the order of the machine
!     precision.
!
!   wa is a work array of length m.
!
!   subprograms called:
!     user-supplied ...... fcn
!     minpack-supplied ... dpmpar
!     fortran-supplied ... dabs,dmax1,dsqrt
!
!   argonne national laboratory. minpack project. march 1980.
!   burton s. garbow, kenneth e. hillstrom, jorge j. more
!
      integer i,j
      real*8 eps,epsmch,h,temp,zero,  dpmpar
      data zero /0.0d0/
!
!     epsmch is the machine precision.
!      epsmch = dpmpar(1)
      epsmch = 2.22044604926d-16
!
      eps = dsqrt(dmax1(epsfcn,epsmch))
      do j = 1, n
        temp = x(j)
        h = eps*dabs(temp)
        if (h .eq. zero) h = eps
        x(j) = temp + h
        call fcn(m,n,x,wa,iflag)
        if (iflag .lt. 0) return
        x(j) = temp
        do i = 1, m
          fjac(i,j) = (wa(i) - fvec(i))/h
        enddo
      enddo
      return
!
      end subroutine fdjac2
!      
!-----------------------------------------------------------------------
!
      real*8 function enorm(n,x)
      implicit none
      integer n
      real*8 x(n)
!
!   Given an n-vector x, this function calculates the euclidean norm of x.
!
!   The euclidean norm is computed by accumulating the sum of squares in three 
!   different sums. The sums of squares for the small and large components are 
!   scaled so that no overflows occur. Non-destructive underflows are permitted. 
!   Underflows and overflows do not occur in the computation of the unscaled sum
!   of squares for the intermediate components. The definitions of small, 
!   intermediate and large components depend on two constants, rdwarf and rgiant. 
!   The main restrictions on these constants are that rdwarf**2 not underflow 
!   and rgiant**2 not overflow. The constants given here are suitable for every 
!   known computer.
!
!   The function statement is
!
!      double precision function enorm(n,x)
!
!    where
!
!    n is a positive integer input variable.
!    x is an input array of length n.
!
!    argonne national laboratory. minpack project. march 1980.
!    burton s. garbow, kenneth e. hillstrom, jorge j. more
!
      integer i
      real*8 agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      s1 = zero;   s2 = zero;    s3 = zero
      x1max = zero;  x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do i = 1, n
        xabs = dabs(x(i))
        if (.not.(xabs .gt. rdwarf .and. xabs .lt. agiant)) then
          if (xabs .gt. rdwarf) then
!
!  sum for large components.
            if (xabs .gt. x1max) then
              s1 = one + s1*(x1max/xabs)**2
              x1max = xabs
            else
              s1 = s1 + (xabs/x1max)**2
            endif
          else
!
!  sum for small components.
            if (xabs .gt. x3max) then
              s3 = one + s3*(x3max/xabs)**2
              x3max = xabs
            else
             if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
            endif
          endif
        else
!
!  sum for intermediate components.
        s2 = s2 + xabs**2
        endif
      enddo  ! i
!
!  calculation of norm.
      if (s1 .ne. zero) then
        enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
      else
        if (s2 .ne. zero) then
          if (s2 .ge. x3max) enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
          if (s2 .lt. x3max) enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
        else
          enorm = x3max*dsqrt(s3)
        endif
      endif
!
      return
!
      end function enorm
!      
!-----------------------------------------------------------------------
!
      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
      IMPLICIT NONE
      integer n,ldr, ipvt(n)
      real*8 r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa(n)
!
!  Given an m by n matrix a, an n by n diagonal matrix d, and an m-vector b, 
!  the problem is to determine an x which solves the system
!
!           a*x = b ,     d*x = 0 ,
!
!  in the least squares sense.
!
!  This subroutine completes the solution of the problem if it is provided 
!  with the necessary information from the qr factorization, with column pivoting,
!  of a. that is, if a*p = q*r, where p is a permutation matrix, q has orthogonal
!  columns, and r is an upper triangular matrix with diagonal elements of nonincreasing 
!  magnitude, then qrsolv expects the full upper triangle of r, the permutation 
!  matrix p, and the first n components of (q transpose)*b. The system a*x = b, 
!  d*x = 0, is then equivalent to
!
!                  t       t
!          r*z = q *b ,  p *d*p*z = 0 ,
!
!  where x = p*z. If this system does not have full rank, then a least squares solution
!  is obtained. On output qrsolv also provides an upper triangular matrix s such that
!
!           t   t               t
!          p *(a *a + d*d)*p = s *s .
!
!  s is computed within qrsolv and may be of separate interest.
!
!  The subroutine statement is
!
!      subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
!
!  n is a positive integer input variable set to the order of r.
!
!  r is an n by n array. on input the full upper triangle
!     must contain the full upper triangle of the matrix r.
!     on output the full upper triangle is unaltered, and the
!     strict lower triangle contains the strict upper triangle
!     (transposed) of the upper triangular matrix s.
!
!  ldr is a positive integer input variable not less than n
!     which specifies the leading dimension of the array r.
!
!  ipvt is an integer input array of length n which defines the
!     permutation matrix p such that a*p = q*r. column j of p
!     is column ipvt(j) of the identity matrix.
!
!  diag is an input array of length n which must contain the
!     diagonal elements of the matrix d.
!
!  qtb is an input array of length n which must contain the first
!     n elements of the vector (q transpose)*b.
!
!  x is an output array of length n which contains the least
!     squares solution of the system a*x = b, d*x = 0.
!
!  sdiag is an output array of length n which contains the
!     diagonal elements of the upper triangular matrix s.
!
!  wa is a work array of length n.
!
!  subprograms called
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
      integer i,j,jp1,k,kp1,l,nsing
      double precision cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero
      data p5,p25,zero /5.0d-1,2.5d-1,0.0d0/
!
!  Copy r and (q transpose)*b to preserve input and initialize s.
!  In particular, save the diagonal elements of r in x.
      do j = 1, n
        do i = j, n
          r(i,j) = r(j,i)
        enddo ! i
        x(j) = r(j,j)
        wa(j) = qtb(j)
      enddo ! j
!
!  Eliminate the diagonal matrix d using a givens rotation.
      do j = 1, n
!
!  Prepare the row of d to be eliminated, locating the diagonal element 
!  using p from the qr factorization.
        l = ipvt(j)
        if (diag(l) .eq. zero) go to 90
        do k = j, n
          sdiag(k) = zero
        enddo ! k
        sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of d modify only a single element 
!  of (q transpose)*b beyond the first n, which is initially zero.
        qtbpj = zero
        do k = j, n
!
!  Determine a givens rotation which eliminates the appropriate element in the current row of d.
          if (sdiag(k) .eq. zero) go to 70
          if (dabs(r(k,k)) .ge. dabs(sdiag(k))) go to 40
            cotan = r(k,k)/sdiag(k)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            go to 50
   40     continue
            tan = sdiag(k)/r(k,k)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
   50     continue
!
!  Compute the modified diagonal element of r and the modified element of ((q transpose)*b,0).
          r(k,k) = cos*r(k,k) + sin*sdiag(k)
          temp = cos*wa(k) + sin*qtbpj
          qtbpj = -sin*wa(k) + cos*qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of s.
          kp1 = k + 1
          if (n .lt. kp1) go to 70
          do i = kp1, n
            temp = cos*r(i,k) + sin*sdiag(i)
            sdiag(i) = -sin*r(i,k) + cos*sdiag(i)
            r(i,k) = temp
          enddo ! i 
 70       continue
        enddo  ! k 
 90     continue
!
!  Store the diagonal element of s and restore the corresponding diagonal element of r.
        sdiag(j) = r(j,j)
        r(j,j) = x(j)
      enddo  ! j
!
!  Solve the triangular system for z. If the system is singular, 
!  then obtain a least squares solution.
      nsing = n
      do j = 1, n
        if (sdiag(j) .eq. zero .and. nsing .eq. n) nsing = j - 1
        if (nsing .lt. n) wa(j) = zero
      enddo ! j
      if (nsing .lt. 1) go to 150
      do k = 1, nsing
        j = nsing - k + 1
        sum = zero
        jp1 = j + 1
        if (nsing .lt. jp1) go to 130
        do i = jp1, nsing
          sum = sum + r(i,j)*wa(i)
        enddo ! i 
 130    continue
        wa(j) = (wa(j) - sum)/sdiag(j)
      enddo ! k
 150  continue
!
!  Permute the components of z back to components of x.
      do j = 1, n
        l = ipvt(j)
        x(l) = wa(j)
      enddo ! j
      return
!
      end subroutine qrsolv
!      
!-----------------------------------------------------------------------
!
      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      IMPLICIT NONE
      integer m,n,lda,lipvt, ipvt(lipvt)
      logical pivot
      REAL*8 a(lda,n),rdiag(n),acnorm(n),wa(n)
!
!   This subroutine uses householder transformations with column pivoting (optional) 
!   to compute a qr factorization of the m by n matrix a. 
!   That is, qrfac determines an orthogonal  matrix q, a permutation matrix p, 
!   and an upper trapezoidal matrix r with diagonal elements of nonincreasing magnitude,
!   such that a*p = q*r. The householder transformation for column k, k = 1,2,...,min(m,n), 
!   is of the form
!
!                         t
!         i - (1/u(k))*u*u
!
!   where u has zeros in the first k-1 positions. The form of this transformation 
!   and the method of pivoting first appeared in the corresponding linpack subroutine.
!
!   The subroutine statement is
!
!   subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
!  where
!
!  m is a positive integer input variable set to the number of rows of a.
!
!  n is a positive integer input variable set to the number of columns of a.
!
!  a is an m by n array. 
!      On input a contains the matrix for which the qr factorization is to be computed. 
!      On output the strict upper trapezoidal part of a contains the strict upper 
!      trapezoidal part of r, and the lower trapezoidal part of a contains a factored 
!      form of q (the non-trivial elements of the u vectors described above).
!
!  lda is a positive integer input variable not less than m which specifies the leading 
!      dimension of the array a.
!
!  pivot is a logical input variable. If pivot is set true, then column pivoting is enforced. 
!      If pivot is set false, then no column pivoting is done.
!
!  ipvt is an integer output array of length lipvt. ipvt defines the permutation matrix p 
!      such that a*p = q*r. Column j of p is column ipvt(j) of the identity matrix.
!      If pivot is false, ipvt is not referenced.
!
!  lipvt is a positive integer input variable. If pivot is false, then lipvt may be as 
!      small as 1. if pivot is true, then lipvt must be at least n.
!
!  rdiag is an output array of length n which contains the diagonal elements of r.
!
!  acnorm is an output array of length n which contains the norms of the 
!      corresponding columns of the input matrix a. If this information is 
!      not needed, then acnorm can coincide with rdiag.
!
!  wa is a work array of length n. if pivot is false, then wa can coincide with rdiag.
!
!  subprograms called
!    minpack-supplied ... dpmpar,enorm
!    fortran-supplied ... dmax1,dsqrt,min0
!
!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more
!
      integer i,j,jp1,k,kmax,minmn
      REAL*8 ajnorm,epsmch,one,p05,sum,temp,zero  ! , dpmpar,enorm
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
!
!     epsmch is the machine precision.
!      epsmch = dpmpar(1)
      epsmch = 2.22044604926d-16
!
!  compute the initial column norms and initialize several arrays.
      do j = 1, n
        acnorm(j) = enorm(m,a(1,j))
        rdiag(j) = acnorm(j)
        wa(j) = rdiag(j)
        if (pivot) ipvt(j) = j
      enddo   
!
!  reduce a to r with householder transformations.
      minmn = min0(m,n)
      do j = 1, minmn
        if (.not.pivot) go to 40
!
!  bring the column of largest norm into the pivot position.
        kmax = j
        do k = j, n
          if (rdiag(k) .gt. rdiag(kmax)) kmax = k
        enddo
        if (kmax .ne. j) then
          do i = 1, m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
          enddo
          rdiag(kmax) = rdiag(j)
          wa(kmax) = wa(j)
          k = ipvt(j)
          ipvt(j) = ipvt(kmax)
          ipvt(kmax) = k
        endif
 40     continue       
!
!  compute the householder transformation to reduce the
!  j-th column of a to a multiple of the j-th unit vector.
        ajnorm = enorm(m-j+1,a(j,j))
        if (ajnorm .eq. zero) go to 100
        if (a(j,j) .lt. zero) ajnorm = -ajnorm
        do i = j, m
          a(i,j) = a(i,j)/ajnorm
        enddo
        a(j,j) = a(j,j) + one
!
!  apply the transformation to the remaining columns and update the norms.
        jp1 = j + 1
        if (n .ge. jp1) then
          do k = jp1, n
            sum = zero
            do i = j, m
              sum = sum + a(i,j)*a(i,k)
            enddo
            temp = sum/a(j,j)
            do i = j, m
              a(i,k) = a(i,k) - temp*a(i,j)
            enddo
            if (.not.pivot .or. rdiag(k) .eq. zero) go to 80
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
            if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80
            rdiag(k) = enorm(m-j,a(jp1,k))
            wa(k) = rdiag(k)
 80         continue
          enddo
        endif
 100    continue
        rdiag(j) = -ajnorm
      enddo ! j
      return
!
      end subroutine qrfac
!      
!-----------------------------------------------------------------------
!
      subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)
      IMPLICIT NONE
      integer n,ldr, ipvt(n)
      real*8 delta,par, r(ldr,n),diag(n),qtb(n),x(n),sdiag(n),wa1(n),wa2(n)
!
!   Given an m by n matrix a, an n by n nonsingular diagonal matrix d, 
!   an m-vector b, and a positive number delta, the problem is to determine 
!   a value for the parameter par such that if x solves the system
!
!          a*x = b ,     sqrt(par)*d*x = 0 ,
!
!   in the least squares sense, and dxnorm is the euclidean norm of d*x, 
!   then either par is zero and
!
!           (dxnorm-delta) .le. 0.1*delta ,
!
!   or par is positive and
!
!           abs(dxnorm-delta) .le. 0.1*delta .
!
!   This subroutine completes the solution of the problem if it is provided 
!   with the necessary information from the qr factorization, with column pivoting,
!   of a. That is, if a*p = q*r, where p is a permutation matrix, q has orthogonal
!   columns, and r is an upper triangular matrix with diagonal elements of nonincreasing 
!   magnitude, then lmpar expects the full upper triangle of r, the permutation matrix p,
!   and the first n components of (q transpose)*b. On output lmpar also provides an upper
!   triangular matrix s such that
!
!            t   t                   t
!           p *(a *a + par*d*d)*p = s *s .
!
!   s is employed within lmpar and may be of separate interest.
!
!   Only a few iterations are generally needed for convergence of the algorithm. 
!   If, however, the limit of 10 iterations is reached, then the output par will 
!   contain the best value obtained so far.
!
!   The subroutine statement is
!
!     subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,wa1,wa2)
!
!  n is a positive integer input variable set to the order of r.
!
!  r is an n by n array. on input the full upper triangle must contain 
!    the full upper triangle of the matrix r. On output the full upper 
!    triangle is unaltered, and the strict lower triangle contains the 
!    strict upper triangle (transposed) of the upper triangular matrix s.
!
!  ldr is a positive integer input variable not less than n which specifies 
!    the leading dimension of the array r.
!
!  ipvt is an integer input array of length n which defines the permutation 
!    matrix p such that a*p = q*r. Column j of p is column ipvt(j) of the 
!    identity matrix.
!
!  diag is an input array of length n which must contain the diagonal elements 
!    of the matrix d.
!
!  qtb is an input array of length n which must contain the first n elements 
!    of the vector (q transpose)*b.
!
!  delta is a positive input variable which specifies an upper bound on 
!    the euclidean norm of d*x.
!
!  par is a nonnegative variable. on input par contains an initial estimate 
!    of the levenberg-marquardt parameter. On output par contains the final estimate.
!
!  x is an output array of length n which contains the least squares solution 
!    of the system a*x = b, sqrt(par)*d*x = 0, for the output par.
! 
!  sdiag is an output array of length n which contains the diagonal elements 
!    of the upper triangular matrix s.
!
!  wa1 and wa2 are work arrays of length n.
!
!  subprograms called
!       minpack-supplied ... dpmpar,enorm,qrsolv
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
      integer i,iter,j,jm1,jp1,k,l,nsing
      real*8 dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001,sum,temp,zero
!      real*8  dpmpar,enorm
      data p1,p001,zero /1.0d-1,1.0d-3,0.0d0/
!
!  dwarf is the smallest positive magnitude.
!      dwarf = dpmpar(2)
      dwarf = 2.22507385852d-308
!
!  compute and store in x the gauss-newton direction. if the
!  jacobian is rank-deficient, obtain a least squares solution.
      nsing = n
      do j = 1, n
        wa1(j) = qtb(j)
        if (r(j,j) .eq. zero .and. nsing .eq. n) nsing = j - 1
        if (nsing .lt. n) wa1(j) = zero
      enddo
      if (nsing .ge. 1) then
        do k = 1, nsing
          j = nsing - k + 1
          wa1(j) = wa1(j)/r(j,j)
          temp = wa1(j)
          jm1 = j - 1
          if (jm1 .ge. 1) then
            do i = 1, jm1
              wa1(i) = wa1(i) - r(i,j)*temp
            enddo ! i
          endif
        enddo  ! k
      endif
      do j = 1, n
        l = ipvt(j)
        x(l) = wa1(j)
      enddo
!
!  initialize the iteration counter.
!  evaluate the function at the origin, and test
!  for acceptance of the gauss-newton direction.
      iter = 0
      do j = 1, n
        wa2(j) = diag(j)*x(j)
      enddo ! j
      dxnorm = enorm(n,wa2)
      fp = dxnorm - delta
      if (fp .le. p1*delta) go to 220
!
!     if the jacobian is not rank deficient, the newton
!     step provides a lower bound, parl, for the zero of
!     the function. otherwise set this bound to zero.
      parl = zero
      if (nsing .ge. n) then
        do j = 1, n
          l = ipvt(j)
          wa1(j) = diag(l)*(wa2(l)/dxnorm)
        enddo ! j 
        do j = 1, n
          sum = zero
          jm1 = j - 1
          if (jm1 .ge. 1) then 
          do i = 1, jm1
            sum = sum + r(i,j)*wa1(i)
          enddo
          endif
          wa1(j) = (wa1(j) - sum)/r(j,j)
        enddo ! j
        temp = enorm(n,wa1)
        parl = ((fp/delta)/temp)/temp
      endif
!
!  calculate an upper bound, paru, for the zero of the function.
      do j = 1, n
        sum = zero
        do i = 1, j
          sum = sum + r(i,j)*qtb(i)
        enddo ! i
        l = ipvt(j)
        wa1(j) = sum/diag(l)
      enddo ! j
      gnorm = enorm(n,wa1)
      paru = gnorm/delta
      if (paru .eq. zero) paru = dwarf/dmin1(delta,p1)
!
!  if the input par lies outside of the interval (parl,paru), set par to the closer endpoint.
      par = dmax1(par,parl)
      par = dmin1(par,paru)
      if (par .eq. zero) par = gnorm/dxnorm
!
!  beginning of an iteration.
  150 continue
        iter = iter + 1
!
!  evaluate the function at the current value of par.
        if (par .eq. zero) par = dmax1(dwarf,p001*paru)
        temp = dsqrt(par)
        do j = 1, n
          wa1(j) = temp*diag(j)
        enddo ! j
        call qrsolv(n,r,ldr,ipvt,wa1,qtb,x,sdiag,wa2)
        do j = 1, n
          wa2(j) = diag(j)*x(j)
        enddo ! j
        dxnorm = enorm(n,wa2)
        temp = fp
        fp = dxnorm - delta
!
!  if the function is small enough, accept the current value of par. 
!  also test for the exceptional cases where parl
!  is zero or the number of iterations has reached 10.
        if (dabs(fp).le.p1*delta .or. parl.eq.zero .and. fp.le.temp .and. temp.lt.zero .or. iter.eq.10) go to 220
!
!  compute the newton correction.
        do j = 1, n
          l = ipvt(j)
          wa1(j) = diag(l)*(wa2(l)/dxnorm)
        enddo ! j
        do j = 1, n
          wa1(j) = wa1(j)/sdiag(j)
          temp = wa1(j)
          jp1 = j + 1
          if (n .ge. jp1) then
            do i = jp1, n
              wa1(i) = wa1(i) - r(i,j)*temp
            enddo ! i
          endif
        enddo ! j
        temp = enorm(n,wa1)
        parc = ((fp/delta)/temp)/temp
!
!  depending on the sign of the function, update parl or paru.
        if (fp .gt. zero) parl = dmax1(parl,par)
        if (fp .lt. zero) paru = dmin1(paru,par)
!
!  compute an improved estimate for par.
        par = dmax1(parl,par+parc)
!
!  end of an iteration.
        go to 150
  220 continue
!
!  termination.
      if (iter .eq. 0) par = zero
      return
!
      end subroutine lmpar
!      
!-----------------------------------------------------------------------
!
END MODULE f_e_fit
!    2165