MODULE f_e_check

USE f_e_data
USE f_e_parameters
!
!  This module makes various checks of the input data for consistency.
!   checkAltp
!   checkAssignments
!   checkBkq (not used?)
!   checkC2y
!   checkDiagOpts
!   checkTimeRev
!   checkESO
!   checkLinkFitP
!   checkPrint
!   checkValidERanges
!

IMPLICIT NONE

! PRIVATE

PUBLIC 

CONTAINS
!-----------------------------------------------------------------------
!
      SUBROUTINE checkAltp()
!
!   This subroutine checks to see if the A(ltp) parameters are consistent
!     Altp, l=2,4,6; t=l-1,l,l+1; p=0,+/-1,+/-2,...+/-t
!     Altp(3,3,8)  ! matrix size: 3x3x8=72. Only 45 values used.
!  Program stopped if a non-zero Altp is found outside allowed values.
!
!   The subroutine also checks to see if the A(ltp) parameters are consistent
!   with the point group being used. A warning is print if it is not.
!   
!
      IMPLICIT NONE
      integer*4 i,j,k,lamda,t,p
      integer*4 Acheck(35,37),Ncheck(8,7)
      logical DEBUG
      real*8 Z
      Z=0.0d+00
!
!  Acheck contains the possible Ckq for particular point groups. 
!  This is essentially Appendix 3 of Goerller-Walrand & Binnemans, HPCRE, 23,(1996).
!  groups marked * are not included (not in the above table)
!    kq
!    10,11
!    20,21,22
!    30,31,32,33
!    40,41,42,43,44
!    50,51,52,53,54,55
!    60,61,62,63,64,65,66
!    70,71,72,73,74,75,76,77  (2+3+4+5+6+7+8)=35
!
!   1. C1    9. C4    17. C3i   25. C6v   33. TdT  
!   2. Ci   10. S4    18. D3    26. D3h   34. OT*
!   3. C2   11. C4h   19. C3v   27. D6h   35. OhT 
!   4. Cs   12. D4    20. D3d   28. T*    36. C2Y*
!   5. C2h  13. C4v   21. C6    29. Th*   37. CsY*
!   6. D2   14. D2d   22. C3h   30. O*
!   7. C2v  15. D4h   23. C6h   31. Td
!   8. D2h  16. C3    24. D6    32. Oh
!
!  iAltp = 0,1,2 for Altp parameters all Real, all Imag, or Complex
!  Note:  Alt0 must be all Real for t odd
!                      all Imag for t even
!
! Ncheck(t,p) identifies where it is in Acheck.  Remember FORTRAN is column major!
      DATA Ncheck/ 1, 2, 0, 0, 0, 0, 0, 0,   3, 4, 5, 0, 0, 0, 0, 0,   6, 7, 8, 9, 0, 0, 0, 0,   &
                  10,11,12,13,14, 0, 0, 0,  15,16,17,18,19,20, 0, 0,  21,22,23,24,25,26,27, 0,   &
                  28,29,30,31,32,33,34,35/ 
!                C1q  C2q    C3q       C4q        C5q          C6q             C7q      
      DATA Acheck/1,1, 1,1,1, 1,1,1,1, 1,1,1,1,1, 1,1,1,1,1,1, 1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,   &
                  0,0, 1,1,1, 0,0,0,0, 1,1,1,1,1, 0,0,0,0,0,0, 1,1,1,1,1,1,1, 0,0,0,0,0,0,0,0,   &
                  1,0, 1,0,1, 1,0,1,0, 1,0,1,0,1, 1,0,1,0,1,0, 1,0,1,0,1,0,1, 1,0,1,0,1,0,1,0,   &
                  0,1, 1,0,1, 0,1,0,1, 1,0,1,0,1, 0,1,0,1,0,1, 1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,   &
                  0,0, 1,0,1, 0,0,0,0, 1,0,1,0,1, 0,0,0,0,0,0, 1,0,1,0,1,0,1, 0,0,0,0,0,0,0,0,   &
                  0,0, 1,0,1, 0,0,1,0, 1,0,1,0,1, 0,0,1,0,1,0, 1,0,1,0,1,0,1, 0,0,1,0,1,0,1,0,   &
                  1,0, 1,0,1, 1,0,1,0, 1,0,1,0,1, 1,0,1,0,1,0, 1,0,1,0,1,0,1, 1,0,1,0,1,0,1,0,   &
                  0,0, 1,0,1, 0,0,0,0, 1,0,1,0,1, 0,0,0,0,0,0, 1,0,1,0,1,0,1, 0,0,0,0,0,0,0,0,   &
                 
                  1,0, 1,0,0, 1,0,0,0, 1,0,0,0,1, 1,0,0,0,1,0, 1,0,0,0,1,0,0, 1,0,0,0,1,0,0,0,   &
                  0,0, 1,0,0, 0,0,1,0, 1,0,0,0,1, 0,0,1,0,1,0, 1,0,0,0,1,0,0, 0,0,1,0,0,0,1,0,   &
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,0,1, 0,0,0,0,0,0, 1,0,0,0,1,0,0, 0,0,0,0,0,0,0,0,   &
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,0,1, 0,0,0,0,1,0, 1,0,0,0,1,0,0, 0,0,0,0,1,0,0,0,   &
                  1,0, 1,0,0, 1,0,0,0, 1,0,0,0,1, 1,0,0,0,1,0, 1,0,0,0,1,0,0, 1,0,0,0,1,0,0,0,   &
                  0,0, 1,0,0, 0,0,1,0, 1,0,0,0,1, 0,0,1,0,0,0, 1,0,0,0,1,0,0, 0,0,1,0,0,0,1,0,   &
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,0,1, 0,0,0,0,0,0, 1,0,0,0,1,0,0, 0,0,0,0,0,0,0,0,   &
                  1,0, 1,0,0, 1,0,0,1, 1,0,0,1,0, 1,0,0,1,0,0, 1,0,0,1,0,0,1, 1,0,0,1,0,0,1,0,   &
                 
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,1,0, 0,0,0,0,0,0, 1,0,0,1,0,0,1, 0,0,0,0,0,0,0,0,   &
                  0,0, 1,0,0, 0,0,0,1, 1,0,0,1,0, 0,0,0,1,0,0, 1,0,0,1,0,0,1, 0,0,0,1,0,0,1,0,   &
                  1,0, 1,0,0, 1,0,0,1, 1,0,0,1,0, 1,0,0,1,0,0, 1,0,0,1,0,0,1, 1,0,0,1,0,0,1,0,   &
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,1,0, 0,0,0,0,0,0, 1,0,0,1,0,0,1, 0,0,0,0,0,0,0,0,   &
                  1,0, 1,0,0, 1,0,0,0, 1,0,0,0,0, 1,0,0,0,0,0, 1,0,0,0,0,0,1, 1,0,0,0,0,0,1,0,   &
                  0,0, 1,0,0, 0,0,0,1, 1,0,0,0,0, 0,0,0,1,0,0, 1,0,0,0,0,0,1, 0,0,0,1,0,0,0,0,   &
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,0,0, 0,0,0,0,0,0, 1,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,   &
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,0,0, 0,0,0,0,0,0, 1,0,0,0,0,0,1, 0,0,0,0,0,0,1,0,   &

                  1,0, 1,0,0, 1,0,0,0, 1,0,0,0,0, 1,0,0,0,0,0, 1,0,0,0,0,0,1, 1,0,0,0,0,0,1,0,   &
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,0,0, 0,0,0,0,0,1, 1,0,0,0,0,0,0, 0,0,0,0,0,1,0,0,   &
                  0,0, 1,0,0, 0,0,0,0, 1,0,0,0,0, 0,0,0,0,0,0, 1,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,   &
                  0,0, 0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,   &
                  0,0, 0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,   &
                  0,0, 0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,   &
                  0,0, 0,0,0, 0,0,1,0, 1,0,0,0,1, 0,0,0,0,0,0, 1,0,0,0,0,0,1, 0,0,1,0,0,0,1,0,   &
                  0,0, 0,0,0, 0,0,0,0, 1,0,0,0,1, 0,0,0,0,0,0, 1,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,   &

                  0,0, 0,0,0, 1,0,0,1, 1,0,0,1,0, 1,0,0,0,0,0, 1,0,0,1,0,0,1, 1,0,0,1,0,0,1,0,   &
                  0,0, 0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,   &
                  0,0, 0,0,0, 0,0,0,0, 1,0,0,1,0, 0,0,0,0,0,0, 1,0,0,1,0,0,1, 0,0,0,0,0,0,0,0,   &
                  0,0, 0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,   &
                  0,0, 0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0/
! TODO check the last 2 lines for C2Y and CsY
!
      DEBUG=.false.; AltpReal=.false.; AltpImag=.false.
      IF (IntN1.eq.1) then
      do i=1,3
        lamda=2*i
        do j=1,3
          t=lamda-2+j
          if (mod(t,2).eq.0) then
            if (altpR(i,j,1).ne.z) then
              write(io,'(" Alt0 must be pure Imaginary for t even. a",2I1,"0R(i,j,1)=",F10.2)') lamda,t,altpR(i,j,1) 
              stop              
            endif  
          else            
            if (altpI(i,j,1).ne.z) then
              write(io,'(" Alt0 must be pure real for t odd. A",2I1,"0I(i,j,1)=",F10.2)') lamda,t,altpI(i,j,1) 
              stop              
            endif           
          endif            
          do k=1,8
            if (altpR(i,j,k).ne.z) AltpReal=.true.
            if (altpI(i,j,k).ne.z) AltpImag=.true.
            if (altpR(i,j,k).ne.z .and. k-1.gt.t) then
              write(io,'("****FATAL: Invalid value for real(A(ltp)) lamda=",I2,"; t=",i2, & 
                         "; p=",i2" = ",E14.4)') lamda,t,k-1,altpR(i,j,k)
              stop
            elseif (altpI(i,j,k).ne.z .and. k-1.gt.t) then
              write(io,'("****FATAL: Invalid value for Imag(A(ltp)) lamda=",I2,"; t=",i2, & 
                         "; p=",i2" = ",E14.4)') lamda,t,k-1,altpI(i,j,k)
              stop
            endif
          enddo ! k
        enddo ! j
      enddo ! i
      endif ! IntN1
      
      IF (IntN1.eq.2) then
        do j=1,3
          if (BlkiI(1,j).ne.z) then
            write(io,'(" B20i must be pure real. Imag(B20i)= ",3F10.2)') (BlkiI(1,i),i=1,3); stop              
          elseif (BlkiI(4,j).ne.z) then
            write(io,'(" B40i must be pure real. Imag(B40i)= ",3F10.2)') (BlkiI(4,i),i=1,3); stop              
          elseif (BlkiI(9,j).ne.z) then
            write(io,'(" B60i must be pure real. Imag(B60i)= ",3F10.2)') (BlkiI(9,i),i=1,3); stop
          endif      
        enddo
      endif ! IntN1
      
      if (IDG.eq.0) return

      if (warnings(2).and.DEBUG) then
        write(io,'(/,"Ncheck")')
        do i=1,7
          write(io,'(8I3)') (Ncheck(j,i),j=1,8)
        enddo
        write(io,'("Acheck, Group=",A4)') Groups(idg)
        write(io,'(" C1q   C2q     C3q       C4q         C5q           C6q             C7q  ")')
        write(io,'(2i2,2x,3i2,2x,4i2,2x,5i2,2x,6i2,2x,7i2,2x,8i2)') (Acheck(j,idg),j=1,35)
      endif

      do i=1,3
        lamda=2*i
        do j=1,3
          t=lamda-2+j  !  1,2,3 or 3,4,5 or 5,6,7
          do k=1,8
            p=k-1
!            write(io,'("lamda,t,p=",3i2,", Ncheck(k,t)=",I2)') lamda,t,p,Ncheck(k,t)
            if (Ncheck(k,t).gt.0) then
              if (altpR(i,j,k).ne.z .or. altpI(i,j,k).ne.z) then
                if (Acheck(Ncheck(k,t),idg).eq.0) then
                  write(io,'("****FATAL: Invalid value for A(ltp) lamda=",I2,"; t=",i2,"; p=",i2,    &
                                              " for point group ",A4)') lamda,t,p,Groups(idg)
                  stop  
                else
                  if (warnings(2).and.DEBUG) write(io,'("A(ltp) lamda=",I2,"; t=",i2,"; p=",i2, & 
                                              " is valid:",2F10.2)') lamda,t,p,altpR(i,j,k),altpI(i,j,k)
                endif
              else 
                if ((Acheck(Ncheck(k,t),idg).ne.0) .and. warnings(2).and.DEBUG) then
                  write(io,'("A(ltp) lamda=",I2,"; t=",i2,"; p=",i2," can be non-zero.")') lamda,t,p
                endif 
              endif ! altpR(i,j,k)#0 or altpI(i,j,k)#0
            endif ! Ncheck(k,t)#0
          enddo ! k
        enddo ! j
      enddo ! i
        
      do i=1,3; do j=1,i+1
        if (ABS(AltpR(1,i,j)).gt.1000.0) write(io,'("****Warning: AltpR(2,",2i1,")=",F10.2)') i,j,AltpR(1,i,j)
        if (ABS(AltpI(1,i,j)).gt.1000.0) write(io,'("****Warning: AltpI(2,",2i1,")=",F10.2)') i,j,AltpI(1,i,j)
      enddo; enddo
      do i=1,3; do j=1,i+3 
        if (ABS(AltpR(2,i,j)).gt.1000.0) write(io,'("****Warning: AltpR(4,",2i1,")=",F10.2)') i,j,AltpR(1,i,j)
        if (ABS(AltpI(2,i,j)).gt.1000.0) write(io,'("****Warning: AltpI(4,",2i1,")=",F10.2)') i,j,AltpI(1,i,j)
      enddo; enddo
      do i=1,3; do j=1,i+5
        if (ABS(AltpR(3,i,j)).gt.1000.0) write(io,'("****Warning: AltpR(6,",2i1,")=",F10.2)') i,j,AltpR(1,i,j)
        if (ABS(AltpI(3,i,j)).gt.1000.0) write(io,'("****Warning: AltpI(6,",2i1,")=",F10.2)') i,j,AltpI(1,i,j)
      enddo; enddo
!      
      end subroutine checkAltp
!
!-----------------------------------------------------------------------
!
      subroutine checkAssignments()  
! Must be called after BlockBasis().
!  check if assignMethod is valid.
      IMPLICIT NONE
      integer*4 i,j 
      if (assignMethod.eq.1) then
        if (.not.block) then
          write(io,'("***FATAL: the BLOC command must used if assigning levels using crystal quantum numbers.")')
          stop
        endif
        do i=1,nexp
          if (NAssign(i,2).lt.1 .or. NAssign(i,2).gt.nblocks) then
            write(io,'("***FATAL: Invalid assignment to block:",i2," of level:",i2)') NAssign(i,2),i
            write(io,'("***nblocks=",i2)') nblocks
            stop
          endif
          if (i.ne.nexp) then
            do j=i+1,nexp
              if (NAssign(i,1).eq.NAssign(j,1) .and. NAssign(i,2).eq.NAssign(j,2)) then
                write(io,'("***FATAL: EXPE levels:",2i3," have the same assignment:",2i3)') i,j,NAssign(i,1),NAssign(i,2)
                stop
              endif !
            enddo ! j
          endif !
        enddo ! i
      else if (assignMethod.eq.2) then
        if (Lvalue.ne.2 .or. zeta.ne.0.0d+00) then
          write(io,'("***FATAL: Assigning levels using spin only valid for d-electrons, with SOC=0.")')
          stop
        endif  
        if (.not.option(2)) then
          write(io,'("***FATAL: spin-projections must be calculated (OPTN(2)) if assigning levels using spin.")')
          stop
        endif 
        do i=1,nexp
          if (Nelectrons.eq.1 .and. NAssign(i,2).ne.2) then
            write(io,'("***FATAL: Invalid assignment to 2S+1=",i2," of level:",i2)') NAssign(i,2),i; stop
          else if (Nelectrons.eq.2 .and. (NAssign(i,2).ne.1.and.NAssign(i,2).ne.3)) then
            write(io,'("***FATAL: Invalid assignment to 2S+1=",i2," of level:",i2)') NAssign(i,2),i; stop
          else if (Nelectrons.eq.3 .and. (NAssign(i,2).ne.2.and.NAssign(i,2).ne.4)) then
            write(io,'("***FATAL: Invalid assignment to 2S+1=",i2," of level:",i2)') NAssign(i,2),i; stop
          else if (Nelectrons.eq.4 .and. (NAssign(i,2).ne.1.and.NAssign(i,2).ne.3.and.NAssign(i,2).ne.5)) then
            write(io,'("***FATAL: Invalid assignment to 2S+1=",i2," of level:",i2)') NAssign(i,2),i; stop
          else if (Nelectrons.eq.5 .and. (NAssign(i,2).ne.2.and.NAssign(i,2).ne.4.and.NAssign(i,2).ne.6)) then
            write(io,'("***FATAL: Invalid assignment to 2S+1=",i2," of level:",i2)') NAssign(i,2),i; stop
          endif
          if (i.ne.nexp) then
            do j=i+1,nexp
              if (NAssign(i,1).eq.NAssign(j,1) .and. NAssign(i,2).eq.NAssign(j,2)) then
                write(io,'("***FATAL: EXPE levels:",2i3," have the same assignment:",2i3)') i,j,NAssign(i,1),NAssign(i,2)
                stop
              endif
            enddo ! j
          endif
        enddo ! i 
      endif  ! assignMethod
!      
      end subroutine checkAssignments
!
!-----------------------------------------------------------------------
!
      SUBROUTINE checkBkq()
!
!   This subroutine checks to see if the Bkq parameters are consistent
!   with the point group being used. A warning is print if it is not.
!   
!
      IMPLICIT NONE
      integer*4 i,Bcheck(37,32)
      real*8 Z,bit
      Z=0.0d+00; bit=1.0d-06
!  This is essentially Veven of Appendix 3 of Goerller-Walrand & Binnemans, HPCRE, 23,(1996).
!  
!   1. C1    9. C4    17. C3i   25. C6v   33. TdT  
!   2. Ci   10. S4    18. D3    26. D3h   34. OT*
!   3. C2   11. C4h   19. C3v   27. D6h   35. OhT 
!   4. Cs   12. D4    20. D3d   28. T*    36. C2Y*
!   5. C2h  13. C4v   21. C6    29. Th*   37. CsY*
!   6. D2   14. D2d   22. C3h   30. O*
!   7. C2v  15. D4h   23. C6h   31. Td
!   8. D2h  16. C3    24. D6    32. Oh
!  
!  Bckeck(i,j)= 0/1 : zero/non-zero
!         i the 36 point groups
!         j the 16+16 real and imaginary Bkq (B00 not included; B20',B40',B60' always zero)
!                B00  B2q    B4q         B6q           B00' B2q'    B4q'        B6q'  
      DATA Bcheck/0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1,  0, 0,1,1,  0,1,1,1,1,  0,1,1,1,1,1,1,   &  
                  0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1,  0, 0,1,1,  0,1,1,1,1,  0,1,1,1,1,1,1,   &
                  0, 1,0,1, 1,0,1,0,1,  1,0,1,0,1,0,1,  0, 0,0,1,  0,0,1,0,1,  0,0,1,0,1,0,1,   &
                  0, 1,0,1, 1,0,1,0,1,  1,0,1,0,1,0,1,  0, 0,0,1,  0,0,1,0,1,  0,0,1,0,1,0,1,   &
                  0, 1,0,1, 1,0,1,0,1,  1,0,1,0,1,0,1,  0, 0,0,1,  0,0,1,0,1,  0,0,1,0,1,0,1,   &
                  0, 1,0,1, 1,0,1,0,1,  1,0,1,0,1,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,1, 1,0,1,0,1,  1,0,1,0,1,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,1, 1,0,1,0,1,  1,0,1,0,1,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                                                      
                  0, 1,0,0, 1,0,0,0,1,  1,0,0,0,1,0,0,  0, 0,0,0,  0,0,0,0,1,  0,0,0,0,1,0,0,   &
                  0, 1,0,0, 1,0,0,0,1,  1,0,0,0,1,0,0,  0, 0,0,0,  0,0,0,0,1,  0,0,0,0,1,0,0,   &
                  0, 1,0,0, 1,0,0,0,1,  1,0,0,0,1,0,0,  0, 0,0,0,  0,0,0,0,1,  0,0,0,0,1,0,0,   &
                  0, 1,0,0, 1,0,0,0,1,  1,0,0,0,1,0,0,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,0,1,  1,0,0,0,1,0,0,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,0,1,  1,0,0,0,1,0,0,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,0,1,  1,0,0,0,1,0,0,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,1,0,  1,0,0,1,0,0,1,  0, 0,0,0,  0,0,0,1,0,  0,0,0,1,0,0,1,   &
                                                      
                  0, 1,0,0, 1,0,0,1,0,  1,0,0,1,0,0,1,  0, 0,0,0,  0,0,0,1,0,  0,0,0,1,0,0,1,   &
                  0, 1,0,1, 1,0,1,0,1,  1,0,1,0,1,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,1,0,  1,0,0,1,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,1,0,  1,0,0,1,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,0,0,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,1,   &
                  0, 1,0,0, 1,0,0,0,0,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,1,   &
                  0, 1,0,0, 1,0,0,0,0,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,1,   &
                  0, 1,0,0, 1,0,0,0,0,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                                                      
                  0, 1,0,0, 1,0,0,0,0,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,0,0,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,0,0, 1,0,0,0,0,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1,  0, 0,1,1,  0,1,1,1,1,  0,1,1,1,1,1,1,   &
                  0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1,  0, 0,1,1,  0,1,1,1,1,  0,1,1,1,1,1,1,   &
                  0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1,  0, 0,1,1,  0,1,1,1,1,  0,1,1,1,1,1,1,   &
                  0, 0,0,0, 1,0,0,0,1,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 0,0,0, 1,0,0,0,1,  1,0,0,0,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                                                      
                  0, 0,0,0, 1,0,0,1,0,  1,0,0,1,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1,  0, 0,1,1,  0,1,1,1,1,  0,1,1,1,1,1,1,   &
                  0, 0,0,0, 1,0,0,1,0,  1,0,0,1,0,0,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0,   &
                  0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1,  0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0/    
!                B00   B2q    B4q         B6q           B00'   B2q'    B4q'        B6q'             
! The last 2 lines are for C2Y and CsY                  
                                                        
!  This is essentially Appendix 3 of Goerller-Walrand & Binnemans, HPCRE, 23,(1996).
!  groups marked * are not included (not int the above)
!                  1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16   
!DATA BkqR_label/"B00 ","B20 ","B21 ","B22 ","B40 ","B41 ","B42 ","B43 ","B44 ","B60 ","B61 ","B62 ","B63 ","B64 ","B65 ","B66 "/
!DATA BkqI_label/"****","****","B21'","B22'","****","B41'","B42'","B43'","B44'","****","B61'","B62'","B63'","B64'","B65'","B66'"/
!
      if (BkqI(1).ne.z) then; write(io,'("****FATAL: Im(B00) must be 0")'); stop; endif
      if (BkqI(2).ne.z) then; write(io,'("****FATAL: Im(B20) must be 0")'); stop; endif
      if (BkqI(5).ne.z) then; write(io,'("****FATAL: Im(B40) must be 0")'); stop; endif
      if (BkqI(10).ne.z) then; write(io,'("****FATAL: Im(B60) must be 0")'); stop; endif
      if (IDG.eq.0) return
      if (warnings(1)) then
        do i=1,16
          if (Abs(BkqR(i)).gt.bit .and. Bcheck(idg,i).eq.0) then
            write(io,'("****The LF parameter: ",A4," should be zero for point group ",A4)') BkqR_label(i),Groups(idg)
          endif
          if (Abs(BkqI(i)).gt.bit .and. Bcheck(idg,i+16).eq.0) then
            write(io,'("****The LF parameter: ",A4," should be zero for point group ",A4)') BkqI_label(i),Groups(idg)
          endif
        enddo ! i
      endif
!      
      end subroutine checkBkq
!
!-----------------------------------------------------------------------
!
      SUBROUTINE checkC2y()
!
!   This subroutine checks to see if the Bkq parameters are consistent with a C2(y) symmetry element. 
!   Important for determining whether the ligand field and therefore matrix is real or complex.
!
      IMPLICIT NONE
      integer*4 i,Bcheck(32),n
      real*8 Z,bit
      Z=0.0d+00; bit=1.0d-06
!  Bckeck(i)= 0/1 : zero/non-zero
!         i the 16+16 real and imaginary Bkq (B00 not included; B20,B40,B60 always zero)
!                B00  B2q    B4q         B6q         B00'   B2q'    B4q'        B6q'  
      DATA Bcheck/0, 1,1,1, 1,1,1,1,1,  1,1,1,1,1,1,1, 0, 0,0,0,  0,0,0,0,0,  0,0,0,0,0,0,0/      ! C2Y             
!
      n=16; if (Lvalue.eq.2) n=9
      C2yExists=.true.
      do i=1,16
        if (Abs(BkqR(i)).gt.bit .and. Bcheck(i).eq.0) C2yExists=.false.
        if (Abs(BkqI(i)).gt.bit .and. Bcheck(i+16).eq.0) C2yExists=.false.
      enddo ! i
      if (.not.C2yExists .and. fastMat(2)) then
        write(io,'("****WARNING: no C2Y symmetry element; matrix will be complex.")')
      endif  
!     
      return 
      end subroutine checkC2y
!
!-----------------------------------------------------------------------
!
      SUBROUTINE checkDiagOpts()
!
!   This subroutine checks the diagonalisation options. 
!
      IMPLICIT NONE
      integer*4 i 
      real*8 Z,bit
      Z=0.0d+00; bit=1.0d-06
!      
      if (lancIt.gt.0) then
      endif
! 
      if (fastMat(1)) then
        write(io,'("****Warning: selective restriction of MJ blocks, not yet implemented")')
        write(io,'("****         FAST(1) reset to false")')
        fastMat(1)=.false.; endif  ! TODO
      if (fastMat(2).and. .not.n_odd) then
        write(io,'("****Warning: cannot block a Kramers doublet unless an odd number of electrons")')
        write(io,'("****         FAST(2) reset to false")')
        fastMat(2)=.false.; endif
      if (fastMat(2) .and. .not.fastMat(3) .and. calcVecs) then
        write(io,'("****Warning: Wavefunction properties may not be as expected as Kramers doublet halved and")')
        write(io,'("****         other component not regenerated: FAST(2)=T, FAST(3)=F")')
      endif
      if (fastMat(3).and. .not.fastMat(2)) then
        write(io,'("****Warning: Do not need to regenerate the Kramers doublet as FAST(2) not true")')
        write(io,'("****         FAST(3) reset to false")')
        fastMat(3)=.false.; endif
      if (fastMat(3).and. .not.n_odd) then
        write(io,'("****Warning: cannot regenerate a Kramers doublet unless an odd number of electrons")')
        write(io,'("****         FAST(3) reset to false")')
        fastMat(3)=.false.; endif
      do i=4,10
        if (fastMat(i)) then
        write(io,'("****Warning: fastMat(",i2,") not yet implemented")') i; fastMat(i)=.false.; endif
      enddo  
!      
      end subroutine checkDiagOpts
!
!-----------------------------------------------------------------------
!
      subroutine checkESO()
      integer f2J(13) 
      DATA f2J/5,8,9,8,5, 0,7,12,15,16,15,12,7/
      if (Nelectrons.eq.0) return   ! CONF command not called yet.
      if (ESO2J.ne.f2J(Nelectrons)) then
        write(io,'("***WARNING: 2J in EXSO command =",I4," does not equal the ground state 2J value (",I2,") for this f^", &
                   i2," electron configuration.",/,  &
                   "            Make sure you know what you''re doing!")') ESO2J, f2J(Nelectrons),Nelectrons
      endif
      return
!      
      end subroutine checkESO
!
!-----------------------------------------------------------------------
!
      subroutine checkLinkFitP()
!  Check the LINK and FIT commands.
! fitType=1   Fit Eng/Int 
! fitType=2   Fit Judd-Ofelt Eng/Int 
! fitType=3   Fit AOM parameters to Crystal Field Coefficients
! fitType=4   Fit calculated MEs to ESO MEs
! fitType=5   Fit g-values
! fitType=6   Fit calculated MEs to one electron LF MEs
      IMPLICIT none
      integer i,j,ifit,n
!      
      if (nlinks.gt.1) then
        do i=1,nlinks
          n=nlink(1,i)
          do j=i+1,nlinks
            if (nlink(1,j).eq.n) then
              write(io,'("***FATAL: Cannot have two LINK commands for parameter ",A9,"(number:",i2,")")') PLab(n),n
              stop
            endif
          enddo
        enddo
      endif
!
!      write(io,'("FitType=",i2,"; nfit=",I4,"; FitN(1)=",I4)') FitType,nfit,FitN(1)
      if (nfit.eq.0) return
!      
      if (FitType.eq.1 .or. FitType.eq.2) then
        if (nfit.eq.1 .and. FitN(1).eq.1) EaveOnly=.true.
        do i=1,nfit
          if (FitN(i).eq.1) then
            if (Option(1)) then
              write(io,'("***FATAL: Cannot set E0=0 (OPTN(1)=T) if varying EAVE in a fit")')
              stop
            endif
          endif  
        enddo  
      endif
!     
      If (Fittype.eq.3) then
        if (LFtype.eq."CF  ") then
          write(io,'("***FATAL: The ligand field cannot be specified by Bkq parameters to fit Bkq parameters.")') 
          write(io,'("          Use AOM or LF1E instead of CF.")'); stop; 
        endif      
      endif

      If (Fittype.eq.6) then
        if (LFtype.eq."LF1E" .or. LFtype.eq."CF  ") then
          write(io,'("***FATAL: The ligand field cannot be specified by one-e LF MEs or Bkq to fit one-e LF MEs.")') 
          write(io,'("          Use AOM instead of LF1E or CF.")'); stop; 
        endif      
      endif

      call loadP() 
      if (grid) then
        grid=.false.
        write(io,'("***WARNING: Cannot calculate a grid of variables if a fit is being made.",/," GRID set to .false.")')
      endif
      if (nvar.ne.0) then
        do i=1,nVar
          if (Nvarvals(i).gt.1) write(io,'("***WARNING: Only the first value of each variable is used for a fit.")')
          Nvarvals(i)=1 
          P(120+i)=allVarVals(i,1)
        enddo
      endif
      
      if (fitType.eq.5) then
        if (nfit.gt.nexpg) then
          write(io,'("***FATAL: Cannot fit more parameters (=",I2,") than experimental values (=",I2,")")') nfit,nexpg
          stop
        endif
        if (nfit.eq.nexpg) then
          write(io,'("***WARNING: Number of parameters fitted is equal to the number of experimental values (", &
                     I2,") Fit unlikely to be unique")') nexpg 
        endif
      else
        if (nfit.gt.nexp .and. fitType.ne.5) then
          write(io,'("***FATAL: Cannot fit more parameters (=",I2,") than experimental values (=",I2,")")') nfit,nexp
          stop
        endif
        if (nfit.eq.nexp) then
          write(io,'("***WARNING: Number of parameters fitted is equal to the number of experimental values (", &
                     I2,") Fit unlikely to be unique")') nexp 
        endif
      endif
      
      DO ifit=1,nfit
        n=FitN(ifit)
        
        If (Fittype.eq.3 .and. n.le.20) then
          write(io,'("***FATAL: Cannot vary parameter:",A4," to fit crystal field coefficients")') PLab(n)
          stop
        endif      

        if (nlinks.gt.0) then
          do i=1,nlinks
            if (nlink(1,i).eq.n) then
              write(io,'("***FATAL: Cannot have parameter ",A9,"(number:",i2,  &
                         ") both varied in a fit and determined by a LINK")') PLab(n),n
              stop
            endif
          enddo
        endif    
!        
        if (LFtype.eq."CF  ") then
          if (n.ge.36 .and. n.le.41) write(io,20) n 
          if (n.eq.44 .or. n.eq.49 .or. n.gt.55 .and. n.le.120) write(io,20) n 
        endif
        if (LFtype.eq."AOM ") then
          if (n.ge. 20+Nligands .and. n.le. 40) write(io,20) n 
          if (n.ge. 40+Nligands .and. n.le. 60) write(io,20) n 
          if (n.ge. 60+Nligands .and. n.le. 80) write(io,20) n 
          if (n.ge. 80+Nligands .and. n.le.100) write(io,20) n 
          if (n.ge.100+Nligands .and. n.le.120) write(io,20) n 
!          If (.not.AOM_Angles .and.(n.ge.61 .and. n.le.120)) then
!            write(io,'("***FATAL: Cannot vary the cartesian coordinates in a fit; cannot vary parameter:",I4)') FitN(i)
!            stop; 
!          endif
        endif        
          
        if (Nelectrons.eq.1) then
          if (n.gt.1  .and.  n.lt.5 ) write(io,10) ifit,Para_label(ifit),Nelectrons
          if (n.gt.5  .and.  n.lt.21) write(io,10) ifit,Para_label(ifit),Nelectrons
        endif
        if (Nelectrons.eq.2) then
          if (n.gt.16 .and. n.lt.21) write(io,10) ifit,Para_label(ifit),Nelectrons
        endif
        if (n.le.20) then
          if (P(n).eq.0.0d+00) write(io,'("****Warning: trying to vary parameter",I2,2X,A4, &
                                          " initially set to zero.")')n,Para_label(n)
          if (P(n).lt.fitL(ifit)) write(io,'("****Warning: initial value of parameter",I2,2X,A4, &
                               " set below lower limit of ",F8.2)')n,Para_label(n),FitL(ifit)
          if (P(n).gt.fitH(ifit)) write(io,'("****Warning: initial value of parameter",I2,2X,A4, &
                               " set above upper limit of ",F8.2)')n,Para_label(n),FitH(ifit)
        else if (n.gt.20 .and. n.le.35) then
          if (P(n).eq.0.0d+00) write(io,'("****Warning: trying to vary parameter",I2,2X,A4, &
                                          " initially set to zero.")')n,BkqR_label(n-20+1)
          if (P(n).lt.fitL(ifit)) write(io,'("****Warning: initial value of parameter",I2,2X,A4, &
                               " set below lower limit of ",F8.2)')n,BkqR_label(n-20+1),FitL(ifit)
          if (P(n).gt.fitH(ifit)) write(io,'("****Warning: initial value of parameter",I2,2X,A4, &
                               " set above upper limit of ",F8.2)')n,BkqR_label(n-20+1),FitH(ifit)
        else if (n.gt.35 .and. n.le.55) then
          if (P(n).eq.0.0d+00) write(io,'("****Warning: trying to vary parameter",I2,2X,A4, &
                                          " initially set to zero.")')n,BkqI_label(n-35+1)
          if (P(n).lt.fitL(ifit)) write(io,'("****Warning: initial value of parameter",I2,2X,A4, &
                               " set below lower limit of ",F8.2)')n,BkqI_label(n-35+1),FitL(ifit)
          if (P(n).gt.fitH(ifit)) write(io,'("****Warning: initial value of parameter",I2,2X,A4, &
                               " set above upper limit of ",F8.2)')n,BkqI_label(n-35+1),FitH(ifit)
        endif

        if (JuddOfelt) then
          if ((n.gt.20 .and. n.lt.120).or.n.gt.304) then
            write(io,'("****FATAL: Cannot vary parameter:",I3," in a Judd-Ofelt fit.")') n
            stop
          endif
        else 
          if (n.eq.301 .or. n.eq.302 .or.n.eq.303) then
            write(io,'("****FATAL: Cannot vary parameter:",I3," if you are not doing a Judd-Ofelt fit.")') n
            stop
          endif        
        endif

      enddo ! ifit

 10   format("****Warning: Cannot vary parameter",i3,2x,A4," for",i2," electrons" )
 20   format("****Warning: Invalid parameter number:",I4)
! 
      return
      end subroutine checkLinkFitP
!
!-----------------------------------------------------------------------
!
      subroutine checkPrint()  
!
! Checks the print options.
! Also OUTP
!
      IMPLICIT none
      integer i
!      
      if (printMat2(3).and. .not.JuddOfelt) then
        write(io,'("****Warning: Cannot have PRT2(3) without JUDO command")'); printMat2(3)=.false.; endif
      if (printMat2(4).and. (IMD.eq.0.and.ngexs.eq.0)) then
        write(io,'("****Warning: Cannot have PRT2(4) without MDIP or GEXS command")'); printMat2(4)=.false.; endif
      if (printMat2(5).and. (IMD.eq.0.and.ngexs.eq.0)) then
        write(io,'("****Warning: Cannot have PRT2(5) without MDIP or GEXS command")'); printMat2(5)=.false.; endif
      if (printMat3(3).and. IED.EQ.0) then
        write(io,'("****Warning: Cannot have PRT3(3) without EDIP command")'); printMat3(3)=.false.; endif
      if (printMat3(4).and. IED.EQ.0) then
        write(io,'("****Warning: Cannot have PRT3(4) without EDIP command")'); printMat3(4)=.false.; endif
      if (printMat3(5).and. IED.EQ.0) then
        write(io,'("****Warning: Cannot have PRT3(5) without EDIP command")'); printMat3(5)=.false.; endif
      if (printMat3(6).and. (IMD.eq.0.and.ngexs.eq.0)) then
        write(io,'("****Warning: Cannot have PRT3(6) without MDIP or GEXS command")'); printMat3(6)=.false.; endif
      if (printMat3(7).and. (IMD.eq.0.and.ngexs.eq.0)) then
        write(io,'("****Warning: Cannot have PRT3(7) without MDIP or GEXS command")'); printMat3(7)=.false.; endif
      if (printMat3(8).and. nPRED.eq.0) then
        write(io,'("****Warning: Cannot have PRT3(8) without PRED command")'); printMat3(8)=.false.; endif
      if (printMat4(1).and. ESO2J.eq.0) then
        write(io,'("****Warning: Cannot have PRT4(1) without EXSO command")'); printMat4(1)=.false.; endif
      if (printMat4(2).and. ESO2J.eq.0) then
        write(io,'("****Warning: Cannot have PRT4(2) without EXSO command")'); printMat4(2)=.false.; endif
!
      if (PLOTout1(5).and. IMD.eq.0) then
        write(io,'("****Warning: Cannot have PLOT(5) without MDIP command")'); PLOTout1(5)=.false.; endif
      if (PLOTout1(6).and. IED.eq.0) then
        write(io,'("****Warning: Cannot have PLOT(6) without EDIP command")'); PLOTout1(6)=.false.; endif
      if (PLOTout1(7).and. Nfit.eq.0) then
        write(io,'("****Warning: Cannot have PLOT(7) without FIT command")'); PLOTout1(7)=.false.; endif
      if (PLOTout1(9).and. ngexs.eq.0) then
        write(io,'("****Warning: Cannot have PLOT(9) without GEXS command")'); PLOTout1(9)=.false.; endif
!
      if (OUTP(1).and. LFtype.ne."AOM ") then
        write(io,'("****Warning: Cannot have OUTP(1) without AOM command")'); OUTP(1)=.false.; endif
      if (OUTP(2).and. LFtype.ne."AOM " .and. LFtype.ne."LF1E") then
        write(io,'("****Warning: OUTP(2) only appropriate for AOM or LF1E command")'); OUTP(2)=.false.; endif
      if (OUTP(3).and. LFtype.ne."AOM ") then
        write(io,'("****Warning: OUTP(3) only appropriate for AOM command")'); OUTP(3)=.false.; endif
      if (OUTP(4).and. IntN1.eq.0) then
        write(io,'("****Warning: OUTP(4) only appropriate if ALTP command given")'); OUTP(4)=.false.; endif
      if (OUTP(9).and. LFtype.ne."AOM ") then
        write(io,'("****Warning: Cannot have OUTP(9) without AOM command")'); OUTP(9)=.false.; endif
!
      return
      end subroutine checkPrint
!
!-----------------------------------------------------------------------
!
      SUBROUTINE checkTimeRev()
!
!   This subroutine checks to see if the Bkq parameters are consistent with 
!   being able to split the matrix due to Time Reversal symmetry.
!
      IMPLICIT NONE
      integer*4 i,Bcheck(32),n
      real*8 Z,bit
      Z=0.0d+00; bit=1.0d-06
!  Bckeck(i)= 0/1 : zero/non-zero
!         i the 16+16 real and imaginary Bkq (B00 not included; B20,B40,B60 always zero)
!                 B00  B2q    B4q         B6q          B00'  B2q'    B4q'        B6q'  
      DATA Bcheck/0,  1,0,1, 1,0,1,0,1,  1,0,1,0,1,0,1,  0, 0,0,1,  0,0,1,0,1,  0,0,1,0,1,0,1/      ! C2Y             
!
      TR_Exists=.true.
      n=16; if (Lvalue.eq.2) n=9
      do i=1,n
        if (Abs(BkqR(i)).gt.bit .and. Bcheck(i).eq.0) TR_Exists=.false.
        if (Abs(BkqI(i)).gt.bit .and. Bcheck(i+16).eq.0) TR_Exists=.false.
      enddo ! i
      if (.not.TR_Exists .and. fastMat(2)) then
        write(io,'("****FATAL: Bkq (q odd) non-zero; cannot do Time-reversal symmetry blocking")')
        STOP
      endif  
!      
      end subroutine checkTimeRev
!
!-----------------------------------------------------------------------
!
      subroutine checkValidERanges(Fn,nExplicitbasis,nExBasisNums)  
! Must be called after getLSBasis(), to have the TS_full values set.
! FitType
!    1   EXPE    Energies/Intensities  NEXP
!    2   EXPE    Judd-Ofelt            NEXP
!    3   EXPB    Bkq data              NEXP
!    4   EXPM    ESO matrix elements    -
!    5   EXPG    g-values              NexpG
!    6   EXP1    one electron LF MEs   1/2*N*(N+1)   N = 2L+1
      IMPLICIT none
      logical Fn
      integer i,j,n,nExplicitbasis,nExBasisNums(max_Jstates),ierr
      real*8 z
      parameter(z=0.0d+00)
!      
!  check experimental energy range; etc for FIT
      if (nexp.gt.0) then
        if (FitType.eq.1 .or. FitType.eq.2) then
          if (assignMethod.eq.0) then
          do i=1,nexp
            do j=1,Nass(i)
              if (NAssign(i,j).lt.1 .or. NAssign(i,j).gt.TS_full(Nelectrons)) then
                write(io,'("***FATAL: assigned level",i3," is not in range [1,",i3,"].")') NAssign(i,j),TS_full(Nelectrons)
                stop
              endif
            enddo ! j
            if (i.ne.nexp) then
              if (Nassign(i,Nass(i)).ge.Nassign(i+1,1)) then
                write(io,'("***FATAL: cannot have overlapping ranges: Nassign(",i3,",",i2,")=",i3," >= Nassign(",i3,", 1)=",I3)') &
                                                                               i,Nass(i),Nassign(i,Nass(i)), i+1,Nassign(i+1,1)
                stop
              endif
            endif            
          enddo ! i
          endif ! assignMethod.eq.0
        else if (FitType.eq.3) then  
          do i=1,nexp
            if (abs(Nassign(i,1)).gt.16) then
               write(io,'("***FATAL: cannot fit to a Crystal Field level with i,NASSIGN(i)=",2i4)') i,Nassign(i,1)
            endif 
          enddo  
        endif  
      else if (FitType.eq.4) then  
        if (ESO2J.eq.0) then
          write(io,'("***FATAL: cannot fit with no values to fit against. You must set the EXSO command.")') 
          stop
        endif  
        if (nPRED.eq.0) then
          write(io,'("***WARNING: Without using the PRED command, your fit may not be valid.")') 
!          stop
        endif  
      endif
!      
!  check the Fit settings:
      if (FitType.eq.1 .or. FitType.eq.2 .or. FitType.eq.3) then  
        if (nexp.eq.0) then
          write(io,'("***FATAL: cannot fit with no values to fit against. Nexp=",i4)') nexp
          stop
         endif  
      endif
      if (FitType.eq.5) then  
        if (.not.N_odd) then
          write(io,'("***WARNING: fitting g-values for even electron systems")') 
        endif
        if (NexpG.eq.0) then
          write(io,'("***FATAL: cannot fit with no values to fit against. NexpG=",i4)') NexpG
          stop
         endif  
      endif
      if (FitType.eq.6) then  
        if (Lvalue.eq.2 .and. nexp.ne.15) then
          write(io,'("***FATAL: to fit real d orbital 1e- LF MEs, you need exactly 15 values; n=",i4)') nexp
          stop
        endif
        if (Lvalue.eq.3 .and. nexp.ne.28) then
          write(io,'("***FATAL: to fit real f orbital 1e- LF MEs, you need exactly 28 values; n=",i4)') nexp
          stop
        endif
      endif
! 
      if (nExplicitBasis.ne.0) then
        if (Lvalue.eq.2) then
          write(io,'("***FATAL: Cannot use EXPL command for d-electrons.")') 
          STOP
        endif
        do i=1,nExplicitBasis
          if (nExBasisNums(i).lt.1 .or. nExBasisNums(i).gt.TSLJ_f(Nelectrons)) then
            write(io,'("***FATAL: Invalid values for basis(",i3,")=",I3," in EXPL command.")') nExplicitBasis, nExBasisNums(i)
            STOP
          endif
        enddo  
      endif      
!      
!  check g-values range
      if (ngexs.ne.0) then
        n=TS_full(Nelectrons)
        do i=1,ngexs
          if (gexs(2*i)-gexs(2*i-1).ne.1) then
            write(io,'("***FATAL: Invalid values for GEXS(i),GEXS(i+1)=",2I4," values must be consecutive.")') &
                        gexs(2*i-1),gexs(2*i); stop
          endif
          if (gexs(2*i-1).le.0 .or. gexs(2*i).le.0) then
            write(io,'("***FATAL: Invalid values for GEXS(i),GEXS(i+1)=",2I4," values must be > 0.")') gexs(2*i-1),gexs(2*i); stop
          endif
          if (gexs(2*i-1).gt.n .or. gexs(2*i).gt.n) then
            write(io,'("***FATAL: Invalid values for GEXS(i),GEXS(i+1)=",2I4," values must be <",I4)') gexs(2*i-1),gexs(2*i),n; stop
          endif
        enddo
      endif
      if (ngexs.ne.0 .and. zeta.eq.z) then; write(io,'("****WARNING: zeta=0; g-values are likely to be nonsense")'); endif
!
!  check EDIP range
      if (IED.gt.0) then
        if (IED.gt.3) then; write(io,'("***FATAL: Invalid value for IED =",I3)') IED; stop; endif
        if (EDunits.lt.1 .or. EDunits.gt.4) then; write(io,'("***FATAL: Invalid value for EDunits =",I2)') EDunits; stop; endif
        if (EDunits.eq.1 .or. EDunits.eq.4) then; write(io,'("***FATAL: EDunits =",I2," is not implemented yet.")') EDunits; 
          stop; endif
        if (ED(1,1).lt.1 .or. ED(2,1).lt.1 .or. ED(1,2).lt.1 .or. ED(2,2).lt.1) THEN
          write(io,'("***FATAL: Invalid ranges in EDIP N,N1-N4=",4I5)') IED,ED(1,1),ED(2,1),ED(1,2),ED(2,2); stop
        ENDIF
        n=TS_full(Nelectrons)
        if (ED(1,1).gt.n .or. ED(2,1).gt.n .or. ED(1,2).gt.n .or. ED(2,2).gt.n) THEN
          write(io,'("***FATAL: Invalid ranges in EDIP N,N1-N4=",4I5)') IED,ED(1,1),ED(2,1),ED(1,2),ED(2,2); stop
        ENDIF
      ENDIF  
!
!  check MDIP range
      if (IMD.gt.0) then
        if (IMD.gt.3) then; write(io,'("***FATAL: Invalid value for IMD =",I3)') IMD; stop; endif 
        if (MDunits.lt.1 .or. MDunits.gt.3) then; write(io,'("***FATAL: Invalid value for MDunits =",I3)') MDunits; stop; endif
        if (MD(1,1).lt.1 .or. MD(2,1).lt.1 .or. MD(1,2).lt.1 .or. MD(2,2).lt.1) THEN
          write(io,'("***FATAL: Invalid ranges in MDIP N,N1-N4=",4I5)') IMD,MD(1,1),MD(2,1),MD(1,2),MD(2,2); stop
        ENDIF
        n=TS_full(Nelectrons)
        if (MD(1,1).gt.n .or. MD(2,1).gt.n .or. MD(1,2).gt.n .or. MD(2,2).gt.n) THEN
          write(io,'("***FATAL: Invalid ranges in MDIP N,N1-N4=",4I5)') IMD,MD(1,1),MD(2,1),MD(1,2),MD(2,2); stop
        ENDIF
        if (.not.FullMJBasis) then
          n=TSLJ_f(Nelectrons)
          if (MD(1,1).gt.n .or. MD(2,1).gt.n .or. MD(1,2).gt.n .or. MD(2,2).gt.n) THEN
            write(io,'("***FATAL: Invalid N1 N2 N3 N4 =",4I3,"  in MDIP; for |SLJ> basis maximum number of states:",I3)') &
              MD(1,1),MD(2,1),MD(1,2),MD(2,2),n
            STOP                                             
          endif
        endif  
      ENDIF  
!
!  check JUDO range 
      if (JuddOfelt) then
        n=TSLJ_f(Nelectrons)
        if (ED(1,1).lt.1 .or. ED(2,1).lt.1 .or. ED(1,2).lt.1 .or. ED(2,2).lt.1) then
          write(io,'("***FATAL: Invalid N1 N2 N3 N4 =",4I3,"  in JUDO")') ED(1,1),ED(2,1),ED(1,2),ED(2,2)
          STOP                                             
        endif
        if ((ED(1,1).gt.n .or. ED(2,1).gt.n .or. ED(1,2).gt.n .or. ED(2,2).gt.n)) then
          write(io,'("***FATAL: N1 N2 N3 N4 =",4I4," greater than number of SLJ states:",I4," in JUDO")')  &
                                               ED(1,1),ED(2,1),ED(1,2),ED(2,2), n
          STOP                                             
        endif
        IF (ED(1,2)-ED(1,1).GT.max_GInt) then
          WRITE(IO,'("***FATAL: Array overflow. Too many ground state levels in JUDO command.",/,  &
                     " MAX(ED(1,2)-ED(1,1))=",I4)') max_GInt
          stop
        endif 
        if (ED(2,2)-ED(2,1).GT.max_EInt) then
          WRITE(IO,'("***FATAL: Array overflow. Too many excited state levels in JUDO command.",/,  &
                     " MAX(ED(2,2)-ED(2,1))=",I4)') max_EInt
          stop
        endif
        if (.not.option(8)) WRITE(IO,'("***WARNING: Judd-Ofelt calcuation to be made, you should set OPTN(8)=true.")')
      ENDIF  
!   Different atomic parameters in the prediagonalisation:      
      if (doAtomic) then
        if (nPRED.eq.0) then
          Write(iO,'("If the ATOM command is used, then the PRED command must also be specified")')  
          STOP
        endif
      Endif
!  check PRED range 
      if (nPRED.gt.0) then
        if (nPRED.gt.TS_full(Nelectrons)) then
          write(io,'("***FATAL: Invalid PRED =",I3)') nPRED
          STOP                                             
        endif
        if (vectLow.gt.nPRED .or. vectHigh.gt.nPRED) then
          write(io,'("***FATAL: Cannot print eigenvectors in range [",I4,",",I4,"] for nPRED=",I4)') vectLow,vectHigh,nPRED; stop
        endif
      endif  
!  check Offsets 
      if (offsetType.eq.1) then
        if (nPRED.eq.0) then
          write(io,'("***FATAL: Cannot have offsetType=1 in OFFS if PRED is not used.")'); stop
        endif
        do i=1,offsetN
          if (ioffset(i,1).gt.nPRED) then
            write(io,'("***FATAL: Cannot have offsetType=1 and offset level:",i2," if nPRED=",i4)') ioffset(i,1),nPRED; stop
          endif
          if (offsetN.gt.1) then
            do j=i+1,offsetN
              if (ioffset(i,1).gt.ioffset(j,1)) then
                write(io,'("***FATAL: Offsets cannot have overlapping ranges: offset(i,1),offset(i+1,1):",2i4)') ioffset(i,1), &
                                                                                                              ioffset(j,1); Stop
              endif
            enddo ! j
          endif

        enddo  
      endif
      if (offsetType.eq.2) then
        do i=1,offsetN
          if (ioffset(i,1).gt.TS_full(Nelectrons) .or. ioffset(i,2).gt.TS_full(Nelectrons)) then
            write(io,'("***FATAL: Cannot offset levels:",2i4," if basis is",i4)') ioffset(i,1),ioffset(i,2),TS_full(Nelectrons)
            stop
          endif
          if (ioffset(i,1).gt.ioffset(i,2)) then
            write(io,'("***FATAL: Offset(1) must be less or equal to offset(2):",2i4)') ioffset(i,1),ioffset(i,2); Stop
          endif
          if (offsetN.gt.1) then
            do j=i+1,offsetN
              if (ioffset(i,2).gt.ioffset(j,1)) then
                write(io,'("***FATAL: Offsets cannot have overlapping ranges: offset(i,2),offset(i+1,1):",2i4)') ioffset(i,2), &
                                                                                                              ioffset(j,1); Stop
              endif
            enddo ! j
          endif
        enddo  
      endif
!
! Check for mixed F2/F4 and B/C parameters 
      if (Racah) then
        if (Fn) then
          WRITE(IO,'("***FATAL: Cannot specify mixed F2/F4 and B/C parameters for electron repulsion.")'); stop
        endif
        IF (Nelectrons.eq.1) then
          WRITE(IO,'("***Warning: No electron-electron repulsion for 1e-, Racah parameters B,C set to 0.0,0.0")')
          Bracah=z; Cracah=z 
        else
          Para_label(2)="B   "; Para_label(3)="C   "
        endif
      endif
!  Check for unnecessary electron-electron repulsion parameters.      
      IF (Nelectrons.eq.1) then 
        if (Fn .and. .not.racah) then
          IF (F2.ne.z .or. F4.ne.z .or. F6.ne.z) then
            WRITE(IO,'("***Warning: No electron-electron repulsion for 1e-, F2,F4,F6 set to zero.")')
            F2=z; F4=z; F6=z
          endif
        endif
!  Check for unnecessary higher order Mn, Pn parameters.      
        IF (M0.ne.z .or. M2.ne.z .or. M4.ne.z .or. P2.ne.z .or. P4.ne.z .or. P6.ne.z) then
          WRITE(IO,'("***Warning: must have 2e- or greater for 2-electron parameters, all Mn,Pn set to zero.")')
          M0=z; M2=z; M4=z; P2=z; P4=z; P6=z
        endif
        IF (ALPHA.ne.z .or. BETA.ne.z .or. GAMMA.ne.z) then
          WRITE(IO,'("***Warning: must have 2e- or greater for 2-electron parameters, ALPHA, BETA, GAMMA set to zero.")')
          ALPHA=z; BETA=z; GAMMA=z
        endif
      endif        
!  Check for unnecessary Tn parameters.      
      IF (Nelectrons.le.2) then 
        IF (T2.ne.z .or. T3.ne.z .or. T4.ne.z .or. T6.ne.z .or. T7.ne.z .or. T8.ne.z) then
          WRITE(IO,'("***Warning: must have 3e- or greater for 3-electron parameters, all Tn set to zero.")')
          T2=z; T3=z; T4=z; T6=z; T7=z; T8=z
        endif
      endif
!  Check values for printing Eigenvectors:          
      if (vectLow.ne.0) then
        n=TS_full(Nelectrons)
        if (vectLow.le.0 .or. vectLow.gt.n) WRITE(IO,'("***Warning: Invalid value for vectLow=",I4," in VECT command.")') vectLow
        if (vectHigh.le.0.or.vectHigh.gt.n) WRITE(IO,'("***Warning: Invalid value for vectHigh=",I4," in VECT command.")')vectHigh
        if (vectN.le.0 .or. vectN.gt.n) vectN=n
      endif
!  Check CCF 
      if (CCFtype.ne.0 .and. LFtype.ne."CF  ") then
        WRITE(IO,'("***Fatal: Ligand field must be specified by the CF command to use a correlate crystal field")'); stop
        STOP
      endif 
!  Check LF rotation 
      if (OUTP(5).and. abs(RotLF(1))+abs(RotLF(2))+abs(RotLF(3)).eq.0.0d+00) then
        write(io,'("****Warning: OUTP(5) only appropriate if ROTL command given")'); OUTP(5)=.false.
      endif
!  Check Print matrix ranges:
      if (PrintMat1(10)) then
        ierr=0
        if (PrintRng1(1).lt.1 .or. PrintRng1(2).le.PrintRng1(1) .or. PrintRng1(2).gt.TS_n(Nelectrons)  &
                                                                .or. PrintRng1(2)-PrintRng1(1).gt.20) ierr=ierr+1
        if (PrintRng1(3).lt.1 .or. PrintRng1(4).le.PrintRng1(3) .or. PrintRng1(4).gt.TS_n(Nelectrons)  &
                                                                .or. PrintRng1(4)-PrintRng1(3).gt.20) ierr=ierr+1
        if (ierr.gt.0) then
          write(io,'("***Fatal: Invalid range for matrix printing in PRT1; i1,i2,j1,j2=",4i3)') (PrintRng1(i),i=1,4); stop
        endif
      endif  
      if (PrintMat2(10)) then
        ierr=0
        if (PrintRng2(1).lt.1 .or. PrintRng2(2).le.PrintRng2(1) .or. PrintRng2(2).gt.TSLJ_f(Nelectrons)  &
                                                                .or. PrintRng2(2)-PrintRng2(1).gt.20) ierr=ierr+1
        if (PrintRng2(3).lt.1 .or. PrintRng2(4).le.PrintRng2(3) .or. PrintRng2(4).gt.TSLJ_f(Nelectrons)  &
                                                                .or. PrintRng2(4)-PrintRng2(3).gt.20) ierr=ierr+1
        if (ierr.gt.0) then
          write(io,'("***Fatal: Invalid range for matrix printing in PRT2; i1,i2,j1,j2=",4i3)') (PrintRng2(i),i=1,4); stop
        endif
      endif  
      if (PrintMat3(10)) then
        ierr=0
        if (PrintRng3(1).lt.1 .or. PrintRng3(2).le.PrintRng3(1) .or. PrintRng3(2).gt.TS_full(Nelectrons)) ierr=ierr+1
        if (PrintRng3(3).lt.1 .or. PrintRng3(4).le.PrintRng3(3) .or. PrintRng3(4).gt.TS_full(Nelectrons)) ierr=ierr+1
        if (ierr.gt.0) then
          write(io,'("***Fatal: Invalid range for matrix printing in PRT3; i1,i2,j1,j2=",4i3)') (PrintRng3(i),i=1,4); stop
        endif
      endif  
!
      return
      end subroutine checkValidERanges
!
!-----------------------------------------------------------------------

END MODULE f_e_check
!  1035