MODULE f_e_fncrosserrors
USE f_e_data
CONTAINS
!
!--------------------------------------------------------------------------------
!
      SUBROUTINE FnCrossErrors(UseFncrossErrors)
!  Inserts the ME errors of fncross
! UseFncrossErrors(k)  k= 1    2    3    4    5    6      7      8
!                         Un  V11   Fn  abg  MnPn  Tn  wrongT2  fnmpUsed
!
!  An error is if |(MEcorr-MEfncross)/MEcorr| > 0.001  for MEcorr > 1.0e-8
!                     or |(MEcorr-MEfncross)| > 0.001  for MEcorr <= 1.0e-8
!
!  Some MEs that are incorrect in fncross file, but correct in the fnmp files.
!  fnmpUsed = .true. means that the fnmp files have been used.
!             (Only a subset of the incorrect MEs in fncross are incorrect in fnmp.)
!
!  The in correct MEs for Mn abnd Pn are given in Chen, J.Luminescence,128,421,(2008)
!  The MEs in Table 1 are incorrect in both fncross and fnmp files.
!  The MEs in Table 2 are incorrectly all zeros in fncross, but are correct in fnmp. 
!
      IMPLICIT none
      integer i,j,k,Lij,Sij
      logical UseFncrossErrors(8),F8
      real*8 D1
      parameter(D1=1.0d+00)
      F8=.not.UseFncrossErrors(8)
! 100  FORMAT('**** WARNING: Possible errors in matrix elements:',A6,' have NOT been introduced.')
! 110  FORMAT('**** WARNING: Errors HAVE been introduced into',A6,' matrix elements.')
! 120  FORMAT('**** WARNING: There are no errors in the ',A6,' matrix elements.')
!--------------------------------------------------------------------------------
!   Errors in the f^1 configuration
      IF (nelectrons.eq.1) then
        if (UseFncrossErrors(1)) then
          WRITE(IO,120) '  Un  '
        endif
        if (UseFncrossErrors(2)) then
          WRITE(IO,120) ' V11  '
        endif
        if (UseFncrossErrors(3)) then
          WRITE(IO,120) '  Fn  '
        endif
        if (UseFncrossErrors(4)) then
          WRITE(IO,120) ' abg  '
        endif
        if (UseFncrossErrors(5)) then
          WRITE(IO,120) ' Mn/Pn'
        endif
        if (UseFncrossErrors(6)) then
          WRITE(IO,120) '  Tn  '
        endif
! ***Total disagreements for f^1 matrix elements (greater than 0.10E-02) is   0
!
!--------------------------------------------------------------------------------
!   Errors in the f^2 configuration
      ELSE IF (nelectrons.eq.2) then
        if (UseFncrossErrors(1)) then
          WRITE(IO,120) '  Un  '
        endif
        if (UseFncrossErrors(2)) then
          WRITE(IO,120) ' V11  '
        endif
        if (UseFncrossErrors(3)) then
          WRITE(IO,120) '  Fn  '
        endif
        if (UseFncrossErrors(4)) then
          WRITE(IO,120) ' abg  '
        endif
        if (UseFncrossErrors(5)) then
          MnMat(  3,  7,1)=  -41.00461400 !  <3H | M0  |1I > should be  -39.51740123
          MnMat(  3,  7,2)=    0.70820000 !  <3H | M2  |1I > should be  -12.63166198
          MnMat(  3,  7,3)=   -8.03126500 !  <3H | M4  |1I > should be   -0.63616184
          PnMat(  3,  7,1)=   -0.21482000 !  <3H | P2  |1I > should be    0.09914760
          PnMat(  3,  7,2)=    0.16661200 !  <3H | P4  |1I > should be   -0.00409701
          PnMat(  3,  7,3)=   -0.03565700 !  <3H | P6  |1I > should be    0.01030313
          WRITE(IO,110) ' Mn/Pn'
        else
          WRITE(IO,100) ' Mn/Pn' 
        endif
        if (UseFncrossErrors(6)) then
          WRITE(IO,120) '  Tn  '
        endif
! ***Total disagreements for f^2 matrix elements (greater than 0.10E-02) is   6 for fncross & fnmp
!
!--------------------------------------------------------------------------------
!   Errors in the f^3 configuration
      ELSE IF (nelectrons.eq.3) then
        if (UseFncrossErrors(1)) then
          WRITE(IO,120) '  Un  '
        endif
        if (UseFncrossErrors(2)) then
          WRITE(IO,120) ' V11  '
        endif
        if (UseFncrossErrors(3)) then
          WRITE(IO,120) '  Fn  '
        endif
        if (UseFncrossErrors(4)) then
          WRITE(IO,120) ' abg  '
        endif
        if (UseFncrossErrors(5)) then
          if (F8) MnMat(  2, 10,2)=    0.00000000 !  <4D | M2  |2F2> should be  -38.40594791
          if (F8) MnMat(  2, 10,3)=    0.00000000 !  <4D | M4  |2F2> should be  -58.78846572
          MnMat(  5, 14,1)=    0.00000000 !  <4I | M0  |2H2> should be  -61.35506861
          MnMat(  5, 14,2)=    0.00000000 !  <4I | M2  |2H2> should be   13.96123600
          MnMat(  5, 14,3)=    0.00000000 !  <4I | M4  |2H2> should be   17.39510576
          MnMat(  5, 15,1)=  -61.35506500 !  <4I | M0  |2I > should be  -49.29052196
          MnMat(  5, 15,2)=   13.96123000 !  <4I | M2  |2I > should be  -16.68770023
          MnMat(  5, 15,3)=   17.39511300 !  <4I | M4  |2I > should be    4.23027366
          MnMat(  5, 16,1)=  -49.29051900 !  <4I | M0  |2K > should be   70.23009801
          MnMat(  5, 16,2)=  -16.68769900 !  <4I | M2  |2K > should be   17.18196248
          MnMat(  5, 16,3)=    4.23027400 !  <4I | M4  |2K > should be   15.20632241
          if (F8) MnMat(  5, 17,1)=   70.23009600 !  <4I | M0  |2L > should be    0.00000000
          if (F8) MnMat(  5, 17,2)=   17.18195600 !  <4I | M2  |2L > should be    0.00000000
          if (F8) MnMat(  5, 17,3)=   15.20632900 !  <4I | M4  |2L > should be    0.00000000
          PnMat(  5, 14,1)=    0.00000000 !  <4I | P2  |2H2> should be   -0.40903379
          PnMat(  5, 14,2)=    0.00000000 !  <4I | P4  |2H2> should be   -0.03329226
          PnMat(  5, 14,3)=    0.00000000 !  <4I | P6  |2H2> should be    0.16706743
          PnMat(  5, 15,1)=   -0.40903400 !  <4I | P2  |2I > should be   -0.20773783
          PnMat(  5, 15,2)=   -0.03329200 !  <4I | P4  |2I > should be   -0.04057989
          PnMat(  5, 15,3)=    0.16706700 !  <4I | P6  |2I > should be    0.10585958
          PnMat(  5, 16,1)=   -0.20773800 !  <4I | P2  |2K > should be    0.22951012
          PnMat(  5, 16,2)=   -0.04058000 !  <4I | P4  |2K > should be   -0.09246792
          PnMat(  5, 16,3)=    0.10586000 !  <4I | P6  |2K > should be    0.00491030
          if (F8) PnMat(  5, 17,1)=    0.22951000 !  <4I | P2  |2L > should be    0.00000000
          if (F8) PnMat(  5, 17,2)=   -0.09246800 !  <4I | P4  |2L > should be    0.00000000
          if (F8) PnMat(  5, 17,3)=    0.00491000 !  <4I | P6  |2L > should be    0.00000000
          WRITE(IO,110) ' Mn/Pn'
        else
          WRITE(IO,100) ' Mn/Pn'  
        endif
        if (UseFncrossErrors(6)) then
          if (F8) TnMat( 17, 17,5)=   -0.02605300 !  <2L | T7  |2L > should be   -0.02650326  (fncross error only)
          WRITE(IO,110) '  Tn  '
        else 
          WRITE(IO,100) '  Tn  '
        endif
! ***Total disagreements for f^3 matrix elements (greater than 0.10E-02) is  27 for fncross & 19 for fnmp
!
!--------------------------------------------------------------------------------
!   Errors in the f^4 configuration
      ELSE IF (nelectrons.eq.4) then
        if (UseFncrossErrors(1)) then
          WRITE(IO,120)'  Un  '
        endif
        if (UseFncrossErrors(2)) then
          WRITE(IO,120)' V11  '
        endif
        if (UseFncrossErrors(3)) then
          WRITE(IO,120)'  Fn  '
        endif
        if (UseFncrossErrors(4)) then
          WRITE(IO,120)' abg  '
        endif
        if (UseFncrossErrors(5)) then
          if (F8) MnMat(  6, 30,1)=    0.00000000 !  <3P1| M0  |1D1> should be  -14.40510326
          if (F8) MnMat(  6, 30,2)=    0.00000000 !  <3P1| M2  |1D1> should be   -7.83243257
          if (F8) MnMat(  6, 30,3)=    0.00000000 !  <3P1| M4  |1D1> should be   -6.24901645
          if (F8) PnMat(  6, 30,1)=    0.00000000 !  <3P1| P2  |1D1> should be   -0.08824419
          if (F8) PnMat(  6, 30,2)=    0.00000000 !  <3P1| P4  |1D1> should be   -0.27801069
          if (F8) PnMat(  6, 30,3)=    0.00000000 !  <3P1| P6  |1D1> should be   -0.09309581
          WRITE(IO,110) ' Mn/Pn'
        else
          WRITE(IO,100) ' Mn/Pn' 
        endif
        if (UseFncrossErrors(6)) then
          WRITE(IO,120) '  Tn  '
        endif
! ***Total disagreements for f^4 matrix elements (greater than 0.10E-02) is   6 for fncross & fnmp
!
!--------------------------------------------------------------------------------
!   Errors in the f^5 configuration
      ELSE IF (nelectrons.eq.5) then
        if (UseFncrossErrors(1)) then
          WRITE(IO,120) '  Un  '
        endif
        if (UseFncrossErrors(2)) then
          Vmat( 45, 45)   =   -0.02564700 !  <2G2| V11 |2G2> should be   -0.02561204
          WRITE(IO,110) ' V11  '
        else
          WRITE(IO,100) ' V11  ' 
        endif
        if (UseFncrossErrors(3)) then
          EEmat( 55, 56,3)=    0.00029100 !  <2H6| F6  |2H7> should be    0.00029139
          WRITE(IO,110) '  EE  '
        else
          WRITE(IO,100) '  EE  ' 
        endif
        if (UseFncrossErrors(4)) then
          WRITE(IO,120) ' abg  '
        endif
        if (UseFncrossErrors(5)) then
          if (F8) MnMat(  1,  9,2)=    0.00000000 !  <6P | M2  |4D3> should be  -37.71389217
          if (F8) MnMat(  1,  9,3)=    0.00000000 !  <6P | M4  |4D3> should be  -57.72912732
          if (F8) MnMat(  2, 15,1)=    0.00000000 !  <6F | M0  |4G2> should be  -32.59272487
          if (F8) MnMat(  2, 15,2)=    0.00000000 !  <6F | M2  |4G2> should be  -34.15018608
          if (F8) MnMat(  2, 15,3)=    0.00000000 !  <6F | M4  |4G2> should be  -35.88031483
          if (F8) MnMat( 10, 44,1)=    0.00000000 !  <4F1| M0  |2G1> should be  -13.02744356
          if (F8) MnMat( 10, 44,2)=    0.00000000 !  <4F1| M2  |2G1> should be  -15.68118206
          if (F8) MnMat( 10, 44,3)=    0.00000000 !  <4F1| M4  |2G1> should be  -26.53738503
          if (F8) MnMat( 14, 51,1)=    0.00000000 !  <4G1| M0  |2H2> should be  -19.33672384
          if (F8) MnMat( 14, 51,2)=    0.00000000 !  <4G1| M2  |2H2> should be  -25.98670357
          if (F8) MnMat( 14, 51,3)=    0.00000000 !  <4G1| M4  |2H2> should be  -33.32756305
          if (F8) MnMat( 20, 60,1)=    0.00000000 !  <4H3| M0  |2I4> should be  -39.06231110
          if (F8) MnMat( 20, 60,2)=    0.00000000 !  <4H3| M2  |2I4> should be  -34.72075317
          if (F8) MnMat( 20, 60,3)=    0.00000000 !  <4H3| M4  |2I4> should be  -41.00669620
          if (F8) MnMat( 28, 33,2)=    0.00000000 !  <2P1| M2  |2D2> should be  -18.85694608
          if (F8) MnMat( 28, 33,3)=    0.00000000 !  <2P1| M4  |2D2> should be  -28.86456366
          if (F8) MnMat( 40, 48,1)=    0.00000000 !  <2F4| M0  |2G5> should be  -11.03456954
          if (F8) MnMat( 40, 48,2)=    0.00000000 !  <2F4| M2  |2G5> should be  -19.81345691
          if (F8) MnMat( 40, 48,3)=    0.00000000 !  <2F4| M4  |2G5> should be  -30.89655174
          MnMat( 42, 45,1)=   -4.80846400 !  <2F6| M0  |2G2> should be   -6.40502035
          MnMat( 42, 45,2)=   13.63966400 !  <2F6| M2  |2G2> should be    6.47935328
          MnMat( 42, 45,3)=   15.11492200 !  <2F6| M4  |2G2> should be    7.17613059
          MnMat( 42, 47,1)=   10.34052400 !  <2F6| M0  |2G4> should be    8.01236869
          MnMat( 42, 47,2)=    8.95767700 !  <2F6| M2  |2G4> should be   -9.56379237
          MnMat( 42, 47,3)=  -17.61396700 !  <2F6| M4  |2G4> should be  -23.97342772
          MnMat( 43, 45,1)=    7.77385100 !  <2F7| M0  |2G2> should be    9.37040840
          MnMat( 43, 45,2)=    0.27276700 !  <2F7| M2  |2G2> should be    7.43307860
          MnMat( 43, 45,3)=    0.30996200 !  <2F7| M4  |2G2> should be    8.24876006
          MnMat( 43, 47,1)=  -14.18679000 !  <2F7| M0  |2G4> should be  -11.85863534
          MnMat( 43, 47,2)=   -6.00368500 !  <2F7| M2  |2G4> should be   12.51778247
          MnMat( 43, 47,3)=    2.29904500 !  <2F7| M4  |2G4> should be    8.65850543
          if (F8) MnMat( 44, 45,1)=    0.00000000 !  <2G1| M0  |2G2> should be   -6.66245956
          if (F8) MnMat( 44, 45,2)=    0.00000000 !  <2G1| M2  |2G2> should be  -32.29765748
          if (F8) MnMat( 44, 45,3)=    0.00000000 !  <2G1| M4  |2G2> should be  -40.24695231
          MnMat( 47, 47,1)=  -76.94950200 !  <2G4| M0  |2G4> should be  -72.76885407
          MnMat( 47, 47,2)=   49.04269100 !  <2G4| M2  |2G4> should be   41.62791550
          MnMat( 47, 47,3)=   55.31669400 !  <2G4| M4  |2G4> should be   45.58376248
          MnMat( 47, 48,1)=    4.25047300 !  <2G4| M0  |2G5> should be    0.06982120
          MnMat( 47, 48,2)=   -1.71319200 !  <2G4| M2  |2G5> should be    5.70159451
          MnMat( 47, 48,3)=   -0.42979600 !  <2G4| M4  |2G5> should be    9.30313940
          MnMat( 53, 58,1)=    7.40787400 !  <2H4| M0  |2I2> should be    7.41619849
          MnMat( 53, 58,2)=  -30.11753500 !  <2H4| M2  |2I2> should be  -30.08021003
          MnMat( 53, 58,3)=  -20.57354900 !  <2H4| M4  |2I2> should be  -20.53216462
          if (F8) PnMat(  2, 15,1)=    0.00000000 !  <6F | P2  |4G2> should be   -0.11142812
          if (F8) PnMat(  2, 15,2)=    0.00000000 !  <6F | P4  |4G2> should be   -0.05582917
          if (F8) PnMat(  2, 15,3)=    0.00000000 !  <6F | P6  |4G2> should be   -0.03099161
          if (F8) PnMat( 10, 44,1)=    0.00000000 !  <4F1| P2  |2G1> should be   -0.77199666
          if (F8) PnMat( 10, 44,2)=    0.00000000 !  <4F1| P4  |2G1> should be   -0.33794792
          if (F8) PnMat( 10, 44,3)=    0.00000000 !  <4F1| P6  |2G1> should be   -0.05367904
          if (F8) PnMat( 14, 51,1)=    0.00000000 !  <4G1| P2  |2H2> should be   -0.48466724
          if (F8) PnMat( 14, 51,2)=    0.00000000 !  <4G1| P4  |2H2> should be   -0.40777788
          if (F8) PnMat( 14, 51,3)=    0.00000000 !  <4G1| P6  |2H2> should be   -0.11010705
          if (F8) PnMat( 20, 60,1)=    0.00000000 !  <4H3| P2  |2I4> should be   -0.25726666
          if (F8) PnMat( 20, 60,2)=    0.00000000 !  <4H3| P4  |2I4> should be   -0.06150814
          if (F8) PnMat( 40, 48,1)=    0.00000000 !  <2F4| P2  |2G5> should be   -0.09012418
          if (F8) PnMat( 40, 48,2)=    0.00000000 !  <2F4| P4  |2G5> should be   -0.02471764
          if (F8) PnMat( 40, 48,3)=    0.00000000 !  <2F4| P6  |2G5> should be   -0.02005919
          PnMat( 42, 45,1)=    0.45462700 !  <2F6| P2  |2G2> should be    0.36592906
          PnMat( 42, 45,2)=    0.16892600 !  <2F6| P4  |2G2> should be    0.14693475
          PnMat( 42, 45,3)=   -0.06758900 !  <2F6| P6  |2G2> should be   -0.01824953
          PnMat( 42, 47,1)=   -0.25570700 !  <2F6| P2  |2G4> should be   -0.12839066
          PnMat( 42, 47,2)=    0.15697000 !  <2F6| P4  |2G4> should be    0.04583270
          PnMat( 42, 47,3)=   -0.06123200 !  <2F6| P6  |2G4> should be   -0.00539391
          PnMat( 43, 45,1)=   -0.27731300 !  <2F7| P2  |2G2> should be   -0.18861531
          PnMat( 43, 45,2)=   -0.12398500 !  <2F7| P4  |2G2> should be   -0.10199376
          PnMat( 43, 45,3)=   -0.02781700 !  <2F7| P6  |2G2> should be   -0.07715620
          PnMat( 43, 47,1)=   -0.04431900 !  <2F7| P2  |2G4> should be   -0.17163553
          PnMat( 43, 47,2)=   -0.00913800 !  <2F7| P4  |2G4> should be    0.10199902
          PnMat( 43, 47,3)=    0.03075300 !  <2F7| P6  |2G4> should be   -0.02508541
          if (F8) PnMat( 44, 45,1)=    0.00000000 !  <2G1| P2  |2G2> should be   -0.13856397
          if (F8) PnMat( 44, 45,2)=    0.00000000 !  <2G1| P4  |2G2> should be   -0.16667522
          if (F8) PnMat( 44, 45,3)=    0.00000000 !  <2G1| P6  |2G2> should be   -0.04670155
          PnMat( 47, 47,1)=    0.05372400 !  <2G4| P2  |2G4> should be    0.23736576
          PnMat( 47, 47,2)=    0.17099200 !  <2G4| P4  |2G4> should be    0.13919595
          PnMat( 47, 47,3)=    0.04159600 !  <2G4| P6  |2G4> should be    0.00807678
          PnMat( 47, 48,1)=    0.02495000 !  <2G4| P2  |2G5> should be   -0.15869238
          PnMat( 47, 48,2)=    0.01520400 !  <2G4| P4  |2G5> should be    0.04699983
          PnMat( 47, 48,3)=    0.01585000 !  <2G4| P6  |2G5> should be    0.04936926
          WRITE(IO,110) ' Mn/Pn'
        else
          WRITE(IO,100) ' Mn/Pn' 
        endif
        if (UseFncrossErrors(6)) then
          WRITE(IO,120) '  Tn  '
        endif
! ***Total disagreements for f^5 matrix elements (greater than 0.10E-02) is  80 for fncross & 41 fnmp
!
!--------------------------------------------------------------------------------
!   Errors in the f^6 configuration
      ELSE IF (nelectrons.eq.6) then
        if (UseFncrossErrors(1)) then
          WRITE(IO,120)'  Un  '
        endif
        if (UseFncrossErrors(2)) then
          WRITE(IO,120)' V11  '
        endif
        if (UseFncrossErrors(3)) then
          EEmat( 63, 63,3)=   -0.00010400 !  <3K4| F6  |3K4> should be   -0.00010565
          WRITE(IO,110)'  EE  '
        else
          WRITE(IO,100) '  EE  ' 
        endif
        if (UseFncrossErrors(4)) then
          WRITE(IO,120)' abg  '
        endif
        if (UseFncrossErrors(5)) then
        if (F8) then
          MnMat(  1, 11,1)=    0.00000000 !  <7F | M0  |5G3> should be  -15.39236824
          MnMat(  1, 11,2)=    0.00000000 !  <7F | M2  |5G3> should be  -51.49792337
          MnMat(  1, 11,3)=    0.00000000 !  <7F | M4  |5G3> should be  -52.47066044
          MnMat(  3, 26,1)=    0.00000000 !  <5P | M0  |3D3> should be   -0.09449112
          MnMat(  3, 26,2)=    0.00000000 !  <5P | M2  |3D3> should be   -3.18120098
          MnMat(  3, 26,3)=    0.00000000 !  <5P | M4  |3D3> should be   -3.02085242
          MnMat(  4, 32,1)=    0.00000000 !  <5D1| M0  |3F4> should be  -11.14970150
          MnMat(  4, 32,2)=    0.00000000 !  <5D1| M2  |3F4> should be  -17.87590838
          MnMat(  4, 32,3)=    0.00000000 !  <5D1| M4  |3F4> should be  -27.26894167
          MnMat( 11, 33,1)=    0.00000000 !  <5G3| M0  |3F5> should be   -0.68563058
          MnMat( 11, 33,2)=    0.00000000 !  <5G3| M2  |3F5> should be   -2.29389983
          MnMat( 11, 33,3)=    0.00000000 !  <5G3| M4  |3F5> should be   -2.33722898
          MnMat( 11, 53,1)=    0.00000000 !  <5G3| M0  |3H9> should be  -23.36699952
          MnMat( 11, 53,2)=    0.00000000 !  <5G3| M2  |3H9> should be  -59.86562945
          MnMat( 11, 53,3)=    0.00000000 !  <5G3| M4  |3H9> should be  -68.54152797
          MnMat( 15, 56,1)=    0.00000000 !  <5I2| M0  |3I3> should be  -10.84919474
          MnMat( 15, 56,2)=    0.00000000 !  <5I2| M2  |3I3> should be   -7.32165149
          MnMat( 15, 56,3)=    0.00000000 !  <5I2| M4  |3I3> should be  -13.38758929
          MnMat( 16, 68,1)=    0.00000000 !  <5K | M0  |3L3> should be  -39.35064182
          MnMat( 16, 68,2)=    0.00000000 !  <5K | M2  |3L3> should be  -10.06422476
          MnMat( 16, 68,3)=    0.00000000 !  <5K | M4  |3L3> should be  -17.33296124
          MnMat( 18, 79,1)=    0.00000000 !  <3P1| M0  |1D1> should be   -9.14696671
          MnMat( 18, 79,2)=    0.00000000 !  <3P1| M2  |1D1> should be   -9.87726345
          MnMat( 18, 79,3)=    0.00000000 !  <3P1| M4  |1D1> should be   -7.77600055
          MnMat( 20, 25,1)=    0.00000000 !  <3P3| M0  |3D2> should be  -10.69044968
          MnMat( 20, 25,2)=    0.00000000 !  <3P3| M2  |3D2> should be  -23.90125385
          MnMat( 20, 25,3)=    0.00000000 !  <3P3| M4  |3D2> should be  -26.65691605
          MnMat( 38, 47,1)=    0.00000000 !  <3G1| M0  |3H3> should be  -25.86336009
          MnMat( 38, 47,2)=    0.00000000 !  <3G1| M2  |3H3> should be  -11.57237072
          MnMat( 38, 47,3)=    0.00000000 !  <3G1| M4  |3H3> should be   -7.75979381
          MnMat( 39, 40,1)=    0.00000000 !  <3G2| M0  |3G3> should be   -4.00266115
          MnMat( 39, 40,2)=    0.00000000 !  <3G2| M2  |3G3> should be   -9.98790763
          MnMat( 39, 40,3)=    0.00000000 !  <3G2| M4  |3G3> should be  -12.76735480
          MnMat( 48, 98,1)=    0.00000000 !  <3H4| M0  |1H2> should be   -8.55983226
          MnMat( 48, 98,2)=    0.00000000 !  <3H4| M2  |1H2> should be   -9.98141795
          MnMat( 48, 98,3)=    0.00000000 !  <3H4| M4  |1H2> should be  -11.97096887
          MnMat( 51,107,1)=    0.00000000 !  <3H7| M0  |1I7> should be  -25.18187940
          MnMat( 51,107,2)=    0.00000000 !  <3H7| M2  |1I7> should be  -13.03208337
          MnMat( 51,107,3)=    0.00000000 !  <3H7| M4  |1I7> should be  -18.13067322
          MnMat( 54,102,1)=    0.00000000 !  <3I1| M0  |1I2> should be   -0.78788781
          MnMat( 54,102,2)=    0.00000000 !  <3I1| M2  |1I2> should be   -5.89093786
          MnMat( 54,102,3)=    0.00000000 !  <3I1| M4  |1I2> should be  -13.36534031
          MnMat( 55,103,1)=    0.00000000 !  <3I2| M0  |1I3> should be  -38.78377560
          MnMat( 55,103,2)=    0.00000000 !  <3I2| M2  |1I3> should be  -56.23469392
          MnMat( 55,103,3)=    0.00000000 !  <3I2| M4  |1I3> should be  -63.81942205
          MnMat( 59,110,1)=    0.00000000 !  <3I6| M0  |1K3> should be  -24.56633391
          MnMat( 59,110,2)=    0.00000000 !  <3I6| M2  |1K3> should be  -46.16664058
          MnMat( 59,110,3)=    0.00000000 !  <3I6| M4  |1K3> should be  -57.68860337
          PnMat(  1, 11,1)=    0.00000000 !  <7F | P2  |5G3> should be   -0.03167154
          PnMat(  1, 11,2)=    0.00000000 !  <7F | P4  |5G3> should be   -0.27352693
          PnMat(  1, 11,3)=    0.00000000 !  <7F | P6  |5G3> should be   -0.01490726
          PnMat(  3, 26,1)=    0.00000000 !  <5P | P2  |3D3> should be   -0.03534668
          PnMat(  3, 26,2)=    0.00000000 !  <5P | P4  |3D3> should be   -0.09146868
          PnMat(  3, 26,3)=    0.00000000 !  <5P | P6  |3D3> should be   -0.01070704
          PnMat(  4, 32,1)=    0.00000000 !  <5D1| P2  |3F4> should be   -0.51650754
          PnMat(  4, 32,2)=    0.00000000 !  <5D1| P4  |3F4> should be   -0.42962225
          PnMat(  4, 32,3)=    0.00000000 !  <5D1| P6  |3F4> should be   -0.17167041
          PnMat( 11, 33,1)=    0.00000000 !  <5G3| P2  |3F5> should be   -0.00141076
          PnMat( 11, 33,2)=    0.00000000 !  <5G3| P4  |3F5> should be   -0.01218386
          PnMat( 11, 53,1)=    0.00000000 !  <5G3| P2  |3H9> should be   -0.07220006
          PnMat( 11, 53,2)=    0.00000000 !  <5G3| P4  |3H9> should be   -0.11430063
          PnMat( 11, 53,3)=    0.00000000 !  <5G3| P6  |3H9> should be   -0.01606077
          PnMat( 15, 56,1)=    0.00000000 !  <5I2| P2  |3I3> should be   -0.19222297
          PnMat( 15, 56,2)=    0.00000000 !  <5I2| P4  |3I3> should be   -0.02692576
          PnMat( 15, 56,3)=    0.00000000 !  <5I2| P6  |3I3> should be   -0.02188713
          PnMat( 16, 68,1)=    0.00000000 !  <5K | P2  |3L3> should be   -0.68207779
          PnMat( 16, 68,2)=    0.00000000 !  <5K | P4  |3L3> should be   -0.01565835
          PnMat( 16, 68,3)=    0.00000000 !  <5K | P6  |3L3> should be   -0.00374175
          PnMat( 18, 79,1)=    0.00000000 !  <3P1| P2  |1D1> should be   -0.38036289
          PnMat( 18, 79,2)=    0.00000000 !  <3P1| P4  |1D1> should be   -0.39972682
          PnMat( 18, 79,3)=    0.00000000 !  <3P1| P6  |1D1> should be   -0.09309581
          PnMat( 20, 25,1)=    0.00000000 !  <3P3| P2  |3D2> should be   -0.19242809
          PnMat( 20, 25,2)=    0.00000000 !  <3P3| P4  |3D2> should be   -0.32567095
          PnMat( 20, 25,3)=    0.00000000 !  <3P3| P6  |3D2> should be   -0.09911150
          PnMat( 38, 47,1)=    0.00000000 !  <3G1| P2  |3H3> should be   -0.20868209
          PnMat( 38, 47,2)=    0.00000000 !  <3G1| P4  |3H3> should be   -0.00812898
          PnMat( 38, 47,3)=    0.00000000 !  <3G1| P6  |3H3> should be   -0.04262666
          PnMat( 39, 40,1)=    0.00000000 !  <3G2| P2  |3G3> should be   -0.06549809
          PnMat( 39, 40,2)=    0.00000000 !  <3G2| P4  |3G3> should be   -0.22136764
          PnMat( 39, 40,3)=    0.00000000 !  <3G2| P6  |3G3> should be   -0.08978772
          PnMat( 48, 98,1)=    0.00000000 !  <3H4| P2  |1H2> should be   -0.15192776
          PnMat( 48, 98,2)=    0.00000000 !  <3H4| P4  |1H2> should be   -0.30598902
          PnMat( 48, 98,3)=    0.00000000 !  <3H4| P6  |1H2> should be   -0.14825210
          PnMat( 51,107,1)=    0.00000000 !  <3H7| P2  |1I7> should be   -0.17394255
          PnMat( 51,107,2)=    0.00000000 !  <3H7| P4  |1I7> should be   -0.00417587
          PnMat( 51,107,3)=    0.00000000 !  <3H7| P6  |1I7> should be   -0.04483573
          PnMat( 54,102,1)=    0.00000000 !  <3I1| P2  |1I2> should be   -0.15158777
          PnMat( 54,102,2)=    0.00000000 !  <3I1| P4  |1I2> should be   -0.03785030
          PnMat( 54,102,3)=    0.00000000 !  <3I1| P6  |1I2> should be   -0.05597185
          PnMat( 55,103,1)=    0.00000000 !  <3I2| P2  |1I3> should be   -0.83429371
          PnMat( 55,103,2)=    0.00000000 !  <3I2| P4  |1I3> should be   -0.22993291
          PnMat( 55,103,3)=    0.00000000 !  <3I2| P6  |1I3> should be   -0.06860399
          PnMat( 59,110,1)=    0.00000000 !  <3I6| P2  |1K3> should be   -0.32693773
          PnMat( 59,110,2)=    0.00000000 !  <3I6| P4  |1K3> should be   -0.07857496
          PnMat( 59,110,3)=    0.00000000 !  <3I6| P6  |1K3> should be   -0.00223096
          WRITE(IO,110) ' Mn/Pn'
        endif  
        else
          WRITE(IO,100) ' Mn/Pn' 
        endif
        if (UseFncrossErrors(6)) then
          TnMat( 74, 76,2)=   -5.73709700 !  <1S1| T3  |1S3> should be    5.73709732
          TnMat( 75, 76,1)=    0.44824900 !  <1S2| T2  |1S3> should be   -0.44824953*
          TnMat( 75, 76,4)=   -3.55841700 !  <1S2| T6  |1S3> should be    3.55841786*
          TnMat( 76, 77,3)=   -0.29277000 !  <1S3| T4  |1S4> should be    0.29277002*
          TnMat( 76, 77,5)=    2.53546300 !  <1S3| T7  |1S4> should be   -2.53546276*
          TnMat(119,119,3)=    0.00000000 !  <1Q | T4  |1Q > should be   -0.85689307
          WRITE(IO,110) '  Tn  '
        else
          WRITE(IO,100) '  Tn  ' 
        endif
! ***Total disagreements for f^6 matrix elements (greater than 0.10E-02) is 102 for fncross & 7 fnmp
!
!--------------------------------------------------------------------------------
!   Errors in the f^7 configuration
      ELSE IF (nelectrons.eq.7) then
        if (UseFncrossErrors(1)) then
          Umat(112,117,4) =   -0.04524600 !  <2M1| U4  |2N2> should be   -0.04519385
          WRITE(IO,110) '  Un  '
        else
          WRITE(IO,100) '  Un  ' 
        endif
        if (UseFncrossErrors(2)) then
          WRITE(IO,120) ' V11  '
        endif
        if (UseFncrossErrors(3)) then
          EEmat(105,105,3)=    0.00063800 !  <2K6| F6  |2K6> should be    0.00063392
          WRITE(IO,110) '  EE  '
        else
          WRITE(IO,100) '  EE  ' 
        endif
        if (UseFncrossErrors(4)) then
          WRITE(IO,120) ' abg  '
        endif
        if (UseFncrossErrors(5)) then
        if (F8) then
          MnMat(  3, 22,1)=    0.00000000 !  <6D | M0  |4F5> should be   -9.68848113
          MnMat(  3, 22,2)=    0.00000000 !  <6D | M2  |4F5> should be  -70.60847614
          MnMat(  3, 22,3)=    0.00000000 !  <6D | M4  |4F5> should be  -89.74522813
        endif  
          MnMat(  5, 34,1)=    2.59229600 !  <6G | M0  |4H5> should be    0.00000000
          MnMat(  5, 34,2)=   69.75631700 !  <6G | M2  |4H5> should be   93.00844833
          MnMat(  5, 34,3)=  101.87080400 !  <6G | M4  |4H5> should be   88.98074832
        if (F8) then
          MnMat( 10, 13,1)=    0.00000000 !  <4P1| M0  |4D2> should be   -1.79417045
          MnMat( 10, 13,2)=    0.00000000 !  <4P1| M2  |4D2> should be   -0.79740909
          MnMat( 10, 13,3)=    0.00000000 !  <4P1| M4  |4D2> should be   -2.80905475
          MnMat( 11, 53,1)=    0.00000000 !  <4P2| M0  |2P4> should be   -7.22361098
          MnMat( 11, 53,2)=    0.00000000 !  <4P2| M2  |2P4> should be  -11.60928281
          MnMat( 11, 53,3)=    0.00000000 !  <4P2| M4  |2P4> should be  -16.25215673
          MnMat( 18, 72,1)=    0.00000000 !  <4F1| M0  |2G1> should be   -7.23746864
          MnMat( 18, 72,2)=    0.00000000 !  <4F1| M2  |2G1> should be  -19.15370490
          MnMat( 18, 72,3)=    0.00000000 !  <4F1| M4  |2G1> should be  -30.21925613
          MnMat( 20, 77,1)=    0.00000000 !  <4F3| M0  |2G6> should be   -8.19662871
          MnMat( 20, 77,2)=    0.00000000 !  <4F3| M2  |2G6> should be  -18.91373093
          MnMat( 20, 77,3)=    0.00000000 !  <4F3| M4  |2G6> should be  -34.26943609
          MnMat( 23, 83,1)=    0.00000000 !  <4G1| M0  |2H2> should be  -12.59135506
          MnMat( 23, 83,2)=    0.00000000 !  <4G1| M2  |2H2> should be  -31.00598472
          MnMat( 23, 83,3)=    0.00000000 !  <4G1| M4  |2H2> should be  -39.38745206
          MnMat( 30, 36,1)=    0.00000000 !  <4H1| M0  |4I2> should be  -19.48662092
          MnMat( 30, 36,2)=    0.00000000 !  <4H1| M2  |4I2> should be   -3.81556214
          MnMat( 30, 36,3)=    0.00000000 !  <4H1| M4  |4I2> should be   -0.06194094
          MnMat( 30, 84,1)=    0.00000000 !  <4H1| M0  |2H3> should be  -30.40239391
          MnMat( 30, 84,2)=    0.00000000 !  <4H1| M2  |2H3> should be  -12.75281317
          MnMat( 30, 84,3)=    0.00000000 !  <4H1| M4  |2H3> should be   -2.47452492
          MnMat( 32, 94,1)=    0.00000000 !  <4H3| M0  |2I4> should be  -28.17767701
          MnMat( 32, 94,2)=    0.00000000 !  <4H3| M2  |2I4> should be  -27.84686376
          MnMat( 32, 94,3)=    0.00000000 !  <4H3| M4  |2I4> should be  -30.35307249
          MnMat( 40,101,1)=    0.00000000 !  <4K1| M0  |2K2> should be  -39.19186902
          MnMat( 40,101,2)=    0.00000000 !  <4K1| M2  |2K2> should be  -12.12548480
          MnMat( 40,101,3)=    0.00000000 !  <4K1| M4  |2K2> should be   -3.15132027
          MnMat( 40,102,1)=    0.00000000 !  <4K1| M0  |2K3> should be   -0.82134230
          MnMat( 40,102,2)=    0.00000000 !  <4K1| M2  |2K3> should be  -26.03903984
          MnMat( 40,102,3)=    0.00000000 !  <4K1| M4  |2K3> should be  -53.37611033
          MnMat( 46,112,1)=    0.00000000 !  <4M | M0  |2M1> should be   -4.50966041
          MnMat( 46,112,2)=    0.00000000 !  <4M | M2  |2M1> should be  -22.33125958
          MnMat( 46,112,3)=    0.00000000 !  <4M | M4  |2M1> should be  -35.19309238
          MnMat( 65, 76,1)=    0.00000000 !  <2F4| M0  |2G5> should be   -4.81616596
          MnMat( 65, 76,2)=    0.00000000 !  <2F4| M2  |2G5> should be  -13.71517220
          MnMat( 65, 76,3)=    0.00000000 !  <2F4| M4  |2G5> should be  -22.63530843
          MnMat( 72, 73,1)=    0.00000000 !  <2G1| M0  |2G2> should be   -2.96109314
          MnMat( 72, 73,2)=    0.00000000 !  <2G1| M2  |2G2> should be  -33.18633054
          MnMat( 72, 73,3)=    0.00000000 !  <2G1| M4  |2G2> should be  -41.99213991
        endif
          MnMat( 85, 86,1)=  -24.32701100 !  <2H4| M0  |2H5> should be  -14.71006710
          MnMat( 85, 86,2)= -103.51234400 !  <2H4| M2  |2H5> should be  -60.38179621
          MnMat( 85, 86,3)= -112.60659800 !  <2H4| M4  |2H5> should be  -64.78680523
          MnMat( 85, 87,1)=   -9.30105100 !  <2H4| M0  |2H6> should be   -9.60417519
          MnMat( 85, 87,2)=   -2.94706600 !  <2H4| M2  |2H6> should be   -4.30652235
          MnMat( 85, 87,3)=   37.35589600 !  <2H4| M4  |2H6> should be   35.84865452
          MnMat( 85, 92,1)=    7.41786100 !  <2H4| M0  |2I2> should be    7.40121627
          MnMat( 85, 92,2)=  -22.78819300 !  <2H4| M2  |2I2> should be  -22.86286645
          MnMat( 85, 92,3)=   -8.02838800 !  <2H4| M4  |2I2> should be   -8.11116201
        if (F8) then
          MnMat(103,109,1)=    0.00000000 !  <2K4| M0  |2L3> should be   -9.81963428
          MnMat(103,109,2)=    0.00000000 !  <2K4| M2  |2L3> should be   -0.11454168
          MnMat(103,109,3)=    0.00000000 !  <2K4| M4  |2L3> should be   -2.96101589
        endif  
          PnMat(  5, 34,1)=   -0.14401700 !  <6G | P2  |4H5> should be   -0.69127901
          PnMat(  5, 34,2)=   -0.03570700 !  <6G | P4  |4H5> should be    0.26184811
          PnMat(  5, 34,3)=    0.08011100 !  <6G | P6  |4H5> should be    0.00000000
        if (F8) then
          PnMat( 10, 13,1)=    0.00000000 !  <4P1| P2  |4D2> should be   -0.12523412
          PnMat( 10, 13,2)=    0.00000000 !  <4P1| P4  |4D2> should be   -0.29275509
          PnMat( 10, 13,3)=    0.00000000 !  <4P1| P6  |4D2> should be   -0.20330214
          PnMat( 11, 53,1)=    0.00000000 !  <4P2| P2  |2P4> should be   -0.28988870
          PnMat( 11, 53,2)=    0.00000000 !  <4P2| P4  |2P4> should be   -0.21772673
          PnMat( 11, 53,3)=    0.00000000 !  <4P2| P6  |2P4> should be   -0.02929321
          PnMat( 18, 72,1)=    0.00000000 !  <4F1| P2  |2G1> should be   -1.06685649
          PnMat( 18, 72,2)=    0.00000000 !  <4F1| P4  |2G1> should be   -0.45990159
          PnMat( 18, 72,3)=    0.00000000 !  <4F1| P6  |2G1> should be   -0.05069687
          PnMat( 20, 77,1)=    0.00000000 !  <4F3| P2  |2G6> should be   -0.39266181
          PnMat( 20, 77,2)=    0.00000000 !  <4F3| P4  |2G6> should be   -0.07187428
          PnMat( 20, 77,3)=    0.00000000 !  <4F3| P6  |2G6> should be   -0.05717707
          PnMat( 23, 83,1)=    0.00000000 !  <4G1| P2  |2H2> should be   -0.79278902
          PnMat( 23, 83,2)=    0.00000000 !  <4G1| P4  |2H2> should be   -0.53475681
          PnMat( 23, 83,3)=    0.00000000 !  <4G1| P6  |2H2> should be   -0.10547472
          PnMat( 30, 36,1)=    0.00000000 !  <4H1| P2  |4I2> should be   -0.64955403
          PnMat( 30, 36,2)=    0.00000000 !  <4H1| P4  |4I2> should be   -0.18995223
          PnMat( 30, 36,3)=    0.00000000 !  <4H1| P6  |4I2> should be   -0.08338204
          PnMat( 30, 84,1)=    0.00000000 !  <4H1| P2  |2H3> should be   -0.59576408
          PnMat( 30, 84,2)=    0.00000000 !  <4H1| P4  |2H3> should be   -0.26338533
          PnMat( 30, 84,3)=    0.00000000 !  <4H1| P6  |2H3> should be   -0.02847087
          PnMat( 32, 94,1)=    0.00000000 !  <4H3| P2  |2I4> should be   -0.75462503
          PnMat( 32, 94,2)=    0.00000000 !  <4H3| P4  |2I4> should be   -0.28624942
          PnMat( 32, 94,3)=    0.00000000 !  <4H3| P6  |2I4> should be   -0.05810351
          PnMat( 40,101,1)=    0.00000000 !  <4K1| P2  |2K2> should be   -0.62044827
          PnMat( 40,101,2)=    0.00000000 !  <4K1| P4  |2K2> should be   -0.18990221
          PnMat( 40,101,3)=    0.00000000 !  <4K1| P6  |2K2> should be   -0.01910174
          PnMat( 40,102,1)=    0.00000000 !  <4K1| P2  |2K3> should be   -0.43622402
          PnMat( 40,102,2)=    0.00000000 !  <4K1| P4  |2K3> should be   -0.03230563
          PnMat( 40,102,3)=    0.00000000 !  <4K1| P6  |2K3> should be   -0.06521299
          PnMat( 46,112,1)=    0.00000000 !  <4M | P2  |2M1> should be   -0.18569190
          PnMat( 46,112,2)=    0.00000000 !  <4M | P4  |2M1> should be   -0.00852580
          PnMat( 46,112,3)=    0.00000000 !  <4M | P6  |2M1> should be   -0.09648278
          PnMat( 65, 76,1)=    0.00000000 !  <2F4| P2  |2G5> should be   -0.31116792
          PnMat( 65, 76,2)=    0.00000000 !  <2F4| P4  |2G5> should be   -0.13135726
          PnMat( 65, 76,3)=    0.00000000 !  <2F4| P6  |2G5> should be   -0.06798686
          PnMat( 72, 73,1)=    0.00000000 !  <2G1| P2  |2G2> should be   -0.23030724
          PnMat( 72, 73,2)=    0.00000000 !  <2G1| P4  |2G2> should be   -0.20863816
          PnMat( 72, 73,3)=    0.00000000 !  <2G1| P6  |2G2> should be   -0.05901993
        endif  
          PnMat( 85, 86,1)=   -0.61068500 !  <2H4| P2  |2H5> should be   -0.07640936
          PnMat( 85, 86,2)=   -0.07754100 !  <2H4| P4  |2H5> should be    0.05492363
          PnMat( 85, 86,3)=    0.42180700 !  <2H4| P6  |2H5> should be    0.12460944
          PnMat( 85, 87,1)=    0.60450000 !  <2H4| P2  |2H6> should be    0.58765946
          PnMat( 85, 87,2)=   -0.18440200 !  <2H4| P4  |2H6> should be   -0.18857775
          PnMat( 85, 87,3)=    0.00210100 !  <2H4| P6  |2H6> should be    0.01146807
        if (F8) then
          PnMat(103,109,1)=    0.00000000 !  <2K4| P2  |2L3> should be   -0.37021335
          PnMat(103,109,2)=    0.00000000 !  <2K4| P4  |2L3> should be   -0.02577373
          PnMat(103,109,3)=    0.00000000 !  <2K4| P6  |2L3> should be   -0.04949658
        endif  
          WRITE(IO,110) ' Mn/Pn'
        else
          WRITE(IO,100) ' Mn/Pn' 
        endif
        if (UseFncrossErrors(6)) then
          TnMat( 62, 69,2)=   -3.23969500 !  <2F1| T3  |2F8> should be   -4.53557368*
          TnMat( 63, 64,1)=    0.19675600 !  <2F2| T2  |2F3> should be    0.23596995*
          TnMat( 63, 64,4)=   -1.14142100 !  <2F2| T6  |2F3> should be   -1.45271801*
          TnMat( 63, 69,1)=   -0.35970100 !  <2F2| T2  |2F8> should be   -0.41032590*
          TnMat( 64, 65,1)=   -0.81351800 !  <2F3| T2  |2F4> should be   -0.79303945*
          TnMat( 64, 67,1)=   -0.63509500 !  <2F3| T2  |2F6> should be   -0.69985421*
          TnMat( 64, 70,3)=    0.49089700 !  <2F3| T4  |2F9> should be    0.47809144*
          TnMat( 64, 70,5)=   -0.11090300 !  <2F3| T7  |2F9> should be    0.00000000*
          TnMat( 64, 71,3)=   -0.32015000 !  <2F3| T4  |2F0> should be   -0.35856858*
          TnMat( 64, 71,5)=    1.39245300 !  <2F3| T7  |2F0> should be    1.72516390*
          TnMat( 65, 69,1)=    0.27643800 !  <2F4| T2  |2F8> should be    0.25000000*
          TnMat( 65, 69,4)=   -0.20987600 !  <2F4| T6  |2F8> should be    0.00000000*
          TnMat( 66, 69,3)=   -0.35771600 !  <2F5| T4  |2F8> should be   -0.37796447*
          TnMat( 66, 69,5)=   -0.17535400 !  <2F5| T7  |2F8> should be    0.00000000*
          TnMat( 67, 69,1)=   -0.08360400 !  <2F6| T2  |2F8> should be    0.00000000*
          TnMat( 67, 69,4)=   -2.06185300 !  <2F6| T6  |2F8> should be   -2.72554058*
          TnMat( 68, 69,3)=    0.03507100 !  <2F7| T4  |2F8> should be    0.00000000*
          TnMat( 68, 69,5)=    1.58610000 !  <2F7| T7  |2F8> should be    1.88982237*
          WRITE(IO,110) '  Tn  '
        else
          WRITE(IO,100) '  Tn  ' 
        endif
! ***Total disagreements for f^7 matrix elements (greater than 0.10E-02) is 128 for fncross & 41 for fnmp
      endif
      DO i=2,TS_n(Nelectrons)
        DO j=1,i-1
          Sij=(TS_f_Bases(1,i,Nelectrons)-TS_f_Bases(1,j,Nelectrons))/2  !  (S - S')  always integer
          Lij=(TS_f_Bases(2,i,Nelectrons)-TS_f_Bases(2,j,Nelectrons))    !  (L - L')  always integer
          VMat(i,j)=(-1)**(Lij-Sij)*VMat(j,i)
          do k=1,3
            Umat(i,j,k)=(-D1)**Lij*Umat(j,i,k)
            Umat(i,j,k+3)=(-D1)**Lij*Umat(j,i,k+3)
            EEMat(i,j,k)=EEMat(j,i,k)
            MnMat(i,j,k)=(-D1)**(Lij-Sij)*MnMat(j,i,k)
            PnMat(i,j,k)=(-D1)**(Lij-Sij)*PnMat(j,i,k)
            TnMat(i,j,k)=TnMat(j,i,k)
            TnMat(i,j,k+3)=TnMat(j,i,k+3)
          enddo ! k
        enddo ! j
      Enddo ! i
 100  FORMAT('**** WARNING: Possible errors in matrix elements:',A6,' have NOT been introduced.')
 110  FORMAT('**** WARNING: Errors HAVE been introduced into',A6,' matrix elements.')
 120  FORMAT('**** WARNING: There are no errors in the ',A6,' matrix elements.')
      return
      END SUBROUTINE FnCrossErrors
!--------------------------------------------------------------------------------
!
END MODULE f_e_fncrosserrors
!  606 lines