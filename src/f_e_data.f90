MODULE f_e_data  

IMPLICIT NONE

PUBLIC 

CHARACTER(LEN=6)  :: version="1.53.5"
!                                    12345678901234567890 
CHARACTER(LEN=20) :: CompileDate   =" 15th Sept 2022     "
CHARACTER(LEN=100):: F_E_PATH
CHARACTER(LEN=80) :: INFILE,OUTFIL

INTEGER, PARAMETER :: max_LSstates=119, max_Jstates=327  ! 
INTEGER, PARAMETER :: max_N=3432   ! for f^7  ! f^5 2002  ! 
INTEGER, PARAMETER :: max_exp=1000  ! maximum number of experimental energies that can be read
INTEGER, PARAMETER :: max_MM=7100 !36000 ! maximum number of Lx,Ly,Lz,Sx,Sy,Sz lower triangle matrix elements.
                                   ! biggest will be 7061(35,912) for f7 Lx,Ly with(without) tLS blocking.
INTEGER, PARAMETER :: max_GInt=50  ! The maximum number of ground states in intensity calculations.
INTEGER, PARAMETER :: max_EInt=100 ! The maximum number of excited states in intensity calculations.

INTEGER, PARAMETER :: in=2, io=3, iccf=7, iMatrix=5, IG=9, idata=10,IP1=11,IP2=12,Ilic=13,ifit=14,itest=15
INTEGER idebug ! set to 4 if one of CHCK is True
! 
LOGICAL TestOutput
logical option(10)/10*.false./,    check1(10)/10*.false./,    check2(10)/10*.false./ 
LOGICAL printMat1(10)/10*.false./, printMat2(10)/10*.false./, printMat3(10)/10*.false./, printMat4(10)/10*.false./
LOGICAL fastMat(10)/10*.false./
INTEGER printRng1(4)/4*0/,         printRng2(4)/4*0/,         printRng3(4)/4*0/,         printRng4(4)/4*0/ 
LOGICAL OUTP(10)/10*.false./,      PLOTout1(10)/10*.false./,  PLOTout2(10)/10*.false./,  warnings(10)/10*.false./
real*8  Yplot2(10)

! default values, can be changed with the INTR command
REAL*8 CFbit/1.0D-03/    ! The limit for a CF element to be considered = 0. 
REAL*8 Edegen/0.01D+00/  ! The threshold for two energy levels to be considered degenerate. 
REAL*8 FImax/0.05d+00/   ! The fraction of FreeIon character to be printed in the WF decomposition. See OPTN(3,4) 
REAL*8 MJmax/0.05d+00/   ! The fraction of MJ character to be printed in the WF decomposition. See OPTN(6,7)
REAL*8 PFactor/1.0d+00/  ! The factor by which the output to "matrix.dat" is multiplied.
REAL*8 TestCrit/1.0d+00/ ! The criteria that will be used to flag when a test is not met.
INTEGER*4 Ndecpts/1/     ! The number of decimal points used to print the energies.
!
!   nelectrons: the electron configuration.
!   n_states: the number of SL basis function for this configuration.
!   Lvalue: =2 for d and =3 for f electrons.
!   n_full: is the size of the full basis.
!   n_matrix: is the actual size of the |SLJMJ> basis. If not truncated, n_matrix=n_full.
!   n_jmatrix: is the size of the |SLJ> basis. Used in EDIP calculations.
      integer nelectrons/0/,n_states,Lvalue,n_matrix,n_jmatrix,n_full
      logical MatComplex          ! If true, the matrix to be diagonalised is complex.
      logical CFcomplex           ! If true, the ligand field and therefore the matrix to be diagonalised is complex
      logical MFcomplex           ! If true, the magnetic field has a non-zero y-component and therefore the matrix is complex
      logical n_odd               ! If true, there is an odd number of electrons
      logical calcVecs/.false./   ! If true, caclulate eigenvectors in diagonalisation.
      LOGICAL Complementary/.false./ ! If true, complementary number of electrons. 5>n>10 for d and 7>n>14 for f.
      logical FullMJBasis/.true./    ! If true the full |SLJMJ> basis is used, otherwise a |SLJ> basis is used.
      logical OrbitalsMl/.false./    ! TODO calculate the population of the ml  orbitals 
      logical OrbitalsReal/.false./  ! TODO calculate the population of the real d/f orbitals 
      LOGICAL wordy/.false./      ! If true, print extra output.
      LOGICAL EaveOnly/.false./   ! True if only EAVE is being varied in a fit.
      character*1 e_lab ! "d" or "f"
      integer vectLow/0/,vectHigh/0/,vectN/0/   ! printing eignevector, using the VECT command
!
      INTEGER lancIt/0/    ! Number of Lanczos iterations. If non-zero, Lanczos routine used for diagonalisation.  
      real*8 RLB,RUB       ! The Energy range used to search of levels in Lanczos routine.
!
      real*8 Vmat(max_LSstates,max_LSstates),    Umat(max_LSstates,max_LSstates,6)
      REAL*8 EEmat(max_LSstates,max_LSstates,4) 
! EEmat stores the e-e repulsion RMEs of in terms of the f2,f4,f6 coefficients to F^2, F^4, F^6. 
! The fourth value is e3, the coefficient to E^3 and is used in testing the Tn MREs.
      real*8 ABGMat(max_LSstates,3)
      real*8 MnMat(max_LSstates,max_LSstates,3),PnMat(max_LSstates,max_LSstates,3)
      real*8 SnMat(max_LSstates,max_LSstates,3),TnMat(max_LSstates,max_LSstates,6)
      real*8 save_EAVE, save_engs(max_N),  MAT1(max_N,max_N)  ! Mat1 is workspace (fit/makeFitMat, magnetics/calcMu)
      real*8 OrbEngs(7)  ! d- or f- orbital energies
      
!  PREDiagonalisation calculation. 
!  nPRED=0 unless PRED command used.
! Highest J = 14 (f7, 2Q) -> 2J+1 = 29
      INTEGER nPRED/0/,nPRmat(29)/29*0/, nPRmult(29,max_N)  
      REAL*8 PReng(29,max_N) 
!REAL*8 LFMAT(max_N,max_N),PRMAT(max_N,max_N)   ! 6 x 3432 x 3432 x 8 = 565,374,000 ~0.56G  (for Nblocks=6)
!COMPLEX*16 CLFMAT(max_N,max_N)                 ! 6 x 3432 x 3432 x 8 = 565,374,000 ~0.56G x2
REAL*8, dimension(:,:,:), allocatable :: PRvec  !PRMAT holds PRED eigenveectors.
REAL*8, dimension(:,:,:), allocatable :: LFMAT  !LFMAT holds the ligand field only matrix evaluated in the original basis (for PRED option).
COMPLEX*16, dimension(:,:,:), allocatable :: CLFMAT !CLFMAT holds the ligand field if complex (for PRED option).
!
! CCF
      INTEGER CCFtype/0/,N_CCF/0/, CCFindex(60,3)
      REAL*8 CCF(60),CCFtest(max_LSstates,max_LSstates,3)
      real*8 CCFmat(max_LSstates,max_LSstates,12,6) ! reduced matrix elements  
!
! matrices for the magnetic moments
      real*8    MM(2,3,max_MM)          
      integer*4 ImuIndex(2,3,max_MM),JmuIndex(2,3,max_MM),nMu(2,3)
!
! Eigenvector matrices 
      REAL*8 AR(max_N,max_N),AI(max_N,max_N)
      REAL*8 VR(max_N,max_N),VI(max_N,max_N)  ! Eignevectors in DIAGCH & workspace in GROUP & workspace buildMatrix to insert mag fld MEs.
      real*8 MAT(max_N,max_N),engs(max_N)     ! Mat contains the real eigenvectors, engs the energies
      COMPLEX*16 CMAT(max_N,max_N)            ! Contains eigenvectors after diagonalisation; 
                                              ! Added to when symmetry blocking is done
                                       
! Experimental energies:  (Each experimental energy can be assign to up to 20 calculated energies)    
      integer*4 NAssign(max_exp,20),NAss(max_exp),Nexp/0/,NAssignG(max_exp),NexpG/0/
      real*8 EXP_E(max_exp),EXP_I(max_exp),EXP_G(3,max_exp),WGT_E(max_exp),WGT_I(max_exp),WGT_G(3,max_exp)
      CHARACTER*1 WGTchar(max_exp) ! indicates whether weighting is different from 1
      logical wgtNote/.false./
!
      integer NrangeL/0/,NrangeH/0/
      real*8 PrangeL/0.0d+00/,PrangeH/0.0d+00/,MAGF(3)/0.0d+00,0.0d+00,0.0d+00/
!
!  g-values
      integer*4, PARAMETER :: maxNgexs=200  ! The maximum number of g-values that can be calculated. 
      real*8 RK(3)/3*1.0d+00/     ! Orbital reduction parameters.
      real*8 GVAL(3,maxNgexs),g_prin(3,3,maxNgexs)    ! Calculated g-values and principal directions for up to full f7 pseudo- S=1/2 pairs of states.
      integer*4 gexs(maxNgexs)/maxNgexs*0/,ngexs 
      REAL*8, PARAMETER :: BOHR=0.4668604D+00  ! The Bohr magneton in cm-1/Telsa 
      REAL*8, PARAMETER :: GE=2.0023D+00       ! The free electron g-value  
      logical g_axes/.false./       ! If true, the axes of g-tensor calculated
!
!  transitions:
      real*8 RENGG(max_GInt),RENGE(max_EInt)
!
!  magnetic dipole transitions
      integer*4 IMD/0/,MDunits/1/,MD(2,2),MDG,MDE,IDEGG(max_GInt),IDEGE(max_EInt)  ! MDunits=1,2 BM, x10^-7D
      real*8 RMD(max_GInt,max_EInt,5),MDconst
      COMPLEX*16 CMD(max_GInt,max_EInt,3)
!  electric dipole transitions
      integer*4 EDG,EDE,EngInt/0/,assignMethod/0/ ! assignMethod=0,1,2: energy order, CQN, Spin
      real*8 RED(max_GInt,max_EInt,5),EDconst
      COMPLEX*16 CED(max_GInt,max_EInt,3)
!  spin-allowed transitions
      integer SpAllow(2)
      real*8 SpinAllowed(max_N)

!REAL*8 ULpq(3,5,max_N,max_N)  ! 3 x 5 x 3432 x 3432 x 8 = 1,413,434,880 ~1.4G
!REAL*8 EDtrans(3,max_N,max_N)      ! 3 x 3432 x 3432 x 8 = 1,413,434,880 ~1.4G x2
COMPLEX*16, dimension(:,:,:), allocatable :: EDtrans
REAL*8, dimension(:,:,:), allocatable :: Ukq
integer, dimension(:,:,:), allocatable :: Ikq,Jkq
integer Nkq(3,0:6)

!  Judd-Ofelt
     logical JuddOfelt/.false./
     character*6 ED_IntMech/"      "/
     integer IED/0/,EDunits/1/,ED(2,2), wfg(max_GInt),wfe(max_EInt)  !IED=0(none)
     real*8 REDJO(max_GInt,max_EInt,5),EngJO(max_Jstates,max_Jstates),IntJO(max_Jstates,max_Jstates)
     real*8 SOtest(max_Jstates,max_Jstates)
     real*8 Dielectric/1.0d+0/,JO_omega(3),JO_U(max_Jstates,max_Jstates,3) ! <||U2||>, <||U4||>, <||U6||>
!
     real*8 MD_J(max_Jstates,max_Jstates)  ! magnetic dipole transition moments in |SLJ> basis.
!      
!  Numbers
CHARACTER(LEN=1), PARAMETER :: digits(0:15)=                           &
    (/"0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"/)
!  Characters
CHARACTER(LEN=1), PARAMETER ::                                         &
    arrow="^", space=" ", dot=".", plus="+", minus="-", uci="I",       &
    lci="i", slash="/", star="*", shriek="!", bra="(", ket=")"

CHARACTER(LEN=1), PARAMETER :: orbital(0:12)=                          &
     (/'S','P','D','F','G','H','I','K','L','M','N','O','Q'/)   !  note no J, no 2nd P!
!       0   1   2   3   4   5   6   7   8   9  10  11  12
!
! Q.N. in BASIS are 2S, L, 2J, 2MJ, K, seniority where K is the appropriate position in the 
! L,S reduced matrix elements in the tabulation of Nielson & Koster.
      INTEGER FullBasis(6,Max_N),JBASIS(5,max_Jstates)
! Q.N. in TS_Bases are 2S, L, seniority and implicit position in the "Standard order" of Nielson & Koster
      integer TS_d_Bases(3,16,5),TS_f_Bases(3,119,7),TS_f_QNs(3,119,7)
      INTEGER TS_Bases(3,max_LSstates,7)
      integer, PARAMETER ::  TSL_d(5)=(/ 1, 5,  8, 16, 16 /)          !  total SL free ions states for d electrons
      integer, PARAMETER ::  TSL_f(7)=(/ 1, 7, 17, 47, 73,119,119/)   !  total SL free ions states for f electrons
      integer, PARAMETER :: TSLJ_f(7)=(/ 2,13, 41,107,198,295,327/)   !  total SLJ free ions states for f electrons
      integer, PARAMETER :: TMsS_d(5)=(/ 5,25,50,100,100/)            !  total SLMlMs=S states for d electrons
      integer, PARAMETER :: TMsS_f(7)=(/ 7,49,147,441,735,1225,1225/) !  total SLMlMs=S states for f electrons
      integer TS_n(7)     ! total LS free ions states (TSL_d or TSL_f copied into this)
      integer TS_full(7)  ! total basis size
      integer S_mult(4),Nspin  ! Nspin is the number of different spin multiplicities contained in S_mult
      integer, PARAMETER :: TJp1_d(5)=(/ 6,  9, 12, 13, 14 /)         !  largest 2J+1 d electrons
      integer, PARAMETER :: TJp1_f(7)=(/ 8, 13, 18, 21, 24, 25, 26/)  !  largest 2J+1 f electrons
      integer TJp1(7)   ! largest 2J+1 (TJp1_d or TJp1_f copied into this)
!                                              
      character*3 TS_d_labels(16,5),TS_f_labels(119,7),TS_labels(max_LSstates,7)
! The sum of "L" degeneracies for each S
! 2S+1=    1    2    3    4    5    6    7    8     Total L
!-------------------------------------------------------------------
!   d0     1                                         1
!   d1          5                                    5 
!   d2    15        10                              25  = 5x5
!   d3         40        10                         50  = 2x 25
!   d4    50        45         5                   100  = 2x2x 25
!   d5         75        24         1              100  = 2x2x 25
!   d6    50        45         5
!   d7         40        10 
!   d8    15        10      
!   d9          5        
!  d10     1
!---------------------------------------------------------------------
!   f0     1                                         1
!   f1          7                                    7  
!   f2    28        21                              49  = 7x7
!   f3        112        35                        147  = 3x 49
!   f4   196       210        35                   441  = 3x3 49
!   f5        490       224        21              735  = 3x5x 49
!   f6   490       588       140        7         1225  = 5x5x 49 
!   f7        784       392        48        1    1225  = 5x5x 49
!
! Term Symbols "Standard order" of Nielson & Koster "coefficients..", 
!      highest spin mutiplicity first, then low to high orbital mutiplicity
      DATA TS_d_labels/                                                                                   &
        "2D ",15*"XXX",                                                                                   &
        "3P ","3F ","1S ","1D ","1G ",11*"XXX",                                                           &
        "4P ","4F ","2P ","2D1","2D2","2F ","2G ","2H ",8*"XXX",                                          &
        "5D ","3P1","3P2","3D ","3F1","3F2","3G ","3H ","1S1","1S2","1D1","1D2","1F ","1G1","1G2","1I ",  &
        "6S ","4P ","4D ","4F ","4G ","2S ","2P ","2D1","2D2","2D3","2F1","2F2","2G1","2G2","2H ","2I "/
      DATA TS_f_labels/                                                                                   &
        "2F ",118*"XXX" ,                                                                                 &
        "3P ","3F ","3H ","1S ","1D ","1G ","1I ",112*"XXX",                                              &
        "4S ","4D ","4F ","4G ","4I ","2P ","2D1","2D2","2F1","2F2","2G1","2G2","2H1","2H2","2I ","2K ",  &
              "2L ",102*"XXX",                                                                            &
        "5S ","5D ","5F ","5G ","5I ","3P1","3P2","3P3","3D1","3D2","3F1","3F2","3F3","3F4","3G1","3G2",  &
              "3G3","3H1","3H2","3H3","3H4","3I1","3I2","3K1","3K2","3L ","3M ","1S1","1S2","1D1","1D2",  &
              "1D3","1D4","1F ","1G1","1G2","1G3","1G4","1H1","1H2","1I1","1I2","1I3","1K ","1L1","1L2",  &
              "1N ",72*"XXX",                                                                             &
        "6P ","6F ","6H ","4S ","4P1","4P2","4D1","4D2","4D3","4F1","4F2","4F3","4F4","4G1","4G2","4G3",  &
              "4G4","4H1","4H2","4H3","4I1","4I2","4I3","4K1","4K2","4L ","4M ","2P1","2P2","2P3","2P4",  &
              "2D1","2D2","2D3","2D4","2D5","2F1","2F2","2F3","2F4","2F5","2F6","2F7","2G1","2G2","2G3",  &
              "2G4","2G5","2G6","2H1","2H2","2H3","2H4","2H5","2H6","2H7","2I1","2I2","2I3","2I4","2I5",  &
              "2K1","2K2","2K3","2K4","2K5","2L1","2L2","2L3","2M1","2M2","2N ","2O ",46*"XXX",           &
        "7F ","5S ","5P ","5D1","5D2","5D3","5F1","5F2","5G1","5G2","5G3","5H1","5H2","5I1","5I2","5K ",  &
              "5L ","3P1","3P2","3P3","3P4","3P5","3P6","3D1","3D2","3D3","3D4","3D5","3F1","3F2","3F3",  &
              "3F4","3F5","3F6","3F7","3F8","3F9","3G1","3G2","3G3","3G4","3G5","3G6","3G7","3H1","3H2",  &
              "3H3","3H4","3H5","3H6","3H7","3H8","3H9","3I1","3I2","3I3","3I4","3I5","3I6","3K1","3K2",  &
              "3K3","3K4","3K5","3K6","3L1","3L2","3L3","3M1","3M2","3M3","3N ","3O ","1S1","1S2","1S3",  &
              "1S4","1P ","1D1","1D2","1D3","1D4","1D5","1D6","1F1","1F2","1F3","1F4","1G1","1G2","1G3",  &
              "1G4","1G5","1G6","1G7","1G8","1H1","1H2","1H3","1H4","1I1","1I2","1I3","1I4","1I5","1I6",  &
              "1I7","1K1","1K2","1K3","1L1","1L2","1L3","1L4","1M1","1M2","1N1","1N2","1Q ",              &
        "8S ","6P ","6D ","6F ","6G ","6H ","6I ","4S1","4S2","4P1","4P2","4D1","4D2","4D3","4D4","4D5",  &
              "4D6","4F1","4F2","4F3","4F4","4F5","4G1","4G2","4G3","4G4","4G5","4G6","4G7","4H1","4H2",  &
              "4H3","4H4","4H5","4I1","4I2","4I3","4I4","4I5","4K1","4K2","4K3","4L1","4L2","4L3","4M ",  &
              "4N ","2S1","2S2","2P1","2P2","2P3","2P4","2P5","2D1","2D2","2D3","2D4","2D5","2D6","2D7",  &
              "2F1","2F2","2F3","2F4","2F5","2F6","2F7","2F8","2F9","2F0","2G1","2G2","2G3","2G4","2G5",  &
              "2G6","2G7","2G8","2G9","2G0","2H1","2H2","2H3","2H4","2H5","2H6","2H7","2H8","2H9","2I1",  &
              "2I2","2I3","2I4","2I5","2I6","2I7","2I8","2I9","2K1","2K2","2K3","2K4","2K5","2K6","2K7",  &
              "2L1","2L2","2L3","2L4","2L5","2M1","2M2","2M3","2M4","2N1","2N2","2O ","2Q "/
!           0   1   2   3   4   5   6   7   8   9  10  11  12
!          'S','P','D','F','G','H','I','K','L','M','N','O','Q'   !  note no J, no 2nd P

!       Three quantum numbers + implicit position in Neilson&Koster standard order:
!            (2S,L,seniority), 1:16 terms,  1:5 d^n config
!            (2S,L,seniority), 1:119 terms, 1:7 f^n config  
      DATA TS_d_Bases/                                                                   &
         1, 2, 1, 45*-1,                                                                 &
     
         2, 1, 2,  2, 3, 2,  0, 0, 0,  0, 2, 2,  0, 4, 2, 33*-1,                         &
     
         3, 1, 3,  3, 3, 3,  1, 1, 3,  1, 2, 1,  1, 2, 3,  1, 3, 3,  1, 4, 3,  1, 5, 3,24*-1,  &
      
         4, 2, 4,  2, 1, 2,  2, 1, 4,  2, 2, 4,  2, 3, 2,  2, 3, 4,  2, 4, 4,  2, 5, 4,  &
         0, 0, 0,  0, 0, 4,  0, 2, 2,  0, 2, 4,  0, 3, 4,  0, 4, 2,  0, 4, 4,  0, 6, 4,  &
     
         5, 0, 5,  3, 1, 3,  3, 2, 5,  3, 3, 3,  3, 4, 5,  1, 0, 5,  1, 1, 3,  1, 2, 1,  &
         1, 2, 3,  1, 2, 5,  1, 3, 3,  1, 3, 5,  1, 4, 3,  1, 4, 5,  1, 5, 3,  1, 6, 5/
      DATA TS_f_Bases/                                                                   &
         1, 3, 1, 354*-1,                                                                &
     
         2, 1, 2,  2, 3, 2,  2, 5, 2,  0, 0, 0,  0, 2, 2,  0, 4, 2,  0, 6, 2,  336*-1,   &
     
         3, 0, 3,  3, 2, 3,  3, 3, 3,  3, 4, 3,  3, 6, 3,  1, 1, 3,  1, 2, 3,  1, 2, 3,  &
         1, 3, 1,  1, 3, 3,  1, 4, 3,  1, 4, 3,  1, 5, 3,  1, 5, 3,  1, 6, 3,  1, 7, 3,  &
         1, 8, 3, 306*-1,                                                                &

         4, 0, 4,  4, 2, 4,  4, 3, 4,  4, 4, 4,  4, 6, 4,  2, 1, 2,  2, 1, 4,  2, 1, 4,  &
         2, 2, 4,  2, 2, 4,  2, 3, 2,  2, 3, 4,  2, 3, 4,  2, 3, 4,  2, 4, 4,  2, 4, 4,  &
         2, 4, 4,  2, 5, 2,  2, 5, 4,  2, 5, 4,  2, 5, 4,  2, 6, 4,  2, 6, 4,  2, 7, 4,  &
         2, 7, 4,  2, 8, 4,  2, 9, 4,  0, 0, 0,  0, 0, 4,  0, 2, 2,  0, 2, 4,  0, 2, 4,  &
         0, 2, 4,  0, 3, 4,  0, 4, 2,  0, 4, 4,  0, 4, 4,  0, 4, 4,  0, 5, 4,  0, 5, 4,  &
         0, 6, 2,  0, 6, 4,  0, 6, 4,  0, 7, 4,  0, 8, 4,  0, 8, 4,  0,10, 4,   216*-1,  &
     
         5, 1, 5,  5, 3, 5,  5, 5, 5,  3, 0, 3,  3, 1, 5,  3, 1, 5,  3, 2, 3,  3, 2, 5,  &
         3, 2, 5,  3, 3, 3,  3, 3, 5,  3, 3, 5,  3, 3, 5,  3, 4, 3,  3, 4, 5,  3, 4, 5,  &
         3, 4, 5,  3, 5, 5,  3, 5, 5,  3, 5, 5,  3, 6, 3,  3, 6, 5,  3, 6, 5,  3, 7, 5,  &
         3, 7, 5,  3, 8, 5,  3, 9, 5,  1, 1, 3,  1, 1, 5,  1, 1, 5,  1, 1, 5,  1, 2, 3,  &
         1, 2, 3,  1, 2, 5,  1, 2, 5,  1, 2, 5,  1, 3, 1,  1, 3, 3,  1, 3, 5,  1, 3, 5,  &
         1, 3, 5,  1, 3, 5,  1, 3, 5,  1, 4, 3,  1, 4, 3,  1, 4, 5,  1, 4, 5,  1, 4, 5,  &
         1, 4, 5,  1, 5, 3,  1, 5, 3,  1, 5, 5,  1, 5, 5,  1, 5, 5,  1, 5, 5,  1, 5, 5,  &
         1, 6, 3,  1, 6, 5,  1, 6, 5,  1, 6, 5,  1, 6, 5,  1, 7, 3,  1, 7, 5,  1, 7, 5,  &
         1, 7, 5,  1, 7, 5,  1, 8, 3,  1, 8, 5,  1, 8, 5,  1, 9, 5,  1, 9, 5,  1,10, 5,  &
         1,11, 5,  138*-1,                                                               &
         
         6, 3, 6,  4, 0, 4,  4, 1, 6,  4, 2, 4,  4, 2, 6,  4, 2, 6,  4, 3, 4,  4, 3, 6,  &
         4, 4, 4,  4, 4, 6,  4, 4, 6,  4, 5, 6,  4, 5, 6,  4, 6, 4,  4, 6, 6,  4, 7, 6,  &
         4, 8, 6,  2, 1, 2,  2, 1, 4,  2, 1, 4,  2, 1, 6,  2, 1, 6,  2, 1, 6,  2, 2, 4,  &  
         2, 2, 4,  2, 2, 6,  2, 2, 6,  2, 2, 6,  2, 3, 2,  2, 3, 4,  2, 3, 4,  2, 3, 4,  &
         2, 3, 6,  2, 3, 6,  2, 3, 6,  2, 3, 6,  2, 3, 6,  2, 4, 4,  2, 4, 4,  2, 4, 4,  &
         2, 4, 6,  2, 4, 6,  2, 4, 6,  2, 4, 6,  2, 5, 2,  2, 5, 4,  2, 5, 4,  2, 5, 4,  &
         2, 5, 6,  2, 5, 6,  2, 5, 6,  2, 5, 6,  2, 5, 6,  2, 6, 4,  2, 6, 4,  2, 6, 6,  &
         2, 6, 6,  2, 6, 6,  2, 6, 6,  2, 7, 4,  2, 7, 4,  2, 7, 6,  2, 7, 6,  2, 7, 6,  &
         2, 7, 6,  2, 8, 4,  2, 8, 6,  2, 8, 6,  2, 9, 4,  2, 9, 6,  2, 9, 6,  2,10, 6,  &
         2,11, 6,  0, 0, 0,  0, 0, 4,  0, 0, 6,  0, 0, 6,  0, 1, 6,  0, 2, 2,  0, 2, 4,  &
         0, 2, 4,  0, 2, 4,  0, 2, 6,  0, 2, 6,  0, 3, 4,  0, 3, 6,  0, 3, 6,  0, 3, 6,  &
         0, 4, 2,  0, 4, 4,  0, 4, 4,  0, 4, 4,  0, 4, 6,  0, 4, 6,  0, 4, 6,  0, 4, 6,  &
         0, 5, 4,  0, 5, 4,  0, 5, 6,  0, 5, 6,  0, 6, 2,  0, 6, 4,  0, 6, 4,  0, 6, 6,  &
         0, 6, 6,  0, 6, 6,  0, 6, 6,  0, 7, 4,  0, 7, 6,  0, 7, 6,  0, 8, 4,  0, 8, 4,  &
         0, 8, 6,  0, 8, 6,  0, 9, 6,  0, 9, 6,  0,10, 4,  0,10, 6,  0,12, 6,            &
         
         7, 0, 7,  5, 1, 5,  5, 2, 7,  5, 3, 5,  5, 4, 7,  5, 5, 5,  5, 6, 7,  3, 0, 3,  &
         3, 0, 7,  3, 1, 5,  3, 1, 5,  3, 2, 3,  3, 2, 5,  3, 2, 5,  3, 2, 7,  3, 2, 7,  &
         3, 2, 7,  3, 3, 3,  3, 3, 5,  3, 3, 5,  3, 3, 5,  3, 3, 7,  3, 4, 3,  3, 4, 5,  &
         3, 4, 5,  3, 4, 5,  3, 4, 7,  3, 4, 7,  3, 4, 5,  3, 5, 5,  3, 5, 5,  3, 5, 7,  &
         3, 5, 7,  3, 5, 7,  3, 6, 3,  3, 6, 5,  3, 6, 5,  3, 6, 7,  3, 6, 7,  3, 7, 5,  &
         3, 7, 5,  3, 7, 7,  3, 8, 5,  3, 8, 7,  3, 8, 7,  3, 9, 5,  3,10, 7,  1, 0, 7,  &
         1, 0, 7,  1, 1, 3,  1, 1, 5,  1, 1, 5,  1, 1, 5,  1, 1, 7,  1, 2, 3,  1, 2, 3,  &
         1, 2, 5,  1, 2, 5,  1, 2, 5,  1, 2, 7,  1, 2, 7,  1, 3, 1,  1, 3, 3,  1, 3, 5,  &
         1, 3, 5,  1, 3, 5,  1, 3, 5,  1, 3, 5,  1, 3, 7,  1, 3, 7,  1, 3, 7,  1, 4, 3,  &
         1, 4, 3,  1, 4, 5,  1, 4, 5,  1, 4, 5,  1, 4, 5,  1, 4, 7,  1, 4, 7,  1, 4, 7,  &
         1, 4, 7,  1, 5, 3,  1, 5, 3,  1, 5, 5,  1, 5, 5,  1, 5, 5,  1, 5, 5,  1, 5, 5,  &
         1, 5, 7,  1, 5, 7,  1, 6, 3,  1, 6, 5,  1, 6, 5,  1, 6, 5,  1, 6, 5,  1, 6, 7,  &
         1, 6, 7,  1, 6, 7,  1, 6, 7,  1, 7, 3,  1, 7, 5,  1, 7, 5,  1, 7, 5,  1, 7, 5,  &
         1, 7, 7,  1, 7, 7,  1, 8, 3,  1, 8, 5,  1, 8, 5,  1, 8, 7,  1, 8, 7,  1, 9, 5,  &
         1, 9, 5,  1, 9, 7,  1, 9, 7,  1,10, 5,  1,10, 7,  1,11, 5,  1,12, 7/

!  Three more quantum numbers (in Neilson&Koster standard order):
!    ((w1,w2,w3),(u1,u2),AB),  1:119 terms, 1:7 f^n configurations  
!    (w1,w2,w3) and (u1,u2) are stored as 3 and 2 digit integers  w1,w2,w3 and u1,u2.
!    AB is zero unless the extra label A(1) or B(2) are required

      DATA TS_f_QNs/                                                                      &
         100,10,1, 354*-1,                                                                &
     
         110,11,0, 110,10,0, 110,11,0, 000,00,0, 200,20,0, 200,20,0, 200,20,0,  336*-1,   &
     
         111,00,0, 111,20,0, 111,10,0, 111,20,0, 111,20,0, 210,11,0, 210,20,0, 210,21,0,  &
         100,10,0, 210,21,0, 210,20,0, 210,21,0, 210,11,0, 210,21,0, 210,20,0, 210,21,0,  &
         210,21,0,   306*-1,                                                              &

         111,00,0, 111,20,0, 111,10,0, 111,20,0, 111,20,0, 110,11,0, 211,11,0, 211,30,0,  &
         211,20,0, 211,21,0, 110,10,0, 211,10,0, 211,21,0, 211,30,0, 211,20,0, 211,21,0,  &
         211,30,0, 110,11,0, 211,11,0, 211,21,0, 211,30,0, 211,20,0, 211,30,0, 211,21,0,  &
         211,30,0, 211,21,0, 211,30,0, 000,00,0, 220,22,0, 200,20,0, 220,20,0, 220,21,0,  &
         220,22,0, 220,21,0, 200,20,0, 220,20,0, 220,21,0, 220,22,0, 220,21,0, 220,22,0,  &
         200,20,0, 220,20,0, 220,22,0, 220,21,0, 220,21,0, 220,22,0, 220,22,0,   216*-1,  &

         110,11,0, 110,10,0, 110,11,0, 111,00,0, 211,11,0, 211,30,0, 111,20,0, 211,20,0,  &
         211,21,0, 111,10,0, 211,10,0, 211,21,0, 211,30,0, 111,20,0, 211,20,0, 211,21,0,  &
         211,30,0, 211,11,0, 211,21,0, 211,30,0, 111,20,0, 211,20,0, 211,30,0, 211,21,0,  &
         211,30,0, 211,21,0, 211,30,0, 210,11,0, 221,11,0, 221,30,0, 221,31,0, 210,20,0,  &
         210,21,0, 221,20,0, 221,21,0, 221,31,0, 100,10,0, 210,21,0, 221,10,0, 221,21,0,  &
         221,30,0, 221,31,1, 221,31,2, 210,20,0, 210,21,0, 221,20,0, 221,21,0, 221,30,0,  &
         221,31,0, 210,11,0, 210,21,0, 221,11,0, 221,21,0, 221,30,0, 221,31,1, 221,31,2,  &          
         210,20,0, 221,20,0, 221,30,0, 221,31,1, 221,31,2, 210,21,0, 221,21,0, 221,30,0,  &
         221,31,1, 221,31,2, 210,21,0, 221,21,0, 221,31,0, 221,30,0, 221,31,0, 221,31,0,  &
         221,31,0,  138*-1,                                                               &

         100,10,0, 111,00,0, 210,11,0, 111,20,0, 210,20,0, 210,21,0, 111,10,0, 210,21,0,  &
         111,20,0, 210,20,0, 210,21,0, 210,11,0, 210,21,0, 111,20,0, 210,20,0, 210,21,0,  &
         210,21,0, 110,11,0, 211,11,0, 211,30,0, 221,11,0, 221,30,0, 221,31,0, 211,20,0,  &
         211,21,0, 221,20,0, 221,21,0, 221,31,0, 110,10,0, 211,10,0, 211,21,0, 211,30,0,  &
         221,10,0, 221,21,0, 221,30,0, 221,31,1, 221,31,2, 211,20,0, 211,21,0, 211,30,0,  &
         221,20,0, 221,21,0, 221,30,0, 221,31,0, 110,11,0, 211,11,0, 211,21,0, 211,30,0,  &
         221,11,0, 221,21,0, 221,30,0, 221,31,1, 221,31,2, 211,20,0, 211,30,0, 221,20,0,  &
         221,30,0, 221,31,1, 221,31,2, 211,21,0, 211,30,0, 221,21,0, 221,30,0, 221,31,1,  &
         221,31,2, 211,21,0, 221,21,0, 221,31,0, 211,30,0, 221,30,0, 221,31,0, 221,31,0,  &
         221,31,0, 000,00,0, 220,22,0, 222,00,0, 222,40,0, 222,30,0, 200,20,0, 220,20,0,  &
         220,21,0, 220,22,0, 222,20,0, 222,40,0, 220,21,0, 222,10,0, 222,30,0, 222,40,0,  &
         200,20,0, 220,20,0, 220,21,0, 220,22,0, 222,20,0, 222,30,0, 222,40,1, 222,40,2,  &
         220,21,0, 220,22,0, 222,30,0, 222,40,0, 200,20,0, 220,20,0, 220,22,0, 222,20,0,  &
         222,30,0, 222,40,1, 222,40,2, 220,21,0, 222,30,0, 222,40,0, 220,21,0, 220,22,0,  &
         222,40,1, 222,40,2, 222,30,0, 222,40,0, 220,22,0, 222,40,0, 222,40,0,            &
                                                                                     
         000,00,0, 110,11,0, 200,20,0, 110,10,0, 200,20,0, 110,11,0, 200,20,0, 111,00,0,  &
         220,22,0, 211,11,0, 211,30,0, 111,20,0, 211,20,0, 211,21,0, 220,20,0, 220,21,0,  &
         220,22,0, 111,10,0, 211,10,0, 211,21,0, 211,30,0, 220,21,0, 111,20,0, 211,20,0,  &
         211,21,0, 211,30,0, 220,20,0, 220,21,0, 220,22,0, 211,11,0, 211,21,0, 211,30,0,  &
         220,21,0, 220,22,0, 111,20,0, 211,20,0, 211,30,0, 220,20,0, 220,22,0, 211,21,0,  &
         211,30,0, 220,21,0, 211,21,0, 220,21,0, 220,22,0, 211,30,0, 220,22,0, 222,00,0,  &
         222,40,0, 210,11,0, 221,11,0, 221,30,0, 221,31,0, 222,30,0, 210,20,0, 210,21,0,  &
         221,20,0, 221,21,0, 221,31,0, 222,20,0, 222,40,0, 100,10,0, 210,21,0, 221,10,0,  &
         221,21,0, 221,30,0, 221,31,1, 221,31,2, 222,10,0, 222,30,0, 222,40,0, 210,20,0,  &
         210,21,0, 221,20,0, 221,21,0, 221,30,0, 221,31,0, 222,20,0, 222,30,0, 222,40,1,  &
         222,40,2, 210,11,0, 210,21,0, 221,11,0, 221,21,0, 221,30,0, 221,31,1, 221,31,2,  &
         222,30,0, 222,40,0, 210,20,0, 221,20,0, 221,30,0, 221,31,1, 221,31,2, 222,20,0,  &
         222,30,0, 222,40,1, 222,40,2, 210,21,0, 221,21,0, 221,30,0, 221,31,1, 221,31,2,  &
         222,30,0, 222,40,0, 210,21,0, 221,21,0, 221,31,0, 222,40,1, 222,40,2, 221,30,0,  &
         221,31,0, 222,30,0, 222,40,0, 221,31,0, 222,40,0, 221,31,0, 222,40,0/            
         
!  WF properties, SPIN maximum number of different spin multiplicities for f electrons = 4.
REAL*8 SPIN(max_N,4),binMJ(1:max_N,-30:30)
character cFreeIon(max_N)*80,cAllMJ(max_N)*80
!  Symmetry blocking
logical block/.false./,printCQN/.false./,C2yExists/.false./,TR_Exists/.false./
integer nBlocks/1/,Nblock(26)/26*0/,Bbasis(26,7,Max_N),iCQN(Max_N)
character*5 CQN(26)/26*"     "/, allCQN(Max_N)
!  Symmetry Group Data
integer IDG/0/,ISO/1/,IGP/0/,NCLASS,NROT,NSINGLE
integer, parameter :: NGRP=24
INTEGER ISYMEL(NGRP,4)
COMPLEX*16 GROUPT(NGRP,NGRP),CSYM(max_N,NGRP)
REAL*8 RSYMEL(NGRP,3)
CHARACTER*4 GROUP,GROUPO(NGRP),GRPR1(NGRP),GRPR2(NGRP)
REAL*8 REPS(max_N,NGRP)    
CHARACTER*4 IRREP(max_N,3)
CHARACTER(LEN=4), PARAMETER :: Groups(37)=                           &
    (/"C1  ","Ci  ","C2  ","Cs  ","C2h ","D2  ","C2v ","D2h ","C4  ","S4  ","C4h ","D4  ","C4v ","D2d ","D4h ","C3  ", &
      "C3i ","D3  ","C3v ","D3d ","C6  ","C3h ","C6h ","D6  ","C6v ","D3h ","D6h ","T   ","Th  ","Td  ","O   ","Oh  ", &
      "TdT ","OT  ","OhT ","C2Y ","CsY "/)
! Intensity parameters (Real and Imaginary parts)
logical AltpReal,AltpImag ! true if real / imag Altp parameters non-zero
integer IntN1,IntN2  ! IntN1=1/2 for Altp/Blki; IntN2=1(pure real)/2(pure imag)/3(complex)
REAL*8 AltpR(3,3,8),AltpI(3,3,8),BlkiR(15,3),BlkiI(15,3)
! timing info?
logical iTime/.false./

END MODULE f_e_data
!  414