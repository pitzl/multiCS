
! electron multiple Coulomb scattering Monte Carlo
! in pure elements

! based on mcsda.f
! F. Salvat, J.D. Martinez, R. Mayol, J. Parellada
! COMP. PHYS. COMMUN. 42 (1986) 93
! AALC

! Fortran90 by Daniel Pitzl, DESY, Sep 2019

! ELASTIC SCATTERING EVENTS ARE DESCRIBED BY USING AN ANALYTICAL FIT
! TO THE RELATIVISTIC MOTT CROSS-SECTION.
! THE MEAN RATE OF ENERGY LOSS PER UNIT PATH LENGTH IS DERIVED FROM
! THE RELATIVISTIC BETHE STOPPING POWER, I.E. FROM THE CONTINUOUS SLOWING
! DOWN APPROXIMATION (CSDA).
! THE STRAGGLING EFFECT IS SIMULATED BY MEANS OF THE LENZ INELASTIC (TOTAL)
! CROSS-SECTION AND THE CLASSICAL 1/W**2 DISTRIBUTION OF ENERGY LOSSES.
! ANGULAR DEFLECTIONS OF THE ELECTRON PATH DUE TO INELASTIC COLLISIONS ARE
! COMPUTED FROM THE CLASSICAL THEORY OF BINARY COLLISIONS.

! REFERENCES:
! (1) F. SALVAT ET AL., J. PHYS. D: APPL. PHYS. 17 (1984) 185
! (2) F. SALVAT ET AL., J. PHYS. D: APPL. PHYS. 17 (1984) 1545
! (3) R. MAYOL ET AL., ANALES DE FISICA A80 (1984) 130
! (4) M.J. BERGER AND S.M. SELTZER, NBSIR 82-2550-A (1982)
!
! *****    ATOMIC DATA
!
! FREE ATOMS
! --------------------------------------------------------------
! ATOM  Z   AW         A         ALPHA1    ALPHA2     I
! --------------------------------------------------------------
! HE    2   4.0026     1.00000   2.20564   2.20564   41.8 (GAS)
! LI    3   6.941      0.60453   2.81743   0.66246   40.0
! BE    4   9.01218    0.32784   4.54284   0.98509   63.7
! B     5   10.81      0.23268   5.99011   1.21345   76.0
! ! 6   12.011     0.15364   8.04166   1.49140   78.0
! N     7   14.0067    0.09957   10.8127   1.76869   82.0 (GAS)
! O     8   15.9994    0.06245   14.8346   2.04052   95.0 (GAS)
! F     9   18.9984    0.03678   21.4093   2.30606   115  (GAS)
! NE   10   20.179     0.01801   35.0596   2.56643   137  (GAS)
! NA   11   22.98977   0.74440   4.12051   0.87184   149
! MG   12   24.305     0.64234   4.72661   1.00249   156
! AL   13   26.98154   0.60013   5.14064   1.01539   166
! SI   14   28.0855    0.51591   5.84957   1.17330   173
! P    15   30.97376   0.43866   6.67081   1.34106   173
! S    16   32.06      0.37130   7.60432   1.50676   180
! CL   17   35.453     0.31314   8.67077   1.66877   174  (GAS)
! AR   18   39.948     0.26299   9.90260   1.82787   188  (GAS)
! K    19   39.0983    0.62732   5.87615   0.98798   190
! CA   20   40.08      0.58800   6.32192   1.00940   191
! SC   21   44.9559    0.55431   6.63300   1.10236   216
! TI   22   47.88      0.53593   6.89381   1.18766   233
! V    23   50.9415    0.52216   7.11955   1.26815   245
! CR   24   51.996     0.50646   7.31623   1.42531   257
! MN   25   54.9380    0.50356   7.49501   1.41970   272
! FE   26   55.847     0.49731   7.65419   1.49227   286
! CO   27   58.9332    0.49262   7.79761   1.56304   297
! NI   28   58.69      0.48921   7.92752   1.63224   311
! CU   29   63.546     0.46600   8.20246   1.82639   322
! ZN   30   65.38      0.48531   8.15505   1.76715   330
! GA   31   69.72      0.57810   7.45294   1.47922   334
! GE   32   72.59      0.55193   7.81505   1.52669   350
! AS   33   74.9216    0.51813   8.27460   1.60648   347
! SE   34   78.96      0.48352   8.78654   1.69690   348
! BR   35   79.904     0.45033   9.33526   1.79016   357
! KR   36   83.80      0.41895   9.91986   1.88404   352  (GAS)
! RB   37   85.4678    0.66732   7.38694   1.13156   363
! SR   38   87.62      0.65417   7.60914   1.10652   366
! Y    39   88.9059    0.62196   7.98456   1.19369   379
! ZR   40   91.22      0.59406   8.34255   1.27785   393
! NB   41   92.9064    0.55591   8.81359   1.43082   417
! MO   42   95.94      0.53108   9.19190   1.51654   424
! TC   43   97.907     0.52988   9.33183   1.50495   428
! RU   44   101.07     0.48905   9.92161   1.67634   441
! RH   45   102.9055   0.47090   10.2776   1.75229   449
! PD   46   106.42     0.36860   12.2194   2.09718   470
! AG   47   107.868    0.43878   10.9799   1.18991   470
! CD   48   112.41     0.45464   10.8499   1.84361   469
! IN   49   114.82     0.52418   9.95823   1.62066   488
! SN   50   118.69     0.51613   10.1882   1.62396   488
! SB   51   121.75     0.50244   10.5023   1.64840   487
! TE   52   127.60     0.48399   10.9022   1.69203   485
! I    53   126.9045   0.46451   11.3435   1.74292   491
! XE   54   131.29     0.44503   11.8150   1.79736   482  (GAS)
! CS   55   132.9054   0.63099   9.27951   1.21804   488
! BA   56   137.33     0.63153   9.37847   1.16983   491
! LA   57   138.9055   0.60457   9.78322   1.23954   501
! CE   58   140.12     0.60176   9.89662   1.26964   523
! PR   59   140.9077   0.59961   10.0029   1.29867   535
! ND   60   144.24     0.61818   9.84443   1.29305   546
! PM   61   144.913    0.61655   9.94325   1.32107   560
! SM   62   150.36     0.61524   10.0389   1.34861   574
! EU   63   151.96     0.61408   10.0134   1.37597   580
! GD   64   157.25     0.59470   10.4721   1.43349   591
! TB   65   158.9254   0.61231   10.3185   1.42998   614
! DY   66   162.5      0.61163   10.4097   1.45675   628
! HO   67   164.9304   0.61105   10.5006   1.48342   650
! ER   68   167.26     0.61054   10.5913   1.51003   658
! TM   69   168.9342   0.61006   10.6828   1.53672   674
! YB   70   173.04     0.60962   10.7746   1.56339   684
! LU   71   174.967    0.59607   11.0553   1.60432   694
! HF   72   178.49     0.57730   11.4181   1.66925   705
! TA   73   180.9479   0.55928   11.7897   1.73462   718
! W    74   183.85     0.54238   12.1611   1.79841   727
! RE   75   186.207    0.52646   12.5334   1.85940   736
! OS   76   190.2      0.51112   12.9161   1.92009   746
! IR   77   192.22     0.49632   13.3110   1.98039   757
! PT   78   195.08     0.46317   14.0902   2.12491   790
! AU   79   196.9665   0.44855   14.5490   2.18774   790
! HG   80   200.59     0.45535   14.5388   2.15774   800
! TL   81   204.383    0.51102   13.5079   1.94984   810
! PB   82   207.2      0.51405   13.5898   1.92284   823
! BI   83   208.9804   0.52372   13.5429   1.87088   823
! PO   84   208.982    0.51902   13.7725   1.87223   830
! AT   85   209.987    0.50987   14.0959   1.89216   825
! RN   86   222.018    0.49903   14.4661   1.92060   794  (GAS)
! FR   87   223.020    0.63069   12.2398   1.42884   827
! RA   88   226.0254   0.63838   12.2521   1.36251   826
! AC   89   227.0278   0.50967   14.1001   1.89316   841
! TH   90   232.0381   0.60184   13.0792   1.46078   847
! PA   91   231.0359   0.60180   13.1874   1.48806   878
! U    92   238.0289   0.59396   13.4427   1.52547   890
!
!
! WIGNER-SEITZ ATOMS
! ----------------------------------------------
! ATOM   Z    RWS      A      ALPHA1    ALPHA2
! ----------------------------------------------
! LI    3   3.25   0.37731   3.19082   1.28284
! NA   11   3.94   0.49868   5.08096   1.51901
! AL   13   2.99   0.31543   7.28752   1.77451
! SI   14   3.18   0.31186   7.76730   1.71824
! P    15   3.57   0.31007   8.17973   1.67147
! K    19   4.87   0.42851   7.49488   1.42250
! CA   20   4.11   0.40255   8.05513   1.42620
! V    23   2.81   0.22974   12.2348   1.96899
! CR   24   2.68   0.19633   13.6659   2.14626
! MN   25   2.70   0.21946   12.9263   2.10354
! FE   26   2.67   0.21948   13.0520   2.16581
! NI   28   2.60   0.21892   13.2947   2.30276
! CU   29   2.65   0.21374   13.4817   2.42997
! GE   32   3.34   0.35634   10.2486   2.09847
! RB   37   5.21   0.53376   8.61480   1.52016
! SR   38   4.48   0.52902   8.81678   1.48175
! NB   41   3.07   0.37795   11.4676   1.93243
! MO   42   2.93   0.35710   12.0620   1.99956
! RH   45   2.81   0.32521   13.2919   2.13494
! PD   46   2.87   0.31116   13.8007   2.21950
! AG   47   3.02   0.32982   13.4086   2.17450
! CS   55   5.64   0.52081   10.6875   1.54308
! BA   56   4.72   0.51198   10.9485   1.53695
! CE   58   3.87   0.45633   12.1097   1.72503
! EU   63   4.26   0.47914   12.1247   1.81620
! YB   70   4.05   0.47457   12.8767   2.03486
! TA   73   3.08   0.41974   14.4589   2.22411
! W    74   2.94   0.40817   14.8914   2.26673
! IR   77   2.84   0.39894   15.5877   2.30243
! PT   78   2.90   0.39375   15.8661   2.33737
! AU   79   2.98   0.39511   15.9768   2.34208
! PB   82   3.65   0.42298   15.6563   2.24990
! TH   90   3.76   0.49879   15.0981   1.84375

! Z ................... ATOMIC NUMBER
! AW .................. ATOMIC WEIGHT
! A, ALPHA1,ALPHA2 .... DENSITY PARAMETERS DERIVED FROM SELF-CONSISTENT
! DHFS CALCULATIONS (REFERENCE 3)
! I ................... MEAN EXCITATION POTENTIAL FOR SINGLE-ELEMENT
! MATERIALS (EV) RECOMMENDED IN REFERENCE 4
! RWS ................. RADIUS OF THE WIGNER-SEITZ SPHERE (A.U.)

! *****    NUMERICAL AND PHYSICAL CONSTANTS IN THE SOURCE LISTING
! PI              = 3.1415926535897931D+00
! 36/PI           = 1.1459155902616464D+01
! 2*PI            = 6.2831853071795862D+00
! 2**31-1         = 2.1474836470000000D+09
! 1/(2**31-1)     = 4.6566128752457969D-10
! 2*PI/(2**31-1)  = 2.9258361598967678D-09
! (2*E**3)**1/2   = 6.3380654656113586D+00   E=DEXP(1.)
! 2*(2/E)**3/4    = 1.5888453635426958D+00
! 2/3*(E/2)**3/4  = 8.3918382740555711D-01
! RADIAN (RAD)    = 5.7295779513082323D+01 DEG
! RAD/2           = 2.8647889756541161D+01
! LOG(2)          = 6.9314718055994533D-01

! atomic units: hbar = m = e = 1

! AVOGADRO'S NUMBER (AVOG)  = 6.022094D23MOL**-1
! SPEED OF LIGHT    (SL)    = 137.036 ATOMIC UNITS
! BOHR RADIUS       (A0B)   = 5.291771D-9 CM
! HARTREE ENERGY    (HRKEV) = 2.721160D-2 KEV

program mcsmce

  IMPLICIT none

  ! CROSS-SECTION PARAMETERS.

  REAL*8 Z1(4),Z2(4),AL12(4),AL22(4),AA1(4),AA2(4),AA3(4),AA4(4), &
       SE0(4),SE1(4),SE2(4),SET(4),EIFP(4),AB1(4),AB2(4),AB3(4)
  integer is(4)

  ! COUNTERS
  REAL*8 DDF(51),TLM(3),SAEM(3),SAIM(3),EM(2),PSIM(2),DLAT(2)
  integer N(3),NE(2,51),NTL(3,51),NAE(2,50,18),NANG(91)

  ! VARIANCES
  REAL*8 VDDF(51),VTLM(3),VSAEM(3),VSAIM(3),VEM(2),VPSIM(2), &
       VDLAT(2), vscat(2), ntet(2)

  real*8 rd(51)

  integer ir, iw
  real*8 AVOG, A0B, HRKEV
  real*8 SL2, SL22, SL22I
  integer iz, i
  real*8 aw, Z
  real*8 aa, al1, al2, fnt, a1a, expot, fie, d, rho, vmol, expot2
  real*8 e0, e0kev, e0mev, eabs
  real*8 chic2, chia2, f, xnu, vartet, rmstet, abin
  real*8 eloss, debin, zt, c1
  integer ntotal, ntot, itest, nscat, nloss
  real*8 x, y, tl, w, df, u, v, e, einv, fa, rb, beta, sfact, qmi, qfact, qme
  real*8 q02, ft, xifp, fb, tifp, afp, csi, wa, wm, w0, w1, s, rnd
  integer ior, ien, j, ian
  real*8 de, cdt, fcor

  ! INPUT-OUTPUT UNITS
  DATA IR,IW/5,6/

  ! PHYSICAL CONSTANTS: Avogadro, Bohr, Hartree
  DATA AVOG, A0B, HRKEV / 6.022094D23, 5.291771D-9, 2.721160D-2 /

  ! POWERS OF THE SPEED OF LIGHT
  DATA SL2, SL22, SL22I / 1.877887D4, 3.755773D4, 2.662568D-5 /

  ! COUNTERS INITIALIZATION

  DATA DDF,TLM,SAEM,SAIM,EM,PSIM,DLAT/66*0D0/
  DATA NTOT,N,NE,NTL,NAE,NANG/2150*0/
  DATA VDDF,VTLM,VSAEM,VSAIM,VEM,VPSIM,VDLAT/66*0D0/
  data vscat / 2*0d0 /
  data ntet / 2*0 /

  WRITE( IW,2003 )
2003 FORMAT( //6X,'input data:' )

  ! tracks:

  ntotal = 10*1000
  WRITE( IW, 2018 ) NTOTAL
2018 FORMAT( /6X, 'NUMBER OF TRAJECTORIES DESIRED .......  ', I5 )

  ! THICKNESS:

  D = 0.01   ! [cm]
  WRITE( IW, 2012 ) D*1e4
2012 FORMAT( /6X, 'THICKNESS ............................. ', f7.1, ' um' )

  e0mev = 1000d0
  WRITE( IW, 2009 ) e0mev
2009 FORMAT( /6X, 'INITIAL ELECTRON ENERGY ( E0 ) ........ ', F6.3, ' MeV' )

  ! ABSORPTION ENERGY

  eabs = 10d0               ! keV
  WRITE( IW, 2010 ) EABS
2010 FORMAT( /6X, 'LOWER ELECTRON ENERGY ( EABS ) ........ ', F7.3, ' keV' )

  ! ELEMENT DATA:

  iz = 13                   ! Al
  aw = 26.984d0             ! g/mol
  rho = 2.7020d0            ! g/cm3

  WRITE( IW, 2004 ) IZ, AW, rho
2004 FORMAT( &
          /6X, 'ATOMIC NUMBER = ', I3 &
          /6X, 'ATOMIC WEIGHT = ', F7.3, ' g/MOL' &
          /6X, 'DENSITY       = ', F8.4, ' G/CM3' )

  ! SCREENING PARAMETERS:

  AA = 0.31543d0
  AL1 = 7.28752d0
  AL2 = 1.77451d0

  print 2005, AA, AL1, AL2
2005 FORMAT( &
       /6X, 'SCREENING CONSTANTS:' &
       /6X, 'A   = ', F8.5 &
       /6X, 'ALPHA1 = ', F8.5 &
       /6X, 'ALPHA2 = ', F8.5 )

  ! CROSS-SECTION PARAMETERS:

  ZT = 6.283185307179586D0 * Z ! 2pi Z
  I = 1
  Z = DFLOAT( IZ )
  Z1( I ) = 2D0 * 6.283185307179586D0 * Z ! 4pi*Z
  Z2( I ) = Z1( I )*Z
  A1A = 1D0 - AA
  AA1( I ) = AA*AA

  IF( DABS( A1A ) .gT. 1D-3 ) then

     IS( I ) = 1
     AL12( I ) = AL1*AL1
     AL22( I ) = AL2*AL2
     AA2( I ) = A1A*A1A
     AA3( I ) = 2D0*AA*A1A/( AL12( I )-AL22( I ) )
     AA4( I ) = DLOG( AL22( I )/AL12( I ) )
     AB2( I ) = 2D0*AA/AL12( I ) - AA3( I )
     AB3( I ) = 2D0*A1A/AL22( I ) + AA3( I )

  else

     IS( I ) = 0
     AL12( I ) = AL1*AL1
     AB2( I ) = 2D0/AL12( I )

  endif

  print *, "elastic screening mode IS ", IS(1)

  ! MEAN EXCITATION POTENTIAL ( SEMIEMPIRICAL FORMULA )

  IF( IZ .LT. 11 ) then
     EXPOT = 1.35D-2*Z
  else
     EXPOT = Z*( 9.76D-3 + 5.88D-2 / Z**1.19D0 )
  endif

  WRITE( IW, 2007 ) EXPOT
2007 FORMAT( /6X, 'MEAN EXCITATION POTENTIAL ............. ', F7.3,' KEV' )

  EXPOT = EXPOT/HRKEV

  ! FIRST IONIZATION ENERGY:

  fie = 6d0                 ! eV
  IF( FIE .LE. 0.5D0 ) FIE = 5D0
  WRITE( IW, 2006 ) FIE
2006 FORMAT( 17X, 'FIRST IONIZATION ENERGY = ', F7.3, ' EV' )

  FIE = FIE*1D-3/HRKEV      ! Hartree
  AB1( I ) = 0.125D0*FIE*FIE

  ! NUMBER OF 'MOLECULES' PER UNIT VOLUME
  VMOL = ( AVOG*RHO/AW )*A0B**3 ! a0b = Bohr

  e0kev = e0mev*1d3
  e0 = e0kev

  E0 = E0/HRKEV             ! atomic units
  EABS = EABS/HRKEV

  EXPOT = 1D0/EXPOT
  EXPOT2 = EXPOT*EXPOT

  ! Moliere angles:

  chic2 = 0.157 * Z*(Z+1) * d*rho / Aw / e0mev / e0mev
  chia2 = 2.007d-5 * Z**(2d0/3d0) * ( 1 + 3.34*Z/137d0*Z/137d0 ) / e0mev / e0mev

  print *, "Moliere chic ", sqrt(chic2)
  print *, "Moliere chia ", sqrt(chia2)
  print *, "Moliere ratio ", chic2/chia2

  ! Lynch-Dahl rms scattering:

  f = 0.98                  ! central fraction
  xnu = 0.5*chic2/chia2/(1-f)
  vartet = chic2/(1+f*f) * ( (1+xnu)/xnu*log(1+xnu) - 1 )
  rmstet = sqrt(vartet) * sqrt(2d0) ! 2-D
  print *, "Lynch-Dahl RMS scattering ", rmstet*1d3, " mrad"

  abin = 0.1d0 * rmstet

  eloss = 300d-3 * d*1d4     ! MIP [keV]
  dEbin = 5*eloss/50d0

  D = D / A0B               ! atomic units

  ! INITIAL DIRECTION

  c1 = 1d0                  ! forward

  ! tracks:

200 NTOT = NTOT+1

  ! ---------------------------- INITIAL POSITION

  X = 0D0
  Y = 0D0
  Z = 0d0

  ITEST = 0
  TL = 0D0
  NSCAT = 0
  NLOSS = 0

  W = C1                    ! direction cosines
  df = 1d0                  ! phi
  U = DCOS( DF )
  V = DSIN( DF )

  E = E0                    ! atomic units

  ! ------------------------- BORN CROSS-SECTIONS
  ! ( IN UNITS OF BETA*SL2 )

301 EINV = 1D0/E
  FA = E + SL2 ! E+mc2
  RB = E + SL22  ! E+2mc2
  BETA = E*RB/( FA*FA )
  SFACT = BETA*SL2
  QMI = ( E+E )*RB*SL22I
  QFACT = QMI + QMI
  QME = QFACT + QFACT

  ! print *, "  E ", e

  I = 1

  Q02 = AB1( I )*EINV

  IF( IS( I ) .EQ. 0 ) then !  WENTZEL MODEL

     EIFP( I ) = Z2( I )*QME/( AL12( I )*( AL12( I )+QME ) )
     FA = AL12( I ) + Q02
     FT = AL12( I ) + QMI
     XIFP  =  Z1( I ) * ( AB2( I )*DLOG( QMI*FA/( FT*Q02 ) ) + AA1( I )*( 1D0/FT-1D0/FA )  )

  else

     ! ELASTIC CROSS-SECTION ( MOTT'S FORMULA )

     FA = AL12( I ) + QME
     FB = AL22( I ) + QME
     SE0( I ) = DLOG( FB/FA )-AA4( I )
     SE1( I ) = AA1( I )*( QME/( AL12( I )*FA ) )
     SE2( I ) = AA2( I )*( QME/( AL22( I )*FB ) )
     SET( I ) = SE1( I )+SE2( I )+AA3( I )*SE0( I )
     EIFP( I ) = Z2( I )*SET( I )

     ! INELASTIC CROSS-SECTION ( LENZ'S FORMULA )

     FA = AL12( I ) + Q02
     FB = AL22( I ) + Q02
     FT = AL12( I ) + QMI
     FNT = AL22( I ) + QMI
     XIFP = Z1( I ) * ( &
          AB2( I )*DLOG( QMI*FA/( FT*Q02 ) ) + &
          AB3( I )*DLOG( QMI*FB/( FNT*Q02 ) ) + &
          AA1( I )*( 1D0/FT-1D0/FA ) + &
          AA2( I )*( 1D0/FNT-1D0/FB )  )
  endif

  ! ------------------------------ MEAN FREE PATH

  TIFP = eIFP(i) + XIFP
  AFP = SFACT / ( VMOL*TIFP )

  ! print *, "  path ", eifp(i), xifp, afp, afp*a0b*1d4

  ! --------------------- STRAGGLING DISTRIBUTION
  ! AVERAGE ENERGY LOSS PER INELASTIC COLLISION.

  CSI = E*EXPOT
  IF( CSI .lt. 6.338065465611359D0 ) then
     WA = ZT*EINV * 3.1776907270853916D0 * DSQRT( CSI ) / XIFP
  else
     FA = 1D0 - BETA
     FB = DSQRT( FA )
     FT = 1D0 - FB
     WA  =  2D0*ZT*EINV * ( DLOG( E*EXPOT2*SFACT/FA ) - &
          ( FB+FB+BETA ) * 6.931471805599453D-1 + &
          FA + 0.125D0*FT*FT ) / XIFP
  endif

  ! print *, "expot ", expot, ", csi", csi, ", wavg ", wa

  ! MINIMUM ENERGY LOSS. HALLEY'S METHOD.

  WM = WA**( 1D0 + 3.3D0 / ( 3.3D0 - DLOG( WA ) ) )
  FT = DLOG( WM )
  ! print *, "wm ", wm, ", ft ", ft

  W0 = WA*( 1D0-WM ) + WM*FT
  W1 = 1D0 - WA + FT
  WM = WM - 2D0*W0*W1 / ( 2D0*W1*W1 - W0/WM )
  FT = DLOG( WM )
  ! print *, "wm ", wm, ", ft ", ft

  W0 = WA*( 1D0 - WM ) + WM*FT
  W1 = 1D0 - WA + FT
  WM = WM - 2D0*W0*W1 / ( 2D0*W1*W1 - W0/WM )
  ! print *, "  wmin ", wm

  ! ----------------------------------- FREE PATH

401 call random_number( rnd ) ! uniform 0..1
  S = -DLOG( 1d0-RND )*AFP

  TL = TL + S

  ! ---------------------------------------------------- NEW POSITION

  X = X + U*S
  Y = Y + V*S
  Z = Z + W*S

  ! print *, "  z ", z*a0b*1d4 ! [um]

  ! ----------------------------------------------------- TEST BOUNDS

  IF( Z .LT. D ) GO TO 610    ! inside
  ITEST = 1
  GO TO 500

610 CONTINUE
  IF( Z .GT. 0D0 ) GO TO 620  ! inside
  ITEST = 2
  GO TO 500

620 CONTINUE

  ! --------------------------------------------------- NEW COLLISION

  call random_number( rnd )
  FA = RND*TIFP
  I = 1
  FNT = EIFP( I )
  IF( FNT .LT. FA ) GO TO 402

  IOR = I
  GO TO 403

402 CONTINUE

  ! ------------------------- INELASTIC COLLISION

  call random_number( ft )
  DE = 0.5D0 * E*WM / ( 1D0 - FT * ( 1D0 - WM ) )
  NLOSS = NLOSS + 1
  ! NEW ENERGY AND TEST
  E = E-DE

  IF( E .GT. EABS ) GO TO 640

  ITEST = 3                 ! absorbed
  GO TO 500

640 CONTINUE

  ! SCATTERING ANGLES

  CDT = DSQRT( ( 1D0 - DE*EINV ) * RB / ( RB-DE ) )
  call random_number( rnd )
  DF = RND * 6.283185307179586D0 ! 0..2pi

  CALL DIRECT( CDT, DF, U, V, W )

  ! print *, "  inelastic ", cdt, de

  GO TO 301

  ! --------------------------- ELASTIC COLLISION

403 call random_number( fa )

  IF( IS( IOR ) .EQ. 0 ) GO TO 405

  ! SELECTION OF THE PARTIAL DISTRIBUTION

  call random_number( rnd )
  FB = RND*SET( IOR )
  FNT = SE1( IOR )
  IF( FNT .GT. FB ) GO TO 405

  FNT = FNT + SE2( IOR )
  IF( FNT .GT. FB ) GO TO 404

  ! SCATTERING ANGLES
  FT = DEXP( FA*SE0( IOR )+AA4( IOR ) )
  CDT = ( AL12( IOR )*FT-AL22( IOR ) )/( 1D0-FT )
  GO TO 406

404 CDT = FA*AL22( IOR )*QME/( AL22( IOR )+( 1D0-FA )*QME )
  GO TO 406

405 CDT = FA*AL12( IOR )*QME/( AL12( IOR )+( 1D0-FA )*QME )

406 CDT = 1D0-CDT/QFACT
  call random_number( rnd )
  DF = RND*6.283185307179586D0 ! 0..2pi
  ! NEW DIRECTION
  CALL DIRECT( CDT, DF, U, V, W )
  NSCAT = NSCAT+1

  ! print *, "    elastic ", cdt

  GO TO 401

  ! ---------------------------------------------- INCREMENT COUNTERS

500 CONTINUE

  N( ITEST ) = N( ITEST ) + 1

  IF( ITEST .EQ. 3 ) GO TO 501 ! absorbed

  IEN = IDINT( ( e0kev - E*hrkev ) / DEbin ) + 1
  IF( IEN .GT. 50 ) IEN = 50
  NE( ITEST, IEN ) = NE( ITEST, IEN ) + 1

  EM( ITEST ) = EM( ITEST ) + E ! exit energy [atomic]
  VEM( ITEST ) = VEM( ITEST ) + E*E

  FA = DACOS( W )           ! exit angle, 0..pi
  J = IDINT( FA/abin )+1
  IF( j .GT. 50 ) j = 50
  NANG( J ) = NANG( J )+1

  PSIM( ITEST ) = PSIM( ITEST ) + FA
  VPSIM( ITEST ) = VPSIM( ITEST ) + FA*FA
  if( fa < 3*rmstet ) then
     vscat( ITEST ) = vscat( ITEST ) + FA*FA
     ntet( ITEST ) = ntet( ITEST ) + 1
  endif

  FB = DSQRT( X*X + Y*Y )*a0b*1e4 ! displacement [um]
  DLAT( ITEST ) = DLAT( ITEST ) + FB
  VDLAT( ITEST ) = VDLAT( ITEST ) + FB*FB

  IAN = IDINT( FA * 11.45915D0 ) + 1
  IF( ITEST .EQ. 2 ) IAN = IAN-18
  IF( IAN .GT. 18 ) IAN = 18
  NAE( ITEST, IEN, IAN ) = NAE( ITEST, IEN, IAN )+1

  if( ntot < 100 .or. mod(ntot,100) == 0 ) then
     print *, ntot, ": scat ", nscat, ", loss ", nloss, &
          ", angle ", fa*1d3, " mrad, offset ", fb, " um"
  endif

501 TLM( ITEST ) = TLM( ITEST ) + TL
  VTLM( ITEST ) = VTLM( ITEST ) + TL*TL

  SAEM( ITEST ) = SAEM( ITEST ) + DFLOAT( NSCAT )
  VSAEM( ITEST ) = VSAEM( ITEST ) + DFLOAT( NSCAT )**2

  SAIM( ITEST ) = SAIM( ITEST ) + DFLOAT( NLOSS )
  VSAIM( ITEST ) = VSAIM( ITEST )+DFLOAT( NLOSS )**2

  IF( NTOT .LT. NTOTAL ) GO TO 200
  !
  ! --------------------------------------------- COMPUTE MEAN VALUES
  ! ITEST:  1- TRANSMITTED,   2- BACKSCATTERED,   3- ABSORBED

701 FNT = DFLOAT( NTOT )

  Z1( 1 ) = DFLOAT( N( 1 ) )
  Z1( 2 ) = DFLOAT( N( 2 ) )
  Z1( 3 ) = DFLOAT( N( 3 ) )

  FT = Z1( 1 ) / FNT
  FB = Z1( 2 ) / FNT
  FA = Z1( 3 ) / FNT

  FCOR = 1D4*A0B        ! um

  DO 702 I = 1, 3

     IF( Z1( I ).EQ.0D0 ) Z1( I ) = 1D0
     DF = 1D0 / Z1( I ) ! 1/N

     VTLM( I ) = DF*FCOR * DSQRT( VTLM( I ) - DF * TLM( I )**2 )
     TLM( I ) = TLM( I ) * FCOR*DF ! track length

     VSAEM( I ) = DF * DSQRT( VSAEM( I ) - DF * SAEM( I )**2 )
     SAEM( I ) = SAEM( I )*DF ! elastic scatterings

     VSAIM( I ) = DF * DSQRT( VSAIM( I ) - DF * SAIM( I )**2 )
     SAIM( I ) = SAIM( I )*DF ! inelastic scatterings

     IF( I .EQ. 3 ) GO TO 702 ! absorbed

     VEM( I ) = DF*HRKEV * DSQRT( VEM( I ) - DF * EM( I )**2 )
     EM( I ) = HRKEV*EM( I )*DF

     VPSIM( I ) = DSQRT( VPSIM( I ) * DF ) ! rms scat ang
     vscat( I ) = DSQRT( vscat( I ) / ntet(i) ) ! truncated

     VDLAT( I ) = DF * DSQRT( VDLAT( I ) - DF * DLAT( I )**2 )
     DLAT( I ) = DLAT( I )*DF

702 enddo ! cases

  ! **************************************************************
  ! *      ---------------   OUTPUT      *
  ! **************************************************************

  WRITE( IW, 2100 )
2100 FORMAT( /6X, '---------- RESULTS ----------' )

  WRITE( IW, 2101 ) NTOT
2101 FORMAT( //6X, 'NUMBER OF COMPUTED TRAJECTORIES ...... ', I7 )

  WRITE( IW, 2103 ) N( 1 )
2103 FORMAT( /6X, 'NUMBER OF TRANSMITTED ELECTRONS ...... ', I7 )

  WRITE( IW, 2104 ) N( 2 )
2104 FORMAT( 6X, 'NUMBER OF BACKSCATTERED ELECTRONS .... ', I7 )

  WRITE( IW, 2105 ) N( 3 )
2105 FORMAT( 6X, 'NUMBER OF ABSORBED ELECTRONS ......... ', I7 )

  WRITE( IW, 2112 ) SAEM( 1 ), VSAEM( 1 )
2112 FORMAT( /6X, 'MEAN NUMBER OF ELASTIC SCATTERING ACTS:'/13X, &
          'TRANSMITTED ELECTRONS ........... ', F8.1, ' +-', D5.2 )

  WRITE( IW, 2115 ) SAIM( 1 ), VSAIM( 1 )
2115 FORMAT( /6X, 'MEAN NUMBER OF INELASTIC SCATTERING ACTS:'/13X, &
          'TRANSMITTED ELECTRONS ........... ', F8.1, ' +-', D5.2 )

  WRITE( IW, 2116 ) EM( 1 )*1d-3, VEM( 1 )*1d-3
2116 FORMAT( /6X, 'MEAN ENERGY:'/13X, &
          'TRANSMITTED ELECTRONS ......... ', F9.3, ' +-', &
          F9.3, ' MeV' )

  WRITE( IW, 2117 ) E0kev-EM( 1 ), VEM( 1 )
2117 FORMAT( /6X, 'MEAN ENERGY loss:'/13X, &
          'TRANSMITTED ELECTRONS ......... ', F9.3, ' +-', &
          F9.3, ' KEV' )

  WRITE( IW, 2118 ) VPSIM( 1 )*1d3
2118 FORMAT( /6X, 'RMS THETA:' / &
          13X, 'TRANSMITTED ELECTRONS ......... ', F12.3, ' mrad' )

  WRITE( IW, 2119 ) vscat( 1 )*1d3
2119 FORMAT( /6X, 'truncated RMS THETA:' / &
          13X, 'TRANSMITTED ELECTRONS ......... ', F12.3, ' mrad' )

  WRITE( IW, 2120 ) DLAT( 1 ), VDLAT( 1 )
2120 FORMAT( /6X, 'MEAN LATERAL DISPLACEMENT:'/13X, &
          'TRANSMITTED ELECTRONS......... ', f11.1, ' +-', &
          f8.3, ' mum' )

  ! scattering angle:

  IF( N( 1 ) .LE. 0 ) GO TO 705

  DO 704 I = 1, 50
     RD( I ) = DFLOAT( I )*abin*1d3
     DDF( I ) = DFLOAT( NANG( I ) )
     VDDF( I ) = DSQRT( DDF( I ) - DDF( I )*DDF( I ) / FNT )
704 enddo
  WRITE( IW, 2202 )
2202 FORMAT( //15X &
          'ANGLE OF TRANSMITTED ELECTRONS: Y(N), X(mrad)'/ )
  CALL GRAF( RD, DDF, VDDF, 1, 50, IW )

705 CONTINUE

  ! energy loss distribution:

  DO 708 I = 1, 51
     RD( I ) = DFLOAT( I )*debin
     DDF( I ) = DFLOAT( NE( 1, I ) )
     VDDF( I ) = DSQRT( DDF( I )-DDF( I )*DDF( I )/FNT )
708 enddo
  DDF( 50 ) = DDF( 50 )+DDF( 51 )
  VDDF( 50 ) = VDDF( 50 )+VDDF( 51 )

  IF( N( 1 ) .LE. 0 ) GO TO 709

  WRITE( IW, 2204 )
2204 FORMAT( //15X, &
          'ENERGY LOSS OF TRANSMITTED ELECTRONS: Y(N), X(keV)'/ )
  CALL GRAF( RD, DDF, VDDF, 1, 50, IW )

709 CONTINUE

END program mcsmce

! ******************************************************************
SUBROUTINE DIRECT( CDT, DF, U, V, W )

  ! update DIRECTION COSINES AFTER A DEFLECTION DT, DF IN THE
  ! PARTICLE REFERENCE SYSTEM.

  ! INPUT:
  ! U, V, W ..... INITIAL DIRECTION COSINES
  ! CDT ....... COSINE OF THE POLAR SCATTERING ANGLE
  ! DF ........ AZIMUTHAL SCATTERING ANGLE ( RAD )

  ! OUTPUT:
  ! U, V, W ..... NEW DIRECTION COSINES
  ! (CDT AND DF REMAIN UNCHANGED)

  IMPLICIT none
  real*8 cdt, df, u, v, w
  real*8 din(3), sdt, cz, sz, f, sf, cf, aa

  sdt = DSQRT( 1D0 - CDT*CDT )

  DIN(1) = SDT*COS(DF)
  DIN(2) = SDT*SIN(DF)
  DIN(3) = CDT

  cz = w                    ! old direction
  sz = sqrt( 1.0 - cz*cz )

  f = atan2( v, u )
  sf = sin(f)
  cf = cos(f)
  U = cz*cf * din(1) - sf * din(2) + sz*cf * din(3)
  V = cz*sf * din(1) + cf * din(2) + sz*sf * din(3)
  W =-sz    * din(1)               + cz    * din(3)

END SUBROUTINE DIRECT

! ******************************************************************
SUBROUTINE GRAF( X, Y, DY, NF, NL, IW )

  ! THIS SUBROUTINE PLOTES THE DISCRETE POINTS ( X, Y ) OF A DOUBLE
  ! ARRAY. THE VALUES OF THE X VARIABLE HAVE TO BE EQUALY SPACED.
  ! IF DY ( THE ABSOLUTE ERROR IN THE Y VALUES ) IS NOT ZERO,  THE
  ! ERROR BARS ARE ALSO PLOTTED.
  ! X  ............  ORDERED ABCISES.
  ! Y  ............  CORRESPONDING ORDINATES.
  ! DY ............  ABSOLUTE ERROR IN Y.
  ! NF,  NL ........  FIRST AND LAST PLOTTED POINTS.
  ! IW ............  OUTPUT UNIT.

  IMPLICIT REAL*8 ( A-H, O-Z )
  REAL*8 X( 1 ), Y( 1 ), DY( 1 ), S( 9 )
  character*1 L(91)
  character*1 l1
  character*1 l2
  character*1 l3
  character*1 l4
  character*1 l5
  character*1 l6
  character*1 l7
  DATA L1/'+'/, L2/' '/, L3/'I'/, L4/'*'/, L5/'-'/, L6/'('/
  data L7/')'/

  YMAX = 0D0
  YMIN = 0D0

  DO 1 I = NF, NL
     E = DABS( DY( I ) )
     YMAX = DMAX1( Y( I )+E, Y( I )-E, YMAX )
     YMIN = DMIN1( Y( I )+E, Y( I )-E, YMIN )
1 enddo
  F = YMAX-YMIN
  IF( F.LE.0D0 ) RETURN
  IF( YMIN.LT.0D0 ) GO TO 2
  IZERO = 1
  D = YMAX/1.8D1
  GO TO 5
2 CONTINUE
  IF( YMAX.GT.0D0 ) GO TO 3
  IZERO = 19
  D = -YMIN/1.8D1
  GO TO 5
3 K = 19
4 K = K-1
  D = F/DFLOAT( K )
  IZERO = IDINT( -YMIN/D )+1
  I = IDINT( YMAX/D )+1
  IF( IZERO+I.GT.18 ) GO TO 4
  IZERO = IZERO+1
5 F = 1D1
  IF( D.GT.1D1 ) F = 0.1D0
  D = D/F
  DO 6 I = 1, 90
     D = D*F
     IF( D.LT.1D1.OR.D.GE.1D2 ) GO TO 6
     D = DFLOAT( IDINT( D )+1 )/F**( I-1 )
     GO TO 7
6 enddo
7 DO 8 I = 1, 9
     S( I ) = DFLOAT( I+I-IZERO )*D
8 enddo
  D = 0.2D0*D
  IZERO = 5*( IZERO-1 )+1
  YMIN = ( DFLOAT( 1-IZERO )-0.5D0 )*D
  WRITE( IW, 101 ) ( S( I ), I = 1, 9 )
101 FORMAT( 7X, 'Y :', 3X, 1PD9.2, 1X, 8( D9.2, 1X ) )
  WRITE( IW, 102 )
102 FORMAT( ' X :', 9X, 9( '+----I----' ), '+', 4X, 'Y +- DY :' )

  DO 12 I = NF, NL
     L( 1 ) = L1
     L( 91 ) = L1
     DO 9 J = 2, 90
        L( J ) = L2
9    enddo
     L( IZERO ) = L3
     E = DABS( DY( I ) )
     K = IDINT( ( Y( I )-YMIN )/D )+1
     K1 = IDINT( ( Y( I )-E-YMIN )/D )+1
     K2 = IDINT( ( Y( I )+E-YMIN )/D )+1
     IF( K1.GE.K2 ) GO TO 11
     DO 10 J = K1, K2
        L( J ) = L5
10   enddo
     L( K1 ) = L6
     L( K2 ) = L7
11   L( K ) = L4
     WRITE( IW, 103 ) X( I ), ( L( K ), K = 1, 91 ), Y( I ), E
12 enddo

103 FORMAT( 1PD10.3, 3X, 91A1, 3X, D10.3, ' +-', D8.1 )

  WRITE( IW, 104 )
  WRITE( IW, 105 ) ( S( I ), I = 1, 9 )
104 FORMAT( 13X, 9( '+----I----' ), '+' )
105 FORMAT( 13X, 1PD9.2, 1X, 8( D9.2, 1X ) )

END SUBROUTINE GRAF
