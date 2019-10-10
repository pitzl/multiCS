
// Multiple Coulomb scattering Monte Carlo for electron in pure elements

// based on mcsda.f
// F. Salvat, J.D. Martinez, R. Mayol, J. Parellada
// COMP. PHYS. COMMUN. 42 (1986) 93 AALC

// C++ with ROOT by Daniel Pitzl, DESY, Sep 2019

// ELASTIC SCATTERING EVENTS ARE DESCRIBED BY USING AN ANALYTICAL FIT
// TO THE RELATIVISTIC MOTT CROSS-SECTION.
// THE MEAN RATE OF ENERGY LOSS PER UNIT PATH LENGTH IS DERIVED FROM
// THE RELATIVISTIC BETHE STOPPING POWER, I.E. FROM THE CONTINUOUS SLOWING
// DOWN APPROXIMATION (CSDA).
// THE STRAGGLING EFFECT IS SIMULATED BY MEANS OF THE LENZ INELASTIC (TOTAL)
// CROSS-SECTION AND THE CLASSICAL 1/W**2 DISTRIBUTION OF ENERGY LOSSES.
// ANGULAR DEFLECTIONS OF THE ELECTRON PATH DUE TO INELASTIC COLLISIONS ARE
// COMPUTED FROM THE CLASSICAL THEORY OF BINARY COLLISIONS.

// REFERENCES:
// (1) F. SALVAT ET AL., J. PHYS. D: APPL. PHYS. 17 (1984) 185
// (2) F. SALVAT ET AL., J. PHYS. D: APPL. PHYS. 17 (1984) 1545
// (3) R. MAYOL ET AL., ANALES DE FISICA A80 (1984) 130
// (4) M.J. BERGER AND S.M. SELTZER, NBSIR 82-2550-A (1982)
//
// *****    ATOMIC DATA
//
// FREE ATOMS
// --------------------------------------------------------------
// ATOM  Z   AW         A         ALPHA1    ALPHA2     I
// --------------------------------------------------------------
// He    2   4.0026     1.00000   2.20564   2.20564   41.8 (GAS)
// Li    3   6.941      0.60453   2.81743   0.66246   40.0
// Be    4   9.01218    0.32784   4.54284   0.98509   63.7
// B     5   10.81      0.23268   5.99011   1.21345   76.0
// C     6   12.011     0.15364   8.04166   1.49140   78.0
// N     7   14.0067    0.09957   10.8127   1.76869   82.0 (GAS)
// O     8   15.9994    0.06245   14.8346   2.04052   95.0 (GAS)
// F     9   18.9984    0.03678   21.4093   2.30606   115  (GAS)
// Ne   10   20.179     0.01801   35.0596   2.56643   137  (GAS)
// Na   11   22.98977   0.74440   4.12051   0.87184   149
// Mg   12   24.305     0.64234   4.72661   1.00249   156
// Al   13   26.98154   0.60013   5.14064   1.01539   166
// Si   14   28.0855    0.51591   5.84957   1.17330   173
// P    15   30.97376   0.43866   6.67081   1.34106   173
// S    16   32.06      0.37130   7.60432   1.50676   180
// Cl   17   35.453     0.31314   8.67077   1.66877   174  (GAS)
// Ar   18   39.948     0.26299   9.90260   1.82787   188  (GAS)
// K    19   39.0983    0.62732   5.87615   0.98798   190
// Ca   20   40.08      0.58800   6.32192   1.00940   191
// Sc   21   44.9559    0.55431   6.63300   1.10236   216
// Ti   22   47.88      0.53593   6.89381   1.18766   233
// V    23   50.9415    0.52216   7.11955   1.26815   245
// Cr   24   51.996     0.50646   7.31623   1.42531   257
// Mn   25   54.9380    0.50356   7.49501   1.41970   272
// Fe   26   55.847     0.49731   7.65419   1.49227   286
// Co   27   58.9332    0.49262   7.79761   1.56304   297
// Ni   28   58.69      0.48921   7.92752   1.63224   311
// Cu   29   63.546     0.46600   8.20246   1.82639   322
// Zn   30   65.38      0.48531   8.15505   1.76715   330
// Ga   31   69.72      0.57810   7.45294   1.47922   334
// Ge   32   72.59      0.55193   7.81505   1.52669   350
// As   33   74.9216    0.51813   8.27460   1.60648   347
// Se   34   78.96      0.48352   8.78654   1.69690   348
// Br   35   79.904     0.45033   9.33526   1.79016   357
// Kr   36   83.80      0.41895   9.91986   1.88404   352  (GAS)
// Rb   37   85.4678    0.66732   7.38694   1.13156   363
// Sr   38   87.62      0.65417   7.60914   1.10652   366
// Y    39   88.9059    0.62196   7.98456   1.19369   379
// Zr   40   91.22      0.59406   8.34255   1.27785   393
// Nb   41   92.9064    0.55591   8.81359   1.43082   417
// Mo   42   95.94      0.53108   9.19190   1.51654   424
// Tc   43   97.907     0.52988   9.33183   1.50495   428
// Ru   44   101.07     0.48905   9.92161   1.67634   441
// Rh   45   102.9055   0.47090   10.2776   1.75229   449
// Pd   46   106.42     0.36860   12.2194   2.09718   470
// Ag   47   107.868    0.43878   10.9799   1.18991   470
// Cd   48   112.41     0.45464   10.8499   1.84361   469
// In   49   114.82     0.52418   9.95823   1.62066   488
// Sn   50   118.69     0.51613   10.1882   1.62396   488
// Sb   51   121.75     0.50244   10.5023   1.64840   487
// Te   52   127.60     0.48399   10.9022   1.69203   485
// I    53   126.9045   0.46451   11.3435   1.74292   491
// Xe   54   131.29     0.44503   11.8150   1.79736   482  (GAS)
// Cs   55   132.9054   0.63099   9.27951   1.21804   488
// Ba   56   137.33     0.63153   9.37847   1.16983   491
// La   57   138.9055   0.60457   9.78322   1.23954   501
// Ce   58   140.12     0.60176   9.89662   1.26964   523
// Pr   59   140.9077   0.59961   10.0029   1.29867   535
// Nd   60   144.24     0.61818   9.84443   1.29305   546
// Pm   61   144.913    0.61655   9.94325   1.32107   560
// Sm   62   150.36     0.61524   10.0389   1.34861   574
// Eu   63   151.96     0.61408   10.0134   1.37597   580
// Gd   64   157.25     0.59470   10.4721   1.43349   591
// Tb   65   158.9254   0.61231   10.3185   1.42998   614
// Dy   66   162.5      0.61163   10.4097   1.45675   628
// Ho   67   164.9304   0.61105   10.5006   1.48342   650
// Er   68   167.26     0.61054   10.5913   1.51003   658
// Tm   69   168.9342   0.61006   10.6828   1.53672   674
// Yb   70   173.04     0.60962   10.7746   1.56339   684
// Lu   71   174.967    0.59607   11.0553   1.60432   694
// Hf   72   178.49     0.57730   11.4181   1.66925   705
// Ta   73   180.9479   0.55928   11.7897   1.73462   718
// W    74   183.85     0.54238   12.1611   1.79841   727
// Re   75   186.207    0.52646   12.5334   1.85940   736
// Os   76   190.2      0.51112   12.9161   1.92009   746
// Ir   77   192.22     0.49632   13.3110   1.98039   757
// Pt   78   195.08     0.46317   14.0902   2.12491   790
// Au   79   196.9665   0.44855   14.5490   2.18774   790
// Hg   80   200.59     0.45535   14.5388   2.15774   800
// Tl   81   204.383    0.51102   13.5079   1.94984   810
// Pb   82   207.2      0.51405   13.5898   1.92284   823
// Bi   83   208.9804   0.52372   13.5429   1.87088   823
// Po   84   208.982    0.51902   13.7725   1.87223   830
// At   85   209.987    0.50987   14.0959   1.89216   825
// Rn   86   222.018    0.49903   14.4661   1.92060   794  (GAS)
// Fr   87   223.020    0.63069   12.2398   1.42884   827
// Ra   88   226.0254   0.63838   12.2521   1.36251   826
// Ac   89   227.0278   0.50967   14.1001   1.89316   841
// Th   90   232.0381   0.60184   13.0792   1.46078   847
// Pa   91   231.0359   0.60180   13.1874   1.48806   878
// U    92   238.0289   0.59396   13.4427   1.52547   890

// WIGNER-SEITZ ATOMs:
// ----------------------------------------------
// ATOM  Z   RWS      A       ALPHA1    ALPHA2
// ----------------------------------------------
// Li    3   3.25   0.37731   3.19082   1.28284
// Na   11   3.94   0.49868   5.08096   1.51901
// Al   13   2.99   0.31543   7.28752   1.77451
// Si   14   3.18   0.31186   7.76730   1.71824
// P    15   3.57   0.31007   8.17973   1.67147
// K    19   4.87   0.42851   7.49488   1.42250
// Ca   20   4.11   0.40255   8.05513   1.42620
// V    23   2.81   0.22974   12.2348   1.96899
// Cr   24   2.68   0.19633   13.6659   2.14626
// Mn   25   2.70   0.21946   12.9263   2.10354
// Fe   26   2.67   0.21948   13.0520   2.16581
// Ni   28   2.60   0.21892   13.2947   2.30276
// Cu   29   2.65   0.21374   13.4817   2.42997
// Ge   32   3.34   0.35634   10.2486   2.09847
// Rb   37   5.21   0.53376   8.61480   1.52016
// Sr   38   4.48   0.52902   8.81678   1.48175
// Nb   41   3.07   0.37795   11.4676   1.93243
// Mo   42   2.93   0.35710   12.0620   1.99956
// Rh   45   2.81   0.32521   13.2919   2.13494
// Pd   46   2.87   0.31116   13.8007   2.21950
// Ag   47   3.02   0.32982   13.4086   2.17450
// Cs   55   5.64   0.52081   10.6875   1.54308
// Ba   56   4.72   0.51198   10.9485   1.53695
// Ce   58   3.87   0.45633   12.1097   1.72503
// Eu   63   4.26   0.47914   12.1247   1.81620
// Yb   70   4.05   0.47457   12.8767   2.03486
// Ta   73   3.08   0.41974   14.4589   2.22411
// W    74   2.94   0.40817   14.8914   2.26673
// Ir   77   2.84   0.39894   15.5877   2.30243
// Pt   78   2.90   0.39375   15.8661   2.33737
// Au   79   2.98   0.39511   15.9768   2.34208
// Pb   82   3.65   0.42298   15.6563   2.24990
// Th   90   3.76   0.49879   15.0981   1.84375

// Z ................... ATOMIC NUMBER
// AW .................. ATOMIC WEIGHT
// A, ALPHA1,ALPHA2 .... DENSITY PARAMETERS DERIVED FROM SELF-CONSISTENT
// DHFS CALCULATIONS (REFERENCE 3)
// I ................... MEAN EXCITATION POTENTIAL FOR SINGLE-ELEMENT
// MATERIALS (EV) RECOMMENDED IN REFERENCE 4
// RWS ................. RADIUS OF THE WIGNER-SEITZ SPHERE (A.U.)

// *****    NUMERICAL AND PHYSICAL CONSTANTS IN THE SOURCE LISTING
// PI              = 3.1415926535897931D+00
// 36/PI           = 1.1459155902616464D+01
// 2*PI            = 6.2831853071795862D+00
// 2**31-1         = 2.1474836470000000D+09
// 1/(2**31-1)     = 4.6566128752457969D-10
// 2*PI/(2**31-1)  = 2.9258361598967678D-09
// (2*E**3)**1/2   = 6.3380654656113586D+00   E=DEXP(1.)
// 2*(2/E)**3/4    = 1.5888453635426958D+00
// 2/3*(E/2)**3/4  = 8.3918382740555711D-01
// RADIAN (RAD)    = 5.7295779513082323D+01 DEG
// RAD/2           = 2.8647889756541161D+01
// LOG(2)          = 6.9314718055994533D-01

#include <iostream> // cout
#include <cmath> // sqrt, sin, cos, log
#include <random>
#include <ctime>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>

using namespace std;

//------------------------------------------------------------------------------
void DIRECT( double CDT, double DF, double & u, double & v, double & w )
{
  // update DIRECTION COSINES AFTER A DEFLECTION DT, DF IN THE
  // PARTICLE REFERENCE SYSTEM.

  // INPUT:
  // U, V, W ..... INITIAL DIRECTION COSINES
  // CDT ....... COSINE OF THE POLAR SCATTERING ANGLE
  // DF ........ AZIMUTHAL SCATTERING ANGLE ( RAD )

  // OUTPUT:
  // U, V, W ..... NEW DIRECTION COSINES
  // (CDT AND DF REMAIN UNCHANGED)

  double din[3];

  double sdt = sqrt( 1 - CDT*CDT );

  din[0] = sdt*cos(DF);
  din[1] = sdt*sin(DF);
  din[2] = CDT;

  double cz = w;                    // old direction
  double sz = sqrt( 1.0 - cz*cz );

  double f = atan2( v, u );
  double sf = sin(f);
  double cf = cos(f);
  u = cz*cf * din[0] - sf * din[1] + sz*cf * din[2];
  v = cz*sf * din[0] + cf * din[1] + sz*sf * din[2];
  w =-sz    * din[0]               + cz    * din[2];

} // direct

//------------------------------------------------------------------------------
int main()
{
  // tracks:

  int ntotal = 10*1000;
  cout << "  NUMBER OF TRAJECTORIES DESIRED .. " << ntotal << endl;

  // THICKNESS:

  double thck = 0.0150; // [cm]
  cout << "  THICKNESS ....................... " << thck*10 << " mm" << endl;

  // incident energy:

  double e0mev = 5000;
  cout << "  INITIAL ELECTRON ENERGY ( E0 ) .. " << e0mev << " MeV" << endl;

  // ABSORPTION ENERGY:

  double eabs = 10;               // keV
  cout << "  LOWER ELECTRON ENERGY ( EABS ) .. " << eabs << " keV" << endl;

  // ELEMENT DATA:
  /*
  double Z = 13;       // Al
  double A = 26.984;   // g/mol
  double rho = 2.702;  // g/cm3
  */
  double Z = 14;       // Si
  double A = 28.0855;  // g/mol
  double rho = 2.329;  // g/cm3

  cout << "  ATOMIC NUMBER " << Z << endl
       << "  ATOMIC WEIGHT " << A << endl
       << "  DENSITY       " << rho << " g/cm3" << endl;

  // SCREENING PARAMETERS:
  /*
  double AA = 0.31543; // Al Wigner-Seitz
  double AL1 = 7.28752;
  double AL2 = 1.77451;
  */
  /*
  double AA = 0.51591; // Si free
  double AL1 = 5.84957;
  double AL2 = 1.17330;
  */
  double AA = 0.31186; // Si Wigner-Seitz
  double AL1 = 7.76730;
  double AL2 = 1.71824 ;

  cout << "  SCREENING CONSTANTS:" << endl
       << "  AA     = " << AA << endl
       << "  alpha1 = " << AL1 << endl
       << "  alpha2 = " << AL2 << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // atomic units: hbar = m = e = 1

  // AVOGADRO'S NUMBER (AVOG)  = 6.022094D23 MOL**-1
  // SPEED OF LIGHT    (SL)    = 137.036 ATOMIC UNITS
  // BOHR RADIUS       (a0b)   = 5.291771D-9 CM
  // HARTREE ENERGY    (HRKEV) = 2.721160D-2 KEV

  double AVOG = 6.022094e23; // [atoms/Mol]
  double a0b = 5.291771e-9; // Bohr radius [cm]
  double HRKEV  = 2.721160e-2; // Hartree aomic units

  // POWERS OF THE SPEED OF LIGHT in atomic units:

  double SL2 = 1.877887e4; // c^2
  double SL22 = 3.755773e4; // 2c^2
  double SL22I = 2.662568e-5; // 1/2c^2

  double pi = 4*atan(1);
  double twopi = 2*pi;

  // CROSS-SECTION PARAMETERS:

  double ZT = twopi * Z;
  double Z4 = 2*twopi*Z;
  double ZZ4 = Z*Z4;
  double A1A = 1 - AA;
  double AA1 = AA*AA;
  double AL12 = AL1*AL1;
  double AL22 = AL2*AL2;
  double AA2 = A1A*A1A;
  double AA3 = 2*AA*A1A/( AL12-AL22 );
  double AA4 = log( AL22/AL12 );
  double AB2 = 2*AA/AL12 - AA3;
  double AB3 = 2*A1A/AL22 + AA3;

  // MEAN EXCITATION POTENTIAL ( SEMIEMPIRICAL FORMULA )

  double expot = 1.35e-2 * Z;
  if( Z > 11.5 )
    expot = Z * ( 9.76e-3 + 5.88e-2 / pow( Z, 1.19 ) );

  cout << "  MEAN EXCITATION POTENTIAL ....... " << expot*1e3 << " eV" << endl;

  expot = expot/HRKEV;
  expot = 1/expot;
  double expot2 = expot*expot;

  // FIRST IONIZATION ENERGY:

  double fie = 6; // eV

  cout << "  FIRST IONIZATION ENERGY           " << fie << " eV" << endl;

  fie = fie*1e-3 / HRKEV;      // Hartree
  double AB1 = 0.125*fie*fie;

  double e0kev = e0mev*1e3;
  double E0 = e0kev/HRKEV;             // atomic units
  eabs = eabs/HRKEV;

  double eloss = 0.3 * thck*1e4;     // MIP [keV]

  // NUMBER OF 'MOLECULES' PER UNIT VOLUME
  double VMOL = AVOG*rho/A * pow( a0b, 3 ); // a0b = Bohr

  // Moliere angles:

  double chic2 = 0.157 * Z*(Z+1) * thck*rho / A / e0mev / e0mev;
  double chia2 = 2.007e-5 * pow( Z, 2.0/3.0 ) * ( 1 + 3.34*Z/137*Z/137 ) / e0mev / e0mev;

  cout << endl;
  cout << "Moliere chic " << sqrt(chic2) << endl;
  cout << "Moliere chia " << sqrt(chia2) << endl;
  cout << "Moliere ratio " << chic2/chia2 << endl;

  // Lynch-Dahl rms scattering:

  double f = 0.98;                  // central fraction
  double xnu = 0.5 * chic2 / chia2 / (1-f);
  double vartet = chic2 / (1+f*f) * ( (1+xnu) / xnu * log(1+xnu) - 1 );
  double rmstet = sqrt(vartet) * sqrt(2); // 2-D
  cout << "Lynch-Dahl RMS scattering " << rmstet*1e3 << " mrad in 2D" << endl;
  cout << endl;

  double dz = thck / a0b;  // atomic units

  // INITIAL DIRECTION:

  double c0 = 1;                  // forward

  double s0 = sqrt(1-c0*c0);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile * histoFile = new TFile( "multics.root", "RECREATE" );

  // book histos:

  TH1I hwm( "wmin", "min single energy loss;min energy loss [eV];e", 200, 0, 20 );
  TH1I hde00( "de00", "single energy loss;energy loss [eV];inelastics", 200, 0, 20 );
  TH1I hde0( "de0", "single energy loss;energy loss [eV];inelastics", 200, 0, 200 );
  TH1I hde1( "de1", "single energy loss;energy loss [keV];inelastics", 200, 0, 2 );
  TH1I hde2( "de2", "single energy loss;energy loss [keV];inelastics", 200, 0, 200 );

  TH1I hloss( "loss", "total energy loss;energy loss [keV];tracks", 200, 0, 5*eloss );

  double trng = 10*rmstet*1e3;
  TH1I htet( "tet", "exit angle;exit angle [mrad];tracks", 200, 0, trng );
  TH1I htx( "tx", "exit angle x;exit angle x [mrad];tracks", 200, -trng, trng );
  TH1I hty( "ty", "exit angle y;exit angle y [mrad];tracks", 200, -trng, trng );

  double xrng = 5*thck*rmstet*1e4;
  TH2I * h2xy = new
    TH2I( "xy",";x [#mum];y [#mum]",
	  200, -xrng, xrng, 200, -xrng, xrng );

  TH2I * h2rz = new
    TH2I( "rz",";z [#mum];r [#mum]",
	  max(100, int(thck*1e3)), 0, thck*1e4, 200, 0, xrng );

  TProfile rvsz( "rvsz",";z [#mum];<r> [#mum]",
	  max(100, int(thck*1e3)), 0, thck*1e4 );

  xrng = 1*thck*rmstet*1e4;
  TH3I * h3xyz = new
    TH3I( "xyz",";x [#mum];y [#mum];z [#mum]",
	  200, -xrng, xrng, 200, -xrng, xrng, 200, 0, thck*1e4 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ranlux24 rgen;
  rgen.seed( time(NULL) ); // seconds since 1.1.1970
  uniform_real_distribution <double> unirnd( 0, 1 );

  int N[3];
  N[0] = 0;
  N[1] = 0;
  N[2] = 0;

  double EM = 0;
  double VEM = 0;
  double vang = 0;
  double vtet = 0;
  int ntet = 0;
  double disp = 0;
  double vdisp = 0;
  unsigned long mscat = 0;
  unsigned long vmscat = 0;
  unsigned long mloss = 0;
  unsigned long vmloss = 0;

  int ntot = 0;

  bool lout = 1;

  while( ntot < ntotal ) { // tracks

    ntot += 1;

    double x = 0;
    double y = 0;
    double z = 0;

    double TL = 0; // track length
    unsigned long nscat = 0;
    unsigned long nloss = 0;

    double W = c0;                    // direction cosines
    double phi = twopi*unirnd(rgen);
    double U = s0*cos(phi);
    double V = s0*sin(phi);

    double E = E0;                    // atomic units

    bool newE = 1;

    double EMFP = dz;
    double XMFP = dz;
    double TMFP = dz;
    double AFP = dz;
    double WM = fie;
    int iexit = 0;
    double EINV = 1/E;
    double RB = E + SL22;  // E+2mc2
    double QFACT = E;
    double QME = E;
    double SE0 = 0;
    double SE1 = 0;
    double SE2 = 0;
    double SET = 0;

    while(1) { // steps

      // energy changed, update:

      if( newE ) {

	EINV = 1/E;
	double FA = E + SL2; // E+mc2
	RB = E + SL22;  // E+2mc2
	double BETA2 = E*RB/( FA*FA ); // beta^2
	double SFACT = BETA2*SL2; // beta^2 mc^2
	double QMI = ( E+E)*RB*SL22I;
	QFACT = QMI + QMI;
	QME = QFACT + QFACT;

	// cout << "  E " <<  E*HRKEV << endl;

	// ELASTIC CROSS-SECTION ( MOTT'S FORMULA )

	FA = AL12 + QME;
	double FB = AL22 + QME;
	SE0 = log( FB/FA ) - AA4;
	SE1 = AA1 * ( QME / ( AL12*FA ) );
	SE2 = AA2 * ( QME / ( AL22*FB ) );
	SET = SE1 + SE2 + AA3*SE0;
	EMFP = ZZ4 * SET;

	// INELASTIC CROSS-SECTION ( LENZ'S FORMULA )
	// F. Lenz, Z. Naturf. 3A 1954 78

	double Q02 = AB1*EINV;
	FA = AL12 + Q02;
	FB = AL22 + Q02;
	double FT = AL12 + QMI;
	double FNT = AL22 + QMI;
	XMFP = Z4 * (
		     AB2*log( QMI*FA/( FT*Q02 ) ) +
		     AB3*log( QMI*FB/( FNT*Q02 ) ) +
		     AA1*( 1/FT-1/FA ) +
		     AA2*( 1/FNT-1/FB ) );

	TMFP = EMFP + XMFP; // SFACT = beta^2 mc^2
	AFP = SFACT / ( VMOL*TMFP ); // atomic units

	if( lout )
	  cout << "  path: " << endl
	       << "    elastic " << SFACT / ( VMOL*EMFP )*a0b*1e4 << " um" << endl
	       << "  inelastic " << SFACT / ( VMOL*XMFP )*a0b*1e4 << " um" << endl
	       << "  total     " << AFP*a0b*1e4 << " um" << endl;

	// AVERAGE ENERGY LOSS PER INELASTIC COLLISION:

	double WA = fie;
	double CSI = E*expot;
	if( CSI < 6.338065465611359 )
	  WA = ZT*EINV / XMFP * 3.1776907270853916 * sqrt( CSI );
	else {
	  FA = 1 - BETA2; // 1/gamma2
	  FB = sqrt( FA ); // 1/gamma
	  FT = 1 - FB;
	  WA = 2*ZT*EINV / XMFP *           // Bethe-Bloch, without 1/beta^2
	    ( log( E*expot2*SFACT/FA ) -
	      ( FB+FB+BETA2 ) * 6.931471805599453e-1 +
	      FA + 0.125*FT*FT );
	}
	if( lout )
	  cout << "  dEavg " << WA*E*HRKEV*1e3 << " eV/step" << endl;

	// delta rays sampled from 1/dE^2 (free Landau)
	// set relative Wmin such that the integral is WA
	// solve transcendental equation by Halley's method, see
	// F.B. Hildebrand: Introduction to Numerical Analysis, McGraw-Hill 1956 412

	WM = pow( WA, 1 + 3.3 / ( 3.3 - log( WA ) ) );
	FT = log( WM );
	if( lout )
	  cout << "  dEmin " << WM*E*HRKEV*1e3 << " eV" << endl;

	double W0 = WA*( 1-WM ) + WM*FT;
	double W1 = 1 - WA + FT;
	WM = WM - 2*W0*W1 / ( 2*W1*W1 - W0/WM );
	FT = log( WM );
	if( lout )
	  cout << "  dEmin " << WM*E*HRKEV*1e3 << " eV" << endl;

	W0 = WA*( 1 - WM ) + WM*FT;
	W1 = 1 - WA + FT;
	WM = WM - 2*W0*W1 / ( 2*W1*W1 - W0/WM );
	if( lout )
	  cout << "  dEmin " << WM*E*HRKEV*1e3 << " eV" << endl;

	hwm.Fill( WM*E*HRKEV*1e3 ); // [eV], around 5 eV for Si

	newE = 0;

	lout = 0;

      } // newE

      double rnd = unirnd(rgen); // uniform 0..1
      double s = -log( 1-rnd )*AFP;

      TL = TL + s;

      x = x + U*s;
      y = y + V*s;
      z = z + W*s;
      double r = sqrt( x*x + y*y );
      rvsz.Fill( z*a0b*1e4, r*a0b*1e4 );
      h2rz->Fill( z*a0b*1e4, r*a0b*1e4 );
      h2xy->Fill( x*a0b*1e4, y*a0b*1e4 );
      h3xyz->Fill( x*a0b*1e4, y*a0b*1e4, z*a0b*1e4 );

      // cout << "  z " << z*a0b*1e4; // [um]

      // ----------------------------------------------------- TEST BOUNDS

      if( z > dz ) {
	iexit = 0;
	break;
      }
      else if( z < 0 ) {
	iexit = 1;
	break;
      }

      // --------------------------------------------------- NEW COLLISION

      rnd = unirnd(rgen);

      if( rnd*TMFP < XMFP ) { // INELASTIC COLLISION

	double rnd = unirnd(rgen);
	double dE = 0.5*E*WM / ( 1 - rnd*( 1-WM ) ); // sample from 1/dE^2 above WM
	nloss += 1;
	hde00.Fill( dE*HRKEV*1e3 );
	hde0.Fill( dE*HRKEV*1e3 );
	hde1.Fill( dE*HRKEV );
	hde2.Fill( dE*HRKEV );

	E = E-dE;

	if( E <= eabs ) {
	  iexit = 2;
	  break;
	}

	newE = 1;

	// SCATTERING ANGLES

	double CDT = sqrt( ( 1 - dE*EINV ) * RB / ( RB-dE ) ); // RB = E+2mc2
	double DF = twopi*unirnd(rgen);
	DIRECT( CDT, DF, U, V, W );

	// cout << "  inelastic ", cdt, de

      } // inelastic

      else { // ELASTIC COLLISION see Baro et al 1994

	double rnd = unirnd(rgen);
	double FB = rnd*SET;
	double FNT = SE1;
	double fa = unirnd(rgen);
	double CDT = 1;

	if( FNT > FB )
	  CDT = fa * AL12 * QME / ( AL12 + ( 1-fa )*QME );

	else {
	  FNT += SE2;
	  if( FNT > FB )
	    CDT = fa * AL22 * QME / ( AL22 + ( 1-fa )*QME );

	  else {
	    double FT = exp( fa*SE0 + AA4 );
	    CDT = ( AL12*FT - AL22 ) / ( 1 - FT );
	  }
	}

	CDT = 1 - CDT / QFACT;
	double DF = twopi*unirnd(rgen);

	DIRECT( CDT, DF, U, V, W );
	nscat += 1;

	// cout << "    elastic ", cdt

      }

    } // do steps

    // ---------------------------------------------- INCREMENT COUNTERS

    N[iexit] += 1;

    if( iexit == 0 ) { // transmitted

      EM += E; // exit energy [atomic]
      VEM += E*E;

      hloss.Fill( e0kev - E*HRKEV );

      double fa = acos( W );           // exit angle, 0..pi
      vang += fa*fa;
      if( fa < 3*rmstet ) {
	vtet += fa*fa;
	ntet += 1;
      }
      htet.Fill( fa*1e3 );
      double phi = atan2( V, U );
      htx.Fill( fa*cos(phi)*1e3 );
      hty.Fill( fa*sin(phi)*1e3 );

      double dd = sqrt( x*x + y*y )*a0b*1e4; // displacement [um]
      disp += dd;
      vdisp += dd*dd;

      if( ntot < 100 || ntot%100 == 0 )
	cout << ntot << ": scat " << nscat << ", loss " << nloss
	     << ", angle " << fa*1e3 << " mrad, offset " << dd << " um"
	     << endl;

      mscat += nscat;
      vmscat += nscat*nscat;

      mloss += nloss;
      vmloss += nloss*nloss;

    } // exit

  } // tracks

  cout << endl;

  cout << "element    " << Z << endl;
  cout << "thickness  " << dz*a0b*1e4 << " um" << endl;
  cout << "energy     " << e0mev << " MeV" << endl;
  cout << endl;
  cout << "tracks         " << ntot << endl;
  cout << "transmitted    " << N[0] << endl;
  cout << "backscattered  " << N[1] << endl;
  cout << "absorbed       " << N[2] << endl;

  if( N[0] > 0 ) {

    cout << endl << "transmitted:" << endl;

    double n0 = 1.0 / N[0];

    double uscat = n0 * sqrt( vmscat - n0 * mscat*mscat ); // uncertainty
    cout << "  elastics  " << n0*mscat << " +- " << uscat << endl;

    double uloss = n0 * sqrt( vmloss - n0 * mloss*mloss );
    cout << "inelastics  " << n0*mloss << " +- " << uloss << endl;

    double uem = n0 * sqrt( VEM - n0 * EM*EM ) * HRKEV;
    EM = n0 * EM * HRKEV;
    cout << endl << "energy  " << EM*1e-3 << " +- " << uem*1e-3 << " MeV"
	 << endl << "loss    " << e0kev-EM << " +- " << uem << " keV"
	 << endl;

    vang = sqrt( vang * n0 ); // rms scat ang
    cout << endl << "RMS angle  " << vang*1e3 << " mrad" << endl;
    if( ntet )
      cout << "limited    " << sqrt( vtet / ntet )*1E3 << " mrad" << endl;

    double udisp = n0 * sqrt( vdisp - n0 * disp*disp );
    cout << endl << "displacement  " << n0*disp << " +- " << udisp << " um" << endl;

  } // trans

  cout << endl
       << histoFile->GetName() << endl
       << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;

} // main
