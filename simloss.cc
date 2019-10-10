
// Daniel Pitzl, Sep 2019
// energy loss a la Landau and Salvat

// make simloss
// simloss -e 5000 -t 150 -n 10100
// creates simloss.root

#include <cstdlib> // atoi
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
int main( int argc, char* argv[] )
{
  // USER steering:

  int ntr = 10*1000;

  double thck = 150e-4; // [cm]

  //double step = 0.06435e-4; // [cm] inelastic MFP
  double step = 0.1e-4; // [cm]

  double e0mev = 1;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-n" ) )
      ntr = atoi( argv[++i] );

    if( !strcmp( argv[i], "-t" ) )
      thck = atof( argv[++i] )*1e-4; // [um]

    if( !strcmp( argv[i], "-e" ) )
      e0mev = atof( argv[++i] ); // [MeV]

  } // argc

  // ELEMENT DATA:

  double Z = 14;       // Si
  double A = 28.0855;  // g/mol
  double rho = 2.329;  // g/cm3

  cout << endl
    << "  ATOMIC NUMBER " << Z << endl
    << "  ATOMIC WEIGHT " << A << endl
    << "  DENSITY       " << rho << " g/cm3" << endl;

  cout << "  thickness  " << thck*1e4 << " um" << endl;
  cout << "  step size  " << step*1e4 << " um" << endl;
  cout << "  steps      " << thck/step << endl;
  cout << "  energy     " << e0mev << " MeV" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  double AVOG = 6.022094e23; // [atoms/Mol]
  double a0b = 5.291771e-9; // Bohr radius [cm]
  double HRKEV  = 2.721160e-2; // Hartree aomic units

  // atomic units: hbar = m = e = 1

  // double SL = 137.036; // ATOMIC UNITS
  double SL2 = 1.877887e4; // c^2
  double SL22 = 3.755773e4; // 2c^2

  // MEAN EXCITATION POTENTIAL ( SEMIEMPIRICAL FORMULA )

  double expot = 1.35e-2 * Z; // [keV]
  if( Z > 11.5 )
    expot = Z * ( 9.76e-3 + 5.88e-2 / pow( Z, 1.19 ) );

  cout << "  MEAN EXCITATION POTENTIAL " << expot*1e3 << " eV" << endl;

  expot = expot/HRKEV;
  expot = 1/expot;
  double expot2 = expot*expot;

  double e0kev = e0mev*1e3;
  double E0 = e0kev/HRKEV;             // atomic units

  double VMOL = AVOG*rho/A * pow( a0b, 3 ); // atoms / volume

  double eloss = 0.3 * thck*1e4;     // MIP [keV]

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile * histoFile = new TFile( "simloss.root", "RECREATE" );

  // book histos:

  TH1I hde00( "de00", "single energy loss;energy loss [eV];inelastics", 200, 0, 20 );
  TH1I hde0( "de0", "single energy loss;energy loss [eV];inelastics", 200, 0, 200 );
  TH1I hde1( "de1", "single energy loss;energy loss [keV];inelastics", 200, 0, 2 );
  TH1I hde2( "de2", "single energy loss;energy loss [keV];inelastics", 200, 0, 200 );

  TH1I htde( "tde", "total energy loss;energy loss [keV];tracks", 200, 0, 5*eloss );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // randoms:

  ranlux24 rgen;
  rgen.seed( time(NULL) ); // seconds since 1.1.1970
  uniform_real_distribution <double> unirnd( 0, 1 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  double E = E0;                    // atomic units

  double EINV = 1/E;
  double FA = E + SL2; // E+mc2
  double RB = E + SL22;  // E+2mc2
  double BETA2 = E*RB/( FA*FA ); // beta^2
  double SFACT = BETA2 * SL2;

  // AVERAGE ENERGY LOSS PER INELASTIC COLLISION:

  double AFP = step/a0b; // free path, atomic units
  double XMFP = SFACT / ( VMOL*AFP ); // cross section

  double twopi = 8*atan(1);
  double ZT = twopi * Z;
  FA = 1 - BETA2; // 1/gamma2
  double FB = sqrt( FA ); // 1/gamma
  double FT = 1 - FB;

  double WA = 2*ZT * EINV / XMFP *           // Bethe-Bloch per step
    ( log( E * expot2 * SFACT / FA ) -
      ( FB+FB+BETA2 ) * 6.931471805599453e-1 +
      FA + 0.125*FT*FT );

  cout << "  dEavg     " << WA*E*HRKEV*1e3 << " eV/step" << endl;

  // delta rays sampled from 1/dE^2 (free Landau)
  // set relative Wmin such that the integral is WA
  // solve transcendental equation by Halley's method, see
  // F.B. Hildebrand: Introduction to Numerical Analysis, McGraw-Hill 1956 412

  double WM = pow( WA, 1 + 3.3 / ( 3.3 - log( WA ) ) );
  FT = log( WM );
   cout << "  dEmin    " << WM*E*HRKEV*1e3 << " eV" << endl;

  double W0 = WA*( 1-WM ) + WM*FT;
  double W1 = 1 - WA + FT;
  WM = WM - 2*W0*W1 / ( 2*W1*W1 - W0/WM );
  FT = log( WM );
  cout << "  dEmin     " << WM*E*HRKEV*1e3 << " eV" << endl;

  W0 = WA*( 1 - WM ) + WM*FT;
  W1 = 1 - WA + FT;
  WM = WM - 2*W0*W1 / ( 2*W1*W1 - W0/WM );
  cout << "  dEmin     " << WM*E*HRKEV*1e3 << " eV" << endl;

  double allde = 0;
  double varde = 0;

  for( int itr = 0; itr < ntr; ++itr ) {

    if( itr%1000 == 0 ) cout << "  " << itr << flush;

    double sumde = 0;

    for( double z = 0; z < thck; z += step ) {

      double rnd = unirnd(rgen);
      double dE = 0.5*E*WM / ( 1 - rnd*( 1-WM ) ); // sample from 1/dE^2 above WM
      sumde += dE;
      hde00.Fill( dE*HRKEV*1e3 );
      hde0.Fill( dE*HRKEV*1e3 );
      hde1.Fill( dE*HRKEV );
      hde2.Fill( dE*HRKEV );

    } // steps

    htde.Fill( sumde*HRKEV );
    allde += sumde;
    varde += sumde*sumde;

  } // tracks

  cout << endl;

  double avgde = allde/ntr;
  double rmsde = sqrt( varde - avgde*avgde*ntr ) / ntr;
  cout << endl;
  cout << "  tracks     " << ntr << endl;
  cout << "  mean loss  " << avgde*HRKEV << " +- " << rmsde*HRKEV << " keV" << endl;
  cout << "  restricted " << htde.GetMean() << " keV" << endl;
  cout << "  Bethe      " << WA*E0*HRKEV/rho*thck/step << " keV" << endl;

  // PDG:
  //double K = 0.307; // MeV cm2 / Mol
  //double dEdx = 2*K*rho*Z/A/BETA2 * ( log( e0kev * expot2 * SFACT / FA ) -2*BEA ); // Wmax=E0/2

  cout << endl
       << histoFile->GetName() << endl
       << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;

} // main
