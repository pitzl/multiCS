
// from multics.cc
// elastic scattering distribution

// root -l scat.C

#include <iostream> // cout
#include <random>
#include <cmath> // log

#include <TFile.h>
#include "TH1I.h"

using namespace std;

//------------------------------------------------------------------------------
void scat()
{
  int nev = 10*1000;

  // incident energy:

  double Ek = 100; // [MeV]

  double me = 0.511; // [MeV]
  double E = Ek+me;
  double p = sqrt( E*E - me*me );

  cout
    << endl << "Ekin  " << Ek << " MeV"
    << endl << "E     " << E << " MeV"
    << endl << "p     " << p << " MeV/c"
    << endl << "beta  " << p/E
    << endl << "gamma " << E/me
    << endl;

    // Hartree aomic units:

  double HRKEV = 2.721160e-2;

  double SL2 = 1.877887e4; // c^2
  double SL22 = 3.755773e4; // 2c^2
  double SL22I = 2.662568e-5; // 1/2c^2

  double Ekin = Ek*1e3/HRKEV;             // atomic units
  double RB = Ekin + SL22;  // E + 2mc^2
  double QMI = 2*Ekin*RB*SL22I;
  double QFACT = 2*QMI;
  double QME = 4*QMI;

  cout
    << endl << "QME  " << QME
    << endl;

  // ELEMENT DATA:

  // SCREENING PARAMETERS:

  double AA = 0.31186; // Si Wigner-Seitz
  double AL1 = 7.76730;
  double AL2 = 1.71824 ;

  double Z = 14;
  double gn = 1e-6 * 2 * 2.61 * pow( Z, 2.0/3.0 ) / p/p; // Mazziotta, corrected, see Yamazaki 1974
  double tetmin = 2.66 * pow( Z, 1.0/3.0 ) / p; // [mrad] Fruehwirth, Regler, from Jackson

  cout
    << endl << "SCREENING:"
    << endl << "  AA     " << AA
    << endl << "  alpha1 " << AL1
    << endl << "  alpha2 " << AL2
    << endl << "  gn     " << gn
    << endl << "  tetmin " << tetmin << " mrad"
    << endl;

  // CROSS-SECTION PARAMETERS:

  double A1A = 1 - AA;
  double AA1 = AA*AA;
  double AL12 = AL1*AL1;
  double AL22 = AL2*AL2;
  double AA2 = A1A*A1A;
  double AA3 = 2*AA*A1A/( AL12-AL22 );
  double AA4 = log( AL22/AL12 );

  // ELASTIC CROSS-SECTION ( MOTT'S FORMULA )

  double FA = AL12 + QME;
  double FB = AL22 + QME;
  double SE0 = log( FB/FA ) - AA4;
  double SE1 = AA1 * ( QME / ( AL12*FA ) );
  double SE2 = AA2 * ( QME / ( AL22*FB ) );
  double SET = SE1 + SE2 + AA3*SE0;

  cout
    << endl << "elastic cross sections"
    << endl << "  0 " << AA3*SE0/SET
    << endl << "  1 " << SE1/SET
    << endl << "  2 " << SE2/SET
    << endl;

  TFile * histoFile = new TFile( "scat.root", "RECREATE" );

  TH1I hcdt( "cdt", "cos tet;cos(theta);tracks", 100, 0.99, 1 );
  TH1I htet( "tet", "theta Salvat;theta [mrad];tracks", 100, 0, 15*tetmin );
  TH1I hchi( "chi", "theta Chaoui;theta [mrad];tracks", 100, 0, 15*tetmin );
  TH1I hrnd( "rnd", "randoms;unirnd;draws", 100, 0, 1 );

  ranlux24 rgen;
  rgen.seed( time(NULL) ); // seconds since 1.1.1970
  uniform_real_distribution <double> unirnd( 0, 1 );

  for( int iev = 0; iev < nev; ++iev ) {

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

    hcdt.Fill( CDT );
    double tet = acos(CDT);
    htet.Fill( tet*1e3 );

    //double cost = 1 - 2*gn*rnd / ( 2 + gn - 2*rnd ); // Mazziotta, Chaoui
    //double cost = 1 - gn*rnd / (1-rnd + 0.5*gn ); // 1-rnd = rnd'
    double cost = 1 - (1-rnd) / ( 0.5 + rnd/gn ); // Liljequist, Salvat 1989
    double chi = acos(cost);
    hchi.Fill( chi*1e3 );
    hrnd.Fill(rnd); // flat

  } // iev

  cout << endl;
  histoFile->Write();
  histoFile->ls();
  cout << endl << histoFile->GetName() << endl << endl;

}

