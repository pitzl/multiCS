
// from multics.cc and ionizer.cc
// elastic scattering distribution

// root -l elas.C

#include <iostream> // cout
#include <random>
#include <cmath> // log

#include <TFile.h>
#include "TH1I.h"

using namespace std;

//------------------------------------------------------------------------------
void elas()
{
  int nev = 10*1000;

  // incident energy:

  //double Ekin = 100; // [MeV]
  //double Ekin = 0.62; // for p = 1
  //double Ekin = 0.3;
  double Ekin = 0.1;

  double me = 0.511; // [MeV/c2]
  double E = Ekin + me;
  double p = sqrt( E*E - me*me ); // [MeV/c]
  double beta = p/E;
  cout
    << endl << "Ekin  " << Ekin << " MeV"
    << endl << "E     " << E << " MeV"
    << endl << "p     " << p << " MeV/c"
    << endl << "beta  " << beta
    << endl << "gamma " << E/me
    << endl;

    // Hartree aomic units:

  double HRKEV = 2.721160e-2;

  double SL2 = 1.877887e4; // c^2
  double SL22 = 3.755773e4; // 2c^2
  double SL22I = 2.662568e-5; // 1/2c^2

  double Ek = Ekin*1e3/HRKEV;             // atomic units
  double RB = Ek + SL22;  // E + 2mc^2
  double QMI = 2*Ek*RB*SL22I;
  double QFACT = 2*QMI;
  double QME = 4*QMI;

  cout
    << endl << "QME  " << QME
    << endl;

  // Salvat SCREENING PARAMETERS:
  /*
  double AA = 0.51591; // Si free
  double AL1 = 5.84957;
  double AL2 = 1.17330;
  */
  double AA = 0.31186; // Si Wigner-Seitz
  double AL1 = 7.76730;
  double AL2 = 1.71824 ;

  double Z = 14;

  double gm = 2 * 2.61 * pow( Z, 2.0/3.0 ) / (Ekin*1e6); // Mazziotta
  double gn = 2 * 2.61 * pow( Z, 2.0/3.0 ) / (p*p) *1e-6; // Mazziotta, corrected, p [MeV/c]
  // see Yamazaki 1974

  double tetmin = 2.66 * pow( Z, 1.0/3.0 ) / p; // [mrad] Fruehwirth, Regler, from Jackson

  // Boschini et al, 2014
  // cites
  // J.M. Fernandez-Vera et al., Nucl. Instr. and Meth. in Phys. Res. B 73 (1993) 447–473
  // cites
  // Berger, Wang 1988 in Jenkins
  // cites
  // Moliere 1947
  // who did the real work

  double pi = 4*atan(1);
  double CTF = 0.5*pow(0.75*pi,2.0/3.0); // 0.88534
  double a0 = 0.53e-8; // [cm] Bohr radius
  double aTF = CTF * a0 / pow( Z, 1.0/3.0 );
  // double AsM = pow( 0.5*hbar/aTF/p, 2 ); // single scattering Moliere
  double hbar = 197*1e-13; // MeV cm
  double alfa = 1.0/137;
  double AsM =
    0.25 * pow( hbar/p, 2 ) * pow( CTF*a0, -2 ) * pow( Z, 2.0/3.0 ) * ( 1.13 + 3.76 * pow( alfa*Z/beta, 2 ) );
  // 4.4e-6 * pow( Z, 2.0/3.0 )/p/p

  cout
    << endl << "SCREENING:"
    << endl << "  AA     " << AA
    << endl << "  alpha1 " << AL1
    << endl << "  alpha2 " << AL2
    << endl << "  gm     " << 2 * 2.61 * pow( Z, 2.0/3.0 ) / (Ekin*1e6) // Mazziotta
    << endl << "  gn     " << gn // corrected
    << endl << "  tetmin " << tetmin << " mrad"
    << endl << "  AsM    " << AsM
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
    << endl << "elastic cross section fractions:"
    << endl << "  0 " << AA3*SE0/SET
    << endl << "  1 " << SE1/SET
    << endl << "  2 " << SE2/SET
    << endl;

  TFile * histoFile = new TFile( "elas.root", "RECREATE" );

  TH1I hrnd( "rnd", "randoms;unirnd;draws", 100, 0, 1 );
  TH1I htet( "tet", "theta Chaoui;theta [mrad];tracks", 100, 0, 15*tetmin );
  TH1I hpsi( "psi", "theta Pitzl;theta [mrad];tracks", 100, 0, 15*tetmin );
  TH1I hdlt( "dlt", "theta Bethe;theta [mrad];tracks", 100, 0, 15*tetmin );
  TH1I hcdt( "cdt", "cos tet;cos(theta);tracks", 100, 0.99, 1 );
  TH1I hchi( "chi", "theta Salvat;theta [mrad];tracks", 100, 0, 15*tetmin );

  ranlux24 rgen;
  rgen.seed( time(NULL) ); // seconds since 1.1.1970
  uniform_real_distribution <double> unirnd( 0, 1 );

  for( int iev = 0; iev < nev; ++iev ) {

    double rnd = unirnd(rgen);
    hrnd.Fill(rnd); // flat

    //double cost = 1 - 2*gm*rnd / ( 2 + gm - 2*rnd ); // Mazziotta, Chaoui
    //double cost = 1 - gm*rnd / (1-rnd + 0.5*gm ); // 1-rnd = rnd'
    double cost = 1 - (1-rnd) / ( 0.5 + rnd/gm ); // Liljequist, Salvat 1989
    double tet = acos(cost);
    htet.Fill( tet*1e3 );

    double cosp = 1 - (1-rnd) / ( 0.5 + rnd/gn ); // Liljequist, Salvat 1989
    double psi = acos(cosp);
    hpsi.Fill( psi*1e3 );

    double cosd = 1 - (1-rnd) / ( 0.5 + rnd/AsM ); // Moliere
    double dlt = acos(cosd);
    hdlt.Fill( dlt*1e3 );

    // Salvat:

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
    double chi = acos(CDT);
    hchi.Fill( chi*1e3 );

  } // iev

  cout << endl;
  histoFile->Write();
  histoFile->ls();
  cout << endl << histoFile->GetName() << endl << endl;

}

