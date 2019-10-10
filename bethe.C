
// Daniel Pitzl, DESY, Oct 2019
// plot various forms of Bethe-Bloch

// root -l bethe.C

#include <iostream> // cout
#include <cmath> // log

#include <TFile.h>
#include "TProfile.h"

using namespace std;

//------------------------------------------------------------------------------
void bethe()
{
  // ELEMENT DATA:

  double Z = 14;      // Si
  double A = 28.0855; // g/mol
  double rho = 2.329; // g/cm3

  cout << endl << "element: Z " << Z << ", A " << A
       << ", rho " << rho << " g/cm3" << endl;

  // MEAN EXCITATION POTENTIAL ( SEMIEMPIRICAL FORMULA ):

  double IkeV = 1.35e-2 * Z;
  if( Z > 11.5 )
    IkeV = Z * ( 9.76e-3 + 5.88e-2 / pow( Z, 1.19 ) );

  cout << endl
       << "mean excitation potential " << IkeV*1e3 << " eV" << endl;

  double IMeV = IkeV*1e-3; // [MeV]
  double I = IkeV*1e3; // [eV]

  double plasma = 31.05; // [eV] Kolanoski, Wermes
  double CD = 2*log(plasma/I) - 1;
  cout << "plasma CD " << CD << endl;
  double a = 0.149;
  double k = 3.25;
  double d0 = 0.14;
  double zeta0 = 0.202;
  double zeta1 = 2.872;

  double log10 = log(10);
  double pi = 4*atan(1);
  double twopi = 2*pi;

  double NA = 6.022094e23; // [atoms/mol] Avogadro
  double re = 2.818e-13; // [cm] e2 / ( 4 pi eps0 me c2 )

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TFile * histoFile = new TFile( "bethe.root", "RECREATE" );

  TProfile gamvsbg( "gamvsbg",
		    "gamma;log_{10}(#beta#gamma);gamma",
		    4*25+1, -1.02, 3.02 ); // 25 steps/decade
  TProfile betvsbg( "betvsbg",
		    "beta;log_{10}(#beta#gamma);beta",
		    4*25+1, -1.02, 3.02 ); // 25 steps/decade
  TProfile Ekinvsbg( "Ekinvsbg",
		     "Ekin;log_{10}(#beta#gamma);Ekin [MeV]",
		     4*25+1, -1.02, 3.02 ); // 25 steps/decade
  TProfile tmaxvsbg( "tmaxvsbg",
		     "maximum energy loss;log_{10}(#beta#gamma);maximum energy loss [MeV]",
		     4*25+1, -1.02, 3.02 ); // 25 steps/decade
  TProfile deltavsbg( "deltavsbg",
		      "delta correction;log_{10}(#beta#gamma);delta correction",
		      4*25+1, -1.02, 3.02 ); // 25 steps/decade
  TProfile shellvsbg( "shellvsbg",
		      "shell correction;log_{10}(#beta#gamma);shell correction",
		      4*25+1, -1.02, 3.02 ); // 25 steps/decade
  TProfile dEdxvsbg( "dEdxvsbg",
		     "mean energy loss;log_{10}(#beta#gamma);mean energy loss [keV/#mum]",
		     4*25+1, -1.02, 3.02 ); // 25 steps/decade
  TProfile dEdxrvsbg( "dEdxrvsbg",
		     "restricted energy loss;log_{10}(#beta#gamma);restricted energy loss [keV/#mum]",
		     4*25+1, -1.02, 3.02 ); // 25 steps/decade
  TProfile dEdxevsbg( "dEdxevsbg",
		     "electron energy loss;log_{10}(#beta#gamma);electron energy loss [keV/#mum]",
		     4*25+1, -1.02, 3.02 ); // 25 steps/decade

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  double M = 105.6583; // [MeV] muon
  double m =   0.511; // [MeV] e
  double z = 1;

  double K0 = twopi*re*re*m*NA;
  double K = K0 * Z*rho/A * z*z; // constant in Bethe-Bloch
  cout << endl
       << "K_PDG " << 2*K0 << " MeV cm2 / mol" << endl; // 0.307

  double dEmin = 9e9;
  double bgmin = 1;
  double dEdxr9 = 0;

  // beta*gamma = p/M

  for( double bg = 0.1; bg < 1001; bg *= pow( 10, 0.04 ) ) {

    double bg2 = bg*bg;
    double g2 = bg2 + 1;
    double gam = sqrt(g2); // E/m
    double bet = bg/gam;
    double b2 = bet*bet;

    double p = bg*M;
    double E = gam*M;
    double Ekin = E-M;

    double lnbg = log(bg);
    double logbg = lnbg/log10; // zeta in Kolanoski, Wermes

    gamvsbg.Fill( logbg, gam );
    betvsbg.Fill( logbg, bet );
    Ekinvsbg.Fill( logbg, Ekin );

    double Tmax = 2*m*bg2 / ( 1 + ( 2*gam + m/M ) * m/M ); // [MeV]

    double Tcut = Tmax;

    double Tmin = IMeV*IMeV / ( 2 * m * bg2 ); // Kolanoski-Wermes 3.20

    double delta = 2*lnbg + CD;
    if( logbg < zeta0 )
      delta = d0*pow( 10, 2*(logbg-zeta0) );
    else if( logbg < zeta1 )
      delta = 2*lnbg + CD + a * pow( zeta1 - logbg, k );

    double shell = 0;

    double dEdx = // [MeV/cm]
      K / b2 * ( log( Tcut/Tmin ) - b2 * ( 1 + Tcut/Tmax ) - delta - 2*shell / Z );

    Tcut = 1.3; // [MeV]
    if( Tcut > Tmax )
      Tcut = Tmax;

    double dEdxr = // restricted [MeV/cm]
      K / b2 * ( log( Tcut/Tmin ) - b2 * ( 1 + Tcut/Tmax ) - delta - 2*shell / Z );

    // electrons: Knoll 4th (2.10)

    Tmax = 0.5*Ekin; // identical particles
    Tcut = 1.3; // [MeV]
    if( Tmax > Tcut )
      Tmax = Tcut;
    Tmin = IMeV*IMeV / bg2/m; // Kolanoski-Wermes 3.20

    double dEdxe = // restricted [MeV/cm], small effect a low Ekin
      K / b2 * ( log( Tmax/Tmin ) - log(2) * ( 2 * sqrt(1-b2) + b2-1 ) +
		 1-b2 + 0.125*pow( 1 - sqrt( 1-b2 ), 2 ) - delta - 2*shell / Z );

    tmaxvsbg.Fill( logbg, Tmax );
    deltavsbg.Fill( logbg, delta );
    shellvsbg.Fill( logbg, shell );
    dEdxvsbg.Fill( logbg, dEdx*1e2 );
    dEdxrvsbg.Fill( logbg, dEdxr*1e2 );
    dEdxevsbg.Fill( logbg, dEdxe*1e2 );

    if( dEdx < dEmin ) {
      dEmin = dEdx;
      bgmin = bg;
    }
    dEdxr9 = dEdxr; // remember

  } // bg

  cout << endl
       << "min dEdx " << dEmin*1e2 << " eV/um"
       << " at bg " << bgmin
       << endl;

  cout << "restricted dEdx " << dEdxr9*1e2 << " eV/um at last point" << endl;

  histoFile->Write();
  cout << endl;
  histoFile->ls();
  cout << endl << histoFile->GetName() << endl << endl;

}
