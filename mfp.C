
// from multics.cc
// mean free path of electrons vs Ekin

// root -l mfp.C

// produces mfp.root

#include <iostream> // cout
#include <cmath> // log

#include <TFile.h>
#include "TProfile.h"

using namespace std;

//------------------------------------------------------------------------------
void mfp()
{
  // physical constants:

  double AVOG = 6.022094e23; // [atoms/Mol]
  double a0b = 5.291771e-9; // Bohr radius [cm]
  double HRKEV  = 2.721160e-2; // Hartree aomic units
  double mec2 = 511; // [keV]

  // POWERS OF THE SPEED OF LIGHT in atomic units:

  double SL2 = 1.877887e4; // c^2
  double SL22 = 3.755773e4; // 2c^2
  double SL22I = 2.662568e-5; // 1/2c^2

  double pi = 4*atan(1);
  double twopi = 2*pi;
  double log10 = log(10);

  // ELEMENT DATA:

  double Z = 14;      // Si
  double A = 28.0855; // g/mol
  double rho = 2.329;            // g/cm3

  cout << "  ATOMIC NUMBER " << Z << endl
       << "  ATOMIC WEIGHT " << A << endl
       << "  DENSITY       " << rho << " g/cm3"
       << endl;

  // SCREENING PARAMETERS:

  double AA = 0.51591; // Si free
  double AL1 = 5.84957;
  double AL2 = 1.17330;

  //double AA = 0.31186; // Si Wigner-Seitz
  //double AL1 = 7.76730;
  //double AL2 = 1.71824 ;

  cout << "  SCREENING CONSTANTS:" << endl
       << "  A      = " << A << endl
       << "  alpha1 = " << AL1 << endl
       << "  alpha2 = " << AL2 << endl
    ;

  // CROSS-SECTION PARAMETERS:

  double Z2 = twopi * Z;
  double Z4 = 2*twopi * Z;
  double ZZ4 = Z*Z4; // 4pi Z*Z

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

  double IkeV = 1.35e-2 * Z;
  if( Z > 11.5 )
    IkeV = Z * ( 9.76e-3 + 5.88e-2 / pow( Z, 1.19 ) );

  cout << "MEAN EXCITATION POTENTIAL " << IkeV*1e3 << " eV" << endl;

  //dE/dx = n 2π r_0^2 m/b^2
  // n = electron density

  double eden = AVOG*rho/A*Z; // egs5 p50
  double re = 2.818e-13; // [cm] e^2/(4 pi eps0 me c^2)

  double expot = IkeV/HRKEV;
  double expotinv = 1/expot;
  double expotinv2 = expotinv*expotinv;

  // FIRST IONIZATION ENERGY:

  double fie = 6; // eV
  //double fie = 12; // Si M-shell
  //double fie = 1.2; // Si band gap
  cout << "FIRST IONIZATION ENERGY = " << fie << " eV" << endl;

  fie = fie*1e-3 / HRKEV;      // Hartree
  double AB1 = 0.125*fie*fie;

  // NUMBER OF 'MOLECULES' PER UNIT VOLUME
  double VMOL = AVOG*rho/A * pow( a0b, 3 ); // a0b = Bohr

  TFile * histoFile = new TFile( "mfp.root", "RECREATE" );

  TProfile gamvse( "gamvse",
		   "gamma;log_{10}(Ekin [keV]);gamma",
		   6*25+1, -0.02, 6.02 ); // 25 steps/decade

  TProfile betavse( "betavse",
		    "beta = p/E;log_{10}(Ekin [keV]);beta = p/E",
		    6*25+1, -0.02, 6.02 ); // 25 steps/decade

  TProfile se0vse( "se0vse",
		   "elastic 0;log_{10}(Ekin [keV]);elastic 0 / elastic total",
		    6*25+1, -0.02, 6.02 ); // 25 steps/decade
  TProfile se1vse( "se1vse",
		   "elastic 1;log_{10}(Ekin [keV]);elastic 1 / elastic total",
		    6*25+1, -0.02, 6.02 ); // 25 steps/decade
  TProfile se2vse( "se2vse",
		   "elastic 2;log_{10}(Ekin [keV]);elastic 2 / elastic total",
		    6*25+1, -0.02, 6.02 ); // 25 steps/decade
  TProfile emfpvse( "emfpvse",
		    "elastic mean free path;log_{10}(Ekin [keV]);elastic mean feee path [#mum]",
		    6*25+1, -0.02, 6.02 ); // 25 steps/decade
  TProfile xmfpvse( "xmfpvse",
		    "inelastic mean free path;log_{10}(Ekin [keV]);inelastic mean feee path [#mum]",
		    6*25+1, -0.02, 6.02 ); // 25 steps/decade
  TProfile wavse( "wavse",
		  "average energy loss / step;log_{10}(Ekin [keV]);average energy loss per step [keV]",
		  6*25+1, -0.02, 6.02 ); // 25 steps/decade
  TProfile wmvse( "wmvse",
		  "minimum energy loss / step;log_{10}(Ekin [keV]);minimum energy loss per step [keV]",
		  6*25+1, -0.02, 6.02 ); // 25 steps/decade

  // incident energy:

  double e0kev = 1;            // [keV]

  while( e0kev < 1.001e6 ) {

    double Ekin = e0kev/HRKEV;             // atomic units
    double E = Ekin + SL2; // Ekin + mc^2
    double gamma = E/SL2;
    gamvse.Fill( log(e0kev)/log10, gamma );

    double RB = Ekin + SL22;  // E + 2mc^2
    double beta2 = Ekin*RB/( E*E ); // beta^2
    betavse.Fill( log(e0kev)/log10, sqrt(beta2) );
    double SFACT = beta2 * SL2; // beta^2 mc^2
    double bg2 = beta2*gamma*gamma; // bg = p/m
    double QMI = 2*Ekin*RB*SL22I; // for inelastic
    double QME = 4*QMI; // for elastic

    // ELASTIC CROSS-SECTION ( MOTT'S FORMULA )

    double FA = AL12 + QME;
    double FB = AL22 + QME;
    double SE0 = log( FB/FA ) - AA4;
    double SE1 = AA1 * ( QME / ( AL12*FA ) );
    double SE2 = AA2 * ( QME / ( AL22*FB ) );
    double SET = SE1 + SE2 + AA3*SE0;
    se0vse.Fill( log(e0kev)/log10, AA3*SE0/SET ); // 0.215
    se1vse.Fill( log(e0kev)/log10, SE1/SET ); // 0.034
    se2vse.Fill( log(e0kev)/log10, SE2/SET ); // 0.751
    double EMFP = ZZ4 * SET;
    double emfp = SFACT / ( VMOL*EMFP )*a0b*1e4; // [um]
    emfpvse.Fill( log(e0kev)/log10, emfp );

    // INELASTIC CROSS-SECTION ( LENZ'S FORMULA )

    double Einv = 1/Ekin;
    double Q02 = AB1*Einv;
    FA = AL12 + Q02;
    FB = AL22 + Q02;
    double FT = AL12 + QMI;
    double FNT = AL22 + QMI;
    double XMFP = Z4 * ( AB2*log( QMI*FA/( FT*Q02 ) ) +
			 AB3*log( QMI*FB/( FNT*Q02 ) ) +
			 AA1*( 1/FT-1/FA ) +
			 AA2*( 1/FNT-1/FB ) );
    double xmfp = SFACT / ( VMOL*XMFP )*a0b*1e4; // [um]
    xmfpvse.Fill( log(e0kev)/log10, xmfp );

    double TMFP = EMFP + XMFP;
    double tmfp = SFACT / ( VMOL*TMFP )*a0b*1e4; // [um]

    cout << e0kev << " keV: " << emfp << "  " << xmfp << "  " << tmfp << " um" << endl;

    // AVERAGE relative ENERGY LOSS PER INELASTIC COLLISION: Bethe-Bloch

    double Tmax = SL2*(gamma-1); // Kolanoski, Wermes (3.19)
    double Tmin = expot*expot/(SL2*bg2);
    double dedx = eden*twopi*re*re*mec2/beta2;// egs5 p74
    double WA = fie;
    double CSI = Ekin*expotinv;
    if( CSI < 6.338065465611359 ) // at 1 keV
      WA = Z2 * Einv / XMFP * 3.1776907270853916 * sqrt( CSI );
    else { // above few keV
      FA = 1 - beta2; // 1/gamma2
      FB = sqrt( FA ); // 1/gamma
      FT = 1 - FB;
      WA = 2*Z2*Einv / XMFP * // relative Bethe-Bloch without 1/beta2
	( log( Ekin*expotinv2*SFACT/FA ) -
	  ( FB + FB + beta2 ) * 6.931471805599453e-1 +
	  FA + 0.125*FT*FT );
    }

    wavse.Fill( log(e0kev)/log10, WA*Ekin*HRKEV/beta2 ); // Bethe-Bloch [keV/?]

    // MINIMUM ENERGY LOSS. HALLEY'S METHOD:

    double WM = pow( WA, 1 + 3.3 / ( 3.3 - log( WA ) ) );
    FT = log( WM );

    double W0 = WA*( 1-WM ) + WM*FT;
    double W1 = 1 - WA + FT;
    WM = WM - 2*W0*W1 / ( 2*W1*W1 - W0/WM );
    FT = log( WM );

    W0 = WA*( 1 - WM ) + WM*FT;
    W1 = 1 - WA + FT;
    WM = WM - 2*W0*W1 / ( 2*W1*W1 - W0/WM );

    wmvse.Fill( log(e0kev)/log10, WM*Ekin*HRKEV );

    cout << "  Ekin " <<  Ekin*HRKEV
	 << ", elas " << SFACT / ( VMOL*EMFP )*a0b*1e4
	 << ", inel " << SFACT / ( VMOL*XMFP )*a0b*1e4 << " um"
	 << ", csi " << CSI
	 << ", wavg " << WA*Ekin*HRKEV << " keV"
	 << " = " << WA*Ekin*HRKEV / ( SFACT / ( VMOL*XMFP )*a0b*1e4 ) << " keV/um"
	 << ", wmin " << WM*Ekin*HRKEV << " keV"
	 << endl;

    e0kev *= pow(10,0.04);

  } // while


  histoFile->Write();
  histoFile->ls();
  cout << endl << histoFile->GetName() << endl << endl;

}
