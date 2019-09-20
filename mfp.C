
// from multics.C
// mean free path of electrons vs Ekin

{
  // physical constants:

  double AVOG = 6.022094e23; // [atoms/Mol]
  double a0b = 5.291771e-9; // Bohr radius [cm]
  double HRKEV  = 2.721160e-2; // Hartree aomic units

  // POWERS OF THE SPEED OF LIGHT in atomic units:

  double SL2 = 1.877887e4; // c^2
  double SL22 = 3.755773e4; // 2c^2
  double SL22I = 2.662568e-5; // 1/2c^2

  double pi = 4*atan(1);
  double twopi = 2*pi;

  cout << "input data:" << endl;

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

  cout << "MEAN EXCITATION POTENTIAL " << expot*1e3 << " eV" << endl;

  expot = expot/HRKEV;
  expot = 1/expot;
  double expot2 = expot*expot;

  // FIRST IONIZATION ENERGY:

  double fie = 6; // eV
  //double fie = 12; // Si M-shell
  //double fie = 1.2; // Si band gap
  cout << "FIRST IONIZATION ENERGY = " << fie << " eV" << endl;

  fie = fie*1e-3 / HRKEV;      // Hartree
  double AB1 = 0.125*fie*fie;

  // NUMBER OF 'MOLECULES' PER UNIT VOLUME
  double VMOL = AVOG*rho/A * pow( a0b, 3 ); // a0b = Bohr

  // incident energy:

  double e0kev = 1;            // [keV]

  while( e0kev < 1e6 ) {

    double E = e0kev/HRKEV;             // atomic units

    double ZT = twopi * Z;

    double EINV = 1/E;
    double FA = E + SL2; // E + mc^2
    double RB = E + SL22;  // E + 2mc^2
    double BETA = E*RB/( FA*FA ); // beta^2
    double SFACT = BETA*SL2; // beta^2 mc^2
    double QMI = ( E+E )*RB*SL22I;
    double QFACT = QMI + QMI;
    double QME = QFACT + QFACT;

    // ELASTIC CROSS-SECTION ( MOTT'S FORMULA )

    FA = AL12 + QME;
    double FB = AL22 + QME;
    double SE0 = log( FB/FA ) - AA4;
    double SE1 = AA1 * ( QME / ( AL12*FA ) );
    double SE2 = AA2 * ( QME / ( AL22*FB ) );
    double SET = SE1 + SE2 + AA3*SE0;
    double EMFP = ZZ4 * SET;

    // INELASTIC CROSS-SECTION ( LENZ'S FORMULA )

    double Q02 = AB1*EINV;
    FA = AL12 + Q02;
    FB = AL22 + Q02;
    double FT = AL12 + QMI;
    double FNT = AL22 + QMI;
    double XMFP = Z4 * ( 
			AB2*log( QMI*FA/( FT*Q02 ) ) + 
			AB3*log( QMI*FB/( FNT*Q02 ) ) + 
			AA1*( 1/FT-1/FA ) + 
			AA2*( 1/FNT-1/FB ) );

    double TMFP = EMFP + XMFP;
    double AFP = SFACT / ( VMOL*TMFP ); // atomic units

    // AVERAGE ENERGY LOSS PER INELASTIC COLLISION: Bethe-Bloch

    double WA = fie;
    double CSI = E*expot;
    if( CSI < 6.338065465611359 ) // at 1 keV
      WA = ZT * EINV * 3.1776907270853916 * sqrt( CSI ) / XMFP;
    else { // above few keV
      FA = 1 - BETA; // 1/gamma2
      FB = sqrt( FA ); // 1/gamma
      FT = 1 - FB;
      WA = 2*ZT*EINV / XMFP *
	( log( E*expot2*SFACT/FA ) -
	  ( FB + FB + BETA ) * 6.931471805599453e-1 +
	  FA + 0.125*FT*FT );
    }

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

    cout << "  E " <<  E*HRKEV
	 << ", elas " << SFACT / ( VMOL*EMFP )*a0b*1e4
	 << ", inel " << SFACT / ( VMOL*XMFP )*a0b*1e4 << " um"
	 << ", csi " << CSI
	 << ", wavg " << WA*E*HRKEV << " keV"
	 << " = " << WA*E*HRKEV / ( SFACT / ( VMOL*XMFP )*a0b*1e4 ) << " keV/um"
	 << ", wmin " << WM*E*HRKEV << " keV"
	 << endl;

    e0kev *= 1.2;

  } // while

}
