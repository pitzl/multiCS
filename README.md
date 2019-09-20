
# MultiCS

## Multiple Coulomb Scattering of electrons in matter

F. Salvat, J.D. Martinez, R. Mayol, J. Parellada:  
A Monte Carlo program to simulate the penetration and energy loss of keV
electrons through matter  
Computer Physics Communications 42 (1986) 93—104

Code mcsda.f stored as AALC at the CPC program library in Belfast

translated into Fortran 90 and C++
- multics.f90  Fortran90
- multicspp.cpp  plain C++
- multics.cc with ROOT, vertical through a thin infinite layer
- edgesc.cc with ROOT, edge-on through a thin layer

Successfully applied to the multiple scattering of GeV electrons in thin
( 0.1 mm to 5 mm ) layers of silicon.

The energy loss calculation uses the original Landau theory (d sigma/d dE ~ 1/dE^2)
without shell corrections and only approximates the width and most probable value
of the distribton, while the mean value is quite OK (Bethe-Bloch).

## bugs fixed from mcsda.f

- CALL RMBRG( rngd, EC, E, D, 1D-6, 10, 6 ) ! 10 gets overwritten!

  SUBROUTINE RMBRG( fct, XL, XU, Y, EPS, NPART, IPRNT )      

  NPART = -NPART ! flag low accuracy

- NaN: bug in DIRECT
  U = AW*DSQRT( 1D0-V*V-W*W ) ! NaN

- DIRECT is wrong at small angles
  completely changed with proper rotation matrix
  multiple scattering correct down to mrad for GeV electrons

- Wentzel model needs AA1 but it is not assigned in case of IS = 0
