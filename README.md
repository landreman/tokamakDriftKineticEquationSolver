tokamakDriftKineticEquationSolver
=================================

This repository contains several variants of a Fokker-Planck drift-kinetic code for neoclassical calculations in tokamaks.

The codes in this repository solve kinetic equations with 3 independent variables: poloidal angle in real space, and speed and pitch angle in velocity space. These codes were predecessors to the more sophisticated 4D codes [PERFECT](https://github.com/landreman/perfect) and [SFINCS](https://github.com/landreman/sfincs).  All of these codes are time-independent.

The codes in this repository implement the algorithms discussed in the following paper: [Landreman and Ernst, J. Comp. Phys. 243, 130 (2013)](http://dx.doi.org/10.1016/j.jcp.2013.02.041).  Versions of these codes were also used to generate figures 2-7 in [Landreman and Ernst, Plasma Phys. Controlled Fusion 54, 115006 (2012)](http://dx.doi.org/10.1088/0741-3335/54/11/115006).  (However, later figures in this paper were made using a different time-dependent 4D code.)

The versions in this repository that use a Fourier modal discretization in theta demonstrate somewhat better convergence in this coordinate than the versions that use finite differencing in theta. This advantage can be important when the collisionality is low so very high resolution is required.  However, the versions that use finite differencing in theta are shorter and easier to understand.

Several options are available for the magnetic geometry in the Matlab codes:
  1. A circular concentric flux surface model
  2. Miller geometry
  3. Reading an eqdsk/EFIT g file
The third option is not yet available in the fortran versions.

All codes in this repository implement the full linearized Fokker-Planck collision operator.

The fortran versions rely on the [PETSc library](http://www.mcs.anl.gov/petsc/). The makefiles may need to be adapted to link to PETSc libraries on your computing system.

At the moment, the Matlab versions compute both ion and electron quantities, while the fortran versions compute only the ion quantities. It is straightforward to add the electron quantities to the fortran versions - I just haven't got around to it.  If you would like such a version, send me an email at matt dot landreman at gmail dot com.

I have other versions of these codes available that use different discretizations in the speed and pitch-angle coordinates. If you are interested in these versions, send me an email.

Versions that allow for multiple ion species will be uploaded soon...