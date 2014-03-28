tokamakDriftKineticEquationSolver
=================================

This repository contains several variants of a Fokker-Planck drift-kinetic code for neoclassical calculations in tokamaks.

The codes in this repository solve kinetic equations with 3 independent variables: poloidal angle in real space, and speed and pitch angle in velocity space. These codes were predecessors to the more sophisticated 4D codes [PERFECT](https://github.com/landreman/perfect) and [SFINCS](https://github.com/landreman/sfincs).  All of these codes are time-independent.

The codes in this repository implement the algorithms discussed in the following paper: [Landreman and Ernst, J. Comp. Phys. 243, 130 (2013)](http://dx.doi.org/10.1016/j.jcp.2013.02.041).  Versions of these codes were also used to generate figures 2-7 in [Landreman and Ernst, Plasma Phys. Controlled Fusion 54, 115006 (2012)](http://dx.doi.org/10.1088/0741-3335/54/11/115006).  (However, later figures in this paper were made using a different time-dependent 4D code.)

The versions in this repository that use a Fourier modal discretization in theta demonstrate somewhat better convergence in this coordinate than the versions that use finite differencing in theta. This advantage can be important when the collisionality is low so very high resolution is required.  However, the versions that use finite differencing in theta are shorter and easier to understand.

Several options are available for the magnetic geometry:
  1. A circular concentric flux surface model
  2. Miller geometry
  3. Reading an eqdsk/EFIT g file

All codes in this repository implement the full linearized Fokker-Planck collision operator.
