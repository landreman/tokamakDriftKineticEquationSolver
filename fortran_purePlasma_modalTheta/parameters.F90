  module parameters
    implicit none
#include <finclude/petscsysdef.h>

    PetscScalar :: one = 1., oneHalf = 0.5d+0
    PetscScalar :: zero = 0., two = 2., three = 3., four = 4., five = 5.
    
    PetscScalar, parameter :: pi = 3.14159265358979d+0
    PetscScalar, parameter :: sqrtpi = 1.77245385090552d+0

    ! Enumerate options for 'programMode':
    integer, parameter :: SINGLE_SOLVE = 1
    integer, parameter :: CONVERGENCE_SCAN = 2
    
    
    ! Enumerate options for geometry:
    integer, parameter :: GEOMETRY_CIRCULAR_CONCENTRIC = 0
    integer, parameter :: GEOMETRY_MILLER = 1
    integer, parameter :: GEOMETRY_CIRCULAR_CONCENTRIC_BOOZER = 2
    
    ! Options for program flow control
    integer :: programMode = SINGLE_SOLVE
    logical :: saveMatlabOutput = .false.
    character(len=200) :: MatlabOutputFilename = "dke.m"
    character(len=200) :: convergenceScanFilename = "dkeConvergence.dat"

    ! Physics input parameters:
    PetscScalar :: epsil = 0.3
    PetscScalar :: nuPrime = 0.2
    ! No default value is assigned to nuStar - it is calculated from nuPrime.
    PetscScalar :: nuStar
    integer :: geometryToUse = GEOMETRY_CIRCULAR_CONCENTRIC_BOOZER

    ! Miller parameters: (these are only used when geometryToUse = 1.)
    PetscScalar :: Miller_kappa = 1.66
    PetscScalar :: Miller_delta = 0.416
    PetscScalar :: Miller_s_delta = 1.37
    PetscScalar :: Miller_s_kappa = 0.70
    PetscScalar :: Miller_dRdr = -0.354
    PetscScalar :: Miller_q = 3.03
    ! The inverse aspect ratio epsil is also used for Miller geometry.
    
    ! Numerical input parameters:
    integer :: Ntheta = 10
    PetscScalar :: NthetaMinFactor = 0.8d+0, NthetaMaxFactor=2d+0
    integer :: NthetaNumRuns = 5
    
    integer :: Nxi = 13
    PetscScalar :: NxiMinFactor = 0.8d+0, NxiMaxFactor=2d+0
    integer :: NxiNumRuns = 5
    
    integer :: NL = 4
    PetscScalar :: NLMinFactor = 0.8d+0, NLMaxFactor=2d+0
    integer :: NLNumRuns = 5
    
    integer :: Nx = 7
    PetscScalar :: NxMinFactor = 0.8d+0, NxMaxFactor=2d+0
    integer :: NxNumRuns = 5

    PetscScalar  :: NxPotentialsPerVth = 2d+0
    PetscScalar :: NxPotentialsPerVthMinFactor = 0.8d+0, NxPotentialsPerVthMaxFactor=2d+0
    integer :: NxPotentialsPerVthNumRuns = 5
    
    PetscScalar :: xMax = 5.
    PetscScalar :: xMaxMinFactor = 0.8d+0, xMaxMaxFactor=2d+0
    integer :: xMaxNumRuns = 5

    integer :: whereToPutbDotGradTheta = 0
    ! 0 = my traditional form of the DKE, with b dot grad theta on the
    ! collision term but not streaming or mirror terms: this makes the
    ! streaming term diagonal in M but makes the collision term dense in M.

    ! 1 = put the b dot grad theta on the streaming and mirror terms. This
    ! makes the collision term diagonal in M but makes the streaming term dense
    ! in M.

    logical :: useIterativeSolver = .false.
    
    PetscScalar :: multipleOfNuPrimeToAddOnDiagonalOfPreconditioner = 1d+0

    ! Outputs:
    PetscScalar :: q, particleFlux, k
    PetscLogDouble :: elapsedTime
    integer :: didItConverge

  end module

