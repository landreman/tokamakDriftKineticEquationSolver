! For PETSc versions 3.4 and later, the line below should be un-commented:
!#define PetscGetTime PetscTime

! Main program
  program tokamakDriftKineticEquationSolver
    use parameters
    use geometry
    use writeOutput
    use solveDKEModule
    use readInput
    use convergenceScan
    use petscsysdef

    implicit none

#include <finclude/petscsysdef.h>
    PetscErrorCode ierr
    PetscLogDouble :: time1, time2
    integer :: runNum

    call printGreeting()
    
    ! Default values for all input parameters are set in parameters.F90.
    ! Defaults are replaced by any values given in the input namelist file:
    call readNamelistInput("input.namelist")
    
    nuStar = nuPrime / (epsil * sqrt(epsil))
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    CHKERRQ(ierr)

    ! In case of Miller geometry or a geometry loaded from a file, do any required initialization:
    call initializeGeometry()

    select case (programMode)
      case (SINGLE_SOLVE)
        call printInputs()
        call solveDKE()
        call printOutputs()
      case (CONVERGENCE_SCAN)
        print *,"Beginning convergence scan."
        call processConvergenceScanParameters()
        call openConvergenceScanOutputFile()
        call PetscGetTime(time1, ierr)

        do runNum = 1,numRunsForConvergenceScan
          print *,"--------------------------------------------------------------"
          print *,"Run", runNum, "of", numRunsForConvergenceScan
          Ntheta = NthetasForConvergenceScan(runNum)
          Nxi = NxisForConvergenceScan(runNum)
          NL = NLsForConvergenceScan(runNum)
          Nx = NxsForConvergenceScan(runNum)
          NxPotentialsPerVth = NxPotentialsPerVthsForConvergenceScan(runNum)
          xMax = xMaxsForConvergenceScan(runNum)
          call printInputs()
          call solveDKE()
          call printOutputs()
          call addLineToConvergenceScanFile()
        end do
        call PetscGetTime(time2, ierr)
        call closeConvergenceScanOutputFile()
        print *,"--------------------------------------------------------------"
        print *,"Total time for convergence scan: ", time2-time1, "seconds."
      case default
        print *,"Error: invalid programMode"
        stop
    end select
    
    
    call PetscFinalize(ierr)

  end program
