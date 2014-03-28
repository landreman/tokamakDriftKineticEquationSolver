! For PETSc versions 3.4 and later, the line below should be un-commented:
!#define PetscGetTime PetscTime

! Main program

! To run with 1 processor, you can just use
! ./tokamakDriftKineticEquationSolver
! To run with more processors than runs in a parameter scan, use
! mpiexec -n ## ./tokamakDriftKineticEquationSolver -pc_factor_mat_solver_package superlu_dist
! or
! mpiexec -n ## ./tokamakDriftKineticEquationSolver -pc_factor_mat_solver_package mumps
! depending on whether superLU_dist or Mumps is installed.

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

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  CHKERRQ(ierr)

  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
  ! This code uses 1-based indices to number the MPI processes, consistent with the Fortran
  ! convention for array indices, not the 0-based indexing used by C and MPI.
  myRank = myRank + 1
  masterProc = (myRank==1)

  call printGreeting()
  
  ! Default values for all input parameters are set in parameters.F90.
  ! Defaults are replaced by any values given in the input namelist file:
  call readNamelistInput("input.namelist")

  nuStar = nuPrime / (epsil * sqrt(epsil))

  call chooseParallelDirectSolver()

  ! In case of Miller geometry or a geometry loaded from a file, do any required initialization:
  call initializeGeometry()

  select case (programMode)
  case (SINGLE_SOLVE)
     ! There is no parallelization over runs, so set the sub-commmunicator 
     ! to be equal to the global communicator:
     MPIComm = PETSC_COMM_WORLD
     myRankInSubComm = myRank
     numProcsInSubComm = numProcs
     masterProcInSubComm = masterProc
     myCommunicatorIndex = 1

     call printInputs()
     call solveDKE()
     call printOutputs()
  case (CONVERGENCE_SCAN)
     call processConvergenceScanParameters()
     call allocateArraysForResults()
     if (masterProc) then
        print *,"Beginning a convergence scan involving ",numRunsForConvergenceScan," runs."
        print *,"The numbers in brackets below indicate the MPI communicator."
     end if
     call setMPICommunicatorsForScan()
     call PetscGetTime(time1, ierr)

     do runNum = minScanUnit,maxScanUnit
        if (masterProcInSubComm) then
           print *,"[",myCommunicatorIndex,"] --------------------------------------------------------------"
           print *,"[",myCommunicatorIndex,"] Run", runNum, "of", numRunsForConvergenceScan
        end if
        Ntheta = NthetasForConvergenceScan(runNum)
        Nxi = NxisForConvergenceScan(runNum)
        NL = NLsForConvergenceScan(runNum)
        Nx = NxsForConvergenceScan(runNum)
        NxPotentialsPerVth = NxPotentialsPerVthsForConvergenceScan(runNum)
        xMax = xMaxsForConvergenceScan(runNum)
        call printInputs()
        call solveDKE()
        call printOutputs()
        qScan(runNum) = q
        kScan(runNum) = k
        particleFluxScan(runNum) = particleFlux
        didItConvergeScan(runNum) = didItConverge
        elapsedTimeScan(runNum) = elapsedTime
     end do
     call PetscGetTime(time2, ierr)
     if (masterProcInSubComm) then
        print *,"[",myCommunicatorIndex,"] --------------------------------------------------------------"
        print *,"[",myCommunicatorIndex,"] Total time for convergence scan on this communicator: ", &
             time2-time1, "seconds."
     end if
     call writeConvergenceScanFile()
  case default
     if (masterProc) then
        print *,"Error: invalid programMode"
     end if
     stop
  end select


  call PetscFinalize(ierr)

end program tokamakDriftKineticEquationSolver
