module writeOutput

  use parameters
  use convergenceScan
  use petscsysdef

  implicit none

#include <finclude/petscsysdef.h>

  public printInputs
  integer, private :: convergenceScanFileUnit

contains

  ! -----------------------------------------------------------------------------------

  subroutine printGreeting()
    
    implicit none

    if (masterProc) then
       print *,"****************************************************************************"
       print *,"Local tokamak drift-kinetic equation solver."
       print *,"Collocation grid in theta, Legendre modal in xi, polynomial spectral collocation in x."
#if defined(PETSC_USE_REAL_SINGLE)
       print *,"Using single precision."
#else
       print *,"Using double precision."
#endif
       print *,numProcs," processes detected."
    end if
  end subroutine printGreeting

  ! -----------------------------------------------------------------------------------

  subroutine printInputs()

    implicit none

    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] ---- Physics parameters: ----"
       print *,"[",myCommunicatorIndex,"] epsilon = ", epsil
       print *,"[",myCommunicatorIndex,"] nuStar  = ", nuStar
       print *,"[",myCommunicatorIndex,"] nuPrime = ", nuPrime
       print *,"[",myCommunicatorIndex,"] Geometry scheme: ", geometryToUse
       print *,"[",myCommunicatorIndex,"] ---- Numerical parameters: ----"
       print *,"[",myCommunicatorIndex,"] Ntheta  = ", Ntheta
       print *,"[",myCommunicatorIndex,"] Nxi     = ", Nxi
       print *,"[",myCommunicatorIndex,"] NL      = ", NL
       print *,"[",myCommunicatorIndex,"] Nx      = ", Nx
       print *,"[",myCommunicatorIndex,"] NxPotentialsPerVth = ", NxPotentialsPerVth
       print *,"[",myCommunicatorIndex,"] xMax = ",xMax
       select case (thetaDerivativeScheme)
       case (0)
          print *,"[",myCommunicatorIndex,"] Theta derivative: spectral collocation"
       case (1)
          print *,"[",myCommunicatorIndex,"] Theta derivative: 2nd order finite difference"
       case (2)
          print *,"[",myCommunicatorIndex,"] Theta derivative: 4th order finite difference"
       case default
          print *,"[",myCommunicatorIndex,"] Error! Invalid setting for thetaDerivativeScheme"
          stop
       end select
       if (useIterativeSolver) then
          print *,"[",myCommunicatorIndex,"] Using iterative solver"
       else
          print *,"[",myCommunicatorIndex,"] Using direct solver"
       end if
    end if
  end subroutine printInputs

  ! -----------------------------------------------------------------------------------

  subroutine printOutputs()

    implicit none

    if (masterProcInSubComm) then
       print *,"[",myCommunicatorIndex,"] Total elapsed time: ", elapsedTime, " seconds."
       print *,"[",myCommunicatorIndex,"] q = ", q
       print *,"[",myCommunicatorIndex,"] particleFlux = ", particleFlux
       print *,"[",myCommunicatorIndex,"] k = ", k
    end if

  end subroutine printOutputs

  ! -----------------------------------------------------------------------------------

  subroutine writeConvergenceScanFile()

    implicit none

    integer :: didFileAccessWork, i, size, minUnitOther, maxUnitOther
    PetscErrorCode :: ierr
    integer :: MPIStatus(MPI_STATUS_SIZE)

    if (.not. masterProc) then
       ! Send scan results to masterProc

       if (masterProcInSubComm) then
          size = maxScanUnit - minScanUnit + 1
          call MPI_Send(qScan(minScanUnit:maxScanUnit), size, MPIU_SCALAR, 0, 0, PETSC_COMM_WORLD, ierr)
          call MPI_Send(kScan(minScanUnit:maxScanUnit), size, MPIU_SCALAR, 0, 0, PETSC_COMM_WORLD, ierr)
          call MPI_Send(particleFluxScan(minScanUnit:maxScanUnit), size, MPIU_SCALAR, 0, 0, PETSC_COMM_WORLD, ierr)
          call MPI_Send(elapsedTimeScan(minScanUnit:maxScanUnit), size, MPI_DOUBLE, 0, 0, PETSC_COMM_WORLD, ierr)
          call MPI_Send(didItConvergeScan(minScanUnit:maxScanUnit), size, MPI_INT, 0, 0, PETSC_COMM_WORLD, ierr)
          CHKERRQ(ierr)
       end if
    else
       ! Collect scan results from other processes
       do i=2,numCommunicators
          maxUnitOther = maxUnits(commMinProcs(i))
          minUnitOther = minUnits(commMinProcs(i))
          size = maxUnitOther - minUnitOther + 1

          call MPI_Recv(qScan(minUnitOther:maxUnitOther), size, MPIU_SCALAR, commMinProcs(i)-1, 0, &
               PETSC_COMM_WORLD, MPIStatus, ierr)
          call MPI_Recv(kScan(minUnitOther:maxUnitOther), size, MPIU_SCALAR, commMinProcs(i)-1, 0, &
               PETSC_COMM_WORLD, MPIStatus, ierr)
          call MPI_Recv(particleFluxScan(minUnitOther:maxUnitOther), size, MPIU_SCALAR, commMinProcs(i)-1, 0, &
               PETSC_COMM_WORLD, MPIStatus, ierr)
          call MPI_Recv(elapsedTimeScan(minUnitOther:maxUnitOther), size, MPI_DOUBLE, commMinProcs(i)-1, 0, &
               PETSC_COMM_WORLD, MPIStatus, ierr)
          call MPI_Recv(didItConvergeScan(minUnitOther:maxUnitOther), size, MPI_INT, commMinProcs(i)-1, 0, &
               PETSC_COMM_WORLD, MPIStatus, ierr)
          CHKERRQ(ierr)
       end do

       ! Open output file
       open(unit=convergenceScanFileUnit, file=convergenceScanFilename,    &
            action="write", status="replace", iostat=didFileAccessWork)
       if (didFileAccessWork /= 0) then
          print *,"Error opening ", convergenceScanFilename, "."
          stop
       end if

       ! Write header in output file
       write (unit=convergenceScanFileUnit, fmt="(a)") &
            "Ntheta,Nxi,NL,Nx,NxPotentialsPerVth,xMax,q,particleFlux,k,didItConverge,elapsedTime"

       ! Write results to output file
       do i=1,numRunsForConvergenceScan
          write (unit=convergenceScanFileUnit, &
               fmt="(i5, a, i5, a, i5, a, i5, a, f18.10, a, f18.10, a, f18.10, a, f18.10, a, f18.10, a, i5, a, f18.10)") &
               NthetasForConvergenceScan(i), ",", NxisForConvergenceScan(i), ",", NLsForConvergenceScan(i), ",", &
               NxsForConvergenceScan(i), ",", NxPotentialsPerVthsForConvergenceScan(i), ",", &
               xMaxsForConvergenceScan(i),&
               ",", qScan(i), ",", particleFluxScan(i), &
               ",", kScan(i), ",", didItConvergeScan(i), ",", elapsedTimeScan(i)
       end do

       close(unit=convergenceScanFileUnit)

    end if

  end subroutine writeConvergenceScanFile


end module writeOutput

