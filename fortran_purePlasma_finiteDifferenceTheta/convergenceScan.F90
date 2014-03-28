module convergenceScan

  use parameters
  use petscsysdef

  implicit none

#include <finclude/petscsysdef.h>

  integer, allocatable :: NthetasForConvergenceScan(:), NxisForConvergenceScan(:), &
       NLsForConvergenceScan(:), NxsForConvergenceScan(:)
  PetscScalar, allocatable :: NxPotentialsPerVthsForConvergenceScan(:), xMaxsForConvergenceScan(:)
  integer :: numRunsForConvergenceScan, minScanUnit, maxScanUnit

contains

  subroutine processConvergenceScanParameters()

    implicit none

    integer, allocatable :: Nthetas(:), Nxis(:), NLs(:), Nxs(:)
    PetscScalar, allocatable :: NxPotentialsPerVths(:), xMaxs(:), tempArray(:)
    integer :: i, j, k, currentIndex
    PetscScalar :: verySmall = 1e-10

    allocate(tempArray(NthetaNumRuns))
    allocate(NThetas(NthetaNumRuns))
    call logspace(Ntheta*NthetaMinFactor, Ntheta*NthetaMaxFactor, NthetaNumRuns, tempArray)
    Nthetas = nint(tempArray)
    deallocate(tempArray)

    allocate(tempArray(NxiNumRuns))
    allocate(Nxis(NxiNumRuns))
    call logspace(Nxi*NxiMinFactor, Nxi*NxiMaxFactor, NxiNumRuns, tempArray)
    Nxis = nint(tempArray)
    deallocate(tempArray)

    allocate(tempArray(NLNumRuns))
    allocate(NLs(NLNumRuns))
    call logspace(NL*NLMinFactor, NL*NLMaxFactor, NLNumRuns, tempArray)
    NLs = nint(tempArray)
    deallocate(tempArray)

    allocate(tempArray(NxNumRuns))
    allocate(Nxs(NxNumRuns))
    call logspace(Nx*NxMinFactor, Nx*NxMaxFactor, NxNumRuns, tempArray)
    Nxs = nint(tempArray)
    deallocate(tempArray)

    allocate(NxPotentialsPerVths(NxPotentialsPerVthNumRuns))
    call logspace(NxPotentialsPerVth*NxMinFactor, NxPotentialsPerVth*NxMaxFactor, &
         NxPotentialsPerVthNumRuns, NxPotentialsPerVths)

    allocate(xMaxs(xMaxNumRuns))
    call logspace(xMax*xMaxMinFactor, xMax*xMaxMaxFactor, xMaxNumRuns, xMaxs)

    numRunsForConvergenceScan = 1 + NthetaNumRuns + NxiNumRuns + NLNumRuns + NxNumRuns + NxPotentialsPerVthNumRuns + xMaxNumRuns
    allocate(NthetasForConvergenceScan(numRunsForConvergenceScan))
    allocate(NxisForConvergenceScan(numRunsForConvergenceScan))
    allocate(NLsForConvergenceScan(numRunsForConvergenceScan))
    allocate(NxsForConvergenceScan(numRunsForConvergenceScan))
    allocate(NxPotentialsPerVthsForConvergenceScan(numRunsForConvergenceScan))
    allocate(xMaxsForConvergenceScan(numRunsForConvergenceScan))

    NthetasForConvergenceScan(1) = Ntheta
    NxisForConvergenceScan(1) = Nxi
    NLsForConvergenceScan(1) = NL
    NxsForConvergenceScan(1) = Nx
    NxPotentialsPerVthsForConvergenceScan(1) = NxPotentialsPerVth
    xMaxsForConvergenceScan(1) = xMax

    currentIndex=2

    do i=1, NThetaNumRuns
       NthetasForConvergenceScan(currentIndex) = Nthetas(i)
       NxisForConvergenceScan(currentIndex) = Nxi
       NLsForConvergenceScan(currentIndex) = NL
       NxsForConvergenceScan(currentIndex) = Nx
       NxPotentialsPerVthsForConvergenceScan(currentIndex) = NxPotentialsPerVth
       xMaxsForConvergenceScan(currentIndex) = xMax
       currentIndex = currentIndex + 1
    end do

    do i=1, NxiNumRuns
       NthetasForConvergenceScan(currentIndex) = Ntheta
       NxisForConvergenceScan(currentIndex) = Nxis(i)
       NLsForConvergenceScan(currentIndex) = NL
       NxsForConvergenceScan(currentIndex) = Nx
       NxPotentialsPerVthsForConvergenceScan(currentIndex) = NxPotentialsPerVth
       xMaxsForConvergenceScan(currentIndex) = xMax
       currentIndex = currentIndex + 1
    end do

    do i=1, NLNumRuns
       NthetasForConvergenceScan(currentIndex) = Ntheta
       NxisForConvergenceScan(currentIndex) = Nxi
       NLsForConvergenceScan(currentIndex) = NLs(i)
       NxsForConvergenceScan(currentIndex) = Nx
       NxPotentialsPerVthsForConvergenceScan(currentIndex) = NxPotentialsPerVth
       xMaxsForConvergenceScan(currentIndex) = xMax
       currentIndex = currentIndex + 1
    end do

    do i=1, NxNumRuns
       NthetasForConvergenceScan(currentIndex) = Ntheta
       NxisForConvergenceScan(currentIndex) = Nxi
       NLsForConvergenceScan(currentIndex) = NL
       NxsForConvergenceScan(currentIndex) = Nxs(i)
       NxPotentialsPerVthsForConvergenceScan(currentIndex) = NxPotentialsPerVth
       xMaxsForConvergenceScan(currentIndex) = xMax
       currentIndex = currentIndex + 1
    end do

    do i=1, NxPotentialsPerVthNumRuns
       NthetasForConvergenceScan(currentIndex) = Ntheta
       NxisForConvergenceScan(currentIndex) = Nxi
       NLsForConvergenceScan(currentIndex) = NL
       NxsForConvergenceScan(currentIndex) = Nx
       NxPotentialsPerVthsForConvergenceScan(currentIndex) = NxPotentialsPerVths(i)
       xMaxsForConvergenceScan(currentIndex) = xMax
       currentIndex = currentIndex + 1
    end do

    do i=1, xMaxNumRuns
       NthetasForConvergenceScan(currentIndex) = Ntheta
       NxisForConvergenceScan(currentIndex) = Nxi
       NLsForConvergenceScan(currentIndex) = NL
       NxsForConvergenceScan(currentIndex) = Nx
       NxPotentialsPerVthsForConvergenceScan(currentIndex) = NxPotentialsPerVth
       xMaxsForConvergenceScan(currentIndex) = xMaxs(i)
       currentIndex = currentIndex + 1
    end do

    if (currentIndex /= numRunsForConvergenceScan+1) then
       print *,"Error - something went wrong:"
       print *,"  currentIndex:",currentIndex
       print *,"  numRunsForConvergenceScan:",numRunsForConvergenceScan
       stop
    end if

    ! Now eliminate any duplicates:
    do i=1,numRunsForConvergenceScan-1
       do j=i+1, numRunsForConvergenceScan
          if ((abs(NthetasForConvergenceScan(i)-NthetasForConvergenceScan(j)) < verySmall) .and. &
               (abs(NxisForConvergenceScan(i)-NxisForConvergenceScan(j)) < verySmall) .and. &
               (abs(NLsForConvergenceScan(i)-NLsForConvergenceScan(j)) < verySmall) .and. &
               (abs(NxsForConvergenceScan(i)-NxsForConvergenceScan(j)) < verySmall) .and. &
               (abs(NxPotentialsPerVthsForConvergenceScan(i)-NxPotentialsPerVthsForConvergenceScan(j)) < verySmall) .and. &
               (abs(xMaxsForConvergenceScan(i)-xMaxsForConvergenceScan(j)) < verySmall) ) then

             ! Item j is a duplicate, so remove it:
             do k=j+1, numRunsForConvergenceScan
                NthetasForConvergenceScan(k-1) = NthetasForConvergenceScan(k)
                NxisForConvergenceScan(k-1) = NxisForConvergenceScan(k)
                NLsForConvergenceScan(k-1) = NLsForConvergenceScan(k)
                NxsForConvergenceScan(k-1) = NxsForConvergenceScan(k)
                NxPotentialsPerVthsForConvergenceScan(k-1) = NxPotentialsPerVthsForConvergenceScan(k)
                xMaxsForConvergenceScan(k-1) = xMaxsForConvergenceScan(k)
             end do
             numRunsForConvergenceScan = numRunsForConvergenceScan - 1
          end if
       end do
    end do

  end subroutine processConvergenceScanParameters

  !---------------------------------------------------------------------------------------------

  subroutine linspace(min, max, N, array)

    implicit none

    PetscScalar, intent(in) :: min, max
    integer, intent(in) :: N
    PetscScalar :: array(:)
    integer :: i

    if (N < 0) then
       print *,"Error! 'N' must be at least 0."
       stop
    end if

    do i=1,N
       array(i) = (i-one)*(max-min)/(N-one) + min
    end do

  end subroutine linspace

  !---------------------------------------------------------------------------------------------

  subroutine logspace(min, max, N, array)
    ! NOTE: this function works differently than the MATLAB logspace function.

    implicit none

    PetscScalar, intent(in) :: min, max
    integer, intent(in) :: N
    PetscScalar :: array(:)
    integer :: i

    if (min<0) then
       print *,"Error! 'min' must be >0."
       stop
    end if

    if (max<0) then
       print *,"Error! 'max' must be >0."
       stop
    end if

    call linspace(log(min), log(max), N, array)
    do i=1,N
       array(i) = exp(array(i))
    end do

  end subroutine logspace

  !---------------------------------------------------------------------------------------------

  subroutine setMPICommunicatorsForScan()
    ! This subroutine divides the PETSC_COMM_WORLD MPI communicator into smaller communicators.
    ! Each sub-group of processes can then work independently on different runs which are needed
    ! for a parameter scan.

    implicit none

    integer :: i, j, numProcsForThisComm
    integer :: currentCommunicator, currentMinUnit, firstIndex, lastIndex
    integer, dimension(:), allocatable :: commMaxProcs, procsToInclude
    integer :: convergenceScanUnits
    PetscErrorCode ierr
    MPI_Group :: mpi_group_world
    MPI_Group, dimension(:), allocatable :: mpi_groups_scan
    MPI_Comm, dimension(:), allocatable :: mpi_comms_scan

    ! In this code, the rank of each procedure is considered a 1-based index (consistent
    ! with Fortran convention) rather than a 0-based index (used in MPI calls).

    ! For now, each run required for the parameter scan will be considered a "unit".
    ! In the future, I might want to group several runs into each unit for better efficiency:
    ! some units would consist of 1 or 2 expensive runs, while other units would consist
    ! of many inexpensive runs.
    convergenceScanUnits = numRunsForConvergenceScan


    numCommunicators = min(numProcs, convergenceScanUnits)

    allocate(minUnits(numProcs))
    allocate(maxUnits(numProcs))
    allocate(procsToInclude(numProcs))
    allocate(mpi_comms_scan(numCommunicators))
    allocate(mpi_groups_scan(numCommunicators))
    allocate(commMinProcs(numCommunicators))
    allocate(commMaxProcs(numCommunicators))

    do i=1,numProcs
       minUnits(i) = floor(convergenceScanUnits * (1.d+0) * (i-1) / numProcs) + 1
       maxUnits(i) = floor(convergenceScanUnits * (1.d+0) * i / numProcs)
       maxUnits(i) = max(maxUnits(i), minUnits(i))
    end do
    minScanUnit = minUnits(myRank)
    maxScanUnit = maxUnits(myRank)

    if (masterProc) then
       do i=1,numProcs
          print "(a, i4, a, i3, a, i3)", "MPI proc ",i," is responsible for scan units ",&
               minUnits(i)," to ",maxUnits(i)
       end do
    end if

    commMinProcs(1) = 1
    commMaxProcs(numCommunicators) = numProcs
    currentCommunicator = 1
    currentMinUnit = 1
    myCommunicatorIndex = -1
    do i=2,numProcs
       if (minUnits(i) /= currentMinUnit) then
          currentMinUnit = minUnits(i)
          commMaxProcs(currentCommunicator) = i-1
          currentCommunicator = currentCommunicator + 1
          commMinProcs(currentCommunicator) = i
       end if
       if (myRank == i) then
          myCommunicatorIndex = currentCommunicator
       end if
    end do

    if (myRank == 1) then
       myCommunicatorIndex = 1
    end if
    if (myCommunicatorIndex == -1) then
       print "(a, i4, a)","Error! Somehow, myCommunicatorIndex for proc ",myRank," did not get assigned."
       stop
    end if

    if (currentCommunicator /= numCommunicators) then
       if (masterProc) then
          print *,"Error! Something went wrong in assigning processors to communicators."
       end if
       stop
    end if

    if (masterProc) then
       print "(a, i4, a)","Creating ",numCommunicators, &
            " MPI communicators for parallelizing the parameter scan."
       do i=1,numCommunicators
          print "(a, i4, a, i4, a, i4, a, i3, a, i3)", "Communicator ",i," consists of procs ", &
               commMinProcs(i)," through ",commMaxProcs(i), " and will handle scan units ", &
               minUnits(commMinProcs(i))," through ",maxUnits(commMinProcs(i))
       end do
    end if

    call MPI_Comm_group(PETSC_COMM_WORLD, mpi_group_world, ierr)

    do i=1,numCommunicators
       numProcsForThisComm = commMaxProcs(i)-commMinProcs(i) + 1
       do j=1,numProcsForThisComm
          procsToInclude(j) = commMinProcs(i) + j - 2
          ! Above, we needed to subtract an additional 1 to convert from Fortran 1-based indices
          ! to MPI 0-based indices.
       end do
!!$       if (masterProc) then
!!$          print *,"Communicator",i," includes procs ",procsToInclude(1:numProcsForThisComm)
!!$       end if
       call MPI_Group_incl(mpi_group_world, numProcsForThisComm, procsToInclude(1:numProcsForThisComm), mpi_groups_scan(i), ierr)
       call MPI_Comm_create(PETSC_COMM_WORLD, mpi_groups_scan(i), mpi_comms_scan(i), ierr)
    end do

    ! Next, set the MPI communicator that this process will use for the rest of the program execution:
    MPIComm = mpi_comms_scan(myCommunicatorIndex)
    CHKERRQ(ierr)

    call MPI_COMM_SIZE(MPIComm, numProcsInSubComm, ierr)
    call MPI_COMM_RANK(MPIComm, myRankInSubComm, ierr)
    myRankInSubComm = myRankInSubComm + 1
    masterProcInSubComm = (myRankInSubComm == 1)
  end subroutine setMPICommunicatorsForScan

  !---------------------------------------------------------------------------------------------

  subroutine allocateArraysForResults()

    implicit none

    allocate(qScan(numRunsForConvergenceScan))
    allocate(kScan(numRunsForConvergenceScan))
    allocate(particleFluxScan(numRunsForConvergenceScan))
    allocate(didItConvergeScan(numRunsForConvergenceScan))
    allocate(elapsedTimeScan(numRunsForConvergenceScan))

  end subroutine allocateArraysForResults

end module convergenceScan

