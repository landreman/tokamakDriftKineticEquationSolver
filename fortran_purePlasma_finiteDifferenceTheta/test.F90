program test
  use petscvecdef

  implicit none

#include <finclude/petscvecdef.h>

  PetscErrorCode ierr
  integer :: myRank, numProcs, i, j, numProcsForThisComm
  integer :: minUnit, maxUnit, numCommunicators, myCommunicator
  integer :: currentCommunicator, currentMinUnit, firstIndex, lastIndex
  integer, dimension(:), allocatable :: minUnits, maxUnits, commMinProcs, commMaxProcs, procsToInclude
  integer :: convergenceScanUnits = 5
  integer :: Nx = 8
  logical :: proc0

  Vec :: myVec

  MPI_Group :: mpi_group_world
  MPI_Group, dimension(:), allocatable :: mpi_groups_scan
  MPI_Comm :: mpi_comm_scan
  MPI_Comm, dimension(:), allocatable :: mpi_comms_scan

  ! In this code, the rank of each procedure is considered a 1-based index (consistent
  ! with Fortran convention) rather than a 0-based index (used in MPI calls).

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  CHKERRQ(ierr)

  call MPI_COMM_SIZE(PETSC_COMM_WORLD, numProcs, ierr)
  call MPI_COMM_RANK(PETSC_COMM_WORLD, myRank, ierr)
  proc0 = (myRank == 0)
  myRank = myRank + 1  ! Convert to 1-based Fortran indices
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
  minUnit = minUnits(myRank)
  maxUnit = maxUnits(myRank)

  !  print *,"Rank ",myRank," is responsible for units ", minUnit, " to ", maxUnit
  if (proc0) then
     do i=1,numProcs
        print "(a, i4, a, i3, a, i3)", "Proc ",i," is responsible for scan units ",&
             minUnits(i)," to ",maxUnits(i)
     end do
  end if

  commMinProcs(1) = 1
  commMaxProcs(numCommunicators) = numProcs
  currentCommunicator = 1
  currentMinUnit = 1
  myCommunicator = -1
  do i=2,numProcs
     if (minUnits(i) /= currentMinUnit) then
        currentMinUnit = minUnits(i)
        commMaxProcs(currentCommunicator) = i-1
        currentCommunicator = currentCommunicator + 1
        commMinProcs(currentCommunicator) = i
     end if
     if (myRank == i) then
        myCommunicator = currentCommunicator
     end if
  end do

  if (myRank == 1) then
     myCommunicator = 1
  end if
  if (myCommunicator == -1) then
     print "(a, i4, a)","Error! Somehow, myCommunicator for proc ",myRank," did not get assigned."
     stop
  end if

  if (currentCommunicator /= numCommunicators) then
     if (proc0) then
        print *,"Error! Something went wrong in assigning processors to communicators."
     end if
     stop
  end if

  if (proc0) then
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
     if (proc0) then
        print *,"Communicator",i," includes procs ",procsToInclude(1:numProcsForThisComm)
     end if
     call MPI_Group_incl(mpi_group_world, numProcsForThisComm, procsToInclude(1:numProcsForThisComm), mpi_groups_scan(i), ierr)
     call MPI_Comm_create(PETSC_COMM_WORLD, mpi_groups_scan(i), mpi_comms_scan(i), ierr)
  end do

  mpi_comm_scan = mpi_comms_scan(myCommunicator)

  ! Done building MPI communicators.

  CHKERRQ(ierr)
  call VecCreateMPI(mpi_comm_scan, PETSC_DECIDE, Nx, myVec, ierr)
  CHKERRQ(ierr)

  call VecGetOwnershipRange(myVec, firstIndex, lastIndex, ierr)
  CHKERRQ(ierr)
  lastIndex = lastIndex - 1
  print "(a, i4, a, i5, a, i5, a)","Proc ",myRank," owns indices ", firstIndex, " to ", &
       lastIndex," of a Vec."

  call VecDestroy(myVec, ierr)
  CHKERRQ(ierr)

  call PetscFinalize(ierr)

end program test
