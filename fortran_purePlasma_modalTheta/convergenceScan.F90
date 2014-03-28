  module convergenceScan
    use parameters
    implicit none
#include <finclude/petscsysdef.h>
    integer, allocatable :: NthetasForConvergenceScan(:), NxisForConvergenceScan(:), &
      NLsForConvergenceScan(:), NxsForConvergenceScan(:)
    PetscScalar, allocatable :: NxPotentialsPerVthsForConvergenceScan(:), xMaxsForConvergenceScan(:)
    integer :: numRunsForConvergenceScan
      
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
      
    end subroutine
    
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
!      array = [((i-one)*(max-min)/(N-one) + min, i=1,N)]
    end subroutine
    
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
!      logspace = exp(linspace(log(min), log(max), N))
    end subroutine
  
  
  end module

