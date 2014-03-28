  module fftModule
    implicit none
#include <finclude/petscsysdef.h>
    contains
    subroutine fft(N, input, outputRealPart, outputImagPart)
      PetscScalar, intent(in) :: input(:)
      integer, intent(in) :: N
      PetscScalar, intent(out) :: outputImagPart(:), outputRealPart(:)
      PetscScalar, allocatable :: workspace(:), array(:)
      integer :: halfN
      
      allocate(workspace(2*N+15))
      allocate(array(N))

      ! Copy input array so it is not over-written:
      array(1:N) = input(1:N)
      
      ! Call FFTPACK routines:
#if defined(PETSC_USE_REAL_SINGLE)
      ! Use single precision
      ! First initialize the workspace:
      call rffti(N, workspace)
      ! Now carry out the Fourier transform:
      call rfftf(N, array, workspace)
#else
      ! Use double precision
      ! First initialize the workspace:
      call dffti(N, workspace)
      ! Now carry out the Fourier transform:
      call dfftf(N, array, workspace)
#endif      
!       print *,"array:"
!       do i=1,N
!         print *,array(i)
!       end do

      if (mod(N,2)==0) then
        ! N is even
        halfN = (N)/2
        outputRealPart(1) = array(1)
        outputRealPart(2:halfN+1) = array(2:N:2)
        outputRealPart(halfN+2:N) = array(N-2:2:-2)
        outputImagPart(1) = 0d+0
        outputImagPart(halfN+1) = 0d+0
        outputImagPart(2:halfN) = array(3:N-1:2)
        outputImagPart(halfN+2:N) = array(N-1:3:-2)
      else
        ! N is odd
        halfN = (N+1)/2
        outputRealPart(1) = array(1)
        outputRealPart(2:halfN) = array(2:N-1:2)
        outputRealPart(halfN+1:N) = array(N-1:2:-2)
        outputImagPart(1) = 0d+0
        outputImagPart(2:halfN) = array(3:N:2)
        outputImagPart(halfN+1:N) = array(N:3:-2)
      end if
      
    end subroutine
  end module

