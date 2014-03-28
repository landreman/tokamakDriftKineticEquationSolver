    module modalMultiplicationMatrix
    use fftModule
    use parameters
    implicit none
#include <finclude/petscsysdef.h>
#if defined(PETSC_USE_REAL_SINGLE)
#define zero 0.
#else
#define zero 0d+0
#endif

    contains 

    subroutine toeplitz(col, row, toeplitzResult)
      ! This subroutine emulates Matlab's 'toeplitz' function
      PetscScalar, intent(in), dimension(:) :: col, row
      PetscScalar, intent(out) :: toeplitzResult(:,:)
      integer :: numRows, numCols
      integer :: i

      numRows=size(col)
      numCols=size(row)

      do i=1,numRows
        toeplitzResult(i,i+1:numCols) = row(2:numCols-i+1)
        toeplitzResult(i,1:i) = col(i:1:-1)
      end do

    end subroutine

    subroutine generateThetaMultiplicationMatrix(functionOnThetaGrid, threshhold, m)
      PetscScalar, dimension(:), intent(in) :: functionOnThetaGrid
      PetscScalar, intent(in) :: threshhold
      PetscScalar, intent(out) :: m(:,:)
      
      integer :: Ntheta, thetaGridSize, i, j
      PetscScalar, dimension(:), allocatable :: realModeCoefficients, imagModeCoefficients
      PetscScalar, dimension(:), allocatable :: toeplitzGenerator
      PetscScalar, dimension(:,:), allocatable :: topLeftBlockBeforeMirroring, bottomLeftBlockBeforeMirroring
      
      if (modulo(size(functionOnThetaGrid),2)==1) then
        print *,"Error! Input to generateThetaMultiplicationMatrix must have even size."
        stop
      end if
      Ntheta = size(functionOnThetaGrid)/2
      thetaGridSize = Ntheta*2 - 1

      m = zero

      ! Allocate arrays:
      allocate(realModeCoefficients(2*Ntheta))
      allocate(imagModeCoefficients(2*Ntheta))
          

      call fft(2*Ntheta, functionOnThetaGrid, realModeCoefficients, imagModeCoefficients);
      realModeCoefficients = realModeCoefficients/(Ntheta + Ntheta)
      imagModeCoefficients = imagModeCoefficients/(Ntheta + Ntheta)

      ! Make the top-left block of matrix:
      allocate(toeplitzGenerator(2*Ntheta))
      toeplitzGenerator = [zero, realModeCoefficients(Ntheta:1:-1), realModeCoefficients(2:Ntheta)]
      allocate(topLeftBlockBeforeMirroring(Ntheta, 2*Ntheta))
      call toeplitz([(zero, i=1,Ntheta)], toeplitzGenerator, topLeftBlockBeforeMirroring)
      m(1:Ntheta, 1:Ntheta) = topLeftBlockBeforeMirroring(:,(Ntheta+1):Ntheta)
      m(1:Ntheta, 2:Ntheta) = m(1:Ntheta, 2:Ntheta) + topLeftBlockBeforeMirroring(:,Ntheta:2:-1)

      ! Make bottom-right block of matrix:
      m((Ntheta+1):thetaGridSize, (Ntheta+1):thetaGridSize) = &
      & topLeftBlockBeforeMirroring(1:(Ntheta-1), (Ntheta+1):thetaGridSize) &
      &   - topLeftBlockBeforeMirroring(1:(Ntheta-1), (Ntheta-1):1:-1)

      ! Make bottom-left block of matrix:
      toeplitzGenerator = [zero, imagModeCoefficients(Ntheta:2:-1), zero, -imagModeCoefficients(2:Ntheta)]
      allocate(bottomLeftBlockBeforeMirroring(Ntheta, 2*Ntheta))
      call toeplitz([(zero, i=1,Ntheta)],toeplitzGenerator, bottomLeftBlockBeforeMirroring)
      m((Ntheta+1):thetaGridSize, 1:Ntheta) = bottomLeftBlockBeforeMirroring(1:(Ntheta-1),Ntheta:thetaGridSize)
      m((Ntheta+1):thetaGridSize, 2:(Ntheta-1)) = m((Ntheta+1):thetaGridSize, 2:(Ntheta-1)) &
        & + bottomLeftBlockBeforeMirroring(1:(Ntheta-1),(Ntheta-1):2:-1)

      ! Make top-right block of matrix:
      m(1:Ntheta, (Ntheta+1):thetaGridSize) = -bottomLeftBlockBeforeMirroring(:,(Ntheta+2):(thetaGridSize+1))
      m(1:Ntheta, (Ntheta+1):thetaGridSize) = m(1:Ntheta, (Ntheta+1):thetaGridSize) + &
        & bottomLeftBlockBeforeMirroring(:,Ntheta:2:-1)


      where (abs(m)<threshhold)
        m=zero
      end where

    end subroutine
    
    subroutine makeFourierSpectralDifferentiationMatrix(N, xMin, xMax, x, matrix)
      integer :: N, i
      PetscScalar :: xMin, xMax, matrix(:,:), x(:)
      PetscScalar :: Delta, h
      PetscScalar, allocatable :: col(:), row(:)
      
      if (mod(N,2)==1) then
        print *,"Error! N must be even in input to makeFourierSpectralDifferentiationMatrix."
        stop
      endif
      
      Delta = xMax-xMin;
      h = 2*pi/N;
      x = [(i*h/(2*pi)*Delta + xMin, i=1, N)]
      
      allocate(col(N))
      allocate(row(N))
      col(1)=0
      row(1)=0
      do i=1,N-1
        col(i+1) = 0.5*((-1)**i)/tan(i*h/2)*2*pi/Delta
      end do
      do i=2,N
        row(i) = col(N-i+2)
      end do
      call toeplitz(col, row, matrix)
      
    end subroutine

  end module
