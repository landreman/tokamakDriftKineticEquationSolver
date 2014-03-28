! For PETSc versions 3.4 and later, the line below should be un-commented:
!#define PetscGetTime PetscTime

  module solveDKEModule
    use geometry
    use parameters
    use modalMultiplicationMatrix
    use poldif
    use polyInterp
    use makeXGridModule
    use petscksp
    implicit none
#include <finclude/petsckspdef.h>

    ! 'threshhold' is used to zero out small elements in theta multiplication matrices.
    PetscScalar, parameter :: threshhold = 1d-7

    contains

    subroutine solveDKE()

      PetscErrorCode :: ierr
      Vec :: rhs, soln
      Mat :: matrix, preconditionerMatrix
      PetscViewer MatlabOutput
      PetscScalar, allocatable :: thetas(:), thetaWeights(:), bs(:), dbdthetas(:), oneOverqRbDotGradThetas(:)
      integer :: i, j, ix, itheta, L, NxPotentials, thetaGridSize, matrixSize
      integer :: ithetaRow, ithetaCol
      PetscScalar, dimension(:), allocatable :: x, xWeights, xPotentials, xWeightsPotentials
      PetscScalar, dimension(:), allocatable :: x2, xPartOfRHS, thetaPartOfMirrorTermOnThetaGrid
      PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
      PetscScalar, dimension(:,:), allocatable :: regridPolynomialToUniform, regridUniformToPolynomial
      PetscScalar :: VPrime, dtheta, xMaxNotTooSmall
      PetscScalar, dimension(:), allocatable :: thetaPartOfRHSOnThetaGrid, thetaPartOfRHSModeCoefficients
      PetscScalar, dimension(:), allocatable :: realModeCoefficients, imagModeCoefficients
      integer, dimension(:), allocatable :: indices, rowIndices, colIndices
      integer, dimension(:,:), allocatable :: ddtheta
      PetscScalar, allocatable :: thetaPartMatrixForStreamingTerm(:,:), thetaPartMatrixForMirrorTerm(:,:)
      PetscScalar, dimension(:,:), allocatable :: thetaPartMatrixForCollisionTerm
      PetscScalar, dimension(:,:), allocatable :: M11, M21, M32, LaplacianTimesX2WithoutL
      PetscScalar, dimension(:,:), allocatable :: xPartOfCECD, M12IncludingX0, M13IncludingX0
      PetscScalar, dimension(:), allocatable :: erfs, x3, expx2, Psi, nuD, PsiPrime
      PetscScalar, dimension(:,:), allocatable :: KWithoutThetaPart, M22, M33, M12, M13
      PetscScalar, dimension(:), allocatable :: diagonalOfKWithoutThetaPart
      PetscScalar, dimension(:,:), allocatable :: M22BackslashM21, M33BackslashM32
      integer, dimension(:), allocatable :: IPIV  ! Needed by LAPACK
      integer :: LAPACKInfo, predictedNNZForEachRowOfPreconditioner, predictedNNZForEachRowOfTotalMatrix
      PetscScalar, allocatable :: particleFluxBeforeThetaIntegral(:), qBeforeThetaIntegral(:), &
        flowDividedByB(:), density(:), fSlice(:)
      PetscScalar, allocatable :: particleFluxIntegralWeight(:), qIntegralWeight(:), &
        kIntegralWeight(:), densityIntegralWeight(:)
      PetscScalar, dimension(:,:), allocatable :: thetaMatrixForFlow, thetaMatrixForFluxes
      
      character :: trans='n'
      PetscLogDouble :: time1, time2, startTime
      KSP :: KSPInstance
      PC :: preconditionerContext
      KSPConvergedReason :: reason
      PetscScalar, pointer :: solnArray(:)
      
      call PetscGetTime(time1, ierr)
      startTime = time1

      thetaGridSize = 2*Ntheta-1
      matrixSize = thetaGridSize * Nxi * Nx
      allocate(thetas(2*Ntheta))
      allocate(thetaWeights(2*Ntheta))
      thetas = [((i-1)*2*pi/(2*Ntheta), i=1,2*Ntheta)]
      allocate(bs(2*Ntheta))
      allocate(dbdthetas(2*Ntheta))
      allocate(oneOverqRbDotGradThetas(2*Ntheta))
      call computeBs(thetas, bs)
      dtheta = thetas(2)-thetas(1)
      thetaWeights = [(dtheta, i=1, 2*Ntheta)]
      call computedBdthetas(thetas, dbdthetas)
      call computeOneOverqRbDotGradThetas(thetas, oneOverqRbDotGradThetas)

      
      allocate(x(Nx))
      allocate(xWeights(Nx))
      call makeXGrid(Nx, x, xWeights)
      xWeights = xWeights / exp(-x*x)
      xMaxNotTooSmall = max(x(Nx), xMax)
      
      allocate(ddx(Nx,Nx))
      allocate(d2dx2(Nx,Nx))
      call makeXPolynomialDiffMatrices(x,ddx,d2dx2)
      
      
      NxPotentials = ceiling(xMaxNotTooSmall*NxPotentialsPerVth)
      
      
      allocate(xPotentials(NxPotentials))
      allocate(xWeightsPotentials(NxPotentials))
      allocate(ddxPotentials(NxPotentials, NxPotentials))
      allocate(d2dx2Potentials(NxPotentials, NxPotentials))
!!$      call makeUniformHighOrderDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, xPotentials, &
!!$        xWeightsPotentials, ddxPotentials, d2dx2Potentials)
      call uniformDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, 12, xPotentials, &
        xWeightsPotentials, ddxPotentials, d2dx2Potentials)

      allocate(regridPolynomialToUniform(NxPotentials, Nx))
      call makePolynomialInterpolationMatrix(x,xPotentials,exp(-x*x),exp(-xPotentials*xPotentials), regridPolynomialToUniform)
      allocate(regridUniformToPolynomial(Nx,NxPotentials))
      call interpolationMatrix(NxPotentials, Nx, xPotentials, x, regridUniformToPolynomial, -1, 0)

      
      ! *********************************************************
      ! Create right-hand side
      ! *********************************************************
      call VecCreate(PETSC_COMM_WORLD, rhs, ierr);
      CHKERRQ(ierr)
      call VecSetType(rhs, VECSEQ, ierr)
      CHKERRQ(ierr)
      call VecSetSizes(rhs, PETSC_DECIDE, matrixSize, ierr)
      CHKERRQ(ierr)
      
      allocate(x2(Nx))
      allocate(xPartOfRHS(Nx))
      allocate(thetaPartOfRHSOnThetaGrid(2*Ntheta))
      allocate(thetaPartOfRHSModeCoefficients(thetaGridSize))
      x2=x*x
      xPartOfRHS = (x2 - 5/(one+one))*x2*exp(-x2)
      select case (whereToPutbDotGradTheta)
        case (0)
          thetaPartOfRHSOnThetaGrid = dbdthetas / (bs*bs)
        case (1)
          thetaPartOfRHSOnThetaGrid = dbdthetas / (bs*bs*oneOverqRbDotGradThetas)
        case default
          print *,"Error! invalid whereToPutbDotGradTheta"
          stop
      end select
      ! Allocate arrays:
      allocate(realModeCoefficients(2*Ntheta))
      allocate(imagModeCoefficients(2*Ntheta))
      allocate(indices(thetaGridSize))
      call fft(2*Ntheta, thetaPartOfRHSOnThetaGrid, realModeCoefficients, imagModeCoefficients);
      realModeCoefficients = realModeCoefficients/(Ntheta + Ntheta)
      imagModeCoefficients = imagModeCoefficients/(Ntheta + Ntheta)
      thetaPartOfRHSModeCoefficients = (/ realModeCoefficients(1:Ntheta), imagModeCoefficients(2:Ntheta) /)
      do ix=1,Nx
        L=0
        indices = (ix-1)*Nxi*thetaGridSize + L*thetaGridSize + [(j, j=1,thetaGridSize)] - 1
        call VecSetValues(rhs, thetaGridSize, indices, &
          & ((oneHalf * 4)/3)*xPartOfRHS(ix)*thetaPartOfRHSModeCoefficients, ADD_VALUES, ierr)

        L=2
        indices = (ix-1)*Nxi*thetaGridSize + L*thetaGridSize + [(j, j=1,thetaGridSize)] - 1
        call VecSetValues(rhs, thetaGridSize, indices, &
          & oneHalf * (two/3)*xPartOfRHS(ix)*thetaPartOfRHSModeCoefficients, ADD_VALUES, ierr)
      end do 
      CHKERRQ(ierr)
      deallocate(indices)
      
      call VecAssemblyBegin(rhs, ierr)
      call VecAssemblyEnd(rhs, ierr)
      
      
      ! *********************************************************
      ! Build some small matrices and vectors needed for the main matrix:
      ! *********************************************************
      allocate(ddtheta(thetaGridSize, thetaGridSize))
      ddtheta=0
      do i=1,Ntheta-1
        ddtheta(Ntheta+i, i+1)=i
        ddtheta(i+1, Ntheta+i)=-i
      end do
      

      allocate(thetaPartOfMirrorTermOnThetaGrid(2*Ntheta))
      allocate(thetaPartMatrixForCollisionTerm(thetaGridSize, thetaGridSize))
      allocate(thetaPartMatrixForStreamingTerm(thetaGridSize, thetaGridSize))
      select case (whereToPutbDotGradTheta)
        case (0)
          ! The usual way I've done it:
          ! 1/(b dot grad theta) multiplies collision term,
          ! simple streaming term
          thetaPartOfMirrorTermOnThetaGrid = -oneHalf * dbdthetas / bs
          call generateThetaMultiplicationMatrix(oneOverqRbDotGradThetas, threshhold, thetaPartMatrixForCollisionTerm)
          thetaPartMatrixForStreamingTerm = ddtheta
          
          predictedNNZForEachRowOfPreconditioner = 2*thetaGridSize + thetaGridSize
          !predictedNNZForEachRowOfDifferenceMatrix = Nx*thetaGridSize
          predictedNNZForEachRowOfTotalMatrix = 2*thetaGridSize + Nx*thetaGridSize
        case (1)
          ! (b dot grad theta) multiplies streaming and mirror terms,
          ! simple collision term
          thetaPartOfMirrorTermOnThetaGrid = -oneHalf * dbdthetas / (bs * oneOverqRbDotGradThetas)
          
          call generateThetaMultiplicationMatrix(1/(oneOverqRbDotGradThetas), threshhold, thetaPartMatrixForStreamingTerm)
          thetaPartMatrixForStreamingTerm = matmul(thetaPartMatrixForStreamingTerm, ddtheta)
            
          predictedNNZForEachRowOfPreconditioner = 2*thetaGridSize + 1
          predictedNNZForEachRowOfTotalMatrix = 2*thetaGridSize + Nx
        case default
          print *,"Error! invalid whereToPutbDotGradTheta"
          stop
      end select
      allocate(thetaPartMatrixForMirrorTerm(thetaGridSize, thetaGridSize))
      call generateThetaMultiplicationMatrix(thetaPartOfMirrorTermOnThetaGrid, threshhold, thetaPartMatrixForMirrorTerm)
      
      ! PETSc and Fortran are opposite in their logical representation of the rows and columns in matrices,
      ! so we must take the transpose of the Fortran matrices before adding them to the PETSc matrix:
      thetaPartMatrixForMirrorTerm = transpose(thetaPartMatrixForMirrorTerm)
      thetaPartMatrixForStreamingTerm = transpose(thetaPartMatrixForStreamingTerm)
      

      if (useIterativeSolver) then
        call MatCreateSeqAIJ(PETSC_COMM_SELF, matrixSize, matrixSize, predictedNNZForEachRowOfPreconditioner, &
          PETSC_NULL_INTEGER, preconditionerMatrix, ierr)
        call MatCreateSeqAIJ(PETSC_COMM_SELF, matrixSize, matrixSize, predictedNNZForEachRowOfTotalMatrix, &
          PETSC_NULL_INTEGER, matrix, ierr)
      else
        call MatCreateSeqAIJ(PETSC_COMM_SELF, matrixSize, matrixSize, predictedNNZForEachRowOfTotalMatrix, &
          PETSC_NULL_INTEGER, preconditionerMatrix, ierr)
      end if
      CHKERRQ(ierr)
      
      allocate(rowIndices(thetaGridSize))
      allocate(colIndices(thetaGridSize))
      do ix=1,Nx
        do L=0,Nxi-1
          rowIndices = (ix-1)*Nxi*thetaGridSize + L*thetaGridSize + [(i, i=1,thetaGridSize)] - 1
          if (L < Nxi-1) then
            ! Super-diagonal in L:
            colIndices = rowIndices + thetaGridSize
            
            ! Add streaming term:
            call MatSetValues(preconditionerMatrix, thetaGridSize, rowIndices, thetaGridSize, colIndices, &
              thetaPartMatrixForStreamingTerm*x(ix)*(L+1)/(two*L+3), ADD_VALUES, ierr)
              
            ! Add mirror term:
            call MatSetValues(preconditionerMatrix, thetaGridSize, rowIndices, thetaGridSize, colIndices, &
              thetaPartMatrixForMirrorTerm*x(ix)*(L+1)*(L+2)/(two*L+3), ADD_VALUES, ierr)
          end if
           if (L>0) then
             ! Sub-diagonal in L:
             colIndices = rowIndices - thetaGridSize;
             
             ! Streaming term
             call MatSetValues(preconditionerMatrix, thetaGridSize, rowIndices, thetaGridSize, colIndices, &
               thetaPartMatrixForStreamingTerm*x(ix)*L/(two*L-1), ADD_VALUES, ierr)
                         
             ! Mirror term
             call MatSetValues(preconditionerMatrix, thetaGridSize, rowIndices, thetaGridSize, colIndices, &
               -thetaPartMatrixForMirrorTerm*x(ix)*(L-1)*L/(two*L-1), ADD_VALUES, ierr)
             
           end if
        end do
      end do
      deallocate(rowIndices)
      deallocate(colIndices)
      
      allocate(M21(NxPotentials, Nx))
      allocate(M32(NxPotentials, NxPotentials))
      allocate(M22BackslashM21(NxPotentials, Nx))
      allocate(M33BackslashM32(NxPotentials, NxPotentials))
      allocate(LaplacianTimesX2WithoutL(NxPotentials, NxPotentials))
      allocate(diagonalOfKWithoutThetaPart(Nx))
      M32 = zero
      M21 = 4*pi*regridPolynomialToUniform
      do i=2,NxPotentials-1
        M21(i,:) = M21(i,:)*xPotentials(i)*xPotentials(i)
        M32(i,i) = -2*xPotentials(i)*xPotentials(i)
      end do
      M21(1,:)=zero
      M21(NxPotentials,:)=zero
      M32(1,:)=zero
      M32(NxPotentials,:)=zero
      do i=1,NxPotentials
        LaplacianTimesX2WithoutL(i,:) = xPotentials(i)*xPotentials(i)*d2dx2Potentials(i,:) &
          + 2 * xPotentials(i) * ddxPotentials(i,:)
      end do
            
      ! Boundary condition suggested by the referee, both uncorrected and corrected:
      !M32(NxPotentials, NxPotentials) = one
      
      allocate(erfs(Nx))
      allocate(x3(Nx))
      allocate(expx2(Nx))
      allocate(Psi(Nx))
      allocate(nuD(Nx))
      allocate(PsiPrime(Nx))
      ! erf is vectorized in gfortran but not pathscale.
      do i=1,Nx
        erfs(i) = erf(x(i))
      end do
      x3 = x*x2
      expx2 = exp(-x2)
      Psi = (erfs - 2/sqrtpi*x*expx2) / (2*x2)
      nuD = 3*sqrtpi/4*(erfs - Psi) / x3
      PsiPrime = (-erfs + 2/sqrtpi*x*(1+x2)*expx2) / x3
      

      allocate(xPartOfCECD(Nx,Nx))
      allocate(M12IncludingX0(Nx,NxPotentials))
      allocate(M13IncludingX0(Nx,NxPotentials))
      allocate(M12(Nx,NxPotentials))
      allocate(M13(Nx,NxPotentials))
      allocate(M22(NxPotentials,NxPotentials))
      allocate(M33(NxPotentials,NxPotentials))
      M12IncludingX0 = nuPrime * 3/(2*pi)*regridUniformToPolynomial
      M13IncludingX0 = -nuPrime * 3/(2*pi) * matmul(regridUniformToPolynomial, d2dx2Potentials)
      do i=1,Nx
        xPartOfCECD(i,:) = (3*sqrtpi/4)*((Psi(i)/x(i))*d2dx2(i,:)  &
          + (PsiPrime(i)*x(i)  + Psi(i) + 2*Psi(i)*x2(i))/x2(i) * ddx(i,:))
        xPartOfCECD(i,i) = xPartOfCECD(i,i) + (3*sqrtpi/4)*(2*PsiPrime(i) + 4*Psi(i)/x(i)) + 3*expx2(i)
        M12IncludingX0(i,:) = M12IncludingX0(i,:) * expx2(i)
        M13IncludingX0(i,:) = M13IncludingX0(i,:) * x2(i) * expx2(i)
      end do

      
      allocate(M11(Nx,Nx))
      allocate(KWithoutThetaPart(Nx,Nx))
      allocate(IPIV(NxPotentials))
      allocate(rowIndices(Nx))
      allocate(colIndices(Nx))
      allocate(indices(Nx))
      do L=0, Nxi-1
        ! Build M11
        do i=1,Nx
          M11(i,:) = -nuPrime * xPartOfCECD(i,:)
          M11(i,i) = M11(i,i) - nuPrime * (-oneHalf*nuD(i)*L*(L+1))
        end do
        
        if (L < NL) then
          ! Add Rosenbluth potential stuff
          M13 = M13IncludingX0
          M12 = M12IncludingX0
          
          M22 = LaplacianTimesX2WithoutL + zero
          do i=1,NxPotentials
            M22(i,i) = M22(i,i) - L*(L+1)
          end do
          
          ! Add Dirichlet or Neumann boundary condition for potentials at x=0:
          if (L==0) then
            M22(1,:)=ddxPotentials(1,:)
          else
            M22(1,:) = 0
            M22(1,1) = 1
            M12(:,1) = 0
            M13(:,1) = 0
          end if
          M33 = M22;
                    
          ! Add Robin boundary condition for potentials at x=xMax:
          M22(NxPotentials,:) = xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
          M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1
          
          ! My original boundary condition for G:
          M33(NxPotentials,:) = xMaxNotTooSmall*xMaxNotTooSmall*d2dx2Potentials(NxPotentials,:) &
            + (2*L+1)*xMaxNotTooSmall*ddxPotentials(NxPotentials,:)
          M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1)
          
          ! Boundary condition suggested by the referee, corrected by me:
          !M33(NxPotentials,:)=zero
          !M33(NxPotentials, NxPotentials) = -(one-4*L)/(2*xMaxNotTooSmall*xMaxNotTooSmall)
          
          ! Boundary condition suggested by the referee, incorrect:
          !M33(NxPotentials,:)=zero
          !M33(NxPotentials, NxPotentials) = -(L*L-L+one)/(xMaxNotTooSmall*xMaxNotTooSmall)
          
          
          if (L /= 0) then
            M22(NxPotentials,1)=0
            M33(NxPotentials,1)=0
          end if


          
          ! Call LAPACK subroutine DGESV to solve a linear system
          ! Note: this subroutine changes M22 and M33!
          M22BackslashM21 = M21  ! This will be overwritten by LAPACK.
#if defined(PETSC_USE_REAL_SINGLE)
          call SGESV(NxPotentials, Nx, M22, NxPotentials, IPIV, M22BackslashM21, NxPotentials, LAPACKInfo)
#else
          call DGESV(NxPotentials, Nx, M22, NxPotentials, IPIV, M22BackslashM21, NxPotentials, LAPACKInfo)
#endif
          if (LAPACKInfo /= 0) then
            print *, "Error in LAPACK call: info = ", LAPACKInfo
            stop
          end if
          M33BackslashM32 = M32  ! This will be overwritten by LAPACK.
#if defined(PETSC_USE_REAL_SINGLE)
          call SGESV(NxPotentials, NxPotentials, M33, NxPotentials, IPIV, M33BackslashM32, NxPotentials, LAPACKInfo)
#else
          call DGESV(NxPotentials, NxPotentials, M33, NxPotentials, IPIV, M33BackslashM32, NxPotentials, LAPACKInfo)
#endif
          if (LAPACKInfo /= 0) then
            print *, "Error in LAPACK call: info = ", LAPACKInfo
            stop
          end if
        
        
          !KWithoutThetaPart = M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21);
          KWithoutThetaPart = M11 - matmul(M12 - matmul(M13, M33BackslashM32), M22BackslashM21)
        else
          KWithoutThetaPart = M11
        end if
        
        ! PETSc and Fortran use row-major vs column-major:
        KWithoutThetaPart = transpose(KWithoutThetaPart)
        
        if (useIterativeSolver) then
          do i=1,Nx
            diagonalOfKWithoutThetaPart(i) = KWithoutThetaPart(i,i) &
              + nuPrime * multipleOfNuPrimeToAddOnDiagonalOfPreconditioner
            KWithoutThetaPart(i,i) = -nuPrime * multipleOfNuPrimeToAddOnDiagonalOfPreconditioner
          end do
        end if
        
        select case (whereToPutbDotGradTheta)
          case (0)
            ! The usual way I've done it:
            ! 1/(b dot grad theta) multiplies collision term.
            if (useIterativeSolver) then
              do ithetaRow=1,thetaGridSize
                do ithetaCol=1,thetaGridSize
                  if (thetaPartMatrixForCollisionTerm(ithetaRow,ithetaCol) /= zero) then
                    rowIndices = [(i, i=0,Nx-1)]*Nxi*thetaGridSize + L*thetaGridSize + ithetaRow - 1
                    colIndices = [(i, i=0,Nx-1)]*Nxi*thetaGridSize + L*thetaGridSize + ithetaCol - 1
                    call MatSetValues(matrix, Nx, rowIndices, Nx, colIndices, &
                      thetaPartMatrixForCollisionTerm(ithetaRow,ithetaCol) * KWithoutThetaPart, ADD_VALUES, ierr)
                    do ix=1,Nx
                      call MatSetValues(preconditionerMatrix, 1, rowIndices(ix), 1, colIndices(ix), &
                        thetaPartMatrixForCollisionTerm(ithetaRow,ithetaCol) * diagonalOfKWithoutThetaPart(ix), &
                          ADD_VALUES, ierr)
                    end do
                  end if
                end do
              end do
            else
              do ithetaRow=1,thetaGridSize
                do ithetaCol=1,thetaGridSize
                  if (thetaPartMatrixForCollisionTerm(ithetaRow,ithetaCol) /= zero) then
                    rowIndices = [(i, i=0,Nx-1)]*Nxi*thetaGridSize + L*thetaGridSize + ithetaRow - 1
                    colIndices = [(i, i=0,Nx-1)]*Nxi*thetaGridSize + L*thetaGridSize + ithetaCol - 1
                    call MatSetValues(preconditionerMatrix, Nx, rowIndices, Nx, colIndices, &
                      thetaPartMatrixForCollisionTerm(ithetaRow,ithetaCol) * KWithoutThetaPart, ADD_VALUES, ierr)
                    CHKERRQ(ierr)
                  end if
                end do
              end do
            end if
          case (1)
            ! Collision term is diagonal in theta
            if (useIterativeSolver) then
              do itheta=1,thetaGridSize
                indices = [(i, i=0,Nx-1)]*Nxi*thetaGridSize + L*thetaGridSize + itheta - 1
                call MatSetValues(matrix, Nx, indices, Nx, indices, KWithoutThetaPart, ADD_VALUES, ierr)
                do ix=1,Nx
                  call MatSetValues(preconditionerMatrix, 1, indices(ix), 1, indices(ix), &
                    diagonalOfKWithoutThetaPart(ix), ADD_VALUES, ierr)
                end do
                CHKERRQ(ierr)
              end do
            else
              do itheta=1,thetaGridSize
                indices = [(i, i=0,Nx-1)]*Nxi*thetaGridSize + L*thetaGridSize + itheta - 1
                call MatSetValues(preconditionerMatrix, Nx, indices, Nx, indices, KWithoutThetaPart, ADD_VALUES, ierr)
                CHKERRQ(ierr)
              end do
            end if
        end select
        CHKERRQ(ierr)
      end do

      call PetscGetTime(time2, ierr)
      print *,"Time to pre-assemble matrix: ", time2-time1, " seconds."
      call PetscGetTime(time1, ierr)
      
      
      call MatAssemblyBegin(preconditionerMatrix, MAT_FINAL_ASSEMBLY, ierr)
      if (useIterativeSolver) then
        call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
      end if
      CHKERRQ(ierr)
      call MatAssemblyEnd(preconditionerMatrix, MAT_FINAL_ASSEMBLY, ierr)
      CHKERRQ(ierr)
      if (useIterativeSolver) then
        call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
        ! Add preconditioner to matrix and store the result in matrix
!        call MatAXPY(matrix, 1d+0, preconditionerMatrix, SUBSET_NONZERO_PATTERN, ierr)
      end if
      
      call PetscGetTime(time2, ierr)
      print *,"Time to assemble matrices: ", time2-time1, " seconds."
      call PetscGetTime(time1, ierr)
      
      if (useIterativeSolver) then
        call MatAXPY(matrix, one, preconditionerMatrix, DIFFERENT_NONZERO_PATTERN, ierr)
        call PetscGetTime(time2, ierr)
        print *,"Time to add matrices: ", time2-time1, " seconds."
        call PetscGetTime(time1, ierr)
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !! 
      !!  Solve the linear system:
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      
      call KSPCreate(PETSC_COMM_WORLD, KSPInstance, ierr)
      CHKERRQ(ierr)
      
      if (useIterativeSolver) then
        call KSPSetOperators(KSPInstance, matrix, preconditionerMatrix, SAME_PRECONDITIONER, ierr)
        CHKERRQ(ierr)
        call KSPGetPC(KSPInstance, preconditionerContext, ierr)
        CHKERRQ(ierr)
        call PCSetType(preconditionerContext, PCLU, ierr)
        CHKERRQ(ierr)
        call KSPSetType(KSPInstance, KSPGMRES, ierr)
        CHKERRQ(ierr)
!        call KSPSetTolerances(KSPInstance, 1d-5, PETSC_DEFAULT_DOUBLE_PRECISION, &
!          PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, ierr)
        CHKERRQ(ierr)
        call KSPSetFromOptions(KSPInstance, ierr)
        CHKERRQ(ierr)
      else
      ! Direct solver:
        call KSPSetOperators(KSPInstance, preconditionerMatrix, preconditionerMatrix, SAME_PRECONDITIONER, ierr)
        CHKERRQ(ierr)
        call KSPGetPC(KSPInstance, preconditionerContext, ierr)
        CHKERRQ(ierr)
        call PCSetType(preconditionerContext, PCLU, ierr)
        CHKERRQ(ierr)
        call KSPSetType(KSPInstance, KSPPREONLY, ierr)
        CHKERRQ(ierr)
        call KSPSetFromOptions(KSPInstance, ierr)
!        call KSPSetTolerances(KSPInstance, 1d-5, PETSC_DEFAULT_REAL, &
!          PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr)
        CHKERRQ(ierr)
      end if

      call VecDuplicate(rhs, soln, ierr)
      CHKERRQ(ierr)
      print *,"Beginning solve ..."
      call KSPSolve(KSPInstance, rhs, soln, ierr)
      CHKERRQ(ierr)
      
      call PetscGetTime(time2, ierr)
      print *,"Done.  Time to solve: ", time2-time1, " seconds."
      call PetscGetTime(time1, ierr)

      if (useIterativeSolver) then
        call KSPGetConvergedReason(KSPInstance, reason, ierr)
        if (reason>0) then
          print *,"Converged!  KSPConvergedReason = ", reason
          didItConverge = 1
        else
          print *,"Did not converge :(   KSPConvergedReason = ", reason
          didItConverge = -1
        end if
      else
        didItConverge = 1
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      !! 
      !!  Calculate moments of the solution:
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      VPrime = dot_product(oneOverqRbDotGradThetas/bs, thetaWeights)
      
      allocate(particleFluxBeforeThetaIntegral(thetaGridSize))
      allocate(qBeforeThetaIntegral(thetaGridSize))
      allocate(flowDividedByB(thetaGridSize))
      allocate(density(thetaGridSize))
      deallocate(indices)
      allocate(indices(Nx))
      allocate(fSlice(Nx))
      
      allocate(densityIntegralWeight(Nx))
      allocate(kIntegralWeight(Nx))
      allocate(particleFluxIntegralWeight(Nx))
      allocate(qIntegralWeight(Nx))
      
      densityIntegralWeight = 4*pi*x*x
      kIntegralWeight = x*x*x
      particleFluxIntegralWeight = x*kIntegralWeight
      qIntegralWeight = particleFluxIntegralWeight * x * x
      
      call VecGetArrayF90(soln, solnArray, ierr)
      CHKERRQ(ierr)
      
      do itheta=1,thetaGridSize
        L = 0
        indices = [((ix-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta, ix=1,Nx)]
        fSlice = solnArray(indices)
        qBeforeThetaIntegral(itheta) = (8 / three) * dot_product(xWeights, qIntegralWeight*fSlice)
        particleFluxBeforeThetaIntegral(itheta) = (8 / three) * dot_product(xWeights, particleFluxIntegralWeight*fSlice)
        density(itheta) = dot_product(xWeights, densityIntegralWeight*fSlice)
        
        L = 1
        indices = [((ix-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta, ix=1,Nx)]
        fSlice = solnArray(indices)
        flowDividedByB(itheta) = dot_product(xWeights, kIntegralWeight*fSlice)
        
        L = 2
        indices = [((ix-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta, ix=1,Nx)]
        fSlice = solnArray(indices)
        qBeforeThetaIntegral(itheta) = qBeforeThetaIntegral(itheta) + (four / 15) * dot_product(xWeights, qIntegralWeight*fSlice)
        particleFluxBeforeThetaIntegral(itheta) = particleFluxBeforeThetaIntegral(itheta) + &
          (four / 15) * dot_product(xWeights, particleFluxIntegralWeight*fSlice)
      end do
      
      call VecRestoreArrayF90(soln, solnArray, ierr)
      CHKERRQ(ierr)
      
      allocate(thetaMatrixForFlow(thetaGridSize, thetaGridSize))
      allocate(thetaMatrixForFluxes(thetaGridSize, thetaGridSize))
      call generateThetaMultiplicationMatrix(1/bs, threshhold, thetaMatrixForFlow)
      call generateThetaMultiplicationMatrix(dbdthetas/(bs*bs*bs), threshhold, thetaMatrixForFluxes)
      flowDividedByB = (two / 3)*4/sqrtpi * matmul(thetaMatrixForFlow, flowDividedByB)
      qBeforeThetaIntegral = matmul(thetaMatrixForFluxes, qBeforeThetaIntegral)
      particleFluxBeforeThetaIntegral = matmul(thetaMatrixForFluxes, particleFluxBeforeThetaIntegral)
      
      k = flowDividedByB(1)
      q = sqrt(2/(pi*epsil))/(nuPrime*VPrime) * 2*pi * qBeforeThetaIntegral(1);
      particleFlux = sqrt(2/(pi*epsil))/(nuPrime*VPrime) * 2*pi * particleFluxBeforeThetaIntegral(1);

      print *,"Max non-constant coefficient in k: (should be nearly 0):", maxval(abs(flowDividedByB(2:thetaGridSize)))
      
      ! *********************************************************
      ! Create a PETSc viewer to record output
      ! *********************************************************
      if (saveMatlabOutput) then
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD, &
          & MatlabOutputFilename,&
          & MatlabOutput, ierr)
        CHKERRQ(ierr)
        call PetscViewerSetFormat(MatlabOutput, PETSC_VIEWER_ASCII_MATLAB, ierr)
        CHKERRQ(ierr)
        
        call PetscObjectSetName(rhs, "rhs", ierr)
        CHKERRQ(ierr)
        call VecView(rhs, MatlabOutput, ierr)
        CHKERRQ(ierr)
        call PetscObjectSetName(soln, "soln", ierr)
        CHKERRQ(ierr)
        call VecView(soln, MatlabOutput, ierr)
        CHKERRQ(ierr)

        call PetscObjectSetName(preconditionerMatrix, "preconditionerMatrix", ierr)
        CHKERRQ(ierr)
        call MatView(preconditionerMatrix, MatlabOutput, ierr)
        CHKERRQ(ierr)
        if (useIterativeSolver) then
          call PetscObjectSetName(matrix, "matrix", ierr)
          CHKERRQ(ierr)
          call MatView(matrix, MatlabOutput, ierr)
          CHKERRQ(ierr)
        end if
          
        call PetscGetTime(time2, ierr)
        print *,"Time to write output: ", time2-time1, " seconds."
        call PetscGetTime(time1, ierr)

        call PetscViewerDestroy(MatlabOutput, ierr)
        CHKERRQ(ierr)
      end if
      
      call VecDestroy(rhs, ierr)
      call VecDestroy(soln, ierr)
      CHKERRQ(ierr)
      call MatDestroy(preconditionerMatrix, ierr)
      if (useIterativeSolver) then
        call MatDestroy(matrix, ierr)
      end if
      call KSPDestroy(KSPInstance, ierr)
      CHKERRQ(ierr)
      

      call PetscGetTime(time2, ierr)
      elapsedTime = time2 - startTime

    end subroutine
  end module
