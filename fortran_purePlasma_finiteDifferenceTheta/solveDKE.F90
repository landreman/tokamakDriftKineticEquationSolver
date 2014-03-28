! For PETSc versions prior to 3.3, the line below should be un-commented:
!#define MatCreateAIJ MatCreateMPIAIJ

! For PETSc versions 3.4 and later, the line below should be un-commented:
!#define PetscGetTime PetscTime

  module solveDKEModule
    use geometry
    use parameters
    use poldif
    use polyInterp
    use makeXGridModule
    use petscksp
    use petscdmda
    implicit none
#include <finclude/petsckspdef.h>
#include <finclude/petscdmdadef.h>

    contains

      subroutine solveDKE()

        PetscErrorCode :: ierr
        Vec :: rhs, soln, localSoln
        Mat :: matrix, preconditionerMatrix
        PetscViewer MatlabOutput
        PetscScalar, dimension(:), allocatable :: theta, thetaWeights, bs, dbdthetas, oneOverqRbDotGradThetas
        PetscScalar, dimension(:,:), allocatable :: ddtheta, d2dtheta2, thetaPartMatrix, localddtheta
        integer :: i, j, ix, itheta, L, NxPotentials, matrixSize
        integer :: ithetaRow, ithetaCol, scheme
        PetscScalar, dimension(:), allocatable :: x, xWeights, xPotentials, xWeightsPotentials
        PetscScalar, dimension(:), allocatable :: x2, xPartOfRHS, thetaPartOfMirrorTerm
        PetscScalar, dimension(:,:), allocatable :: ddx, d2dx2, ddxPotentials, d2dx2Potentials
        PetscScalar, dimension(:,:), allocatable :: regridPolynomialToUniform, regridUniformToPolynomial
        PetscScalar :: VPrime, dtheta, xMaxNotTooSmall
        PetscScalar, dimension(:), allocatable :: thetaPartOfRHS, localThetaPartOfRHS
        integer, dimension(:), allocatable :: indices, rowIndices, colIndices
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

        character :: trans='n'
        PetscLogDouble :: time1, time2, startTime
        KSP :: KSPInstance
        PC :: preconditionerContext
        KSPConvergedReason :: reason
        PetscScalar, pointer :: solnArray(:)
        DM :: myDM
        integer :: ithetaMin, ithetaMax, localNtheta
        VecScatter :: VecScatterContext

        call PetscGetTime(time1, ierr)
        startTime = time1

        if ((.not. isAParallelDirectSolverInstalled) .and. (numProcsInSubComm > 1)) then
           if (masterProcInSubComm) then
              print *,"Error! Neither mumps nor superlu_dist appears to be installed,"
              print *," yet you have asked for a matrix to be distributed across processsors."
           end if
           stop
        end if

        ! Assign a range of theta indices to each processor.
        ! This is done by creating a PETSc DM that is not actually used for anything else.
        call DMDACreate1d(MPIComm, DMDA_BOUNDARY_NONE, Ntheta, 1, 0, PETSC_NULL_INTEGER, myDM, ierr)
        call DMDAGetCorners(myDM, ithetaMin, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
             localNtheta, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
        ! Switch to 1-based indices:
        ithetaMin = ithetaMin + 1
        ithetaMax = ithetaMin+localNtheta-1
        CHKERRQ(ierr)
        print *,"[",myCommunicatorIndex,"] Processor ",myRank," owns theta indices ",ithetaMin," to ",ithetaMax

        ! Each processor is responsible for building the rows of the matrix and rhs corresponding
        ! to its ithetaMin:ithetaMax, and each processor is resposible for all columns of the matrix.

        matrixSize = Ntheta * Nxi * Nx

        select case (thetaDerivativeScheme)
        case (0)
           scheme = 20
        case (1)
           scheme = 0
        case (2)
           scheme = 10
        case default
           if (masterProcInSubComm) then
              print *,"[",myCommunicatorIndex,"] Error! Invalid setting for thetaDerivativeScheme"
           end if
           stop
        end select

        allocate(theta(Ntheta))
        allocate(thetaWeights(Ntheta))
        allocate(ddtheta(Ntheta,Ntheta))
        allocate(localddtheta(localNtheta,Ntheta))
        allocate(d2dtheta2(Ntheta,Ntheta))
        call uniformDiffMatrices(Ntheta, 0, two*pi, scheme, theta, thetaWeights, ddtheta, d2dtheta2)
        localddtheta = ddtheta(ithetaMin:ithetaMax,:)

        allocate(bs(Ntheta))
        allocate(dbdthetas(Ntheta))
        allocate(oneOverqRbDotGradThetas(Ntheta))
        call computeBs(theta, bs)
        call computedBdthetas(theta, dbdthetas)
        call computeOneOverqRbDotGradThetas(theta, oneOverqRbDotGradThetas)


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
        call uniformDiffMatrices(NxPotentials, zero, xMaxNotTooSmall, 12, xPotentials, &
             xWeightsPotentials, ddxPotentials, d2dx2Potentials)

        allocate(regridPolynomialToUniform(NxPotentials, Nx))
        call makePolynomialInterpolationMatrix(x,xPotentials,exp(-x*x),exp(-xPotentials*xPotentials), regridPolynomialToUniform)
        allocate(regridUniformToPolynomial(Nx,NxPotentials))
        call interpolationMatrix(NxPotentials, Nx, xPotentials, x, regridUniformToPolynomial, -1, 0)


        ! *********************************************************
        ! Create right-hand side
        ! *********************************************************
        call VecCreateMPI(MPIComm, PETSC_DECIDE, matrixSize, rhs, ierr);
        CHKERRQ(ierr)

        allocate(x2(Nx))
        allocate(xPartOfRHS(Nx))
        allocate(thetaPartOfRHS(Ntheta))
        allocate(localThetaPartOfRHS(localNtheta))
        x2=x*x
        xPartOfRHS = (x2 - 5/two)*x2*exp(-x2)
        thetaPartOfRHS = dbdthetas / (bs*bs)
        localThetaPartOfRHS = thetaPartOfRHS(ithetaMin:ithetaMax)
        allocate(indices(localNtheta))

        do ix=1,Nx
           L=0
           indices = (ix-1)*Nxi*Ntheta + L*Ntheta + [( j, j=ithetaMin,ithetaMax )] - 1
           !           indices = (ix-1)*Nxi*Ntheta + L*Ntheta + [( j, j=1,Ntheta )] - 1
           call VecSetValues(rhs, localNtheta, indices, &
                & ((oneHalf * 4)/3)*xPartOfRHS(ix)*localThetaPartOfRHS, ADD_VALUES, ierr)

           L=2
           indices = (ix-1)*Nxi*Ntheta + L*Ntheta + [( j, j=ithetaMin,ithetaMax )] - 1
           call VecSetValues(rhs, localNtheta, indices, &
                & oneHalf * (two/3)*xPartOfRHS(ix)*localThetaPartOfRHS, ADD_VALUES, ierr)
        end do
        CHKERRQ(ierr)
        deallocate(indices)

        call VecAssemblyBegin(rhs, ierr)
        call VecAssemblyEnd(rhs, ierr)

        ! *********************************************************
        ! Done building the right-hand side.
        ! Now build the main matrix:
        ! *********************************************************


        allocate(thetaPartOfMirrorTerm(Ntheta))
        allocate(thetaPartMatrix(localNtheta,Ntheta))
        thetaPartOfMirrorTerm = -oneHalf * dbdthetas / bs

        !predictedNNZForEachRowOfPreconditioner = ceiling((2*Ntheta + 2d+0 + 1) / numProcessors)
        !predictedNNZForEachRowOfTotalMatrix = ceiling((2*Ntheta + 2d+0 + Nx) / numProcessors)
        predictedNNZForEachRowOfPreconditioner = (2*Ntheta + 2 + 1)
        predictedNNZForEachRowOfTotalMatrix = (2*Ntheta + 2 + Nx)


        if (useIterativeSolver) then
           call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
                predictedNNZForEachRowOfPreconditioner, PETSC_NULL_INTEGER, &
                predictedNNZForEachRowOfPreconditioner, PETSC_NULL_INTEGER, &
                preconditionerMatrix, ierr)
           call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
                predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER, &
                predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER, &
                matrix, ierr)
        else
           call MatCreateAIJ(MPIComm, PETSC_DECIDE, PETSC_DECIDE, matrixSize, matrixSize, &
                predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER, &
                predictedNNZForEachRowOfTotalMatrix, PETSC_NULL_INTEGER, &
                preconditionerMatrix, ierr)
        end if
!!$        if (useIterativeSolver) then
!!$           call MatCreateSeqAIJ(PETSC_COMM_SELF, matrixSize, matrixSize, predictedNNZForEachRowOfPreconditioner, &
!!$                PETSC_NULL_INTEGER, preconditionerMatrix, ierr)
!!$           call MatCreateSeqAIJ(PETSC_COMM_SELF, matrixSize, matrixSize, predictedNNZForEachRowOfTotalMatrix, &
!!$                PETSC_NULL_INTEGER, matrix, ierr)
!!$        else
!!$           call MatCreateSeqAIJ(PETSC_COMM_SELF, matrixSize, matrixSize, predictedNNZForEachRowOfTotalMatrix, &
!!$                PETSC_NULL_INTEGER, preconditionerMatrix, ierr)
!!$        end if
        CHKERRQ(ierr)

        allocate(rowIndices(localNtheta))
        allocate(colIndices(Ntheta))
        do ix=1,Nx
           do L=0,Nxi-1
              rowIndices = (ix-1)*Nxi*Ntheta + L*Ntheta + [(i, i=ithetaMin,ithetaMax)] - 1
              if (L < Nxi-1) then
                 ! Super-diagonal in L:
                 colIndices = (ix-1)*Nxi*Ntheta + (L+1)*Ntheta + [(i, i=1,Ntheta)] - 1
                 ! colIndices = rowIndices + Ntheta

                 ! Add streaming term:
                 thetaPartMatrix = x(ix)*(L+1)/(two*L+3)*localddtheta

                 ! Add mirror term:
                 ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                 do itheta=ithetaMin,ithetaMax
                    thetaPartMatrix(itheta-ithetaMin+1,itheta) = x(ix)*(L+1)*(L+2)/(two*L+3)*thetaPartOfMirrorTerm(itheta)
                 end do

                 ! Put values in matrix, noting that Petsc uses a transposed format relative to Fortran
                 call MatSetValues(preconditionerMatrix, localNtheta, rowIndices, Ntheta, colIndices, &
                      transpose(thetaPartMatrix), ADD_VALUES, ierr)
              end if
              if (L>0) then
                 ! Sub-diagonal in L:
                 !colIndices = rowIndices - Ntheta
                 colIndices = (ix-1)*Nxi*Ntheta + (L-1)*Ntheta + [(i, i=1,Ntheta)] - 1

                 ! Streaming term
                 thetaPartMatrix = x(ix)*L/(two*L-1)*localddtheta

                 ! Mirror term
                 ! (I exploit the fact that d/dtheta has zeros on the diagonal)
                 do itheta=ithetaMin,ithetaMax
                    thetaPartMatrix(itheta-ithetaMin+1,itheta) = -x(ix)*L*(L-1)/(two*L-1)*thetaPartOfMirrorTerm(itheta)
                 end do

                 ! Put values in matrix, noting that Petsc uses a transposed format relative to Fortran
                 call MatSetValues(preconditionerMatrix, localNtheta, rowIndices, Ntheta, colIndices, &
                      transpose(thetaPartMatrix), ADD_VALUES, ierr)

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

           if (useIterativeSolver) then
              do itheta=ithetaMin,ithetaMax
                 indices = [(i, i=0,Nx-1)]*Nxi*Ntheta + L*Ntheta + itheta - 1
                 call MatSetValues(matrix, Nx, indices, Nx, indices, &
                      KWithoutThetaPart*oneOverqRbDotGradThetas(itheta), ADD_VALUES, ierr)
                 do ix=1,Nx
                    call MatSetValues(preconditionerMatrix, 1, indices(ix), 1, indices(ix), &
                         diagonalOfKWithoutThetaPart(ix)*oneOverqRbDotGradThetas(itheta), ADD_VALUES, ierr)
                 end do
                 CHKERRQ(ierr)
              end do
           else
              do itheta=ithetaMin,ithetaMax
                 indices = [(i, i=0,Nx-1)]*Nxi*Ntheta + L*Ntheta + itheta - 1
                 call MatSetValues(preconditionerMatrix, Nx, indices, Nx, indices, &
                      KWithoutThetaPart*oneOverqRbDotGradThetas(itheta), ADD_VALUES, ierr)
                 CHKERRQ(ierr)
              end do
           end if
           CHKERRQ(ierr)
        end do

        call PetscGetTime(time2, ierr)
        if (masterProcInSubComm) then
           print *,"[",myCommunicatorIndex,"] Time to pre-assemble matrix: ", time2-time1, " seconds."
        end if
        call PetscGetTime(time1, ierr)


        call MatAssemblyBegin(preconditionerMatrix, MAT_FINAL_ASSEMBLY, ierr)
        CHKERRQ(ierr)
        if (useIterativeSolver) then
           call MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY, ierr)
        end if
        CHKERRQ(ierr)
        CHKERRQ(ierr)
        call MatAssemblyEnd(preconditionerMatrix, MAT_FINAL_ASSEMBLY, ierr)
        CHKERRQ(ierr)
        if (useIterativeSolver) then
           call MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY, ierr)
           ! Add preconditioner to matrix and store the result in matrix
           !        call MatAXPY(matrix, 1d+0, preconditionerMatrix, SUBSET_NONZERO_PATTERN, ierr)
        end if

        call PetscGetTime(time2, ierr)
        if (masterProcInSubComm) then
           print *,"[",myCommunicatorIndex,"] Time to assemble matrices: ", time2-time1, " seconds."
        end if
        call PetscGetTime(time1, ierr)

        if (useIterativeSolver) then
           call MatAXPY(matrix, one, preconditionerMatrix, DIFFERENT_NONZERO_PATTERN, ierr)
           call PetscGetTime(time2, ierr)
           if (masterProcInSubComm) then
              print *,"[",myCommunicatorIndex,"] Time to add matrices: ", time2-time1, " seconds."
           end if
           call PetscGetTime(time1, ierr)
        end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! 
        !!  Solve the linear system:
        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call KSPCreate(MPIComm, KSPInstance, ierr)
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

        if (numProcsInSubComm > 1) then
           select case (whichParallelSolverToFactorPreconditioner)
           case (1)
              call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERMUMPS, ierr)
              if (masterProcInSubComm) then
                 print *,"[",myCommunicatorIndex,"] Using mumps to factorize the preconditioner."
              end if
           case (2)
              call PCFactorSetMatSolverPackage(preconditionerContext, MATSOLVERSUPERLU_DIST, ierr)
              if (masterProcInSubComm) then
                 print *,"[",myCommunicatorIndex,"] Using superlu_dist to factorize the preconditioner."
              end if
           case default
              if (masterProcInSubComm) then
                 print *,"Error: Invalid setting for whichParallelSolverToFactorPreconditioner"
                 stop
              end if
           end select
        else
           if (masterProcInSubComm) then
              print *,"[",myCommunicatorIndex,"] Using PETSc's serial sparse direct solver to factorize the preconditioner."
           end if
        end if

        call VecDuplicate(rhs, soln, ierr)
        CHKERRQ(ierr)
        if (masterProcInSubComm) then
           print *,"[",myCommunicatorIndex,"] Beginning solve ..."
        end if
        call KSPSolve(KSPInstance, rhs, soln, ierr)
        CHKERRQ(ierr)

        call PetscGetTime(time2, ierr)
        if (masterProcInSubComm) then
           print *,"[",myCommunicatorIndex,"] Done.  Time to solve: ", time2-time1, " seconds."
        end if
        call PetscGetTime(time1, ierr)

        if (useIterativeSolver) then
           call KSPGetConvergedReason(KSPInstance, reason, ierr)
           if (reason>0) then
              if (masterProcInSubComm) then
                 print *,"[",myCommunicatorIndex,"] Converged!  KSPConvergedReason = ", reason
              end if
              didItConverge = 1
           else
              if (masterProcInSubComm) then
                 print *,"[",myCommunicatorIndex,"] Did not converge :(   KSPConvergedReason = ", reason
              end if
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

        ! First, give each processor a copy of the entire solution vector:
        call VecScatterCreateToAll(soln, VecScatterContext, localSoln, ierr)
        CHKERRQ(ierr)
        call VecScatterBegin(VecScatterContext, soln, localSoln, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(VecScatterContext, soln, localSoln, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterDestroy(VecScatterContext, ierr)

        VPrime = dot_product(oneOverqRbDotGradThetas/bs, thetaWeights)

        allocate(particleFluxBeforeThetaIntegral(Ntheta))
        allocate(qBeforeThetaIntegral(Ntheta))
        allocate(flowDividedByB(Ntheta))
        allocate(density(Ntheta))
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

        call VecGetArrayF90(localSoln, solnArray, ierr)
        CHKERRQ(ierr)

        do itheta=1,Ntheta
           L = 0
           indices = [((ix-1)*Nxi*Ntheta + L*Ntheta + itheta, ix=1,Nx)]
           fSlice = solnArray(indices)
           qBeforeThetaIntegral(itheta) = (8 / three) * dot_product(xWeights, qIntegralWeight*fSlice)
           particleFluxBeforeThetaIntegral(itheta) = (8 / three) &
                * dot_product(xWeights, particleFluxIntegralWeight*fSlice)
           density(itheta) = dot_product(xWeights, densityIntegralWeight*fSlice)

           L = 1
           indices = [((ix-1)*Nxi*Ntheta + L*Ntheta + itheta, ix=1,Nx)]
           fSlice = solnArray(indices)
           flowDividedByB(itheta) = dot_product(xWeights, kIntegralWeight*fSlice)

           L = 2
           indices = [((ix-1)*Nxi*Ntheta + L*Ntheta + itheta, ix=1,Nx)]
           fSlice = solnArray(indices)
           qBeforeThetaIntegral(itheta) = qBeforeThetaIntegral(itheta) &
                + (four / 15) * dot_product(xWeights, qIntegralWeight*fSlice)
           particleFluxBeforeThetaIntegral(itheta) = particleFluxBeforeThetaIntegral(itheta) + &
                (four / 15) * dot_product(xWeights, particleFluxIntegralWeight*fSlice)
        end do

        call VecRestoreArrayF90(localSoln, solnArray, ierr)
        CHKERRQ(ierr)

        flowDividedByB = (two / 3)*4/sqrtpi * flowDividedByB / bs
        qBeforeThetaIntegral = qBeforeThetaIntegral * dbdthetas / (bs*bs*bs)
        particleFluxBeforeThetaIntegral = particleFluxBeforeThetaIntegral * dbdthetas / (bs*bs*bs)

        k = 1/(2*pi)*dot_product(thetaWeights, flowDividedByB)
        q = sqrt(2/(pi*epsil))/(nuPrime*VPrime) * dot_product(thetaWeights, qBeforeThetaIntegral);
        particleFlux = sqrt(2/(pi*epsil))/(nuPrime*VPrime) &
             * dot_product(thetaWeights, particleFluxBeforeThetaIntegral);

        if (masterProcInSubComm) then
           print *,"[",myCommunicatorIndex,"] Variation in k with theta: (should be nearly 0):", maxval(abs(flowDividedByB-k))
        end if

        ! *********************************************************
        ! Create a PETSc viewer to record output
        ! *********************************************************
        if (saveMatlabOutput) then
           call PetscViewerASCIIOpen(MPIComm, &
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
           if (masterProcInSubComm) then
              print *,"[",myCommunicatorIndex,"] Time to write output: ", time2-time1, " seconds."
           end if
           call PetscGetTime(time1, ierr)

           call PetscViewerDestroy(MatlabOutput, ierr)
           CHKERRQ(ierr)
        end if

        call VecDestroy(rhs, ierr)
        call VecDestroy(soln, ierr)
        call VecDestroy(localSoln, ierr)
        CHKERRQ(ierr)
        call MatDestroy(preconditionerMatrix, ierr)
        if (useIterativeSolver) then
           call MatDestroy(matrix, ierr)
        end if
        call KSPDestroy(KSPInstance, ierr)
        CHKERRQ(ierr)


        call PetscGetTime(time2, ierr)
        elapsedTime = time2 - startTime

      end subroutine solveDKE

      ! ------------------------------------------------------------------------

      subroutine chooseParallelDirectSolver()

        implicit none

        isAParallelDirectSolverInstalled = .false.

        if ((whichParallelSolverToFactorPreconditioner<1) .or. (whichParallelSolverToFactorPreconditioner>2)) then
           print *,"Error! Invalid setting for whichParallelSolverToFactorPreconditioner"
           stop
        end if

#ifdef PETSC_HAVE_MUMPS
        isAParallelDirectSolverInstalled = .true.
        if (masterProc) then
           print *,"mumps detected"
        end if
#else
        whichParallelSolverToFactorPreconditioner = 2
        if (masterProc) then
           print *,"mumps not detected"
        end if
#endif

#ifdef PETSC_HAVE_SUPERLU_DIST
        isAParallelDirectSolverInstalled = .true.
        if (masterProc) then
           print *,"superlu_dist detected"
        end if
#else
        if (masterProc) then
           print *,"superlu_dist not detected"
        end if
        if (whichParallelSolverToFactorPreconditioner==2) then
           whichParallelSolverToFactorPreconditioner = 1
        end if
#endif

      end subroutine chooseParallelDirectSolver

  end module
