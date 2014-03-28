  module geometry
    ! Subroutines and functions related to determining B(theta).
    
    use parameters

    implicit none
#include <finclude/petscsysdef.h>
    logical, private :: initializedYet = .false.
    PetscScalar :: Miller_A, Miller_x, Miller_QQ
    integer, parameter :: NThetaIntegral = 100
    
    contains

    subroutine initializeGeometry()
      implicit none
      integer :: i
      select case (geometryToUse)
        case (GEOMETRY_CIRCULAR_CONCENTRIC)
        case (GEOMETRY_MILLER)
          Miller_x = asin(Miller_delta)
          Miller_A = 1/epsil
          Miller_QQ = 0
          ! Integrate QQIntegrand from 0 to 2*pi.
          do i=1,NThetaIntegral
            Miller_QQ = Miller_QQ + QQIntegrand(2*pi*i/NThetaIntegral)
          end do
          Miller_QQ = Miller_kappa / (2*pi*Miller_A) * (Miller_QQ * 2*pi/NThetaIntegral)
          
        case (GEOMETRY_CIRCULAR_CONCENTRIC_BOOZER)
        case default
          print *,"Error! Invalid geometry."
          stop
      end select

      initializedYet = .true.
    end subroutine


    subroutine computeBs(thetas, bs)
      implicit none
      PetscScalar, intent(in) :: thetas(:)
      PetscScalar, intent(out) :: bs(:)
      integer :: i

      if (.not. initializedYet) then
          print *,"Error!  No geometry has been initialized yet."
          stop
      end if
      
      select case (geometryToUse)
        case (GEOMETRY_CIRCULAR_CONCENTRIC)
          bs = 1 / (1 + epsil * cos(thetas))
        case (GEOMETRY_MILLER)
          do i=1,size(thetas)
            bs(i) = sqrt(BPoloidal(thetas(i))**2 + 1./((RHat(thetas(i)))**2))
          end do
        case (GEOMETRY_CIRCULAR_CONCENTRIC_BOOZER)
          bs = 1 + epsil * cos(thetas)
        case default
          print *,"Error! Invalid geometry."
          stop
      end select
    end subroutine


    subroutine computedBdthetas(thetas, dBdthetas)
      implicit none
      PetscScalar, intent(in) :: thetas(:)
      PetscScalar, intent(out) :: dBdthetas(:)
      PetscScalar, dimension(:), allocatable :: temp
      PetscScalar, allocatable :: spectralDerivative(:,:), thetaFine(:), bs(:), dBdthetaFine(:)
      PetscScalar, allocatable :: d2dtheta2(:,:), weights(:)
      integer :: i, N, multipliedIndex
      integer, parameter :: dBdthetaResolutionMultiplier = 10
    
      if (.not. initializedYet) then
          print *,"Error!  No geometry has been initialized yet."
          stop
      end if
      
      select case (geometryToUse)
        case (GEOMETRY_CIRCULAR_CONCENTRIC)
          allocate(temp(size(thetas)))
          temp = 1 + epsil * cos(thetas)
          dBdthetas = epsil * sin(thetas) / (temp*temp)
          deallocate(temp)
        case (GEOMETRY_MILLER)
          ! It is not worth analytically differentiating the Miller formulae.
          ! Instead, just numerically differentiate b(theta).
          
          N=size(thetas)*dBdthetaResolutionMultiplier
          
          allocate(spectralDerivative(N,N))
          allocate(thetaFine(N))
          allocate(d2dtheta2(N,N))
          allocate(weights(N))
          allocate(bs(N))
          allocate(dBdthetaFine(N))
          call uniformDiffMatrices(N, zero, two*pi, 20, thetaFine, weights, spectralDerivative, d2dtheta2)

          call computeBs(thetaFine, bs)
          dBdthetaFine = matmul(spectralDerivative, bs)
          do i=1,size(thetas)
            multipliedIndex = (i-1)*dBdthetaResolutionMultiplier+1
            dBdthetas(i) = dBdthetaFine(multipliedIndex)
            if (abs(thetas(i) - thetaFine(multipliedIndex)) > 1e-10) then
              print *,"Error! The input theta array to computedBdthetas was not of the expected form."
              print *,"thetas:",thetas
              print *,"thetaFine:",thetaFine
              stop
            end if
          end do
          
          deallocate(spectralDerivative)
          deallocate(thetaFine)
          deallocate(bs)
          deallocate(dBdthetaFine)
        case (GEOMETRY_CIRCULAR_CONCENTRIC_BOOZER)
          dBdthetas = - epsil * sin(thetas)
        case default
          print *,"Error! Invalid geometry."
          stop
      end select
    end subroutine
    
    subroutine computeOneOverqRbDotGradThetas(thetas, oneOverqRbDotGradThetas)
      implicit none
      PetscScalar, intent(in) :: thetas(:)
      PetscScalar, intent(out) :: oneOverqRbDotGradThetas(:)
      integer :: i
      PetscScalar, allocatable :: bs(:)
      
      if (.not. initializedYet) then
          print *,"Error!  No geometry has been initialized yet."
          stop
      end if
      
      select case (geometryToUse)
        case (GEOMETRY_CIRCULAR_CONCENTRIC)
          oneOverqRbDotGradThetas = [(one, i=1,size(thetas))]
        case (GEOMETRY_MILLER)
          allocate(bs(size(thetas)))
          call computeBs(thetas,bs)
          do i=1,size(thetas)
            oneOverqRbDotGradThetas(i) = bs(i) / (Miller_q*BDotGradTheta(thetas(i)));
          end do
          deallocate(bs)
        case (GEOMETRY_CIRCULAR_CONCENTRIC_BOOZER)
          oneOverqRbDotGradThetas = 1 / (1 + epsil*cos(thetas))
        case default
          print *,"Error! Invalid geometry."
          stop
      end select
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Below are a set of functions needed for Miller geometry
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    function RHat(theta)
      implicit none
      PetscScalar :: theta, RHat
      RHat = 1 + (1/Miller_A)*cos(theta + Miller_x*sin(theta))
    end function

    function ZHat(theta)
      implicit none
      PetscScalar :: theta, ZHat
      ZHat = (Miller_kappa/Miller_A)*sin(theta)
    end function

    function QQIntegrand(theta)
      implicit none
      PetscScalar :: theta, QQIntegrand
      QQIntegrand = ((1+Miller_s_kappa)*sin(theta + Miller_x*sin(theta)) * (1+Miller_x*cos(theta)) * sin(theta) &
        + cos(theta) * (Miller_dRdr + cos(theta + Miller_x *sin(theta)) &
        - Miller_s_delta*sin(theta + Miller_x*sin(theta)) * sin(theta))) / RHat(theta)
    end function

    function BPoloidal(theta)
      implicit none
      PetscScalar :: theta, BPoloidal
      BPoloidal = Miller_QQ/(Miller_kappa*Miller_q)*sqrt((sin(theta+Miller_x*sin(theta)) &
        * (1+Miller_x*cos(theta)))**2 + (Miller_kappa*cos(theta))**2) &
        / (RHat(theta) * ( cos(Miller_x*sin(theta)) + Miller_dRdr*cos(theta) + (Miller_s_kappa-Miller_s_delta*cos(theta) &
        + (1+Miller_s_kappa)*Miller_x*cos(theta)) * sin(theta) * sin(theta + Miller_x*sin(theta))))
    end function

    function BDotGradTheta(theta)
      implicit none
      PetscScalar :: theta, BDotGradTheta
      BDotGradTheta = - Miller_A*Miller_QQ/(Miller_kappa*Miller_q*RHat(theta) * &
        ((1+Miller_s_kappa)*sin(theta + Miller_x*sin(theta)) * (1+Miller_x*cos(theta)) * sin(theta) &
        + cos(theta) * (Miller_dRdr + cos(theta + Miller_x *sin(theta)) &
        - Miller_s_delta*sin(theta + Miller_x*sin(theta)) * sin(theta))))
    end function
    
  end module
