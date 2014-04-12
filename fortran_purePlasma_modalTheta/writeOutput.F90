  module writeOutput
    use parameters
    implicit none
    public printInputs
    integer, private :: convergenceScanFileUnit
    
    contains
    
    subroutine printGreeting()
      print *,"****************************************************************************"
      print *,"Local tokamak drift-kinetic equation solver."
      print *,"Fourier modal in theta, Legendre modal in xi, polynomial spectral collocation in x."
#if defined(PETSC_USE_REAL_SINGLE)
      print *,"Using single precision."
#else
      print *,"Using double precision."
#endif
    end subroutine
    
    subroutine printInputs()
      print *,"---- Physics parameters: ----"
      print *,"epsilon = ", epsil
      print *,"nuStar  = ", nuStar
      print *,"nuPrime = ", nuPrime
      print *,"Geometry scheme: ", geometryToUse
      print *,"---- Numerical parameters: ----"
      print *,"Ntheta  = ", Ntheta
      print *,"Nxi     = ", Nxi
      print *,"NL      = ", NL
      print *,"Nx      = ", Nx
      print *,"NxPotentialsPerVth = ", NxPotentialsPerVth
      print *,"xMax = ",xMax
      print *,"whereToPutbDotGradTheta: ", whereToPutbDotGradTheta
      if (useIterativeSolver) then
        print *,"Using iterative solver"
      else
        print *,"Using direct solver"
      end if
    end subroutine

    subroutine printOutputs()
      print *,"Total elapsed time: ", elapsedTime, " seconds."
      print *,"q = ", q
      print *,"particleFlux = ", particleFlux
      print *,"k = ", k
    end subroutine

    subroutine openConvergenceScanOutputFile()
      integer :: didFileAccessWork
      open(unit=convergenceScanFileUnit, file=convergenceScanFilename,    &
        action="write", status="replace", iostat=didFileAccessWork)
      if (didFileAccessWork /= 0) then
        print *,"Error opening ", convergenceScanFilename, "."
        stop
      end if
      write (unit=convergenceScanFileUnit, fmt="(a)") &
        "Ntheta,Nxi,NL,Nx,NxPotentialsPerVth,xMax,q,particleFlux,k,didItConverge,elapsedTime"
    end subroutine
    
    subroutine addLineToConvergenceScanFile()
      write (unit=convergenceScanFileUnit, &
      fmt="(i5, a, i5, a, i5, a, i5, a, f18.10, a, f18.10, a, f18.10, a, f18.10, a, f18.10, a, i5, a, f18.10)") &
        Ntheta, ",", Nxi, ",", NL, ",", Nx, ",", NxPotentialsPerVth, ",", xMax, ",", q, ",", particleFlux, &
        ",", k, ",", didItConverge, ",", elapsedTime
    end subroutine
    
    subroutine closeConvergenceScanOutputFile()
    close(unit = convergenceScanFileUnit)
    end subroutine
  end module

