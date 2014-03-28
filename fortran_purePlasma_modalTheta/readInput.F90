  module readInput
    use parameters
    implicit none
    public readNamelistInput
    
    contains
    
    subroutine readNamelistInput(filename)
      implicit none
      character(len=*), intent(in) :: filename
      integer :: fileUnit, didFileAccessWork
      namelist / parameters / epsil, nuPrime, geometryToUse, whereToPutbDotGradTheta, &
        useIterativeSolver, multipleOfNuPrimeToAddOnDiagonalOfPreconditioner, &
        saveMatlabOutput, programMode, MatlabOutputFilename, convergenceScanFilename, &
        Miller_kappa, Miller_delta, Miller_s_delta, Miller_s_kappa, Miller_dRdr, Miller_q, &
        Ntheta, NthetaMaxFactor, NthetaMinFactor, NthetaNumRuns, &
        Nxi, NxiMaxFactor, NxiMinFactor, NxiNumRuns, &
        NL, NLMaxFactor, NLMinFactor, NLNumRuns, &
        Nx, NxMaxFactor, NxMinFactor, NxNumRuns, &
        xMax, xMaxMaxFactor, xMaxMinFactor, xMaxNumRuns, &
        NxPotentialsPerVth, NxPotentialsPerVthMaxFactor, NxPotentialsPerVthMinFactor, NxPotentialsPerVthNumRuns
      fileUnit=11
      open(unit=fileUnit, file=filename,    action="read", status="old", iostat=didFileAccessWork)
      if (didFileAccessWork /= 0) then
        print *,"Error opening ", filename, ".  Default parameters used instead"
      else
        read(fileUnit, nml=parameters, iostat=didFileAccessWork)
        if (didFileAccessWork /= 0) then
          print *,"I was able to open the file ", filename, " but not read data from it.  Default parameters used instead."
        else
          print *,"Successfully read parameters from ", filename, "."
        end if
      end if
      close(unit = fileUnit)
    end subroutine

  end module

