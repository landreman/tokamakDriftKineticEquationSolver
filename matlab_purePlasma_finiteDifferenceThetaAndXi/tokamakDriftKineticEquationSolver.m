function tokamakDriftKineticEquationSolver()

% Tokamak drift-kinetic equation solver.
% Original version created 2012 by Matt Landreman
% Massachusetts Institute of Technology
% Plasma Science & Fusion Center.

% This program implements the algorithms described in
% Landreman & Ernst, J Comp Phys 243, 130 (2013).

% This version applies to a pure plasma and makes the small orbit width
% approximation: rho_poloidal for the ions is small compared to all
% equilibrium scale lengths.

% This program calculates the bootstrap current coefficients and neoclassical
% conductivity. Specifically, it computes the quantities L_31, L_32, L_34, 
% alpha, and sigma_neo / sigma_Spitzer defined by equations (5)-(6) in
% Sauter, Angioni, and Lin-Liu, Phys Plasmas 6 (1999).
% This program also computes the neoclassical ion conductivity, defined by
% eq (5) in
% Landreman & Ernst, Plasma Phys Controlled Fusion 54, 115006 (2012).
% Notice the ion flow coefficient k_|| in eq (6) of this paper is equal 
% to -alpha from the Sauter paper.

% Discretizations used in this version:
% theta: finite differences or Fourier spectral collocation (see the 'thetaGridMode' input parameter.)
% xi (cosine of pitch angle): finite differences for the distribution
%                             function, Legendre modal expansion for the
%                             Rosenbluth potentials.
% x (speed / thermal speed): spectral collocation

%******************************************************************
% Options for controlling the general program execution:
%******************************************************************

runMode = 0;
% 0 = Single run. This usually takes anywhere from < 1 second to a couple
%     minutes, depending on the resolution used.
% 1 = Do a bunch of solves, scanning the numerical resolution parameters to
%     check convergence.  This scan may take a fraction of a minute to
%     several minutes.
% 2 = Do a bunch of solves, scanning the collisionality and also doing a 
%     convergence scan, and saving the results. This scan will take at 
%     least several minutes.
% 3 = Replot a saved collisionality scan, without doing any solves. This
%     does not take significant time.

% The following parameter specifies the file to plot when runMode = 3:
dataFileToPlot = 'tokamakDriftKineticEquationSolver_2014-04-12_14-05_foobar.mat';

% The switch below determines whether an output file is created when runMode = 2:
saveStuff = true;
%saveStuff = false;

% The string below is appended to the output file name when runMode = 2:
filenameNote = 'foobar';

%******************************************************************
% Parameters for magnetic geometry:
%******************************************************************

geometry = 1;
% 0 = Concentric circular flux surfaces
% 1 = Miller geometry
% 2 = Load an eqdsk / EFIT equilibrium from a file

% The geometry=0 scheme is defined by 
% B \propto 1 / (1 + (Miller_A)^{-1} * cos(theta))
% and
% grad_|| theta = constant.

% For details of the geometry=1 (Miller) scheme, see the appendix of
% Landreman & Ernst, J Comp Phys 243, 130 (2013).

% Parameters for Miller geometry:
Miller_A = 3.17;
Miller_kappa = 1.66;
Miller_delta = 0.416;
Miller_s_delta = 1.37;
Miller_s_kappa = 0.70;
Miller_dRdr = -0.354;
Miller_q = 3.03;

% When geometry=1 is used, set the aspect ratio by setting Miller_A.

% The next parameters matter only when geometry = 2:
EFITFilename = '../equilibria/g_CMod_1120824011.01140';  topCropZ = 0.4;    bottomCropZ = -0.4;

% Specify the value of normalized poloidal flux for the flux surface on
% which we want to do the kinetic computation:
desiredPsi = 0.7;

%******************************************************************
% Phyics parameters:
%******************************************************************

species = 0;
% 0 for ions
% 1 for electrons
% For runMode=2, the value of 'species' here is ignored, and both ion and electron equations are solved.

% Z is the ion charge, which (paradoxically perhaps) matters in the code
% only when species = 1.
Z = 1;

% Next we define the collisionality parameter nuPrime.
% For ions (species = 0), 
% nuPrime = nu_ii * q * R_0 / v_i
% where
% v_i = sqrt(2 * T_i / m_i).
% For electrons (species = 1),
% nuPrime = nu_ee * q * R_0 / v_e
% where
% v_e = sqrt(2 * T_e / m_e).

% For geometry=0, q*R_0 here is taken to be the constant 1 / (grad_|| theta).
% For geometry=1, q = Miller_q is the safety factor and R_0 is the major
%                 radius used to define Miller geometry.
% For geometry=2, q is the safety factor of the flux surface and R_0 is
%                 taken to be the major radius of the magnetic axis in the
%                 eqdsk file.

% The dimensional collision frequencies used to define nuPrime above are:
%
%                  4 * sqrt{2*pi} * n_i * Z^4 * e^4 * ln(Lambda)
% nu_ii = -----------------------------------------------------------   (SI units)
%             3 * (4 * pi * epsilon_0)^2 * sqrt(m_i} * T_i^(3/2)
%
%
%                  4 * sqrt{2*pi} * n_e * e^4 * ln(Lambda)
% nu_ee = -----------------------------------------------------------   (SI units)
%             3 * (4 * pi * epsilon_0)^2 * sqrt(m_e} * T_e^(3/2)
%
% or, equivalently,
%
%                  4 * sqrt{2*pi} * n_i * Z^4 * e^4 * ln(Lambda)
% nu_ii = -----------------------------------------------------------   (Gaussian units)
%                       3 * sqrt(m_i} * T_i^(3/2)
%
%
%                  4 * sqrt{2*pi} * n_e * e^4 * ln(Lambda)
% nu_ee = -----------------------------------------------------------   (Gaussian units)
%                       3 * sqrt(m_e} * T_e^(3/2)
%

nuPrime = 0.3;

nuStar = nuPrime * (Miller_A^1.5);

% Notice that this definition of nu_*i is identical to Sauter's,
% but nu_{*e}^{Sauter} = sqrt(2) * nu_{*e}^{here}.


%******************************************************************
% Numerical resolution parameters:
%******************************************************************
% For each pair of variables below, the first is the value used in a single run.
% The second is the set of values used in a convergence scan.

% Number of grid points in the poloidal angle.
% Memory and time requirements do depend strongly on this parameter.  The
% value of Ntheta required for numerical convergence depends strongly on
% collisionality. At high collisionality, Ntheta = 5 might be plenty. At nu* =
% 10^{-3}, Ntheta > 60 may be required.
NthetaConverged = 17;
Nthetas = floor(linspace(13,31,5));

% Number of Legendre modes retained for the Rosenbluth potentials.
% Unlike the versions of this program that use a Legendre modal discretization
% in xi, for this version in which finite differences in xi are used, the memory
% and time requirements do depend on this parameter.  A value of 2 or 4 
% almost always works well.
NLConverged = 2;
NLs = 2:6;

% Number of grid points in pitch angle retained for the distribution function.
% Memory and time requirements do depend strongly on this parameter.  The
% value of Nxi required for numerical convergence depends strongly on
% collisionality. At high collisionality, Nxi=5 might be plenty. At nu* =
% 10^{-3}, Nxi > 100 may be required.
NxiConverged = 25;
Nxis = floor(linspace(10,40,9));

% Number of grid points in x = v / v_th used for the distribution function.
% Typically, values in the range 5-8 work well, with slightly higher values
% required at higher collisionality.
NxConverged = 6;
Nxs = 6:10;

% Number of grid points in x used for the Rosenbluth potentials.
% Unlike the versions of this program that use a Legendre modal discretization
% in xi, for this version in which finite differences in xi are used, the memory
% and time requirements do depend on this parameter.  A value of 3-10 is
% usually good.
NxPotentialsPerVthConverged = 4;
NxPotentialsPerVths = linspace(2,10,5);

% Tolerance for iterative solver. Ignored when the direct solver is used.
% Memory requirements do not depend on this parameter, and the time
% required by the iterative solver depends somewhat on this parameter.
% Typically, values in the range 5-7 are good.
log10tolConverged = 5.0;
log10tols = 4.5:0.5:5;

%******************************************************************
% Other numerical parameters:
%******************************************************************

% The following parameter should usually be 2.
thetaGridMode = 2;
thetaGridModeForPreconditioner = thetaGridMode;
% 0 = uniform periodic spectral collocation
% 1 = finite difference, 3 point stencil
% 2 = finite difference, 5 point stencil
% 3 = finite difference, 7 point stencil

% This next parameter should be 1 except in rare circumstances.
forceThetaParity = 1;
% 0 = either even or odd Ntheta is fine.
% 1 = force Ntheta to be odd.
% 2 = force Ntheta to be even.

% The code's performance does not seem to depend much on the following two
% parameters:
xiGridMode = 2;
xiGridModeForPreconditioner = 1;
% 1 = finite difference, 3 point stencil
% 2 = finite difference, 5 point stencil

tryIterativeSolvers = true;
%tryIterativeSolvers = false;
% If false, the sparse direct solver will be used.

tryDirectSolverIfIterativeSolversFail = true;
%tryDirectSolverIfIterativeSolversFail = false;

% The setting below determines the order of iterative solvers to try:
orderOfSolversToTry = [1, 1, 3];
% 1 = GMRES
% 2 = BiCGStab
% 3 = BiCGStab(l)
% 4 = TFQMR
% 5 = CGS

% Maximum number of iterations allowed in the iterative solver:
maxIterations = 400;

% The parameter below only matters for the GMRES solver:
restart = maxIterations;

% If the following parameter is true, 2 extra rows are added to the linear
% system to force all the flux-surface-averaged density and pressure to be
% in f_Maxwellian instead of in f_1. Two extra columns are also added
% (associated with Lagrange multipliers) so the linear system remains
% square.  It seems to make little difference whether or not these
% constraints are included.
%includeConstraints = true;
includeConstraints = false;


%******************************************************************
% Options for the preconditioner, controlling how much simpler the 
% preconditioner matrix is than the "actual" matrix.  These options are all
% ignored when the direct solver is used.
%******************************************************************

% The following option controls the pieces of the field term in the 
% Fokker-Planck operator that involve the Rosenbluth potentials.
preconditionerFieldTermOption = 1;
% 0 = Include all available Legendre modes of the potentials.
% 1 = Do not include any of the Legendre modes. This makes the
%     preconditioner much sparser.
% 2 = Include only Legendre modes 0, 1, and 2.
% Typically 1 is a good value for this parameter.

% The following option controls the M21 block, i.e. mapping f_1 from the
% finite difference xi grid to the Legendre modal expansion for the Poisson
% equation for the H Rosenbluth potential.
preconditionerM21Option = 1;
% 0 = Keep M21 in the preconditioner matrix.
% 1 = Do not keep M21 in the preconditioner matrix. This makes the
%     preconditioner much sparser.
% Typically 1 is a good option.

% The following option controls the M32 block, i.e. the source term due to
% the H potential in the Poisson equation for the G potential.
preconditionerM32Option = 1;
% 0 = Keep M32 in the preconditioner matrix.
% 1 = Do not keep M32 in the preconditioner matrix. This makes the
%     preconditioner slightly sparser.
% Typically 1 is a good option.

%******************************************************************
% Switches for plotting:
%******************************************************************

% The number below is added to each figure number. You can change it if you
% want to save windows from a previous run.
figureOffset=30;

plotVelocitySpaceGrid = true;
%plotVelocitySpaceGrid = false;

plotEFITDetails = true;
%plotEFITDetails = false;

%plotPoloidalVariationOfFlow = true;
plotPoloidalVariationOfFlow = false;

makeSpyPlot = true;
%makeSpyPlot = false;

%******************************************************************
%******************************************************************
%
% End of the main input parameters
%
%******************************************************************
%******************************************************************

speedGridFigureHandle = 0;
KrylovFigureHandle = 0;

if geometry==2 && runMode ~= 3
    % Load EFIT data
    
    NPsi=1;
    psi = desiredPsi;
    [thetaData, BData, BDotGradThetaData, qData, as, R0, B0] = getGeometryFromEFITForSeveralFluxSurfaces(EFITFilename, psi, topCropZ, bottomCropZ, plotEFITDetails);
    B0 = abs(B0);
    
    % Filter B(theta)
    NThetaFilter = 100;
    thetaForBFilter = linspace(0,2*pi,NThetaFilter+1);
    thetaForBFilter(end)=[];
    bModes = zeros(NThetaFilter, NPsi);
    for psiIndex = 1:NPsi
        temp = interp1(thetaData{psiIndex}, BData{psiIndex}/B0, thetaForBFilter, 'spline');
        bModes(:,psiIndex) = fft(temp)/NThetaFilter;
    end
    
    numModesToKeep = 5;
    numModesToKeeps = numModesToKeep;
    
    epsilon = as/R0;
    Miller_A = 1/epsilon;
    %fprintf('Epsilon derived from EFIT equilibrium: %g\n',epsilon)
    
end


dBdthetaResolutionMultiplier = 10;
% Compute a few quantities related to the Miller equilibrium:
Miller_x = asin(Miller_delta);
Miller_Q = Miller_kappa / (2*pi*Miller_A) * quadgk(@QIntegrand, 0, 2*pi);


q=0;
k=0;
particleFlux=0;
density=0;
theta=0;
theta1D=linspace(0,2*pi,100);
iteration=0;

if plotVelocitySpaceGrid && (runMode ~= 3)
    figure(figureOffset+7)
    clf
end
epsilon = 1/Miller_A;


if runMode ~= 3
    VPrime = quadgk(@(theta) (oneOverqRbDotGradTheta(theta)./b(theta)), 0, 2*pi);
    FSAB2 = quadgk(@(theta) (oneOverqRbDotGradTheta(theta).*b(theta)), 0, 2*pi) / VPrime;
    
    if geometry==2
        [~,lambdaMax] = fminbnd(@(x) (1./b(x)), pi/2, 3*pi/2);
    else
        lambdaMax = 1/b(pi);
    end
    lambda=0;
    
    fc = 0.75*FSAB2 * quadgk(@outerIntegrand, 0, lambdaMax);
    if abs(imag(fc)/real(fc)) > 1e-2
        error('Something went wrong. fc has an imaginary part that is not negligible.')
    end
    fc=real(fc);
    
    ft=1-fc;
    
    %fprintf('f_t for this magnetic geometry: %g\n',ft)
end

switch runMode
    case 0
        % Single run
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime,nuStar);
        
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        solveDKE();
        
    case 2
        % Make plot for comparison with Sauter
        
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        
        nuStars = 10.^(2:(-0.25):(-3));
        Nthetas =     [ 8, 8, 8, 8, 8, 8, 8, 8, 8, 9,10,11,12,13,14,15,25,32,32,32,32]*2;
        Nxis =        [15,15,15,15,15,15,15,15,15,28,28,35,46,55,60,64,85,100,140,140,140];
        Nxs =         [ 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7];
        %               ^           ^           ^           ^           ^
        
        % Comment out these next 4 lines if you want a higher resolution
        % scan.
        nuStars = nuStars(1:4:end);
        Nthetas = Nthetas(1:4:end);
        Nxis = Nxis(1:4:end);
        Nxs = Nxs(1:4:end);
        
        NthetaMultipliers = [1, 2, 1, 1];
        NxiMultipliers    = [1, 1, 2, 1];
        NxMultipliers     = [1, 1, 1, 2];
        
        numNuStars = numel(nuStars);
        numConvergenceTests = numel(NthetaMultipliers);
        nuPrimes = nuStars*(epsilon^1.5);
        
        NL=2;
        flowCoefficients = zeros(numConvergenceTests,numNuStars);
        conductivities = zeros(numConvergenceTests,numNuStars);
        L31s = zeros(numConvergenceTests,numNuStars);
        L32s = zeros(numConvergenceTests,numNuStars);
        L34s = zeros(numConvergenceTests,numNuStars);
        
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning nu* scan.  Epsilon = %g\n',epsilon)
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        for convergenceScan = 1:numConvergenceTests
            for runNum=1:numNuStars
                nuPrime=nuPrimes(runNum);
                nuStar = nuPrime/(epsilon^1.5);
                if nuStar<0.3
                    tryIterativeSolvers = true;
                else
                    tryIterativeSolvers = false;
                end
                fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
                fprintf('nuStar = %g, nuPrime = %g\n',nuStar,nuPrime);
                fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
                Ntheta=floor(Nthetas(runNum)*NthetaMultipliers(convergenceScan));
                Nxi=floor(Nxis(runNum)*NxiMultipliers(convergenceScan));
                Nx=floor(Nxs(runNum)*NxMultipliers(convergenceScan));
                species=0;
                solveDKE();
                flowCoefficients(convergenceScan,runNum)=k;
                species=1;
                solveDKE();
                conductivities(convergenceScan,runNum)=conductivity;
                L31s(convergenceScan,runNum)=L31;
                L32s(convergenceScan,runNum)=L32;
                L34s(convergenceScan,runNum)=L34;
            end
        end
        
        assignin('base','nuPrimes',nuPrimes)
        assignin('base','nuStars',nuStars)
        assignin('base','flowCoefficients',flowCoefficients)
        assignin('base','conductivities',conductivities)
        assignin('base','L31s',L31s)
        assignin('base','L32s',L32s)
        assignin('base','L34s',L34s)
        
        temp=dbstack;
        nameOfThisProgram=sprintf('%s',temp.file);
        filenameBase=[nameOfThisProgram(1:(end-2)),'_',datestr(now,'yyyy-mm-dd_HH-MM'),'_',filenameNote];
        outputFilename=[filenameBase,'.mat'];
        if saveStuff
            save(outputFilename)
        end
        
        plotFigure()
    case 3
        load(dataFileToPlot)
        plotFigure()
    case 1
        % Convergence scan
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Beginning convergence scans.\n')
        fprintf('nuPrime = %g, nu_* = %g\n',nuPrime, nuStar);
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        
        startTime=clock;
        
        if species==0
            quantitiesToRecord = {'q','k','Particle flux'};
        else
            quantitiesToRecord = {'conductivity','L31','L32','L34'};
        end
        linespecs = {'.-r','.-b','.-g','.-c','.-m','.-k','.-r','.-b'};
        
        parametersToVary = {'N\theta','NL','N\xi','Nx','NxPotentialsPerVth','-log_{10}tol'};
        abscissae = {Nthetas, NLs, Nxis, Nxs, NxPotentialsPerVths, log10tols};
        convergeds = {NthetaConverged, NLConverged, NxiConverged, NxConverged, NxPotentialsPerVthConverged, log10tolConverged};
        
        numQuantities = numel(quantitiesToRecord);
        numParameters = numel(parametersToVary);
        quantities = cell(numParameters,1);
        quantities{1} = zeros(numel(Nthetas), numQuantities);
        quantities{2} = zeros(numel(NLs), numQuantities);
        quantities{3} = zeros(numel(Nxis), numQuantities);
        quantities{4} = zeros(numel(Nxs), numQuantities);
        quantities{5} = zeros(numel(NxPotentialsPerVths), numQuantities);
        quantities{6} = zeros(numel(log10tols), numQuantities);
        thetas = cell(numParameters,1);
        densities = cell(numParameters,1);
        
        % Vary Ntheta, keeping other numerical parameters fixed.
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        parameterScanNum=1;
        for iii = 1:numel(Nthetas)
            Ntheta=Nthetas(iii);
            solveDKE()
            if species==0
                quantities{parameterScanNum}(iii,1)=q;
                quantities{parameterScanNum}(iii,2)=k;
                quantities{parameterScanNum}(iii,3)=particleFlux;
            else
                quantities{parameterScanNum}(iii,1)=conductivity;
                quantities{parameterScanNum}(iii,2)=L31;
                quantities{parameterScanNum}(iii,3)=L32;
                quantities{parameterScanNum}(iii,4)=L34;
            end
        end
        thetas{parameterScanNum} = theta;
        densities{parameterScanNum} = density;
        
        % Vary NL, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        parameterScanNum=2;
        for iii = 1:numel(NLs)
            NL=NLs(iii);
            solveDKE()
            if species==0
                quantities{parameterScanNum}(iii,1)=q;
                quantities{parameterScanNum}(iii,2)=k;
                quantities{parameterScanNum}(iii,3)=particleFlux;
            else
                quantities{parameterScanNum}(iii,1)=conductivity;
                quantities{parameterScanNum}(iii,2)=L31;
                quantities{parameterScanNum}(iii,3)=L32;
                quantities{parameterScanNum}(iii,4)=L34;
            end
        end
        thetas{parameterScanNum} = theta;
        densities{parameterScanNum} = density;
        
        % Vary Nxi, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        parameterScanNum=3;
        for iii = 1:numel(Nxis)
            Nxi=Nxis(iii);
            solveDKE()
            if species==0
                quantities{parameterScanNum}(iii,1)=q;
                quantities{parameterScanNum}(iii,2)=k;
                quantities{parameterScanNum}(iii,3)=particleFlux;
            else
                quantities{parameterScanNum}(iii,1)=conductivity;
                quantities{parameterScanNum}(iii,2)=L31;
                quantities{parameterScanNum}(iii,3)=L32;
                quantities{parameterScanNum}(iii,4)=L34;
            end
        end
        thetas{parameterScanNum} = theta;
        densities{parameterScanNum} = density;
        
        
        % Vary Nx, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        tol = 10^(-log10tolConverged);
        parameterScanNum=4;
        for iii = 1:numel(Nxs)
            Nx = Nxs(iii);
            solveDKE()
            if species==0
                quantities{parameterScanNum}(iii,1)=q;
                quantities{parameterScanNum}(iii,2)=k;
                quantities{parameterScanNum}(iii,3)=particleFlux;
            else
                quantities{parameterScanNum}(iii,1)=conductivity;
                quantities{parameterScanNum}(iii,2)=L31;
                quantities{parameterScanNum}(iii,3)=L32;
                quantities{parameterScanNum}(iii,4)=L34;
            end
        end
        thetas{parameterScanNum} = theta;
        densities{parameterScanNum} = density;
        
        % Vary NxPotentialsPerVth, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        tol = 10^(-log10tolConverged);
        parameterScanNum=5;
        for iii = 1:numel(NxPotentialsPerVths)
            NxPotentialsPerVth = NxPotentialsPerVths(iii);
            solveDKE()
            if species==0
                quantities{parameterScanNum}(iii,1)=q;
                quantities{parameterScanNum}(iii,2)=k;
                quantities{parameterScanNum}(iii,3)=particleFlux;
            else
                quantities{parameterScanNum}(iii,1)=conductivity;
                quantities{parameterScanNum}(iii,2)=L31;
                quantities{parameterScanNum}(iii,3)=L32;
                quantities{parameterScanNum}(iii,4)=L34;
            end
        end
        thetas{parameterScanNum} = theta;
        densities{parameterScanNum} = density;
        
        % Vary tol, keeping other numerical parameters fixed.
        Ntheta=NthetaConverged;
        NL=NLConverged;
        Nxi=NxiConverged;
        Nx=NxConverged;
        NxPotentialsPerVth = NxPotentialsPerVthConverged;
        parameterScanNum=6;
        for iii = 1:numel(log10tols)
            tol = 10^(-log10tols(iii));
            solveDKE()
            if species==0
                quantities{parameterScanNum}(iii,1)=q;
                quantities{parameterScanNum}(iii,2)=k;
                quantities{parameterScanNum}(iii,3)=particleFlux;
            else
                quantities{parameterScanNum}(iii,1)=conductivity;
                quantities{parameterScanNum}(iii,2)=L31;
                quantities{parameterScanNum}(iii,3)=L32;
                quantities{parameterScanNum}(iii,4)=L34;
            end
        end
        thetas{parameterScanNum} = theta;
        densities{parameterScanNum} = density;
        
        maxs=ones(numQuantities,1)*(-1e10);
        mins=ones(numQuantities,1)*(1e10);
        for iParameter = 1:numParameters
            maxs = max([maxs, quantities{iParameter}'],[],2);
            mins = min([mins, quantities{iParameter}'],[],2);
        end
        
        figure(1+figureOffset)
        numRows = numQuantities;
        numCols = numParameters;
        clf
        for iQuantity = 1:numQuantities
            if maxs(iQuantity) <= mins(iQuantity)
                maxs(iQuantity) = mins(iQuantity)+1;
            end
            for iParameter = 1:numParameters
                subplot(numRows, numCols, iParameter  + (iQuantity - 1)*numParameters)
                plot(1./abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
                hold on
                plot(1./[convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
                ylim([mins(iQuantity), maxs(iQuantity)])
                xlabel(['1/',parametersToVary{iParameter}])
                ylabel(quantitiesToRecord{iQuantity})
            end
        end
        switch species
            case 0
                speciesText = 'Ions';
            case 1
                speciesText = 'Electrons';
            otherwise
                error('Program should not get here')
        end
        switch geometry
            case 0
                geometryText = 'concentric circular flux surfaces';
            case 1
                geometryText = 'Miller geometry';
            case 2
                geometryText = 'eqdsk geometry';
            otherwise
                error('Program should not get here')
        end
        stringForTop=sprintf('%s, %s, aspect ratio = %g, nuPrime = %g, nu_* = %g, thetaGridMode = %d, finite difference in xi, polynomial collocation in x',speciesText,geometryText,Miller_A, nuPrime,nuStar,thetaGridMode);
        annotation('textbox',[0 0.95 1 .05],'HorizontalAlignment','center',...
            'Interpreter','none','VerticalAlignment','bottom',...
            'FontSize',10,'LineStyle','none','String',stringForTop);
        
        temp=dbstack;
        nameOfThisProgram=sprintf('%s',temp(1).file);
        stringForBottom = ['Plotted using ',nameOfThisProgram];
        
        annotation('textbox',[0 0 1 .04],'HorizontalAlignment','center',...
            'Interpreter','none','VerticalAlignment','top',...
            'FontSize',10,'LineStyle','none','String', ...
            stringForBottom);

        figure(8+figureOffset)
        numRows = numQuantities;
        numCols = numParameters;
        clf
        for iQuantity = 1:numQuantities
            if maxs(iQuantity) <= mins(iQuantity)
                maxs(iQuantity) = mins(iQuantity)+1;
            end
            for iParameter = 1:numParameters
                subplot(numRows, numCols, iParameter  + (iQuantity - 1)*numParameters)
                plot(abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
                hold on
                plot([convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
                ylim([mins(iQuantity), maxs(iQuantity)])
                xlabel(parametersToVary{iParameter})
                ylabel(quantitiesToRecord{iQuantity})
            end
        end
        
        annotation('textbox',[0 0.95 1 .05],'HorizontalAlignment','center',...
            'Interpreter','none','VerticalAlignment','bottom',...
            'FontSize',10,'LineStyle','none','String',stringForTop);
        
        annotation('textbox',[0 0 1 .04],'HorizontalAlignment','center',...
            'Interpreter','none','VerticalAlignment','top',...
            'FontSize',10,'LineStyle','none','String', ...
            stringForBottom);
        
        
        
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
        fprintf('Total elapsed time: %g seconds.\n',etime(clock,startTime))
        fprintf('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n')
    otherwise
        error('Unrecognized runMode')
end











% **************************************************************
% **************************************************************
%
% Start of main function to solve the kinetic equation
%
% **************************************************************
% **************************************************************


    function solveDKE()
        
        startTimeForThisRun=tic;
        
        sqrtpi=sqrt(pi);
        iteration = iteration+1;
        
        switch forceThetaParity
            case 0
                % Do nothing
            case 1
                % Force Ntheta to be odd
                if mod(Ntheta,2)==0
                    Ntheta=Ntheta+1;
                end
            case 2
                % Force Ntheta to be even
                if mod(Ntheta,2)==1
                    Ntheta=Ntheta+1;
                end
            otherwise
                error('Invalid forceThetaParity')
        end
        
        fprintf('Ntheta = %d,  NL = %d,  Nxi = %d,  Nx = %d, NxPtentialsPerVth = %g, tol = %g\n',Ntheta,NL,Nxi,Nx,NxPotentialsPerVth,tol)
        
        tic
        
        % Generate abscissae, quadrature weights, and derivative matrix for theta grid.
        switch thetaGridMode
            case 0
                % Spectral uniform
                scheme = 20;
            case 1
                % Uniform periodic 2nd order FD
                scheme = 0;
            case 2
                % Uniform periodic 4th order FD
                scheme = 10;
            case 3
                % Uniform periodic FD, 7 point stencil
                scheme = 70;
            otherwise
                error('Error! Invalid thetaGridMode')
        end
        [theta, thetaWeights, ddtheta, ~] = differentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
        theta = theta';
        thetaWeights = thetaWeights';
        
        bs=b(theta);
        oneOverqRbDotGradThetas = oneOverqRbDotGradTheta(theta);
        if geometry==1
            % Miller: too hard to analytically differentiate b(theta) so do
            % it numerically:
            scheme = 20;
            [thetaFine, ~, ddthetaFine, ~] = differentiationMatricesForUniformGrid(Ntheta*dBdthetaResolutionMultiplier, 0, 2*pi, scheme);
            dbdthetaFine = ddthetaFine * b(thetaFine);
            dbdthetas = dbdthetaFine(1:dBdthetaResolutionMultiplier:end)';
        else
            dbdthetas = dbdtheta(theta);
        end
        
        switch thetaGridModeForPreconditioner
            case 0
                % Spectral uniform
                scheme = 20;
            case 1
                % Uniform periodic 2nd order FD
                scheme = 0;
            case 2
                % Uniform periodic 4th order FD
                scheme = 10;
            case 3
                % Uniform periodic FD, 7 point stencil
                scheme = 70;
            otherwise
                error('Error! Invalid thetaGridMode')
        end
        [~, ~, ddthetaForPreconditioner, ~] = differentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);
        
        
        % Generate abscissae, quadrature weights, and derivative matrix for x grid.
        kk=0;
        scale=1;
        pointAtZero=false;
        [x, ddx, d2dx2, xWeights] = spectralNodesWeightsAndDifferentiationMatricesForV(Nx, kk, scale, pointAtZero);
        x = x(:)';
        xWeights = xWeights(:)';
        
        function y=weight(xxx)
            y=exp(-xxx.*xxx);
        end
        
        xMax=max([4.5, max(x)]);
        NxPotentials = ceil(xMax * NxPotentialsPerVth);
        % Uniform, higher order FD
        xMin=0;
        scheme = 12;
        [xPotentials, ~, ddxPotentials, d2dx2Potentials] = differentiationMatricesForUniformGrid(NxPotentials, xMin, xMax, scheme);
        regridPolynomialToUniform = polynomialInterpolationMatrix(x,xPotentials,weight(x),weight(xPotentials));
        regridUniformToPolynomial = makeHighOrderInterpolationMatrix(xPotentials,x,0,'f');
        
        
        % Generate abscissae, quadrature weights, and derivative matrix for xi grid.
        switch xiGridMode
            case 1
                % Finite difference, 3 point stencil
                scheme = 2;
            case 2
                % Finite difference, 5 point stencil
                scheme = 12;
            otherwise
                error('Error! Invalid xiGridMode')
        end
        [xi, xiWeights, ddxi, d2dxi2] = differentiationMatricesForUniformGrid(Nxi, -1, 1, scheme);
        
        switch xiGridModeForPreconditioner
            case 1
                % Finite difference, 3 point stencil
                scheme = 2;
            case 2
                % Finite difference, 5 point stencil
                scheme = 12;
            otherwise
                error('Error! Invalid xiGridMode')
        end
        [~, ~, ddxiPreconditioner, d2dxi2Preconditioner] = differentiationMatricesForUniformGrid(Nxi, -1, 1, scheme);
        xi = xi(:)';
        xiWeights = xiWeights(:)';
        
        % ********************************************
        % Generate Legendre polynomials on the xi grid.
        % ********************************************
        Legendres=ones(NL, Nxi);
        LegendresTimesWeights=zeros(NL, Nxi);
        Legendres(2,:)=xi;
        for L=2:(NL-1)
            Legendres(L+1,:) = ((2*L-1)*xi.*Legendres(L,:) - (L-1)*Legendres(L-1,:))/L;
        end
        for iL=1:NL
            LegendresTimesWeights(iL,:) = Legendres(iL,:) .* xiWeights;
        end        
        
        if plotVelocitySpaceGrid
            if speedGridFigureHandle == 0
                speedGridFigureHandle = figure(figureOffset+7);
            else
                set(0, 'CurrentFigure', speedGridFigureHandle);
            end
            plot(xPotentials,zeros(size(xPotentials))+iteration,'.r')
            hold on
            plot(x, zeros(size(x))+iteration,'o')
            title('Speed grid for distribution function (blue) and Rosenbluth potentials(red)')
            xlabel('x')
            ylabel('Solve number')
        end
        
        matrixSize = Ntheta * Nx * Nxi + 2*NxPotentials*NL*Ntheta;
        if includeConstraints
            matrixSize = matrixSize + 2;
        end
        %estimated_nnz = floor(1.1*(Nx*Nx*Nxi*Ntheta + 2*Nxi*nnz(ddtheta)*Nx+2*Nxi*Nx*Ntheta));
        estimated_nnz = floor(1.1*(Nx*Nx*Nxi*Ntheta ... % energy scattering
            + nnz(d2dxi2)*Nx*Ntheta ... % pitch angle scattering
            + Nxi*nnz(ddtheta)*Nx ... % parallel streaming
            + nnz(d2dxi2)*Nx*Ntheta ... % mirror force
            + 2*NL*nnz(d2dx2Potentials)*Ntheta ... % M22 and M33
            + 2*NL*NxPotentials*Ntheta ... % M21 and M32
            + 2*Ntheta*Nxi*NL*Nx*NxPotentials ... % M12 and M13
            ));
        if includeConstraints
            estimated_nnz = estimated_nnz + 4*Ntheta*Nxi*Nx;
        end
        fprintf('matrixSize: %d.\n',matrixSize)
        
        % Begin timer for matrix construction:
        tic
        
        % *******************************************
        % *******************************************
        %
        % Build the right-hand side(s) for the linear system:
        %
        % *******************************************
        % *******************************************
        
        switch species
            case 0
                rhs=zeros(matrixSize,1);
            case 1
                rhs=zeros(matrixSize,4);
            otherwise
                error('Invalid setting for species');
        end
        
        % Order of columns in the RHS:
        % 1 = L31
        % 2 = L32
        % 3 = L34
        % 4 = conductivity
        
        x2=x.*x;
        expx2=exp(-x2);
        xPartOfRHSForConductivity = x.*expx2;
        xPartOfRHSForP = x2.*expx2;
        xPartOfRHSForT = (x2-5/2).*x2.*expx2;
        
        thetaPartOfRHS = dbdthetas ./ (bs.*bs);
        thetaPartOfL34RHS = dbdthetas;
        thetaPartOfConductivityRHS = bs .* oneOverqRbDotGradThetas;
        
        for ixi = 1:Nxi
            for itheta=1:Ntheta
                indices = getIndices(1,1:Nx,ixi,itheta);
                rhs(indices,1) = 0.5 * (1 + xi(ixi)*xi(ixi)) * thetaPartOfRHS(itheta) * xPartOfRHSForT;
                
                if species ~= 0
                    rhs(indices,2) = 0.5 * (1 + xi(ixi)*xi(ixi)) * thetaPartOfRHS(itheta) * xPartOfRHSForP;
                    rhs(indices,3) = 0.5 * (3*xi(ixi)*xi(ixi) - 1) * thetaPartOfL34RHS(itheta) * xPartOfRHSForP;
                    rhs(indices,4) = xi(ixi) * thetaPartOfConductivityRHS(itheta) * xPartOfRHSForConductivity;
                end
                
            end
        end
        
        % *************************************
        % *************************************
        %
        % Build the main matrix:
        %
        % *************************************
        % *************************************
        
        sparseCreatorIndex=1;
        sparseCreator_i=0;
        sparseCreator_j=0;
        sparseCreator_s=0;
        resetSparseCreator()
        
        if tryIterativeSolvers
            matricesToMake=1:2;
        else
            matricesToMake=1;
        end
        
        for whichMatrixToMake = matricesToMake
            % 1 = main matrix
            % 2 = preconditioner
            
            if whichMatrixToMake==1
                ddthetaForThisMatrix = ddtheta;
                ddxiForThisMatrix = ddxi;
                d2dxi2ForThisMatrix = d2dxi2;
            else
                ddthetaForThisMatrix = ddthetaForPreconditioner;
                ddxiForThisMatrix = ddxiPreconditioner;
                d2dxi2ForThisMatrix = d2dxi2Preconditioner;
            end
            
            matrixStartTime = tic;

            % *********************************************
            % Add streaming term
            % *********************************************
            for ix=1:Nx
                for ixi=1:Nxi
                    indices = getIndices(1,ix,ixi,1:Ntheta);
                    addSparseBlock(indices, indices, x(ix)*xi(ixi)*ddthetaForThisMatrix)
                end
            end
            
            
            % *********************************************
            % Add mirror term
            % *********************************************
            thetaPart = -0.5*dbdthetas./bs;
            mirrorXiOperator = diag(1-xi.*xi)*ddxiForThisMatrix;
            for itheta = 1:Ntheta
                for ix=1:Nx
                    indices = getIndices(1,ix,1:Nxi,itheta);
                    addSparseBlock(indices, indices, x(ix)*thetaPart(itheta)*mirrorXiOperator)
                end
            end
           
            
            % ***************************************
            % Build the collision operator:
            % For a detailed explanation of the implementation of the collision operator, see
            % Landreman & Ernst, Journal of Computational Physics 243, 130 (2013)
            % ***************************************
            
            erfs=erf(x);
            x2 = x.*x;
            x3 = x2.*x;
            expx2 = exp(-x.*x);
            Psi = (erfs - 2/sqrtpi*x .* expx2) ./ (2*x.*x);
            switch species
                case 0
                    % Ions
                    nuD = 3*sqrtpi/4*(erfs - Psi) ./ x3;
                case 1
                    % Electrons
                    nuD = 3*sqrtpi/4*(erfs - Psi + Z) ./ x3;
                otherwise
                    error('Invalid setting for species')
            end
            PsiPrime = (-erfs + 2/sqrtpi*x.*(1+x.*x) .* expx2) ./ x3;
            xPartOfCECD = 3*sqrtpi/4*(diag(Psi./x)*d2dx2   +  diag((PsiPrime.*x  + Psi + 2*Psi.*x2)./x2)*ddx + diag(2*PsiPrime + 4*Psi./x)) + 3*diag(expx2);
            if whichMatrixToMake == 2
                % For preconditioner, drop all terms that are
                % off-diagonal in x.
                xPartOfCECD = diag(diag(xPartOfCECD));
            end
            M12 = nuPrime * 3/(2*pi) * diag(expx2) * regridUniformToPolynomial;
            M13 = -nuPrime * 3/(2*pi) * diag(x2.*expx2) * regridUniformToPolynomial * d2dx2Potentials;
            
            PAS = -diag(xi)*ddxiForThisMatrix + diag((1-xi.^2)*0.5) * d2dxi2ForThisMatrix;

            if whichMatrixToMake==1
                % Main matrix:
                LsToUse = 0:(NL-1);
            else
                % Preconditioner:
                switch preconditionerFieldTermOption
                    case 0
                        LsToUse = 0:(NL-1);
                    case 1
                        LsToUse = [];
                    case 2
                        LsToUse = 0:2;
                    otherwise
                        error('Invalid preconditionerFieldTermOption')
                end
            end
            
            for itheta=1:Ntheta
                % Add pitch-angle scattering:
                for ix = 1:Nx
                    indices = getIndices(1, ix, 1:Nxi, itheta);
                    addSparseBlock(indices, indices, -nuPrime * nuD(ix) * oneOverqRbDotGradThetas(itheta) * PAS);
                end
                
                
                for ixi=1:Nxi
                    % Add energy scattering, as well as the piece of the field term
                    % that involves f rather than the Rosenbluth potentials:
                    indices = getIndices(1, 1:Nx, ixi, itheta);
                    addSparseBlock(indices, indices, -nuPrime * xPartOfCECD * oneOverqRbDotGradThetas(itheta));
                    
                    % Add the part of the field term that depends on the
                    % Rosenbluth potentials
                    %for L = 0:(NL-1)
                    for L = LsToUse
                        colIndices = getIndices(2, 1:NxPotentials, L+1, itheta);
                        addSparseBlock(indices, colIndices, M12 * oneOverqRbDotGradThetas(itheta) * Legendres(L+1,ixi));
                        colIndices = getIndices(3, 1:NxPotentials, L+1, itheta);
                        addSparseBlock(indices, colIndices, M13 * oneOverqRbDotGradThetas(itheta) * Legendres(L+1,ixi));
                    end
                end
                
            end
            
                
%{            
            for itheta = 1:Ntheta
                for ixi = 1:Nxi
                    indices = getIndices(1, 1:Nx, ixi, itheta);
                end
            end
%}
            % *************************************
            % Poisson equation for G
            % *************************************
            LaplacianTimesx2WithoutLTerms = diag(xPotentials.^2) * d2dx2Potentials + diag(2*xPotentials) * ddxPotentials;
            
            for itheta=1:Ntheta
                for L=0:(NL-1)
                    LaplacianTimesx2 = LaplacianTimesx2WithoutLTerms -L*(L+1)*diag(ones(NxPotentials,1));
                    HIndices = getIndices(2,1:NxPotentials,L+1,itheta);
                    GIndices = getIndices(3,1:NxPotentials,L+1,itheta);
                    
                    % M33 terms:
                    addSparseBlock(GIndices(2:end-1), GIndices, LaplacianTimesx2(2:(end-1), :))
                    
                    if whichMatrixToMake == 1 || preconditionerM32Option == 0
                        % M32 terms:
                        addToSparse(GIndices(2:end-1), HIndices(2:end-1), -2*xPotentials(2:end-1).^2)
                    end
                    
                    % Boundary condition at xMax:
                    bc = (2*L+1)*xMax*ddxPotentials(end,:) + xMax*xMax*d2dx2Potentials(end,:);
                    bc(end) = bc(end) + L*L-1;
                    addSparseBlock(GIndices(end), GIndices, bc)

                    % Boundary condition at x=0:
                    if L==0
                        % Neumann boundary condition at x=0 for the L=0 modes:
                        addSparseBlock(GIndices(1),GIndices,ddxPotentials(1,:))
                    else
                        % Dirichlet boundary condition at x=0 for the L~=0 modes:
                        addToSparse(GIndices(1),GIndices(1),1)
                    end
                end
            end
            
            % *************************************
            % Poisson equation for H
            % *************************************
            
            xInterpolationTimesX2 = diag(xPotentials.^2)*regridPolynomialToUniform;
            xInterpolationTimesX2 = xInterpolationTimesX2(2:(end-1), :);
            
            for itheta=1:Ntheta
                for L=0:(NL-1)
                    LaplacianTimesx2 = LaplacianTimesx2WithoutLTerms -L*(L+1)*diag(ones(NxPotentials,1));
                    HIndices = getIndices(2,1:NxPotentials,L+1,itheta);
                    
                    % M22 terms:
                    addSparseBlock(HIndices(2:end-1), HIndices, LaplacianTimesx2(2:(end-1), :))
                    
                    if whichMatrixToMake == 1 || preconditionerM21Option == 0
                        % M21 terms:
                        for ixi=1:Nxi
                            fIndices = getIndices(1,1:Nx,ixi,itheta);
                            addSparseBlock(HIndices(2:end-1), fIndices, 2*pi*(2*L+1)*LegendresTimesWeights(L+1,ixi)*xInterpolationTimesX2)
                        end
                    end
                    
                    % Boundary condition at xMax:
                    bc = xMax * ddxPotentials(end,:);
                    bc(end) = bc(end) + (L+1);
                    addSparseBlock(HIndices(end), HIndices, bc)
                    
                    % Boundary condition at x=0:
                    if L==0
                        % Neumann boundary condition at x=0 for the L=0 modes:
                        addSparseBlock(HIndices(1),HIndices,ddxPotentials(1,:))
                    else
                        % Dirichlet boundary condition at x=0 for the L~=0 modes:
                        addToSparse(HIndices(1),HIndices(1),1)
                    end
                end
            end
            
            
            % *************************************************************
            % *************************************************************
            % Finished adding everything related to the collision operator.
            % *************************************************************
            % *************************************************************
            
            if includeConstraints
                for ix=1:Nx
                    for ixi=1:Nxi
                        rowIndices = getIndices(1,ix,ixi,1:Ntheta);
                        
                        colIndex = getIndices(4,0,0,0);
                        addSparseBlock(rowIndices, colIndex, expx2(ix)*ones(Ntheta,1))
                        colIndex = getIndices(5,0,0,0);
                        addSparseBlock(rowIndices, colIndex, expx2(ix)*x2(ix)*ones(Ntheta,1))
                    end
                end
                
                thetaPartOfConstraint = thetaWeights .* oneOverqRbDotGradThetas ./ bs;
                xPartOfDensityConstraint = xWeights;
                xPartOfPressureConstraint = xWeights .* x2;
                for ix=1:Nx
                    for ixi=1:Nxi
                        colIndices = getIndices(1,ix,ixi,1:Ntheta);
                        
                        rowIndex = getIndices(4,0,0,0);
                        addSparseBlock(rowIndex, colIndices, xPartOfDensityConstraint(ix)*thetaPartOfConstraint)
                        rowIndex = getIndices(5,0,0,0);
                        addSparseBlock(rowIndex, colIndices, xPartOfPressureConstraint(ix)*thetaPartOfConstraint)
                    end
                end
            end
            
            switch whichMatrixToMake
                case 1
                    fprintf('Time to contruct main matrix: %g seconds.\n',toc(matrixStartTime))
                    tic
                    matrix = createSparse();
                    fprintf('Time to sparsify main matrix: %g seconds.\n',toc)
                case 2
                    fprintf('Time to contruct preconditioner: %g seconds.\n',toc(matrixStartTime))
                    tic
                    preconditionerMatrix = createSparse();
                    fprintf('Time to sparsify preconditioner: %g seconds.\n',toc)
            end
        end
        
        
        % *****************************
        % Finalize matrix
        % *****************************
        
        fprintf('Predicted number of nonzeros: %d,  actual: %d,  fraction of nonzero entries: %g\n', estimated_nnz, nnz(matrix), nnz(matrix)/numel(matrix))
        
        if tryIterativeSolvers
            fprintf('LU-factorizing preconditioner...')
            tic
            [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditionerMatrix);
            fprintf('done.  Took %g seconds.\n',toc)
        end
        
        
        function solnVector=preconditioner(rhsVector)
            solnVector = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * rhsVector)));
        end
        
        if makeSpyPlot
            figure(2+figureOffset)
            clf
            if tryIterativeSolvers
                subplot(1,2,1)
            end
            spy(matrix)
            hold on
            index = getIndices(2,1,1,1) - 0.5;
            plot([index, index], [1, matrixSize],'r')
            plot([1, matrixSize], [index, index],'r')
            index = getIndices(3,1,1,1) - 0.5;
            plot([index, index], [1, matrixSize],'r')
            plot([1, matrixSize], [index, index],'r')
            
            if tryIterativeSolvers
                title('Actual matrix')
                subplot(1,2,2)
                spy(preconditionerMatrix)
                hold on
                index = getIndices(2,1,1,1) - 0.5;
                plot([index, index], [1, matrixSize],'r')
                plot([1, matrixSize], [index, index],'r')
                index = getIndices(3,1,1,1) - 0.5;
                plot([index, index], [1, matrixSize],'r')
                plot([1, matrixSize], [index, index],'r')
                title('Preconditioner matrix')
            end
        end
        
        % *******************************************
        % *******************************************
        %
        % Solve!
        %
        % *******************************************
        % *******************************************
        
        if tryIterativeSolvers
            if KrylovFigureHandle == 0
                KrylovFigureHandle = figure(3 + figureOffset);
            else
                set(0, 'CurrentFigure', KrylovFigureHandle);
            end
        end
        
        [soln, didItConverge, residual] ...
            = solveLinearSystem(matrix, rhs, @preconditioner, ...
            tryIterativeSolvers, orderOfSolversToTry, tol, maxIterations, restart, ...
            0, tryDirectSolverIfIterativeSolversFail);

            
        % **********************************************
        % **********************************************
        %
        % Calculate moments of the distribution function
        %
        % **********************************************
        % **********************************************
        
        VPrime = (oneOverqRbDotGradTheta(theta) ./ b(theta))*thetaWeights';
        
        if species==0
            % Ions
            particleFluxBeforeZetaIntegral=zeros(Ntheta,1);
            qBeforeZetaIntegral=zeros(Ntheta,1);
            flowDividedByB = zeros(1,itheta);
            density = zeros(1,itheta);
            pressure = zeros(1,itheta);
            
            particleFluxIntegralWeight = (x.^4);
            qIntegralWeight = (x.^6);
            kIntegralWeight = (x.^3);
            densityIntegralWeight = 4/sqrtpi*(x.*x);
            pressureIntegralWeight = 4/sqrtpi*(x.^4);
            for itheta=1:Ntheta
                for ixi=1:Nxi
                    xiWeight = xiWeights(ixi);
                    indices = getIndices(1,1:Nx,ixi,itheta);
                    fSlice = soln(indices);
                    qBeforeZetaIntegral(itheta) = qBeforeZetaIntegral(itheta) + xiWeight*(1+xi(ixi)*xi(ixi))*xWeights*(qIntegralWeight' .* fSlice);
                    particleFluxBeforeZetaIntegral(itheta) = particleFluxBeforeZetaIntegral(itheta) + xiWeight*(1+xi(ixi)*xi(ixi))*xWeights*(particleFluxIntegralWeight' .* fSlice);
                    density(itheta) = density(itheta) + xiWeight*xWeights*(densityIntegralWeight' .* fSlice);
                    pressure(itheta) = pressure(itheta) + xiWeight*xWeights*(pressureIntegralWeight' .* fSlice);
                    flowDividedByB(itheta) = flowDividedByB(itheta) + xiWeight*xi(ixi)*xWeights*(kIntegralWeight' .* fSlice);
                end
            end
            flowDividedByB = 4/sqrtpi*flowDividedByB./bs;
            %flowDividedByB = 2/3*4/sqrtpi*flowDividedByB./bs;
            qBeforeZetaIntegral = qBeforeZetaIntegral .* (dbdthetas ./ (bs.*bs.*bs))';
            particleFluxBeforeZetaIntegral = particleFluxBeforeZetaIntegral .* (dbdthetas ./ (bs.*bs.*bs))';
            k = 1/(2*pi)*thetaWeights * flowDividedByB';
            avgVParB = (1/VPrime) * thetaWeights * (flowDividedByB .* bs .* oneOverqRbDotGradTheta(theta))';
            q= sqrt(2*Miller_A/pi)/(nuPrime*VPrime) * thetaWeights * qBeforeZetaIntegral;
            particleFlux= sqrt(2*Miller_A/pi)/(nuPrime*VPrime) * thetaWeights * particleFluxBeforeZetaIntegral;
            fprintf('Normalized radial heat flux k_q: %g\n',q)
            fprintf('Normalized radial particle flux: %g  (It should be << 1.)\n',particleFlux)
            fprintf('Parallel flow coefficient k_||: %g\n',k)
            %fprintf('Parallel flow coefficient k_||: %g,   <V_|| B>: %g\n',k, avgVParB)
            fprintf('Poloidal variation of k_||: %g  (It should be << 1.)\n', max(abs(k-flowDividedByB)))
            
            if runMode==0 && plotPoloidalVariationOfFlow
                figure(6+figureOffset)
                clf
                plot(theta,flowDividedByB)
                xlabel('\theta')
                ylabel('parallel flow coefficient k||')
            end
            
        else
            % Electrons
            
            kIntegralWeight = (x.^3);
            
            % Compute L31:
            L31BeforeThetaIntegral = zeros(Ntheta,1);
            for itheta=1:Ntheta
                for ixi=1:Nxi
                    %L=1;
                    %indices = ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    indices = getIndices(1,1:Nx,ixi,itheta);
                    fSlice = soln(indices,2);
                    L31BeforeThetaIntegral(itheta) = L31BeforeThetaIntegral(itheta) + xiWeights(ixi)*xi(ixi)*xWeights*(kIntegralWeight' .* fSlice);
                end
            end
            
            %L31BeforeThetaIntegral = -(4/sqrtpi)*(2/3)/VPrime*oneOverqRbDotGradThetas'.*L31BeforeThetaIntegral;
            L31BeforeThetaIntegral = -(4/sqrtpi)/VPrime*oneOverqRbDotGradThetas'.*L31BeforeThetaIntegral;
            L31 = thetaWeights * L31BeforeThetaIntegral;
            
            % Compute L32:
            L32BeforeThetaIntegral = zeros(Ntheta,1);
            for itheta=1:Ntheta
                for ixi=1:Nxi
                    %L=1;
                    %indices = ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    indices = getIndices(1,1:Nx,ixi,itheta);
                    fSlice = soln(indices,1);
                    L32BeforeThetaIntegral(itheta) = L32BeforeThetaIntegral(itheta) + xiWeights(ixi)*xi(ixi)*xWeights*(kIntegralWeight' .* fSlice);
                end
            end
            
            %L32BeforeThetaIntegral = -(4/sqrtpi)*(2/3)/VPrime*oneOverqRbDotGradThetas'.*L32BeforeThetaIntegral;
            L32BeforeThetaIntegral = -(4/sqrtpi)/VPrime*oneOverqRbDotGradThetas'.*L32BeforeThetaIntegral;
            L32 = thetaWeights * L32BeforeThetaIntegral;
            
            % Compute L34:
            L34BeforeThetaIntegral = zeros(Ntheta,1);
            for itheta=1:Ntheta
                for ixi=1:Nxi
                    %L=1;
                    %indices = ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    indices = getIndices(1,1:Nx,ixi,itheta);
                    fSlice = soln(indices,3);
                    L34BeforeThetaIntegral(itheta) = L34BeforeThetaIntegral(itheta) + xiWeights(ixi)*xi(ixi)*xWeights*(kIntegralWeight' .* fSlice);
                end
            end
            
            %L34BeforeThetaIntegral = (4/sqrtpi)*(2/3)/(FSAB2*VPrime)*oneOverqRbDotGradThetas'.*L34BeforeThetaIntegral;
            L34BeforeThetaIntegral = (4/sqrtpi)/(FSAB2*VPrime)*oneOverqRbDotGradThetas'.*L34BeforeThetaIntegral;
            L34 = thetaWeights * L34BeforeThetaIntegral;
            
            % Compute conductivity:
            ConductivityBeforeThetaIntegral = zeros(Ntheta,1);
            for itheta=1:Ntheta
                for ixi=1:Nxi
                    %L=1;
                    %indices = ((1:Nx)-1)*Nxi*Ntheta + L*Ntheta + itheta;
                    indices = getIndices(1,1:Nx,ixi,itheta);
                    fSlice = soln(indices,4);
                    ConductivityBeforeThetaIntegral(itheta) = ConductivityBeforeThetaIntegral(itheta) + xiWeights(ixi)*xi(ixi)*xWeights*(kIntegralWeight' .* fSlice);
                end
            end
            
            L11=1.9693;
            %ConductivityBeforeThetaIntegral = (nuPrime/L11)*(4/sqrtpi)*(2/3)/(FSAB2 *VPrime) ...
            ConductivityBeforeThetaIntegral = (nuPrime/L11)*(4/sqrtpi)/(FSAB2 *VPrime) ...
                *oneOverqRbDotGradThetas'.*ConductivityBeforeThetaIntegral;
            conductivity = thetaWeights * ConductivityBeforeThetaIntegral;
            
            fprintf('L31 bootstrap current coefficient: %g\n',L31)
            fprintf('L32 bootstrap current coefficient: %g\n',L32)
            fprintf('L34 bootstrap current coefficient: %g\n',L34)
            fprintf('Neoclassical conductivity, normalized to Spitzer: %g\n',conductivity)
        end
        
        
        fprintf('******************** Done ********************\n')

        
        
        
        
        
        
        % *****************************************************************
        % *****************************************************************
        % End of main program.
        % Below are utilities for building sparse matrices:
        % *****************************************************************
        % *****************************************************************
        
        function ii = getIndices(whichBlock, ix, ixi, itheta)
            % Options for whichBlock:
            % 1 = distribution function and drift-kinetic equation
            % 2 = H Rosenbluth potential and Poisson equation for H
            % 3 = G Rosenbluth potential and Poisson equation for G
            % 4 = density source and constraint <n_1> = 0.
            % 5 = heat source and constraint <p_1> = 0.
            
            % Order of rows of the matrix and of the right-hand side:
            % -------------------------------------------------------
            % for ix = 1:Nx
            %   for ixi = 1:Nxi
            %     for itheta = 1:Ntheta
            %       Enforce drift-kinetic equation
            % for ix = 1:NxPotentials
            %   for ixi = 1:NL
            %     for itheta = 1:Ntheta
            %       Enforce Poisson equation for the H Rosenbluth potential
            % for ix = 1:NxPotentials
            %   for ixi = 1:NL
            %     for itheta = 1:Ntheta
            %       Enforce Poisson equation for the G Rosenbluth potential
            
            % Order of columns of the matrix and of the vector of unknowns:
            % -------------------------------------------------------
            % for ix = 1:Nx
            %   for ixi = 1:Nxi
            %     for itheta = 1:Ntheta
            %       distribution function f(ix,ixi,itheta)
            % for ix = 1:NxPotentials
            %   for ixi = 1:NL
            %     for itheta = 1:Ntheta
            %       H Rosenbluth potential: H(ix,ixi,itheta)
            % for ix = 1:NxPotentials
            %   for ixi = 1:NL
            %     for itheta = 1:Ntheta
            %       G Rosenbluth potential: G(ix,ixi,itheta)
            
            % Note that in this getIndices function, ixi is 1-based for the
            % potentials, so the Legendre index L is related to ixi by 
            % ixi = L + 1
            
            if whichBlock == 4
                ii = Nx*Nxi*Ntheta + 2*NxPotentials*NL*Ntheta + 1;
                return
            elseif whichBlock == 5
                ii = Nx*Nxi*Ntheta + 2*NxPotentials*NL*Ntheta + 2;
                return
            end
            
            if numel(ix)>1 && numel(ixi)>1 && numel(itheta) > 1
                error('No more than 1 of ix, ixi, and itheta can have multiple elements')
            end
            if any(ix<1)
                error('ix must be at least 1')
            end
            if any(ixi<1)
                error('ixi must be at least 1')
            end
            if any(itheta<1)
                error('itheta must be at least 1')
            end
            if any(itheta > Ntheta)
                error('itheta cannot exceed Ntheta')
            end
            switch whichBlock
                case 1
                    if any(ix > Nx)
                        error('ix cannot exceed Nx')
                    end
                    if any(ixi > Nxi)
                        error('ixi cannot exceed Nxi')
                    end
                    ii = (ix-1)*Nxi*Ntheta + (ixi-1)*Ntheta + itheta;
                case 2
                    if any(ix > NxPotentials)
                        error('ix cannot exceed NxPotentials')
                    end
                    if any(ixi > NL)
                        error('ixi cannot exceed NL')
                    end
                    ii = Nx*Nxi*Ntheta + (ix-1)*NL*Ntheta + (ixi-1)*Ntheta + itheta;
                case 3
                    if any(ix > NxPotentials)
                        error('ix cannot exceed NxPotentials')
                    end
                    if any(ixi > NL)
                        error('ixi cannot exceed NL')
                    end
                    ii = Nx*Nxi*Ntheta + NxPotentials*NL*Ntheta + (ix-1)*NL*Ntheta + (ixi-1)*Ntheta + itheta;
                otherwise
                    error('Block must be 1, 2, or 3.')
            end
        end
        
        function resetSparseCreator()
            sparseCreatorIndex=1;
            sparseCreator_i=zeros(estimated_nnz,1);
            sparseCreator_j=zeros(estimated_nnz,1);
            sparseCreator_s=zeros(estimated_nnz,1);
        end
        
        function addToSparse(i,j,s)
            n=numel(i);
            if n ~= numel(j)
                error('Error A');
            end
            if n ~= numel(s)
                error('Error B');
            end
            if any(i<1)
                error('Error Q: i<1');
            end
            if any(j<1)
                error('Error Q: j<1');
            end
            sparseCreator_i(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = i;
            sparseCreator_j(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = j;
            sparseCreator_s(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = s;
            sparseCreatorIndex = sparseCreatorIndex+n;
            if sparseCreatorIndex > estimated_nnz
                %fprintf('Error! estimated_nnz is too small.\n')
            end
        end
        
        function addSparseBlock(rowIndices, colIndices, block)
            s=size(block);
            if (s(1) ~= numel(rowIndices)) || (s(2) ~= numel(colIndices))
                s
                size(rowIndices)
                size(colIndices)
                error('Error in addSparseBlock!')
            end
            [rows, cols, values] = find(block);
            addToSparse(rowIndices(rows),colIndices(cols),values)
        end
        
        function sparseMatrix = createSparse()
            sparseMatrix = sparse(sparseCreator_i(1:(sparseCreatorIndex-1)), sparseCreator_j(1:(sparseCreatorIndex-1)), sparseCreator_s(1:(sparseCreatorIndex-1)), matrixSize, matrixSize);
            resetSparseCreator()
        end
        
    end

    function zz = RHat(theta)
        zz = 1 + (1/Miller_A)*cos(theta + Miller_x*sin(theta));
    end

    function zz = ZHat(theta)
        zz = (Miller_kappa/Miller_A)*sin(theta);
    end

    function zz = QIntegrand(theta)
        zz = ((1+Miller_s_kappa)*sin(theta + Miller_x*sin(theta)) .* (1+Miller_x*cos(theta)) .* sin(theta) ...
            + cos(theta) .* (Miller_dRdr + cos(theta + Miller_x *sin(theta)) - Miller_s_delta*sin(theta + Miller_x*sin(theta)) .* sin(theta))) ./ RHat(theta);
    end

    function zz = BPoloidal(theta)
        zz = Miller_Q/(Miller_kappa*Miller_q)*sqrt((sin(theta+Miller_x*sin(theta)) .* (1+Miller_x*cos(theta))).^2 + (Miller_kappa*cos(theta)).^2) ...
            ./ (RHat(theta) .* ( cos(Miller_x*sin(theta)) + Miller_dRdr*cos(theta) + (Miller_s_kappa-Miller_s_delta*cos(theta) + (1+Miller_s_kappa)*Miller_x*cos(theta)) .* sin(theta) .* sin(theta + Miller_x*sin(theta))));
    end

    function zz = BDotGradTheta(theta)
        zz = - Miller_A*Miller_Q./(Miller_kappa*Miller_q*RHat(theta) .* ...
            ((1+Miller_s_kappa)*sin(theta + Miller_x*sin(theta)) .* (1+Miller_x*cos(theta)) .* sin(theta) ...
            + cos(theta) .* (Miller_dRdr + cos(theta + Miller_x *sin(theta)) - Miller_s_delta*sin(theta + Miller_x*sin(theta)) .* sin(theta))));
    end


    function zz = b(theta)
        switch geometry
            case 0
                zz = 1./(1 + cos(theta)/Miller_A);
            case 1
                zz = sqrt(BPoloidal(theta).^2 + 1./((RHat(theta)).^2));
            case 2
                zz = ones(size(theta)) * bModes(1);
                for m=1:numModesToKeep
                    zz = zz + 2*cos(m*theta)*real(bModes(m+1,psiIndex)) - 2*sin(m*theta)*imag(bModes(m+1,psiIndex));
                end
            otherwise
                error('Invalid geometry')
        end
    end

    function zz = dbdtheta(theta)
        switch geometry
            case 0
                zz = sin(theta) ./(Miller_A * (1 + cos(theta)/Miller_A).^2);
            case 1
                error('Program should not get here!')
            case 2
                zz = zeros(size(theta)) * bModes(1);
                for m=1:numModesToKeep
                    zz = zz - 2*m*sin(m*theta)*real(bModes(m+1,psiIndex)) - 2*m*cos(m*theta)*imag(bModes(m+1,psiIndex));
                end
            otherwise
                error('Invalid geometry')
        end
    end

    function zz = oneOverqRbDotGradTheta(theta)
        switch geometry
            case 0
                zz = ones(size(theta));
            case 1
                zz = b(theta) ./ (Miller_q*BDotGradTheta(theta));
            case 2
                zz = interp1(thetaData{psiIndex}, BData{psiIndex} ./ ( R0 * qData(psiIndex) * BDotGradThetaData{psiIndex}), theta, 'spline');
            otherwise
                error('Invalid geometry')
        end
    end

    function [kpar, conductivity, L31, L32, L34] = computeSauterK(nuStarsMe)
        Z = 1;
        
        VPrime = quadgk(@(theta) (oneOverqRbDotGradTheta(theta)./b(theta)), 0, 2*pi);
        FSAB2 = quadgk(@(theta) (oneOverqRbDotGradTheta(theta).*b(theta)), 0, 2*pi) / VPrime;
        
        if geometry==2
            [~,lambdaMax] = fminbnd(@(x) (1./b(x)), pi/2, 3*pi/2);
        else
            lambdaMax = 1/b(pi);
        end
        lambda=0;
        
        fc = 0.75*FSAB2 * quadgk(@outerIntegrand, 0, lambdaMax);
        if abs(imag(fc)/real(fc)) > 1e-2
            error('Something went wrong. fc has an imaginary part that is not negligible.')
        end
        fc=real(fc);
        
        ft=1-fc;
        alpha0 = -1.17*(1-ft)/(1-0.22*ft-0.19*ft*ft);
        nuStarsSauter = nuStarsMe;
        alpha = ((alpha0 + 0.25*(1-ft*ft)*sqrt(nuStarsSauter))./(1+0.5*sqrt(nuStarsSauter)) + 0.315*(nuStarsSauter.^2)*(ft^6)) ...
            ./ (1+0.15*(nuStarsSauter.^2) * (ft^6));
        kpar = -alpha;
        
        conductivity = zeros(size(nuStarsMe));
        L31 = zeros(size(nuStarsMe));
        L32 = zeros(size(nuStarsMe));
        L34 = zeros(size(nuStarsMe));
        
        for iNu = 1:numel(nuStarsMe)
            nuStar = sqrt(2)*nuStarsMe(iNu);
            X = ft/(1 + (0.55-0.1*ft)*sqrt(nuStar) + 0.45*(1-ft)*nuStar);
            conductivity(iNu) = 1 - 1.36*X + 0.59*X*X - 0.23*X*X*X;
            
            X = ft / (1 + (1-0.1*ft)*sqrt(nuStar) + 0.5*(1-ft)*nuStar/Z);
            L31(iNu) = (1 + 1.4/(Z+1))*X - 1.9/(Z+1)*X*X + 0.3/(Z+1)*X*X*X + 0.2/(Z+1)*X*X*X*X;
            
            X = ft / (1 + 0.26*(1-ft)*sqrt(nuStar) + 0.18*(1-0.37*ft)*nuStar/Z);
            Y = ft / (1 + (1+0.6*ft)*sqrt(nuStar) + 0.85 * (1-0.37*ft)*nuStar*(1+Z));
            F32_ee = (0.05+0.62*Z)/(Z*(1+0.44*Z))*(X-X^4) + 1/(1+0.22*Z)*(X*X-X^4 - 1.2*(X^3-X^4)) + 1.2/(1+0.5*Z)*X^4;
            F32_ei = -(0.56+1.93*Z)/(Z*(1+0.44*Z))*(Y-Y^4) + 4.95/(1+2.48*Z)*(Y*Y-Y^4 - 0.55*(Y^3-Y^4)) - 1.2/(1+0.5*Z)*Y^4;
            L32(iNu) = F32_ee + F32_ei;
            
            X = ft / (1 + (1-0.1*ft)*sqrt(nuStar) + 0.5*(1-0.5*ft)*nuStar/Z);
            L34(iNu) = (1 + 1.4/(Z+1))*X - 1.9/(Z+1)*X*X + 0.3/(Z+1)*X*X*X + 0.2/(Z+1)*X*X*X*X;
        end
        
    end
    function yy = outerIntegrand(lambdas)
        yy = zeros(size(lambdas));
        for i=1:numel(lambdas)
            lambda=lambdas(i);
            yy(i) = lambda / quadgk(@innerIntegrand, 0, 2*pi);
        end
    end

    function y=innerIntegrand(thetas)
        y = (sqrt(1-lambda*b(thetas)) .* oneOverqRbDotGradTheta(thetas) ./ b(thetas)) / VPrime;
    end

    function plotFigure()
        figure(1)
        clf
        
        numRows=2;
        numCols=3;
        
        nuStarMin = min(nuStars);
        nuStarMax = max(nuStars);
        
        subplot(numRows,numCols,1)
        nuStarsFine = logspace(log10(nuStars(1)), log10(nuStars(end)));
        [SauterKs, SauterConductivity, SauterL31, SauterL32, SauterL34] = computeSauterK(nuStarsFine);
        semilogx(nuStarsFine, SauterKs,'k')
        hold on
        semilogx(nuStars,flowCoefficients(1,:),'.-')
        semilogx(nuStars,flowCoefficients(2,:),'.--g')
        semilogx(nuStars,flowCoefficients(3,:),'.-.r')
        semilogx(nuStars,flowCoefficients(4,:),'.:m')
        xlabel('\nu_*')
        title('flow coefficient k_{||}')
        xlim([nuStarMin,nuStarMax])
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx')
        
        subplot(numRows,numCols,2)
        semilogx(nuStarsFine, SauterConductivity,'k')
        hold on
        semilogx(nuStars,conductivities(1,:),'.-')
        semilogx(nuStars,conductivities(2,:),'.--g')
        semilogx(nuStars,conductivities(3,:),'.-.r')
        semilogx(nuStars,conductivities(4,:),'.:m')
        xlabel('\nu_*')
        title('conductivity / Spitzer')
        xlim([nuStarMin,nuStarMax])
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx')
        
        subplot(numRows,numCols,3)
        semilogx(nuStarsFine, SauterL31,'k')
        hold on
        semilogx(nuStars,L31s(1,:),'.-')
        semilogx(nuStars,L31s(2,:),'.--g')
        semilogx(nuStars,L31s(3,:),'.-.r')
        semilogx(nuStars,L31s(4,:),'.:m')
        xlabel('\nu_*')
        title('L_{31}')
        xlim([nuStarMin,nuStarMax])
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx')
        
        subplot(numRows,numCols,4)
        semilogx(nuStarsFine, SauterL32,'k')
        hold on
        semilogx(nuStars,L32s(1,:),'.-')
        semilogx(nuStars,L32s(2,:),'.--g')
        semilogx(nuStars,L32s(3,:),'.-.r')
        semilogx(nuStars,L32s(4,:),'.:m')
        xlabel('\nu_*')
        title('L_{32}')
        xlim([nuStarMin,nuStarMax])
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx')
        
        subplot(numRows,numCols,5)
        semilogx(nuStarsFine, SauterL34,'k')
        hold on
        semilogx(nuStars,L34s(1,:),'.-')
        semilogx(nuStars,L34s(2,:),'.--g')
        semilogx(nuStars,L34s(3,:),'.-.r')
        semilogx(nuStars,L34s(4,:),'.:m')
        xlabel('\nu_*')
        title('L_{34}')
        xlim([nuStarMin,nuStarMax])
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx')
        
        
    end
end
