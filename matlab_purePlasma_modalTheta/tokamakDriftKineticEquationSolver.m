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
% This program also computes the neoclassical ion thermal conductivity, 
% i.e. k_q in eq (5) in
% Landreman & Ernst, Plasma Phys Controlled Fusion 54, 115006 (2012).
% Notice the ion flow coefficient k_|| in eq (6) of this paper is equal 
% to -alpha from the Sauter paper.

%******************************************************************
% Options for controlling the general program execution:
%******************************************************************

runMode = 0;
% 0 = Single run.
% 1 = Do a bunch of solves, scanning the numerical resolution parameters to check convergence.
% 2 = Do a bunch of solves, scanning the collisionality, and saving the results.
% 3 = Replot a saved collisionality scan, without doing any solves.

% The following parameter specifies the file to plot when runMode = 3:
dataFileToPlot = 'tokamakDriftKineticEquationSolver_2014-03-27_14-15_foobar.mat';

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
% For runMode = 2, the value of 'species' here is ignored, and both ion and electron equations are solved.

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

% Number of Fourier modes in the poloidal angle.
% Memory and time requirements do depend strongly on this parameter.  The
% value of Ntheta required for numerical convergence depends strongly on
% collisionality. At high collisionality, Ntheta = 5 might be plenty. At nu* =
% 10^{-3}, Ntheta > 30 may be required.
NthetaConverged = 8;
Nthetas = floor(linspace(7,15,5));

% Number of Legendre modes retained for the Rosenbluth potentials.
% The memory and time requirements are not very sensitive to this
% parameter.  A value of 4 almost always works well.
NLConverged = 4;
NLs = 2:2:6;

% Number of Legendre modes retained for the distribution function.
% Memory and time requirements do depend strongly on this parameter.  The
% value of Nxi required for numerical convergence depends strongly on
% collisionality. At high collisionality, Nxi=5 might be plenty. At nu* =
% 10^{-3}, Nxi > 100 may be required.
NxiConverged = 25;
Nxis = floor(linspace(10,40,7));

% Number of grid points in x = v / v_th used for the distribution function.
% Typically, values in the range 5-8 work well, with slightly higher values
% required at higher collisionality.
NxConverged = 7;
Nxs = 6:12;

% Number of grid points in x used for the Rosenbluth potentials.
% Memory and time requirements are not sensitive to this parameter. A value
% of 40 usually works well.
NxPotentialsPerVthConverged = 40;
NxPotentialsPerVths = [40, 81];

% Tolerance for iterative solver. Ignored when the direct solver is used.
% Memory requirements do not depend on this parameter, and the time
% required by the iterative solver depends somewhat on this parameter.
% Typically, values in the range 5-7 are good.
log10tolConverged = 5.0;
log10tols = 4.5:0.5:5;


%******************************************************************
% Other numerical parameters:
%******************************************************************

%tryIterativeSolvers = true;
tryIterativeSolvers = false;
% If false, the sparse direct solver will be used.

tryDirectSolverIfIterativeSolversFail = true;
%tryDirectSolverIfIterativeSolversFail = false;

% The setting below determines the order of iterative solvers to try:
orderOfSolversToTry = [1, 3];
% 1 = GMRES
% 2 = BiCGStab
% 3 = BiCGStab(l)
% 4 = TFQMR
% 5 = CGS

% Maximum number of iterations allowed in the iterative solver:
maxIterations = 100;

% The parameter below only matters for the GMRES solver:
restart = 100;

% The switch below, if set to true, may speed the code up a tiny bit, but
% it must be set to false when the magnetic geometry is not perfectly
% up-down symmetric. Since the gain is minimal from setting this parameter
% to true, I'd recommend keeping it false to be safe.
%zeroModesWithWrongParity = true;
zeroModesWithWrongParity = false;

% The parameter below may have some effect on the memory required to
% factorize the matrix, but in practice it doesn't seem to make much
% difference.
whereToPutbDotGradTheta = 1;
% 0 = my traditional form of the DKE, with b dot grad theta on the
% collision term but not streaming or mirror terms: this makes the
% streaming term diagonal in M but makes the collision term dense in M.

% 1 = put the b dot grad theta on the streaming and mirror terms. This
% makes the collision term diagonal in M but makes the streaming term dense
% in M.

%******************************************************************
% Switches for plotting:
%******************************************************************

% The number below is added to each figure number. You can change it if you
% want to save windows from a previous run.
figureOffset=20;

plotVelocitySpaceGrid = true;
%plotVelocitySpaceGrid = false;

plotEFITDetails = true;
%plotEFITDetails = false;

%plotPoloidalVariationOfFlow = true;
plotPoloidalVariationOfFlow = false;

% ***********************************************
% ***********************************************
%
% End of the main input parameters
%
% ***********************************************
% ***********************************************

speedGridFigureHandle = 0;
KrylovFigureHandle = 0;


if geometry==2 && runMode ~= 3
    % Load EFIT data
    
    zeroModesWithWrongParity = false;
    
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

if plotVelocitySpaceGrid
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
        Nthetas =     [ 8, 8, 8, 8, 8, 8, 8, 8, 8, 9,10,11,12,13,14,15,25,32,32,32,32];
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
        
        NL=4;
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
        stringForTop=sprintf('%s, %s, aspect ratio = %g, nuPrime = %g, nu_* = %g, theta modal, Legendre modal, polynomial colocation in x',speciesText,geometryText,Miller_A, nuPrime,nuStar);
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

    function solveDKE()
        
        
        startTimeForThisRun=tic;
        
        sqrtpi=sqrt(pi);
        iteration = iteration+1;
        
        
        % Order of the rows of the matrix and of the RHS:
        % --------------------------------
        % for x = dx to xMax-dx  (Nx points)
        %   for L = 0:(NL-1)
        %     for M = 0:(Ntheta-1)
        %       Real part of the DKE (equivalent to -M)
        %     for M = 1:(Ntheta-1)
        %       Imag part of the DKE
        
        % Order of the cols of the matrix and of the soln:
        % --------------------------------
        % for x = dx to xMax-dx  (Nx points)
        %   for l = 0:(NL-1)
        %     for m = 0:(Ntheta-1)
        %       Real part of f_(l,m)
        %     for m = 1:(Ntheta-1)
        %       Imag part of f_(l,m)
        
        
        fprintf('Ntheta = %d,  NL = %d,  Nxi = %d,  Nx = %d, NxPtentialsPerVth = %g, tol = %g\n',Ntheta,NL,Nxi,Nx,NxPotentialsPerVth,tol)
        
        tic
        
        % NOTE: Ntheta is the number of MODES in theta.
        % The number of grid points in theta is 2*Ntheta.
        
        theta = linspace(0, 2*pi, 2*Ntheta+1);
        theta(end)=[];
        
        bs=b(theta);
        oneOverqRbDotGradThetas = oneOverqRbDotGradTheta(theta);
        if geometry==1
            % Miller: too hard to analytically differentiate b(theta) so do
            % it numerically:
            thetaMin=0;
            thetaMax = 2*pi;
            [thetaFine, ddthetaFine, NthetaFinal] = makeUniformSpectralDifferentiationMatrix(2*Ntheta*dBdthetaResolutionMultiplier, thetaMin, thetaMax);
            thetaFine = circshift(thetaFine, [0,1]);
            dbdthetaFine = ddthetaFine * b(thetaFine)';
            dbdthetas = dbdthetaFine(1:dBdthetaResolutionMultiplier:end)';
        else
            dbdthetas = dbdtheta(theta);
        end
        
        % Generate abscissae, quadrature weights, and derivative matrix for x grid.
        kk=0;
        scale=1;
        pointAtZero=false;
        [x, ddx, d2dx2, xWeights] = spectralNodesWeightsAndDifferentiationMatricesForV(Nx, kk, scale, pointAtZero);
        x_i = x(:)';
        xWeights = xWeights(:)';
        
        function y=weight(xxx)
            y=exp(-xxx.*xxx);
        end
        
        xMax=max([5, max(x_i)]);
        NxPotentials = ceil(xMax * NxPotentialsPerVth);
        % Uniform finite differences with 5-point stencil.
        xMin=0;
        scheme = 12;
        [xPotentials, ~, ddxPotentials, d2dx2Potentials] = differentiationMatricesForUniformGrid(NxPotentials, xMin, xMax, scheme);
        regridPolynomialToUniform = polynomialInterpolationMatrix(x_i,xPotentials,weight(x_i),weight(xPotentials));
        regridUniformToPolynomial = makeHighOrderInterpolationMatrix(xPotentials,x_i,0,'f');
        
        
        if plotVelocitySpaceGrid
            if speedGridFigureHandle == 0
                speedGridFigureHandle = figure(figureOffset+7);
            else
                set(0, 'CurrentFigure', speedGridFigureHandle);
            end
            plot(xPotentials,zeros(size(xPotentials))+iteration,'.r')
            hold on
            plot(x_i, zeros(size(x_i))+iteration,'o')
            title('Speed grid for distribution function (blue) and Rosenbluth potentials(red)')
            xlabel('x')
            ylabel('Solve number')
        end
        
        thetaGridSize = (Ntheta*2-1);
        matrixSize = thetaGridSize * Nx * Nxi;
        estimated_nnz = ceil(1.1*((nnz(generateThetaMultiplicationMatrix(bs, 1e-7))+nnz(generateThetaMultiplicationMatrix(oneOverqRbDotGradThetas, 1e-7)))*4*Nxi*Nx + thetaGridSize*(NL*Nx*Nx + Nxi*Nx*Nx)));
        fprintf('matrixSize: %d.\n',matrixSize)
        
        isModeAllowedSmall = zeros(1,thetaGridSize*Nxi);
        evens = [ones(1,Ntheta), zeros(1,Ntheta-1)];
        odds = 1-evens;
        for L=0:(Nxi-1)
            indices = L*thetaGridSize + (1:thetaGridSize);
            if mod(L,2)==0
                isModeAllowedSmall(indices) = odds;
            else
                isModeAllowedSmall(indices) = evens;
            end
        end
        isModeAllowedVector = repmat(isModeAllowedSmall,[1,Nx])';
        isModeAllowed = spdiags(isModeAllowedVector, 0, matrixSize, matrixSize);
        thisModeNotAllowed = spdiags(1-isModeAllowedVector, 0, matrixSize, matrixSize);
        
        % Begin timer for matrix construction:
        tic
        
        
        % *************************************
        % *************************************
        %
        % Build the right-hand side of the main linear system
        %
        % *************************************
        % *************************************
        
        switch species
            case 0
                % Ions
                rhs=zeros(matrixSize,1);
            case 1
                % Electrons
                rhs=zeros(matrixSize,4);
            otherwise
                error('Invalid setting for species')
        end
        
        % Order of columns in the RHS:
        % 1 = L31
        % 2 = L32
        % 3 = L34
        % 4 = conductivity
        
        x2=x_i.*x_i;
        expx2=exp(-x2);
        xPartOfRHSForConductivity = x_i.*expx2;
        xPartOfRHSForP = x2.*expx2;
        xPartOfRHSForT = (x2-5/2).*x2.*expx2;
        switch whereToPutbDotGradTheta
            case 0
                thetaPartOfRHSOnThetaGrid = dbdthetas ./ (bs.*bs);
                thetaPartOfL34RHSOnThetaGrid = dbdthetas;
                thetaPartOfConductivityOnThetaGrid = bs .* oneOverqRbDotGradThetas;
            case 1
                thetaPartOfRHSOnThetaGrid = dbdthetas ./ (bs.*bs .* oneOverqRbDotGradThetas);
                thetaPartOfL34RHSOnThetaGrid = dbdthetas ./ oneOverqRbDotGradThetas;
                thetaPartOfConductivityOnThetaGrid = bs;
            otherwise
                error('Invalid whereToPutbDotGradTheta')
        end
        
        complexThetaModeCoefficientsForRHS = fft(thetaPartOfRHSOnThetaGrid)/(2*Ntheta);
        thetaModeCoefficientsForRHS = [real(complexThetaModeCoefficientsForRHS(1:Ntheta)), imag(complexThetaModeCoefficientsForRHS(2:Ntheta))];
        
        complexThetaModeCoefficientsForRHS = fft(thetaPartOfL34RHSOnThetaGrid)/(2*Ntheta);
        thetaModeCoefficientsForL34RHS = [real(complexThetaModeCoefficientsForRHS(1:Ntheta)), imag(complexThetaModeCoefficientsForRHS(2:Ntheta))];
        
        complexThetaModeCoefficientsForRHS = fft(thetaPartOfConductivityOnThetaGrid)/(2*Ntheta);
        thetaModeCoefficientsForConductivityRHS = [real(complexThetaModeCoefficientsForRHS(1:Ntheta)), imag(complexThetaModeCoefficientsForRHS(2:Ntheta))];
        
        for itheta=1:thetaGridSize
            L=0;
            indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
            rhs(indices,1) = 0.5 * (4/3) * thetaModeCoefficientsForRHS(itheta) * xPartOfRHSForT;
            if species ~= 0
                rhs(indices,2) = 0.5 * (4/3) * thetaModeCoefficientsForRHS(itheta) * xPartOfRHSForP;
                
                L=1;
                indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                rhs(indices,4) = thetaModeCoefficientsForConductivityRHS(itheta) * xPartOfRHSForConductivity;
            end
            
            L=2;
            indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
            rhs(indices,1) = 0.5 * (2/3) * thetaModeCoefficientsForRHS(itheta) * xPartOfRHSForT;
            if species ~= 0
                rhs(indices,2) = 0.5 * (2/3) * thetaModeCoefficientsForRHS(itheta) * xPartOfRHSForP;
                rhs(indices,3) = thetaModeCoefficientsForL34RHS(itheta) * xPartOfRHSForP;
            end
        end
        
        % *************************************
        % *************************************
        %
        % Build the main matrix
        %
        % *************************************
        % *************************************
        
        sparseCreatorIndex=1;
        sparseCreator_i=0;
        sparseCreator_j=0;
        sparseCreator_s=0;
        resetSparseCreator()
        
        Ms = 0:(Ntheta-1);
        ddtheta = -diag(Ms, Ntheta-1) + diag(Ms, 1-Ntheta);
        
        if tryIterativeSolvers
            matricesToMake=1:2;
        else
            matricesToMake=1;
        end
        
        for whichMatrixToMake = matricesToMake
            % 1 = main matrix
            % 2 = preconditioner
            
            switch whichMatrixToMake
                case 1
                    threshhold=1e-7;
                case 2
                    threshhold=1e-7;
            end
            switch whereToPutbDotGradTheta
                case 0
                    % 1/(b dot grad theta) multiplies collision term,
                    % simpler streaming term
                    thetaPartOfMirrorTermOnThetaGrid = -0.5*dbdthetas ./ (bs);
                    thetaPartMatrixForCollisionTerm = generateThetaMultiplicationMatrix(oneOverqRbDotGradThetas, threshhold);
                    thetaPartMatrixForStreamingTerm = ddtheta;
                case 1
                    % (b dot grad theta) multiplies streaming and mirror terms,
                    % simpler collision term
                    thetaPartOfMirrorTermOnThetaGrid = -0.5*dbdthetas ./ (bs .* oneOverqRbDotGradThetas);
                    thetaPartOfStreamingTermOnThetaGrid = 1./(oneOverqRbDotGradThetas);
                    thetaPartMatrixForStreamingTerm = generateThetaMultiplicationMatrix(thetaPartOfStreamingTermOnThetaGrid, threshhold) * ddtheta;
                otherwise
                    error('Invalid whereToPutbDotGradTheta')
            end
            thetaPartMatrixForMirrorTerm = generateThetaMultiplicationMatrix(thetaPartOfMirrorTermOnThetaGrid, threshhold);
            
            matrixStartTime = tic;
            
            % ***************************************
            % Add the streaming and mirror terms:
            % ***************************************
            for ix=1:Nx
                
                x=x_i(ix);
                for L=0:(Nxi-1)
                    rowIndices = (ix-1)*Nxi*thetaGridSize + L*thetaGridSize + (1:thetaGridSize);
                    
                    % super-diagonals in L: (streaming & mirror terms)
                    if L<(Nxi-1)
                        colIndices = rowIndices + thetaGridSize;
                        % Streaming term
                        addSparseBlock(rowIndices, colIndices, thetaPartMatrixForStreamingTerm*x*(L+1)/(2*L+3));
                        
                        % Mirror term
                        addSparseBlock(rowIndices, colIndices, thetaPartMatrixForMirrorTerm*x*(L+1)*(L+2)/(2*L+3));
                    end
                    
                    % Sub-diagonals in L: (streaming & mirror terms)
                    if L>0
                        colIndices = rowIndices - thetaGridSize;
                        % Streaming term
                        addSparseBlock(rowIndices, colIndices, thetaPartMatrixForStreamingTerm*x*L/(2*L-1));
                        
                        % Mirror term
                        addSparseBlock(rowIndices, colIndices, -thetaPartMatrixForMirrorTerm*x*(L-1)*L/(2*L-1));
                    end
                end
            end
            
            % ***************************************
            % Build the collision operator:
            % For a detailed explanation of the implementation of the collision operator, see
            % Landreman & Ernst, Journal of Computational Physics (2013)
            % http://dx.doi.org/10.1016/j.jcp.2013.02.041
            % ***************************************
            
            xWith0s = [0; xPotentials(2:(end-1)); 0];
            M21 = 4*pi*diag(xWith0s.^2) * regridPolynomialToUniform;
            M32 = -2*diag(xWith0s.^2);
            LaplacianTimesX2WithoutL = diag(xPotentials.^2)*d2dx2Potentials + 2*diag(xPotentials)*ddxPotentials;
            
            erfs=erf(x_i);
            x2 = x_i.*x_i;
            x3 = x2.*x_i;
            expx2 = exp(-x_i.*x_i);
            Psi = (erfs - 2/sqrtpi*x_i .* expx2) ./ (2*x_i.*x_i);
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
            PsiPrime = (-erfs + 2/sqrtpi*x_i.*(1+x_i.*x_i) .* expx2) ./ x3;
            xPartOfCECD = 3*sqrtpi/4*(diag(Psi./x_i)*d2dx2   +  diag((PsiPrime.*x_i  + Psi + 2*Psi.*x2)./x2)*ddx + diag(2*PsiPrime + 4*Psi./x_i)) + 3*diag(expx2);
            M12IncludingX0 = nuPrime * 3/(2*pi)*diag(expx2)*regridUniformToPolynomial;
            M13IncludingX0 = -nuPrime * 3/(2*pi) * diag(x2.*expx2) * regridUniformToPolynomial* d2dx2Potentials;
            
            for L=0:(Nxi-1)
                M11 = -nuPrime * (-0.5*diag(nuD)*L*(L+1) + xPartOfCECD);
                
                
                if L <= (NL-1)
                    % Add Rosenbluth potential stuff
                    
                    M13 = M13IncludingX0;
                    M12 = M12IncludingX0;
                    
                    M22 = LaplacianTimesX2WithoutL-L*(L+1)*eye(NxPotentials);
                    
                    % Add Dirichlet or Neumann boundary condition for
                    % potentials at x=0:
                    if L==0
                        M22(1,:)=ddxPotentials(1,:);
                    else
                        M22(1,:) = 0;
                        M22(1,1) = 1;
                        M12(:,1) = 0;
                        M13(:,1) = 0;
                    end
                    M33 = M22;
                    
                    % Add Robin boundary condition for potentials at x=xMax:
                    M22(NxPotentials,:) = xMax*ddxPotentials(NxPotentials,:);
                    M22(NxPotentials,NxPotentials) = M22(NxPotentials,NxPotentials) + L+1;
                    
                    M33(NxPotentials,:) = xMax*xMax*d2dx2Potentials(NxPotentials,:) + (2*L+1)*xMax*ddxPotentials(NxPotentials,:);
                    M33(NxPotentials,NxPotentials) = M33(NxPotentials,NxPotentials) + (L*L-1);
                    if L~=0
                        M22(NxPotentials,1)=0;
                        M33(NxPotentials,1)=0;
                    end
                    
                    KWithoutThetaPart = M11 -  (M12 - M13 * (M33 \ M32)) * (M22 \ M21);
                    
                    
                else
                    KWithoutThetaPart = M11;
                end
                if whichMatrixToMake ~= 1
                    KWithoutThetaPart = diag(diag(KWithoutThetaPart));
                end
                
                switch whereToPutbDotGradTheta
                    case 0
                        
                        % 1/(b dot grad theta) multiplies collision term.
                        for ithetaRow=1:thetaGridSize
                            rowIndices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + ithetaRow;
                            for ithetaCol = 1:thetaGridSize
                                colIndices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + ithetaCol;
                                addSparseBlock(rowIndices, colIndices, KWithoutThetaPart*thetaPartMatrixForCollisionTerm(ithetaRow,ithetaCol))
                            end
                        end
                        
                    case 1
                        % Simpler collision term
                        for itheta=1:thetaGridSize
                            indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                            addSparseBlock(indices, indices, KWithoutThetaPart)
                        end
                    otherwise
                        error('Invalid whereToPutbDotGradTheta')
                end
            end
            switch whichMatrixToMake
                case 1
                    fprintf('Time to contruct main matrix: %g seconds.\n',toc(matrixStartTime))
                    tic
                    if zeroModesWithWrongParity
                        matrix = isModeAllowed * createSparse() * isModeAllowed + thisModeNotAllowed;
                    else
                        matrix = createSparse();
                    end
                    fprintf('Time to sparsify main matrix: %g seconds.\n',toc)
                case 2
                    fprintf('Time to contruct preconditioner: %g seconds.\n',toc(matrixStartTime))
                    tic
                    if zeroModesWithWrongParity
                        preconditionerMatrix = isModeAllowed * createSparse() * isModeAllowed + thisModeNotAllowed;
                    else
                        preconditionerMatrix = createSparse();
                    end
                    fprintf('Time to sparsify preconditioner: %g seconds.\n',toc)
            end
        end
        
        
        
        % *****************************
        % Finalize matrix
        % *****************************
        
        fprintf('Predicted nnz: %d.  Actual nnz: %d.  Fraction of nonzero entries: %g\n', estimated_nnz, nnz(matrix), nnz(matrix)/numel(matrix))
        
        if tryIterativeSolvers
            fprintf('LU-factorizing preconditioner...')
            tic
            [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditionerMatrix);
            fprintf('done.  Took %g seconds.\n',toc)
        end
        
        
        function solnVector = preconditioner(rhsVector)
            solnVector = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * rhsVector)));
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
        
        
        % ************************************************
        % ************************************************
        %
        % Calculate moments of the distribution function:
        %
        % ************************************************
        % ************************************************
        
        thetaMatrixForFlows = generateThetaMultiplicationMatrix(1./bs, threshhold);
        thetaMatrixForBootstrapCoefficients = generateThetaMultiplicationMatrix(oneOverqRbDotGradThetas, threshhold);
        
        VPrime = quadgk(@(th) (oneOverqRbDotGradTheta(th) ./ b(th)), 0, 2*pi);
        FSAB2 = quadgk(@(th) (oneOverqRbDotGradTheta(th) .* b(th)), 0, 2*pi) / VPrime;
        fprintf('<B^2> = %g\n',FSAB2)
        
        if species==0
            % Ions
            
            particleFluxBeforeZetaIntegral=zeros(thetaGridSize,1);
            qBeforeZetaIntegral=zeros(thetaGridSize,1);
            flowDividiedByB = zeros(thetaGridSize,1);
            density = zeros(thetaGridSize,1);
            
            particleFluxIntegralWeight = (x_i.^4);
            qIntegralWeight = (x_i.^6);
            kIntegralWeight = (x_i.^3);
            densityIntegralWeight = 4*pi*(x_i.*x_i);
            for itheta=1:thetaGridSize
                L=0;
                indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                fSlice = soln(indices,1);
                qBeforeZetaIntegral(itheta) = 8/3*xWeights*(qIntegralWeight' .* fSlice);
                particleFluxBeforeZetaIntegral(itheta) = 8/3*xWeights*(particleFluxIntegralWeight' .* fSlice);
                density(itheta) = xWeights*(densityIntegralWeight' .* fSlice);
                
                L=1;
                indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                fSlice = soln(indices,1);
                flowDividiedByB(itheta) = xWeights*(kIntegralWeight' .* fSlice);
                
                L=2;
                indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                fSlice = soln(indices,1);
                qBeforeZetaIntegral(itheta) = qBeforeZetaIntegral(itheta) + 4/15*xWeights*(qIntegralWeight' .* fSlice);
                particleFluxBeforeZetaIntegral(itheta) = particleFluxBeforeZetaIntegral(itheta) + 4/15*xWeights*(particleFluxIntegralWeight' .* fSlice);
            end
            
            threshhold = 1e-10;
            flowDividiedByB = 2/3*4/sqrtpi*thetaMatrixForFlows*flowDividiedByB;
            thetaMatrixForFluxes = generateThetaMultiplicationMatrix(dbdthetas ./ (bs.*bs.*bs), threshhold);
            qBeforeZetaIntegral = thetaMatrixForFluxes * qBeforeZetaIntegral;
            particleFluxBeforeZetaIntegral = thetaMatrixForFluxes * particleFluxBeforeZetaIntegral;
            k = flowDividiedByB(1) * FSAB2;
            temp = generateThetaMultiplicationMatrix(bs .* oneOverqRbDotGradTheta(theta), threshhold) * flowDividiedByB;
            avgVParB = (1/VPrime) * 2*pi * temp(1);
            
            q= sqrt(2*Miller_A/pi)/(nuPrime*VPrime) * 2*pi * qBeforeZetaIntegral(1);
            particleFlux= sqrt(2*Miller_A/pi)/(nuPrime*VPrime) * 2*pi * particleFluxBeforeZetaIntegral(1);
            fprintf('Normalized radial heat flux k_q: %g\n',q)
            fprintf('Normalized radial particle flux: %g  (It should be nearly 0.)\n',particleFlux)
            fprintf('Parallel flow coefficient k_||: %g\n',k)
            %fprintf('Parallel flow coefficient k_||: %g,   <V_|| B>: %g\n',k, avgVParB)
            fprintf('max non-constant mode coefficient of k_||: %g  (It should be nearly 0).\n', max(abs(flowDividiedByB(2:end))))
            
            if runMode==0 && plotPoloidalVariationOfFlow
                figure(6+figureOffset)
                clf
                plot(flowDividiedByB)
                xlabel('\theta mode number')
                ylabel('parallel flow coefficient k')
            end

        else
            % Electrons
            
            kIntegralWeight = (x_i.^3);
            
            % Compute L31:
            L31BeforeThetaIntegral = zeros(thetaGridSize,1);
            for itheta=1:thetaGridSize
                L=1;
                indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                fSlice = soln(indices,2);
                L31BeforeThetaIntegral(itheta) = xWeights*(kIntegralWeight' .* fSlice);
            end
            
            L31BeforeThetaIntegral = -(4/sqrtpi)*(2/3)*2*pi/VPrime*thetaMatrixForBootstrapCoefficients*L31BeforeThetaIntegral;
            L31 = L31BeforeThetaIntegral(1);
            
            % Compute L32:
            L32BeforeThetaIntegral = zeros(thetaGridSize,1);
            for itheta=1:thetaGridSize
                L=1;
                indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                fSlice = soln(indices,1);
                L32BeforeThetaIntegral(itheta) = xWeights*(kIntegralWeight' .* fSlice);
            end
            
            L32BeforeThetaIntegral = -(4/sqrtpi)*(2/3)*2*pi/VPrime*thetaMatrixForBootstrapCoefficients*L32BeforeThetaIntegral;
            L32 = L32BeforeThetaIntegral(1);
            
            % Compute L34:
            L34BeforeThetaIntegral = zeros(thetaGridSize,1);
            for itheta=1:thetaGridSize
                L=1;
                indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                fSlice = soln(indices,3);
                L34BeforeThetaIntegral(itheta) = xWeights*(kIntegralWeight' .* fSlice);
            end
            
            L34BeforeThetaIntegral = (4/sqrtpi)*(2/3)*2*pi/(FSAB2*VPrime)*thetaMatrixForBootstrapCoefficients*L34BeforeThetaIntegral;
            L34 = L34BeforeThetaIntegral(1);
            
            % Compute conductivity:
            ConductivityBeforeThetaIntegral = zeros(thetaGridSize,1);
            for itheta=1:thetaGridSize
                L=1;
                indices = ((1:Nx)-1)*Nxi*thetaGridSize + L*thetaGridSize + itheta;
                fSlice = soln(indices,4);
                ConductivityBeforeThetaIntegral(itheta) = xWeights*(kIntegralWeight' .* fSlice);
            end
            
            L11=1.9693;
            ConductivityBeforeThetaIntegral = (nuPrime/L11)*(4/sqrtpi)*(2/3)*2*pi/(FSAB2 *VPrime) ...
                *thetaMatrixForBootstrapCoefficients*ConductivityBeforeThetaIntegral;
            conductivity = ConductivityBeforeThetaIntegral(1);
            
            fprintf('L31 bootstrap current coefficient: %g\n',L31)
            fprintf('L32 bootstrap current coefficient: %g\n',L32)
            fprintf('L34 bootstrap current coefficient: %g\n',L34)
            fprintf('Neoclassical conductivity, normalized to Spitzer: %g\n',conductivity)
        end
        
        
        fprintf('******************** Done ********************\n')
        
        
        
        
        % ************************************
        % Below are utilities for building sparse matrices:
        % ************************************
        
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
                fprintf('Warning! estimated_nnz is too small.\n')
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
        
        function matrix = generateThetaMultiplicationMatrix(functionOnThetaGrid, threshhold)
            % This function generates the matrix that corresponds to
            % multiplication by a function of theta.
            
            complexAModeCoefficients = fft(functionOnThetaGrid)/(2*Ntheta);
            
            matrix = zeros(thetaGridSize);
            
            % Make top-left block of matrix:
            toeplitzGenerator = [0, fliplr(real(complexAModeCoefficients(1:Ntheta))), real(complexAModeCoefficients(2:Ntheta))];
            topLeftBlockBeforeMirroring = toeplitz(zeros(Ntheta,1),toeplitzGenerator);
            matrix(1:Ntheta,1:Ntheta) = topLeftBlockBeforeMirroring(:,(Ntheta+1):end);
            matrix(1:Ntheta, 2:Ntheta) = matrix(1:Ntheta, 2:Ntheta) + fliplr(topLeftBlockBeforeMirroring(:,2:Ntheta));
            
            % Make bottom-right block of matrix:
            matrix((Ntheta+1):thetaGridSize,(Ntheta+1):thetaGridSize) = topLeftBlockBeforeMirroring(1:(Ntheta-1),(Ntheta+1):thetaGridSize) ...
                - fliplr(topLeftBlockBeforeMirroring(1:(Ntheta-1), 1:(Ntheta-1)));
            
            % Make bottom-left block of matrix:
            toeplitzGenerator = [0, fliplr(imag(complexAModeCoefficients(2:Ntheta))), 0, -imag(complexAModeCoefficients(2:Ntheta))];
            bottomLeftBlockBeforeMirroring = toeplitz(zeros(Ntheta,1),toeplitzGenerator);
            matrix((Ntheta+1):thetaGridSize, 1:Ntheta) = bottomLeftBlockBeforeMirroring(1:(Ntheta-1),Ntheta:thetaGridSize);
            matrix((Ntheta+1):thetaGridSize, 2:(Ntheta-1)) = matrix((Ntheta+1):thetaGridSize, 2:(Ntheta-1)) + fliplr(bottomLeftBlockBeforeMirroring(1:(Ntheta-1),2:(Ntheta-1)));
            
            % Make top-right block of matrix:
            matrix(1:Ntheta, (Ntheta+1):thetaGridSize) = -bottomLeftBlockBeforeMirroring(:,(Ntheta+2):end);
            matrix(1:Ntheta, (Ntheta+1):thetaGridSize) = matrix(1:Ntheta, (Ntheta+1):thetaGridSize) + fliplr(bottomLeftBlockBeforeMirroring(:,2:Ntheta));
            
            %threshhold = 1e-7;
            matrix(abs(matrix)<threshhold)=0;
        end
        
    end

% *********************************************************
% *********************************************************
%
% Below are several functions related to magnetic geometry.
%
% *********************************************************
% *********************************************************

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

% *********************************************************
% *********************************************************
%
% Below are several functions related to the Sauter formulae:
%
% *********************************************************
% *********************************************************

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
        semilogx(nuStarsFine, -SauterKs,'k')
        hold on
        semilogx(nuStars,-flowCoefficients(1,:),'.-')
        semilogx(nuStars,-flowCoefficients(2,:),'.--g')
        semilogx(nuStars,-flowCoefficients(3,:),'.-.r')
        semilogx(nuStars,-flowCoefficients(4,:),'.:m')
        xlabel('\nu_*')
        title('Ion coefficient \alpha')
        xlim([nuStarMin,nuStarMax])
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx','Location','northwest')
        
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
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx','Location','northwest')
        
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
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx','Location','northeast')
        
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
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx','Location','southeast')
        
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
        legend('Sauter','me: base case','me: 2x N\theta','me: 2x N\xi','me: 2x Nx','Location','northeast')
        
        
    end
end