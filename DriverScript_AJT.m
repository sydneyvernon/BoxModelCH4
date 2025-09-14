%%% =======================================================================
%%% = DriverScript.m
%%% = Alex Turner
%%% = Originally created on 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES:
%%% =
%%% = This is the driver script for the 4-box model methane inversion.
%%% = There are currently two different inversions implemented: (1) a
%%% = linear or non-linear deterministic inversion following Rodgers (2000)
%%% = and (2) an inversion using the non-linear Covariance Matrix Adaptation
%%% = Evolution Strategy (CMA-ES).  Case (1) requires us to compute gradients
%%% = and only allows for Gaussian errors.  Case (2) is  a stochastic method
%%% = that automatically tunes the proposal distribution to improve the sampling,
%%% = however it does not provide error statistics that are consistent with
%%% = the distributions.  Case (2) also allows one to specify non-analytic
%%% = distributions (e.g., bounded Gaussians or uniform distributions).
%%% =======================================================================


%%
%%% =======================================================================
%%% 1. Initialize
%%% =======================================================================

%%% Clear the MatLab space
clf
clear all
close all
clc

%%% Header
fprintf('\n ***********************************\n')
fprintf(' *** STARTING GLOBAL 4-BOX MODEL ***\n')
fprintf(' ***********************************\n')

%%% Define the directories
baseDir = pwd;
utilDir = sprintf('%s/funcs/', baseDir);
dataDir = sprintf('%s/data/',  baseDir);
outDir  = sprintf('%s/output/',baseDir);

%%% Add the utility functions
addpath(utilDir);
addpath(sprintf('%s/obs',               utilDir));
addpath(sprintf('%s/ems',               utilDir));
addpath(sprintf('%s/model',             utilDir));
addpath(sprintf('%s/util',              utilDir));
addpath(sprintf('%s/plot',              utilDir));
addpath(sprintf('%s/inv',               utilDir));
addpath(sprintf('%s/inv/deterministic', utilDir));
addpath(sprintf('%s/inv/stochastic',    utilDir));

%%% Define the time period
% Are we using a spinup?
spinUp = 1000;
% Time period
sYear = -1000000;
%sYear = -2500000;
%sYear = -500000;
eYear = 0;
sYear = sYear - spinUp;
tRes  = 'kyr';     % Can be 'year', 'month', or 'kyr' (year preferred)
tAvg  = 'kyr';     % Smooth the observations
St    = getTime(sYear,eYear,tRes); % Time vector
pSt   = datevec(St); pSt = double(pSt(:,1))/1000; % Time vector for plotting;
nT    = length(St);

%%% Execute in parallel?
run_parallel = false;
if run_parallel
    nWorkers     = 4;
    setupParallel(run_parallel,nWorkers);
end

%%% What kind of inversions do we want to do?
save_data        = false;
do_deterministic = false;    % Rodgers (2000)
det_linear       = false;    % Use a linear deterministic inversion?
do_cmaes         = true;     % Covariance Matrix Adaptation Evolution Strategy

%%% For reading the observations
% Do we want to reread the raw data?
reread.flag  = true; %false;
% Other flags for re-reading
reread.sYear = sYear;
reread.eYear = eYear;
reread.tRes  = tRes;
reread.tAvg  = tAvg;
reread.dir   = dataDir;

%%% Which case are we doing?
% % Case A: Constant Emissions
% caseName = 'caseAa';   % variable OH
% caseName = 'caseAb';   % fixed OH
% % Case B: Variable Emissions
% caseName = 'caseBa';   % variable OH
% caseName = 'caseBb';   % fixed OH
% % Case C: Variable CH4, CO, and OH Emissions
% caseName = 'caseBa';   % variable OH
% caseName = 'caseBb';   % fixed OH

%%% My choice
caseName = 'caseBa';
add_name = 'newplots'; %[caseName 'newplots'];
add_name_inv = ['cmaes_best_NHdD' ''];

%%% Other options and flags
% General flags
use_strat       = true;     % Use a stratosphere?
use_other_sinks = true;     % Use non-OH sinks?
diagnostics     = true;
% MCF sensitivity test flags
k_co_flag       = true;     % Use k_CO that AJT derived

%%% Set the seed for repeatability
rng('default');

%%% Sort out the case info
cStyle = char(caseName);
if strcmp(cStyle(end),'b')
    interactive_OH = false;
else
    interactive_OH = true;
end

%%% Build the structure with parameters we'll use for the box model
params                 = getParameters(St, true);
params.caseName        = caseName;
params.tRes            = tRes;
params.interactive_OH  = interactive_OH;
params.use_strat       = use_strat;
params.pSt             = pSt;
params.St              = St;
params.spinup          = spinUp;
params.diagnostics     = false;


%%
%%% =======================================================================
%%% 2. Load the obs
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** LOADING THE OBSERVATIONS *** \n');

%%% Load the observations
% Structures with with three fields:
% - "obs":  Observations from each NOAA site (ppb)
% - "tim":  Julian date for the observation
% - "lat":  Latitude of the NOAA site
try % Add a try-catch statement in case the user hasn't downloaded the data
    ch4_NH_obs  = getCH4_NH(dataDir,reread);       % CH4 observations, Greenland (ppb)
    dD_NH_obs  = getdD_NH(dataDir,reread);       % deltaD observations, Greenland (permil)
    ch4_obs  = getCH4(dataDir,reread);          % CH4 observations, Antarctica (ppb)
    d13c_obs = getd13C(dataDir,reread);         % delta13C observations (permil)
    dD_obs   = getdD(dataDir,reread);           % deltaD observations (permil)
    d14c_obs = getd14C(dataDir,reread);         % delta14C observations (permil)
catch % Some data is missing, set the observation structures to NaN
    try % Try just leaving out d14c
        ch4_NH_obs  = getCH4_NH(dataDir,reread);
        dD_NH_obs  = getdD_NH(dataDir,reread);
        ch4_obs  = getCH4(dataDir,reread);      % CH4 observations (ppb)
        d13c_obs = getd13C(dataDir,reread);     % delta13C observations (permil)
        dD_obs   = getdD(dataDir,reread);       % deltaD observations (permil)
        d14c_obs = NaN;
    catch
        try % Try leaving out dD and d14C
            ch4_NH_obs  = getCH4_NH(dataDir,reread); 
            dD_NH_obs  = getdD_NH(dataDir,reread);
            ch4_obs  = getCH4(dataDir,reread);  % CH4 observations (ppb)
            d13c_obs = getd13C(dataDir,reread); % delta13C observations (permil)
            dD_obs   = NaN;
            d14c_obs = NaN;
        catch
            fprintf(' * UNABLE TO READ OBSERVATIONS!\n');
            ch4_NH_obs  = NaN; 
            dD_NH_obs = NaN;
            ch4_obs  = NaN;
            d13c_obs = NaN;
            dD_obs   = NaN;
            d14c_obs = NaN;
        end
    end
end

%%% Make the observation structure
% Structure with 6 fields:
% - NH/SH CH4    obs & err (ppb)
% - NH/SH CH4C13 obs & err (permil)
% - NH/SH CO     obs & err (ppb)
obs = makeObs(St,tAvg,ch4_NH_obs,dD_NH_obs,ch4_obs,d13c_obs,dD_obs,d14c_obs,reread);

%%% Also store these in the params structure
params.obs = obs;


%%
%%% =======================================================================
%%% 3. Set parameters for the 4-box model
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** SET PARAMETERS FOR THE 4-BOX MODEL *** \n');

%%% Specify input parameters for the box model
% Initial conditions
IC = params.IC;

%%% Parameters for the emissions
emsParams.oh_scale = 0.08;       % OH reaction rate scaling factor
% Base emissions
emsParams.base_wet    = 300;     % Baseline tropical wetlands emissions (Tg CH4/yr)
emsParams.base_wetB   = 0;     % Baseline tropical wetlands emissions (Tg CH4/yr)
emsParams.base_fire   = 15;     % Baseline fire emissions (Tg CH4/yr)
emsParams.base_fossil = 6;      % Baseline fossil emissions (Tg CH4/yr)
emsParams.base_ocean  = 54;     % Baseline ocean emissions for CO (Tg CO/yr)
emsParams.base_oh     = 600;    % Baseline OH production (Tg OH/yr) (Murray et al., 2014 estimate is 540-2000 Tg/yr)
emsParams.base_animal = 10;     % Baseline animal emissions (Tg CH4/yr)
emsParams.base_cl     = 1650;   % Chlorine abundance to get a lifetime of 1/302 yr (molec/cm3)
% Other parameters
emsParams.soilOffset = 24;      % Temperature offset for soils
emsParams.fireFacCO  = 10;      % 10x Tg/yr for CO compared to methane
emsParams.base_tauTS = 9.0;     % strat-trop exchange (years)
emsParams.stepChange = datenum(-45000,1,1); % Onset of megafauna extinction (years)
stepFun              = params.step(St,emsParams.stepChange);
% Scale factors for the different sectors
%emsParams.amp_wet          = 70;      % Wetland scaling factor
emsParams.Q10_tropical     = 1.3582;  % From regression (previous work used ~2)
emsParams.Q10_boreal       = 1.3582;  % From regression (previous work used ~2)
emsParams.amp_wet_tropical = 40;      % Wetland scaling factor
emsParams.amp_wet_boreal   = 40;      % Wetland scaling factor
emsParams.amp_fire         = 10;      % Fire temperature scaling factor
emsParams.amp_fossil       = 8;       % Fossil temperature scaling factor
emsParams.amp_oh           = 0;       % Vegetation scaling factor
emsParams.amp_animals      = 0;       % Animal scaling factor

%%% 
% Tropical Wetlands
base_WT   = emsParams.base_wet;
eccentWT  = 0;
obliqWT   = 0;
precessWT = 0;
insolWT   = 0;
iceWT     = 0;
tempWT    = 0;
stepWT    = 0;
% Boreal Wetlands
base_WTB   = emsParams.base_wetB;
eccentWTB  = 0;
obliqWTB   = 0;
precessWTB = 0;
insolWTB   = 0;
iceWTB     = 0;
tempWTB    = 0;
stepWTB    = 0;
% Fossil
base_FF   = emsParams.base_fossil;
eccentFF  = 0;
obliqFF   = 0;
precessFF = 0;
insolFF   = 0;
iceFF     = 0;
tempFF    = 0;
stepFF    = 0;
% Fires
base_BB   = emsParams.base_fire;
eccentBB  = 0;
obliqBB   = 0;
precessBB = 0;
insolBB   = 0;
iceBB     = 0;
tempBB    = 0;
stepBB    = 0;
% OH
base_OH   = emsParams.base_oh;
eccentOH  = 0;
obliqOH   = 0;
precessOH = 0;
insolOH   = 0;
iceOH     = 0;
tempOH    = 0;
stepOH    = 0;
% Strat-Trop exchange
base_tauTS   = emsParams.base_tauTS;
eccenttauTS  = 0;
obliqtauTS   = 0;
precesstauTS = 0;
insoltauTS   = 0;
icetauTS     = 0;
temptauTS    = 0;
steptauTS    = 0;
% Animals
base_AN   = emsParams.base_animal;
eccentAN  = 0;
obliqAN   = 0;
precessAN = 0;
insolAN   = 0;
iceAN     = 0;
tempAN    = 0;
stepAN    = 0;
% Chlorine
base_CL   = emsParams.base_cl;
eccentCL  = 0;
obliqCL   = 0;
precessCL = 0;
insolCL   = 0;
iceCL     = 0;
tempCL    = 0;
stepCL    = 0;

%%% Put them all into structures
emsParams.base_WT      = base_WT;
emsParams.base_WTB     = base_WTB;
emsParams.base_FF      = base_FF;
emsParams.base_BB      = base_BB;
emsParams.base_OH      = base_OH;
emsParams.base_tauTS   = base_tauTS;
emsParams.base_AN      = base_AN;
emsParams.base_CL      = base_CL;

% USE 4 FORCINGS (no Milankovitch)
emsParams.forcings     = [params.insol_NH, params.iceVolN, params.sTempN,   stepFun];
emsParams.WT_scalings  = [insolWT,          iceWT,        tempWT,    stepWT]';
emsParams.WTB_scalings = [insolWTB,         iceWTB,       tempWTB,   stepWTB]';
emsParams.FF_scalings  = [insolFF,          iceFF,        tempFF,    stepFF]';
emsParams.BB_scalings  = [insolBB,          iceBB,        tempBB,    stepBB]';
emsParams.OH_scalings  = [insolOH,          iceOH,        tempOH,    stepOH]';
emsParams.tau_scalings = [insoltauTS,       icetauTS,     temptauTS, steptauTS]';
emsParams.AN_scalings  = [insolAN,          iceAN,        tempAN,    stepAN]';
emsParams.CL_scalings  = [insolCL,          iceCL,        tempCL,    stepCL]';
emsParams.IC           = params.IC;


% % %%% Fitted parameters (fitted without temperature dependence for RxNs)
% % Wetlands
% eccentWT  =  1.1829;
% obliqWT   = -1.2878;
% precessWT = -0.1796;
% insolWT   =  0.0333;
% iceWT     =  0.5131;
% tempWT    =  1.5753;
% % Fossil
% eccentFF  = -0.1070;
% obliqFF   =  0.0900;
% precessFF =  0.0859;
% insolFF   =  0.0311;
% iceFF     =  8.3244;
% tempFF    =  0.9923;
% % Fires
% eccentBB  =  0.2761;
% obliqBB   =  0.0001;
% precessBB =  0.0808;
% insolBB   =  0.0804;
% iceBB     =  3.0602;
% tempBB    = -0.5709;
% % OH
% eccentOH  = -0.3818;
% obliqOH   =  0.9347;
% precessOH = -0.2962;
% insolOH   =  1.1644;
% iceOH     =  1.4802;
% tempOH    =  0.2198;
% % Strat-Trop exchange
% eccenttauTS  = -0.3392;
% obliqtauTS   = -1.0192;
% precesstauTS =  0.8555;
% insoltauTS   =  0.0000;
% icetauTS     =  1.0467;
% temptauTS    =  0.8010;
% steptauTS    =  0;
% % Animals
% eccentAN  = 0;
% obliqAN   = 0;
% precessAN = 0;
% insolAN   = 0;
% iceAN     = 0;
% tempAN    = 0;
% stepAN    = 0;
% % Chlorine
% eccentCL  = 0;
% obliqCL   = 0;
% precessCL = 0;
% insolCL   = 0;
% iceCL     = 0;
% tempCL    = 0;
% stepCL    = 0;
% emsParams.oh_scale = 0.2568;  % alex original
% %emsParams.oh_scale = 0.5; % sv attempt
% emsParams.Q10_tropical     = .2663;
% emsParams.Q10_boreal     = 1.0002;
% emsParams.amp_wet_tropical = 46.8404;
% emsParams.amp_wet_boreal = 45.1881;
% emsParams.stepChange = emsParams.stepChange;
% emsParams.base_WT = 11.1899;
% emsParams.base_WTB = 0;
% emsParams.base_FF = -0.1047;
% emsParams.base_BB = -20.8357;
% emsParams.base_OH = 908.1101;
% emsParams.base_tauTS = 7.6618;
% emsParams.base_AN = 0;
% emsParams.base_CL = 1650;
% emsParams.forcings     = [params.eccent, params.obliq, params.precess, params.insol_NH, params.iceVolN, params.sTempN,   stepFun];
% emsParams.WT_scalings  = [     eccentWT,      obliqWT,      precessWT,         insolWT,          iceWT,        tempWT,    stepWT]';
% emsParams.WTB_scalings = [    eccentWTB,     obliqWTB,     precessWTB,        insolWTB,         iceWTB,       tempWTB,   stepWTB]';
% emsParams.FF_scalings  = [     eccentFF,      obliqFF,      precessFF,         insolFF,          iceFF,        tempFF,    stepFF]';
% emsParams.BB_scalings  = [     eccentBB,      obliqBB,      precessBB,         insolBB,          iceBB,        tempBB,    stepBB]';
% emsParams.OH_scalings  = [     eccentOH,      obliqOH,      precessOH,         insolOH,          iceOH,        tempOH,    stepOH]';
% emsParams.tau_scalings = [  eccenttauTS,   obliqtauTS,   precesstauTS,      insoltauTS,       icetauTS,     temptauTS, steptauTS]';
% emsParams.AN_scalings  = [     eccentAN,      obliqAN,      precessAN,         insolAN,          iceAN,        tempAN,    stepAN]';
% emsParams.CL_scalings  = [     eccentCL,      obliqCL,      precessCL,         insolCL,          iceCL,        tempCL,    stepCL]';
% emsParams.IC           = params.IC;


%%
%%% =======================================================================
%%% 4. Run the 2-box model
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** RUN THE 2-BOX MODEL WITH PRIOR FLUXES *** \n');

%%% Run the box model
% 
% out = boxModel_wrapper(params,St,emsParams);
% plotResults_alt( params, out, obs, add_name, true, false)  % name suffix, plot OH?, plot CO?

%%% Write the prior data
if save_data
    %writeData(pSt, out, outDir, caseName);
end


%%
%%% =======================================================================
%%% 5. Deterministic inversion (Rodgers, 2000)
%%% =======================================================================

if do_deterministic
        
    %%% Diagnostic
    fprintf('\n *** DETERMINISTIC INVERSION *** \n');

    %%% Invert
    [anal_soln]    = invert_methane(St,params,emsParams,IC,obs,det_linear,run_parallel);
    emsParams_anal = anal_soln{1};

    %%% Plot the Jacobians
    %plotJacobian(St,jacobian_ems,tRes,sprintf('%s/%s/jacobian_%%s.%s',outDir,tRes,ftype));

    %%% Try plotting the solution
    out_anal = boxModel_wrapper(params,St,emsParams_anal);
    plotResults( params, out_anal, obs )
    writeData(pSt, out_anal, outDir, caseName);
    
end


%%
%%% =======================================================================
%%% 6. Invert with CMAES
%%% =======================================================================

if do_cmaes
    
    %%% Diagnostics
    fprintf('\n *** STARTING CMA-ES INVERSION *** \n');
    
    %%% Set the function parameters for the box model
    fun_param.run_parallel = run_parallel;
    fun_param.p_prior      = @(emsParams,input_param) define_prior(emsParams,input_param);
    fun_param.p_like       = @(emsParams,input_param) define_likelihood(emsParams,input_param);
    fun_param.use_log      = true;
    fun_param.params       = params;
    fun_param.St           = St;
    fun_param.emsParams    = emsParams;
    fun_param.nT           = nT;
    fun_param.obs          = obs;

    %%% Set the options for CMAES
    CMAES_opts                   = cmaes;
    CMAES_opts.EvalParallel      = 'yes';
    CMAES_opts.SaveFilename      = sprintf('output/%s/cmaes_dat/variablescmaes.mat',tRes);
    CMAES_opts.LogFilenamePrefix = sprintf('output/%s/cmaes_dat/outcmaes',tRes);
    CMAES_opts.Resume            = 'no'; % Default is 'no'
    CMAES_opts.CMA.active        = 2;
    CMAES_opts.DiagonalOnly      = 100;
    % Full
    % CMAES_opts.StopFunEvals      = 5000000;
    % CMAES_opts.MaxIter           = 20000;
    % CMAES_opts.Restarts          = 10;
%     % Medium
%     CMAES_opts.StopFunEvals      = 1000000;
%     CMAES_opts.MaxIter           = 5000;
%     CMAES_opts.Restarts          = 6;
%     % Small
    % CMAES_opts.StopFunEvals      = 10000;
    % CMAES_opts.MaxIter           = 200;
    % CMAES_opts.Restarts          = 2;
    % Very small
    CMAES_opts.StopFunEvals      = 1000;
    CMAES_opts.MaxIter           = 200;
    CMAES_opts.Restarts          = 1;

    %%% Get the starting point and standard deviations
    % Starting point
    if do_deterministic % Use the deterministic inversion as a starting point
        xstart = assembleStateVector(emsParams_anal);
    else
        xstart = assembleStateVector(emsParams);
    end
    % Standard deviations for step size
    divFac                   = 1/20;
    insigma.oh_scale         = divFac * emsParams.oh_scale;
    insigma.Q10_tropical     = divFac * emsParams.Q10_tropical;
    insigma.Q10_boreal       = divFac * emsParams.Q10_boreal;
    insigma.amp_wet_tropical = divFac * emsParams.amp_wet_tropical;
    insigma.amp_wet_boreal   = divFac * emsParams.amp_wet_boreal;
    insigma.stepChange       = divFac * 1000;
    insigma.base_WT          = divFac * emsParams.base_WT;
    insigma.base_FF          = divFac * emsParams.base_FF;
    insigma.base_BB          = divFac * emsParams.base_BB;
    insigma.base_OH          = divFac * emsParams.base_OH;
    insigma.base_tauTS       = divFac * emsParams.base_tauTS;
    insigma.base_AN          = divFac * emsParams.base_AN;
    insigma.base_CL          = divFac * emsParams.base_CL;
    insigma.WT_scalings      = divFac * ones(size(emsParams.WT_scalings));
    insigma.FF_scalings      = divFac * ones(size(emsParams.FF_scalings));
    insigma.BB_scalings      = divFac * ones(size(emsParams.BB_scalings));
    insigma.OH_scalings      = divFac * ones(size(emsParams.OH_scalings));
    insigma.tau_scalings     = divFac * ones(size(emsParams.tau_scalings));
    insigma.AN_scalings      = divFac * ones(size(emsParams.AN_scalings));
    insigma.CL_scalings      = divFac * ones(size(emsParams.CL_scalings));
    insigma.IC               = divFac * emsParams.IC;
    % Reshape into state vector format
    insigma               = abs(assembleStateVector(insigma));
    insigma(insigma <= 0) = 1e-4; insigma_save = insigma;
    
    %%% Use an old estimate?
    oldName = sprintf('./%s',CMAES_opts.SaveFilename);
    use_old = false;
    if use_old
        if (exist(oldName,'file') == 2)
            load(oldName);
            xstart = bestever.x;  % this has now been overwritten with correct size xstart
%             if length(xstart) < length(assembleStateVector(emsParams))
%                 emsParams.oh_scale         = bestever.x(1);
%                 emsParams.Q10_tropical     = bestever.x(2);
%                 emsParams.amp_wet_tropical = bestever.x(3);
%                 emsParams.Q10_boreal       = bestever.x(4);
%                 emsParams.amp_wet_boreal   = bestever.x(5);
%                 emsParams.stepChange       = 0;
%                 emsParams.base_WT          = bestever.x(6);
%                 emsParams.base_FF          = bestever.x(7);
%                 emsParams.base_BB          = bestever.x(8);
%                 emsParams.base_OH          = bestever.x(9);
%                 emsParams.base_tauTS       = bestever.x(10);
%                 emsParams.base_AN          = 0;
%                 emsParams.base_CL          = base_CL;
%                 ii  = 11;
%                 nI  = length(emsParams.WT_scalings)-2;
%                 nIC = length(emsParams.IC)-1;
%                 emsParams.WT_scalings  = [bestever.x(ii:ii+nI); 0]; ii = ii+nI+1;
%                 emsParams.FF_scalings  = [bestever.x(ii:ii+nI); 0]; ii = ii+nI+1;
%                 emsParams.BB_scalings  = [bestever.x(ii:ii+nI); 0]; ii = ii+nI+1;
%                 emsParams.OH_scalings  = [bestever.x(ii:ii+nI); 0]; ii = ii+nI+1;
%                 emsParams.tau_scalings = [bestever.x(ii:ii+nI); 0]; ii = ii+nI+1;
%                 emsParams.AN_scalings(:) = 0;
%                 emsParams.CL_scalings(:) = 0;
%                 emsParams.IC           = bestever.x(ii:ii+nIC);
%                 xstart  = assembleStateVector(emsParams);
%                 insigma = insigma_save;
%             end
        end
    end
    
    % Run the CMA-ES inversion or just plot an old one?
    run_cmaes = true;
    if run_cmaes

        %%% Invert with CMAES
        % INPUTS
        %  - fitfun:  name of objective/fitness function
        %  - xstart:  objective variables initial point, determines N
        %  - insigma: initial coordinate wise standard deviation(s)
        %  - inopts:  options struct
        % OUTPUTS
        %  - xmin:      minimum search point of last iteration
        %  - fmin:      function value of xmin
        %  - counteval: number of function evaluations done
        %  - stopflag:  stop criterion reached
        %  - out:       struct with various histories and solutions
        %  - bestever:  struct containing overall best solution (for convenience)
        
        %%% Use a try catch in case we crash early
        [xmin,fmin,counteval,stopflag,out,bestever] = cmaes('cmaes_fun_eval', xstart, insigma, CMAES_opts, fun_param);

    end

    %%% Plot the best one
    [cmaes_res.emsParams] = disassembleStateVector(bestever.x,emsParams); % this emsParams might have different stepfun forcing
    emsParams_best = cmaes_res.emsParams;
    %emsParams_best.forcings(:,end) = params.step(St,emsParams_best.stepChange); 
        % this compensates for moving this line out of boxmodel_wrapper
    out_best       = boxModel_wrapper(params,St,emsParams_best);
    plotResults_alt( params, out_best, obs, add_name_inv, true, true)
    %load(oldName);close all;plotResults(params,boxModel_wrapper(params,St,disassembleStateVector(bestever.x,emsParams)),obs)
    %writeData(pSt, out_best, sprintf('%s/%s/',outDir,tRes), caseName);
    
    %%% 
    if true
        eval_param = fun_param;
        testParams = emsParams_best;
        sim_ch4  = nan(length(out_best.ch4),size(emsParams.forcings,2)+1);
        sim_d13C = nan(length(out_best.d13c),size(emsParams.forcings,2)+1);
        sim_dD   = nan(length(out_best.dD),size(emsParams.forcings,2)+1);
        sim_d14C = nan(length(out_best.d14c),size(emsParams.forcings,2)+1);
        cov_ch4  = nan(size(emsParams.forcings,2)+1,1);
        cov_d13C = nan(size(emsParams.forcings,2)+1,1);
        cov_dD   = nan(size(emsParams.forcings,2)+1,1);
        cov_d14C = nan(size(emsParams.forcings,2)+1,1);
        for i = 1:length(cov_ch4)
            testParams.forcings = emsParams.forcings;
            if i > 1 % i == 1: All forcings
                testParams.forcings(:,i-1) = mean(emsParams.forcings(:,i-1));
            end
            out_test      = boxModel_wrapper(params,St,testParams);
            sim_ch4(:,i)  = out_test.ch4;
            sim_d13C(:,i) = out_test.d13c;
            sim_dD(:,i)   = out_test.dD;
            sim_d14C(:,i) = out_test.d14c;
            % CH4
            sim        = sim_ch4(:,i);
            meas       = obs.ch4;
            ind        = ~isnan(sim) & ~isnan(meas);
            cov_ch4(i) = corr(sim(ind),meas(ind));   % cov_ch4, etc are 8 x 1
            % d13C
            sim         = sim_d13C(:,i);
            meas        = obs.d13c;
            ind         = ~isnan(sim) & ~isnan(meas);
            cov_d13C(i) = corr(sim(ind),meas(ind));
            % dD
            sim       = sim_dD(:,i);
            meas      = obs.dD;
            ind       = ~isnan(sim) & ~isnan(meas);
            cov_dD(i) = corr(sim(ind),meas(ind));
            % d14C
            sim         = sim_d14C(:,i);
            meas        = obs.d14c;
            ind         = ~isnan(sim) & ~isnan(meas);
            cov_d14C(i) = corr(sim(ind),meas(ind));
        end
        % Get reduction in R2
        X_catN  = {'Insol NH','iceVol','Temp','StepFun'};
        r2_ch4  = cov_ch4.^2;
        r2_d13C = cov_d13C.^2;
        r2_dD   = cov_dD.^2;
        r2_d14C = cov_d14C.^2;
        r2_res  = [r2_ch4(1)  - r2_ch4(2:end), ...
                   r2_d13C(1) - r2_d13C(2:end),...
                   r2_dD(1)   - r2_dD(2:end),  ...
                   r2_d14C(1) - r2_d14C(2:end)];
        % Re-order the variables
        ind_order = [2,3,1,4];  % order is now ice, temp, insol, step.
        X_catN    = X_catN(ind_order);
        r2_res    = r2_res(ind_order,:);
        % Plot
        fN = figure();pos = get(fN,'position');set(gcf,'color','w');
        %set(fN,'position',[0 pos(2) pos(3)*1.7 pos(4)*2.7])
        X_cats = categorical(X_catN);
        X_cats = reordercats(X_cats,X_catN);
        b = bar(X_cats,r2_res);
        for i = 1:length(b)
            b(i).EdgeColor = 'none';
        end
        labels = {'ch4', 'd13C', 'dD', 'd14C'};
        legend(labels, 'Location', 'northeast');
        set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'YMinorTick','on')
        ylabel('Change in R^2','FontSize',20,'FontWeight','bold')
        export_fig(fN,['./output/' add_name_inv 'var_explained_diff.png'],'-png','-m2','-painters','-cmyk');

    end
    
end


%%
%%% Finished simulation
fprintf('\n ***********************************\n')
fprintf(' ***            DONE!            ***\n')
fprintf(' ***********************************\n\n')


%%% =======================================================================
%%% END
%%% =======================================================================