%%% =======================================================================
%%% = 09/14/2025
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% = Load CMAES output files, make relevant plots
%%% =======================================================================

baseDir = pwd;
utilDir = sprintf('%s/funcs/', baseDir);
dataDir = sprintf('%s/data/',  baseDir);
outDir  = sprintf('%s/output/',baseDir);
cmaesDir = sprintf('%s/hermesout/', baseDir);
addpath(utilDir);
addpath(sprintf('%s/obs',               utilDir));
addpath(sprintf('%s/ems',               utilDir));
addpath(sprintf('%s/model',             utilDir));
addpath(sprintf('%s/util',              utilDir));
addpath(sprintf('%s/plot',              utilDir));
addpath(sprintf('%s/inv',               utilDir));
addpath(sprintf('%s/inv/deterministic', utilDir));
addpath(sprintf('%s/inv/stochastic',    utilDir));

% Load CMAES output data
fit = load(sprintf('%s/outcmaesfit.dat', cmaesDir));
    % columns: iteration #, evaluation #, sigma, axis ratio, bestever,
    % best-median-worst fitness function value
xmean = load(sprintf('%s/outcmaesxmean.dat', cmaesDir));
    % columns="iteration #, evaluation #, void, ' ...
    %                           'void, void, xmean"'
axlen = load(sprintf('%s/outcmaesaxlen.dat', cmaesDir));
    % columns="iteration, evaluation, sigma, ' ...
    % 		 'max axis length, min axis length, ' ...
    % 		 'all principal axes lengths (sorted square roots ' ...
    %                   'of eigenvalues of C)"' 
xrecentbest = load(sprintf('%s/outcmaesxrecentbest.dat', cmaesDir));
    % columns="iteration, evaluation, fitness, ' ...
    %                           'void, void, xrecentbest"'

% find a time in our algorithm where things looked good (might change this
% metric)
[~, idx] = min(fit(:,5));
% take the recent best value at that time
xbest = xrecentbest(idx+2,6:end);

% calculate emsParams as in DriverScript to use as inParams so we can disassemble state vector

if true
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
    obs = makeObs(St,tAvg,ch4_NH_obs,dD_NH_obs,ch4_obs,d13c_obs,dD_obs,d14c_obs,reread);
    %%% Also store these in the params structure
    params.obs = obs;
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
end


% plot the best ever
emsParams_best = disassembleStateVector(xbest,emsParams);
out_best       = boxModel_wrapper(params,St,emsParams_best);
plotResults_alt( params, out_best, obs, 'cmaesfull_', true, true)



% 
% fun_param.run_parallel = false;
% 
% plikes = zeros(767);
% ppriors = zeros(767);
% % likeli_params.obs = params.obs;
% % likeli_params.St = params.St;
% % likeli_params.params = params;
% % likeli_params.emsParams = emsParams;
% % likeli_params.use_log = fun_param.use_log;
% for i = 1:767
%     emsi = disassembleStateVector(xmean(i,:),emsParams);
%     ppriors(i) = define_prior(emsi, likeli_params);
%     plikes(i) = define_likelihood(emsi, likeli_params); % pass in function params
%     %plikes(i) = cmaes_fun_eval((xmean(i,:))',fun_param);
% end
% 
% 


 % [xmin,fmin,counteval,stopflag,out,bestever] = cmaes('cmaes_fun_eval', xstart, insigma, CMAES_opts, fun_param);
 % 
 %    end
 % 
 %    %%% Plot the best one
 %    [cmaes_res.emsParams] = disassembleStateVector(bestever.x,emsParams); % this emsParams might have different stepfun forcing
 %    emsParams_best = cmaes_res.emsParams;

%cmaes_fun_eval( stateVector, fun_param )