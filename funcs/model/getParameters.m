%%% =======================================================================
%%% = getParameters.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Computes the parameters that are needed for the box model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St -- Our time vector.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): params -- Structure with the parameters for the box model.
%%% =======================================================================

function [ params ] = getParameters(St, tempdependence)

%%% Turn off an unnecessary warning
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

%%% Flags for the simulation
% General flags
use_other_sinks = true;         % Use non-OH sinks?
k_co_flag       = true;         % Use k_CO that AJT derived

%%% Box model parameters and unit conversions
% Unit conversions
DaysToS = 60 * 60 * 24;         % Days to Seconds
YrToDay = 365.25;               % Years to Days
% Masses
m       = 5.15e21;              % Rough: Total mass of atmosphere in g
n_air   = 2.69e19;              % Rough: number density of air (molec/cm3)
m_air   = 28.8;                 % Average molar mass of the atmosphere in g/mol (mostly N2 and O2)
m_ch4   = 16;                   % Molar mass of CH4
m_oh    = 17.008;               % Molar mass of hydroxyl
m_co    = 28.01;                % Molar mass of carbon monoxide
mConv   = m/m_air/1e12*1e-9;	% Factor to convert mixing ratios to Tg
mm_ch4  = m_ch4 * mConv;        % Convert CH4 mixing ratios (ppb) to Tg
mm_oh   = m_oh  * mConv;        % Convert OH number densities (ppb) to Tg
mm_co   = m_co  * mConv;        % Convert CO mixing ratios (ppb) to Tg
% Isotope standards
R_VPDB  = 0.0112372;            % VPDB standard for 13C/12C (Coplen et al., 2011)
R_VSMOW = 0.00015575;           % VSMOW standard for D/H (DeWit et al., 1980)
R_14CH4 = 1;                    % No standard needed for radiocarbon
% Kinetic isotope effect
KIE_oh    = 1.0039;             % Kinetic isotope effect (k12/k13) from (Saueressig et al., 2001 found 1.0039, recommended by Burkholder et al., 2020)
KIE_strat = 1.013;              % Kinetic isotope effect from (Saueressig et al., 2001 found 1.013) (Bock et al. used 1.022)
KIE_soil  = 1.012;              % Kinetic isotope effect from Bock et al. (2017)
KIE_cl    = 1.066;              % Kinetic isotope effect from (Saueressig et al. 1996 found 1.066) (Bock et al. used 1.060)
% Kinetic isotope effect for for dD
KIE_dD_oh    = 1.294;           % Kinetic isotope effect from (Saueressig et al., 2001 found 1.294) (Bock et al. used 1.231)
KIE_dD_strat = 1.060;           % Kinetic isotope effect from (Saueressig et al., 2001 found 1.060) (Bock et al. used 1.080)
KIE_dD_soil  = 1.160;           % Kinetic isotope effect from Bock et al. (2017)
KIE_dD_cl    = 1.508;           % Kinetic isotope effect from (Saueressig et al., 1996 found 1.508) (Bock et al. used 1.470)
% Gloal mean OH
gmOH     = 0.97d6;              % Global mean OH (molec/cm^3)
gmOHplot = gmOH;                % Save value for plotting
% Reaction rate conversion factor
RxNconv = n_air / 1d9;          % Convert ppb to molec/cm3
% Reaction rates and lifetimes
k_12ch4        = 3.395e-15;     % reaction rate of OH with 12CH4 at 270K (cm3/molec/s)
k_13ch4        = k_12ch4 / KIE_oh; % reaction rate of OH with 13CH4 (cm3/molec/s)
k_14ch4        = k_13ch4 / KIE_oh; % reaction rate of OH with 14CH4 (cm3/molec/s)
k_12ch3d       = k_12ch4 / KIE_dD_oh; % reaction rate of OH with 12CH3D (cm3/molec/s)
t_ch4_strat_nh = 188;           % CH4 lifetime (yr) in the NH (Brown et al., 2013; Table 4)
t_ch4_strat_sh = 200;           % CH4 lifetime (yr) in the SH (Brown et al., 2013; Table 4)
k_co_A         = 1.0133e-12;    % reaction rate of OH with CO at 270K (cm3/molec/s), sum of bi- & termolecular reactions (Burkholder, Table 2.1)
k_co_B         = 2e-13;         % reaction rate (cm3/molec/s) from Prather (1993)
k_co_strat     = 12;            % (yr^-1) Assuming lifetime of 1 month in the stratosphere
k_oh_strat     = 12;            % (yr^-1) Assuming lifetime of 1 month in the stratosphere
k_cl           = 6.3654e-14;    % (cm3/molec/s) from Abbatt, Table 15-10 (k=2.36e-12*(T/298)^1.37*exp(-939/T) )
k_ch4_cl       = 1/302;         % (yr^-1) Tropospheric chlorine (~3% of total loss; Kirschke et al., 2013)
k_ch4_soil     = 1/302;         % (yr^-1) Soil uptake (~3% of total loss; Kirschke et al., 2013)
k_ch4_rad      = 1/8267;        % (yr^-1) Radioactive decay of 14C (Lassey et al., 2007)
k_co_other     = 0;             % (yr^-1) Other CO losses (currently neglecting)
% Temperature dependent reaction rates  %% HERE
k_cl  = @(T) 2.36e-12.*(T./298).^1.37.*exp(-939./T);    % cm3/molec/s
k_ch4 = @(T) 2.80e-14.*(T.^0.667).*exp(-1575./T);       % cm3/molec/s
% Lifetimes
tau_NS       = 1.0;             % Interhemispheric exchange rate (years)
tau_NS_strat = 3.3;             % Interhemispheric exchange rate in the stratosphere (Fabian et al., 1968)
% Are we using the other sinks?
if ~use_other_sinks
    k_ch4_cl   = 0;
    k_ch4_soil = 0;
    k_co_other = 0;
end

%%% Convert so everything has units of days (working with julian days)
tau_NS       = tau_NS * YrToDay;
tau_NS_strat = tau_NS_strat * YrToDay;
RxNconv      = RxNconv * DaysToS;
% Which methyl chloroform reaction rate are we using?
k_co = k_co_A;
if k_co_flag
    k_co = k_co_B;
end

%%% Convert some losses to consistent units 
k_ch4_cl   = k_ch4_cl   / YrToDay;
k_ch4_soil = k_ch4_soil / YrToDay;
k_ch4_rad  = k_ch4_rad  / YrToDay;
k_co_other = k_co_other / YrToDay;
k_co_strat = k_co_strat / YrToDay;
k_oh_strat = k_oh_strat / YrToDay;

%%% Get the atmospheric 14CO2
fspec             = '%f %f %f %f %f';
fileID            = fopen('./data/obs/iceCore/IntCal20.txt','r');
fdat              = textscan(fileID,fspec,'HeaderLines',11,'Delimiter',',');
atmos14C.t        = datenum(1950-fdat{1}(:),1,1);
atmos14C.d14c     = fdat{4}(:);
atmos14C.d14c_sig = fdat{5}(:);
% Put this onto our time domain
d14c     = interp1(atmos14C.t,atmos14C.d14c,St);
d14c_sig = interp1(atmos14C.t,atmos14C.d14c_sig,St);
d14c(isnan(d14c)) = nanmean(d14c);
d14c_sig(isnan(d14c_sig)) = 3*nanmax(d14c_sig);

%%% Guess for the initial conditions for the box model
% Troposphere
nh_12ch4  = 500;                                    % ppb
sh_12ch4  = 500;                                    % ppb
nh_13ch4  = -47.6;                                  % permil
sh_13ch4  = -47.4;                                  % permil
nh_14ch4  = 300;                                    % permil
sh_14ch4  = 300;                                    % permil
nh_12ch3d = -120;                                   % permil
sh_12ch3d = -120;                                   % permil
nh_oh     = gmOHplot * 1e-5;                        % 10^5 molec/cm3;
sh_oh     = gmOHplot * 1e-5;                        % 10^5 molec/cm3
nh_co     = 20;                                     % ppb
sh_co     = 20;                                     % ppb
% Stratosphere
nh_12ch4_S  = nh_12ch4;                             % ppb
sh_12ch4_S  = nh_12ch4;                             % ppb
nh_13ch4_S  = nh_13ch4;                             % permil
sh_13ch4_S  = sh_13ch4;                             % permil
nh_14ch4_S  = nh_14ch4;                             % permil
sh_14ch4_S  = sh_14ch4;                             % permil
nh_12ch3d_S = nh_12ch3d;                            % permil
sh_12ch3d_S = sh_12ch3d;                            % permil
nh_oh_S     = gmOHplot * 1e-5;                      % 10^5 molec/cm3
sh_oh_S     = gmOHplot * 1e-5;                      % 10^5 molec/cm3
nh_co_S     = nh_co;                                % ppb
sh_co_S     = sh_co;                                % ppb

% Assemble the ICs into a vector
IC = [  nh_12ch4,   sh_12ch4,   nh_13ch4,   sh_13ch4,   nh_14ch4,   sh_14ch4,   nh_12ch3d,   sh_12ch3d,    nh_oh,   sh_oh,   nh_co,   sh_co,...
      nh_12ch4_S, sh_12ch4_S, nh_13ch4_S, sh_13ch4_S, nh_14ch4_S, sh_14ch4_S, nh_12ch3d_S, sh_12ch3d_S,  nh_oh_S, sh_oh_S, nh_co_S, sh_co_S]';

%%% ODE45 parameters
Tspan = St;
% opts  = odeset('MaxStep',YrToDay/1,...     % Make sure the max timestep is 1 month
%               'NonNegative',1:length(IC),...% Ensure the result is positive
%               'OutputFcn',@odetpbar); % Progress Bar
% opts  = odeset('NonNegative',1:length(IC),...% Ensure the result is positive
%               'OutputFcn',@odetpbar); % Progress Bar
% clear textprogressbar
opts  = odeset('NonNegative',1:length(IC)); % Ensure the result is positive

%%% Load the climate forcings
% Bitanji et al. (2008) climate reconstruction
fspec               = '%f %f %f %f %f %f %f %f %f';
fileID              = fopen('./data/obs/iceCore/bintanja2008.txt','r');
fdat                = textscan(fileID,fspec,'HeaderLines',109);
reconstruct.t       = datenum(1950-fdat{1}(:)*1000,1,1);
reconstruct.sTemp   = fdat{5}(:); %nanmean(fdat{5}(:)).*ones(size(fdat{5}(:)));  %             % Surface air temperature (relative to present)
reconstruct.sTempN  = (reconstruct.sTemp - min(reconstruct.sTemp)) ./ max(reconstruct.sTemp - min(reconstruct.sTemp)); %ones(size(fdat{5}(:))); 
reconstruct.iceVol  = fdat{7}(:) + fdat{8}(:);   % Northern hemisphere ice volume (m sea level equivalent)
reconstruct.iceVolN = (reconstruct.iceVol - min(reconstruct.iceVol)) ./ max(reconstruct.iceVol - min(reconstruct.iceVol));
reconstruct.veg     = 1 - (reconstruct.iceVol - min(reconstruct.iceVol)) ./ max(reconstruct.iceVol - min(reconstruct.iceVol));
fclose(fileID);
% Read in the Milankovitch forcing
fspec             = '%f %f %f %f %f %f';
fileID            = fopen('./data/obs/iceCore/milankovitch.data.txt','r');
fdat              = textscan(fileID,fspec,'HeaderLines',8);
milank.t          = datenum(fdat{1}(:)*1000,1,1);
milank.eccent     = fdat{2}(:);
milank.obliq      = fdat{3}(:);
milank.peri       = fdat{4}(:);
milank.insol_NH   = fdat{5}(:);
milank.insol_glob = fdat{6}(:);
milank.precess    = milank.eccent.*sin(milank.peri);
fclose(fileID);
% Put them on our time domain
sTemp    = interp1(reconstruct.t,reconstruct.sTemp,St);
sTempN   = interp1(reconstruct.t,reconstruct.sTempN,St);
iceVol   = interp1(reconstruct.t,reconstruct.iceVol,St);
iceVolN  = interp1(reconstruct.t,reconstruct.iceVolN,St);
veg      = interp1(reconstruct.t,reconstruct.veg,St);
eccent   = interp1(milank.t,milank.eccent,St);
obliq    = interp1(milank.t,milank.obliq,St);
insol_NH = interp1(milank.t,milank.insol_NH,St);
precess  = interp1(milank.t,milank.precess,St);
RxNtemp  = sTemp + 270; % K
% A function to make a step function
stepFun = @(St,jump) St > jump;

%%% Compute the reaction rate coefficients with this temperature  %% HERE
k_cl_use  = k_cl(RxNtemp)  * DaysToS;   % cm3/molec/day
k_ch4_use = k_ch4(RxNtemp) * RxNconv;   % Tg/day
%k_ch4_use = k_ch4(mean(RxNtemp)*ones(size(RxNtemp))) * RxNconv;   % Tg/day

%%% Make a structure with the parameters
% Unit conversions
params.mm_ch4   = mm_ch4;
params.mm_oh    = mm_oh;
params.mm_co    = mm_co;
params.YrToDay  = YrToDay;
params.DaysToS  = DaysToS;
params.n_air    = n_air;
params.gmOH     = gmOH;
params.gmOHplot = gmOHplot;
% Rate constants, KIEs, and lifetimes
params.R_VPDB          = R_VPDB;
params.R_VSMOW         = R_VSMOW;
params.R_14CH4         = R_14CH4;
if tempdependence  % temperature dependent reactions
    params.k_ch4           = k_ch4_use; %k_ch4; 
    params.k_cl            = k_cl_use; % k_cl;  % temperature dependent reactions
else               % temperature independent reactions (NOT CORRECT)
    params.k_ch4           = k_12ch4 * ones(size(k_ch4_use));
    params.k_cl            = 6.3654e-14 * ones(size(k_cl_use));    % (cm3/molec/s) from Abbatt, Table 15-10 (k=2.36e-12*(T/298)^1.37*exp(-939/T) )
end
params.k_12ch4         = RxNconv * k_12ch4;   
params.k_ch4_strat_nh  = 1/(t_ch4_strat_nh * YrToDay);
params.k_ch4_strat_sh  = 1/(t_ch4_strat_sh * YrToDay);
params.k_ch4_soil      = k_ch4_soil;
params.k_ch4_cl        = k_ch4_cl;
params.k_14ch4_rad     = k_ch4_rad;
params.k_co            = RxNconv * k_co;
params.k_co_strat      = k_co_strat;
params.k_co_other      = k_co_other;
params.k_oh_strat      = k_oh_strat;
params.KIE_oh          = KIE_oh;
params.KIE_dD_oh       = KIE_dD_oh;
params.KIE_cl          = KIE_cl;
params.KIE_dD_cl       = KIE_dD_cl;
params.KIE_strat       = KIE_strat;
params.KIE_soil        = KIE_soil;
params.KIE_dD_soil     = KIE_dD_soil;
params.KIE_dD_strat    = KIE_dD_strat;
params.tau_NS          = tau_NS;
params.tau_NS_strat    = tau_NS_strat;
params.use_other_sinks = use_other_sinks;
params.k_co_flag       = k_co_flag;
params.RxNtemp         = RxNtemp;
% ODE parameters
params.Tspan   = Tspan;
params.IC      = IC;
params.odeOpts = opts;
% Climate forcings
params.sTemp    = sTemp;
params.sTempN   = sTempN;
params.iceVol   = iceVol;
params.iceVolN  = iceVolN;
params.veg      = veg;
params.eccent   = eccent;
params.obliq    = obliq;
params.insol_NH = insol_NH;
params.precess  = precess;
params.step     = stepFun;
% Atmospheric 14CO2
params.atmos14CO2     = d14c;
params.atmos14CO2_sig = d14c_sig;


end


%%% =======================================================================
%%% = END
%%% =======================================================================
