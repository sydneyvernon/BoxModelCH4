%%% =======================================================================
%%% = boxModel_wrapper.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): A wrapper for the box model.
%%% =  ( 2): Precomputes terms for the box model for speed.
%%% =  ( 3): Runs the box model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St     -- Our time vector.
%%% =  ( 2): ems    -- Structure with emission sources (and OH) for the box model.
%%% =  ( 3): IC     -- Initial conditions for the box model.
%%% =  ( 4): params -- Structure with parameters for the box model.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure with the simulated concentrations.
%%% =======================================================================

function [ out ] = boxModel_wrapper(params,St,emsParams)

%%% ==================
%%% Load the emissions
%%% ==================

%%% Get the initial conditions
IC = emsParams.IC;

%%% Convert isotopes from delta notation to the form needed for the box model
% 13CH4
IC(3)  =  IC(1) * params.R_VPDB * (1 +  IC(3)/1000);
IC(4)  =  IC(2) * params.R_VPDB * (1 +  IC(4)/1000);
% 14CH4
IC(5)  =  IC(1) * params.R_14CH4 * (1 +  IC(5)/1000);
IC(6)  =  IC(2) * params.R_14CH4 * (1 +  IC(6)/1000);
% DCH4
IC(7)  =  IC(1) * params.R_VSMOW * (1 +  IC(7)/1000);
IC(8)  =  IC(2) * params.R_VSMOW * (1 +  IC(8)/1000);
% 13CH4 (stratosphere)
IC(15) = IC(13) * params.R_VPDB * (1 + IC(15)/1000);
IC(16) = IC(14) * params.R_VPDB * (1 + IC(16)/1000);
% 14CH4 (stratosphere)
IC(17) = IC(13) * params.R_14CH4 * (1 + IC(17)/1000);
IC(18) = IC(14) * params.R_14CH4 * (1 + IC(18)/1000);
% DCH4 (stratosphere)
IC(19) = IC(13) * params.R_VSMOW * (1 + IC(19)/1000);
IC(20) = IC(14) * params.R_VSMOW * (1 + IC(20)/1000);

% Convert OH from molec/cm3 to form for box model
IC(9)  =  IC(9) / params.n_air*1d9 * 1d5;
IC(10) = IC(10) / params.n_air*1d9 * 1d5;
IC(21) = IC(21) / params.n_air*1d9 * 1d5;
IC(22) = IC(22) / params.n_air*1d9 * 1d5;

% update stepChange
emsParams.forcings(:,end) = params.step(St,emsParams.stepChange);

%%% Get the emissions
ch4_ems    = getCH4ems(St,params,emsParams);             % CH4 (Tg/yr)
ch4c13_ems = get13CH4ems(St,params,emsParams,ch4_ems);   % d13C (permil)
ch4c14_ems = get14CH4ems(St,params,emsParams,ch4_ems);   % d14C (permil)
ch3D_ems   = get12CH3Dems(St,params,emsParams,ch4_ems);  % dD (permil)
oh_ems     = getOHems(St,params,emsParams);              % OH (Tg/yr)
co_ems     = getCOems(St,params,emsParams);              % CO (Tg/yr)

%%% Other sources/sinks
% Strat-trop exchange
tau_TS  = emsParams.base_tauTS * ones(length(St),1); % years
e_tauTS = zeros(length(St),1);
for i = 1:size(emsParams.forcings,2)
    e_tauTS = e_tauTS + emsParams.forcings(:,i) .* emsParams.tau_scalings(i);
end
tau_TS = tau_TS + e_tauTS; tau_TS(tau_TS < 1e-6) = 1e-6;
if ~params.use_strat; tau_TS(:) = 1e6; end
% Arbitrary reactions with OH
kX_NH = emsParams.oh_scale*ones(length(St),1); % s^-1
kX_SH = emsParams.oh_scale*ones(length(St),1); % s^-1
% Chlorine
cl_conc = emsParams.base_CL*ones(length(St),1);	% molec/cm3
e_cl    = zeros(length(St),1);
for i = 1:size(emsParams.forcings,2)
    e_cl = e_cl + emsParams.forcings(:,i) .* emsParams.CL_scalings(i);
end
cl_conc = cl_conc + e_cl;

%%% Put all this into a single structure with 11 fields:
% - NH CH4 emissions
% - SH CH4 emissions
% - NH CH4C13 composition
% - SH CH4C13 composition
% - NH OH emissions
% - SH OH emissions
% - NH CO emissions
% - SH CO emissions
% - Strat-trop exchange
% - NH arbitrary OH reaction rate
% - SH arbitrary OH reaction rate
ems.nh_ch4   = ch4_ems.nh;
ems.sh_ch4   = ch4_ems.sh;
ems.nh_13ch4 = ch4c13_ems.nh;
ems.sh_13ch4 = ch4c13_ems.sh;
ems.nh_14ch4 = ch4c14_ems.nh;
ems.sh_14ch4 = ch4c14_ems.sh;
ems.nh_ch3d  = ch3D_ems.nh;
ems.sh_ch3d  = ch3D_ems.sh;
ems.nh_oh    = oh_ems.nh;
ems.sh_oh    = oh_ems.sh;
ems.nh_co    = co_ems.nh;
ems.sh_co    = co_ems.sh;
ems.tau_TS   = tau_TS;
ems.kX_NH    = kX_NH;
ems.kX_SH    = kX_SH;
ems.cl_conc  = cl_conc;
% Convert the structure to a matrix
S = assembleEms(ems);


%%% ==================
%%% Structure the emissions
%%% ==================

%%% Set up the emissions for the box model
% Convert CH4, OH, and CO emissions to units of per day
S(:,[1,2,9,10,11,12]) = S(:,[1,2,9,10,11,12]) / params.YrToDay;
% 12CH4
S(:,1)  = 2/params.mm_ch4*S(:,1);                           % NH (factor of two is for mass of one hemisphere)
S(:,2)  = 2/params.mm_ch4*S(:,2);                           % SH
% 13CH4
S(:,3)  = S(:,1) .* params.R_VPDB  .* (1 + S(:,3)/1000);    % NH
S(:,4)  = S(:,2) .* params.R_VPDB  .* (1 + S(:,4)/1000);    % SH
% 14CH4
S(:,5)  = S(:,1) .* params.R_14CH4 .* (1 + S(:,5)/1000);	% NH (precompute S1*S5)
S(:,6)  = S(:,2) .* params.R_14CH4 .* (1 + S(:,6)/1000);	% SH (precompute S2*S6)
% 12CH3D
S(:,7)  = S(:,1) .* params.R_VSMOW .* (1 + S(:,7)/1000);	% NH
S(:,8)  = S(:,2) .* params.R_VSMOW .* (1 + S(:,8)/1000);	% SH
% OH
S(:,9)  = 2/params.mm_oh*S(:,9);                            % NH
S(:,10) = 2/params.mm_oh*S(:,10);                           % SH
% CO
S(:,11) = 2/params.mm_co*S(:,11);                           % NH
S(:,12) = 2/params.mm_co*S(:,12);                           % SH
% Strat/Trop exchange
tau_TS = S(:,13) * params.YrToDay;
% Arbitrary reaction with OH
kX_NH = S(:,14) * (60 * 60 * 24);                           % NH
kX_SH = S(:,15) * (60 * 60 * 24);                           % SH
% Reaction rates
k_cl  = params.k_cl .* S(:,16);                             % CH4 + Cl (1 / day)
k_ch4 = params.k_ch4;                                       % CH4 + OH (Tg / day)


%%% ==================
%%% Run the box model 
%%% ================== 

%%% Make all timesteps positive
t_offset = 0;
if St(1) < 0
    t_offset = 0 - St(1);
end
St = St + t_offset;

%%% Run the box model with ode15s
[T, F] = ode15s(@(t,y) boxModel(t,y,St,S,tau_TS,kX_NH,kX_SH,k_ch4,k_cl,params),St,IC,params.odeOpts);

%%% Convert back to our original time domain
St = St - t_offset;
T  = T  - t_offset;

%%% ================== 
%%% Extract the results
%%% ==================

%%% Isotope conversion
F(:,3)  = ( (F(:, 3)./F(:,1))  ./ params.R_VPDB  - 1) * 1000;
F(:,4)  = ( (F(:, 4)./F(:,2))  ./ params.R_VPDB  - 1) * 1000;
F(:,5)  = ( (F(:, 5)./F(:,1))  ./ params.R_14CH4 - 1) * 1000;
F(:,6)  = ( (F(:, 6)./F(:,2))  ./ params.R_14CH4 - 1) * 1000;
F(:,7)  = ( (F(:, 7)./F(:,1))  ./ params.R_VSMOW - 1) * 1000;
F(:,8)  = ( (F(:, 8)./F(:,2))  ./ params.R_VSMOW - 1) * 1000;
F(:,15) = ( (F(:,15)./F(:,13)) ./ params.R_VPDB  - 1) * 1000;
F(:,16) = ( (F(:,16)./F(:,14)) ./ params.R_VPDB  - 1) * 1000;
F(:,17) = ( (F(:,17)./F(:,13)) ./ params.R_14CH4 - 1) * 1000;
F(:,18) = ( (F(:,18)./F(:,14)) ./ params.R_14CH4 - 1) * 1000;
F(:,19) = ( (F(:,19)./F(:,13)) ./ params.R_VSMOW - 1) * 1000;
F(:,20) = ( (F(:,20)./F(:,14)) ./ params.R_VSMOW - 1) * 1000;

%%% Make the output structure
% Do we need to interpolate?
if any(T ~= St)
    [pindex,index,slope] = lininterp1_ind(T,St);
    F                    = F(pindex,:) * (1 - slope) + slope * F(index,:);
end
% Store the simulated concentrations
out.ch4             = mean([F(:,1),F(:,2)],2);
out.nh_ch4          = F(:,1);
out.sh_ch4          = F(:,2);
out.d13c            = mean([F(:,3),F(:,4)],2);
out.nh_d13c         = F(:,3);
out.sh_d13c         = F(:,4);
out.d14c            = mean([F(:,5),F(:,6)],2);
out.nh_d14c         = F(:,5);
out.sh_d14c         = F(:,6);
out.dD              = mean([F(:,7),F(:,8)],2);
out.nh_dD           = F(:,7);
out.sh_dD           = F(:,8);
out.oh              = mean([F(:,9),F(:,10)],2) * params.n_air/1e9;  % Convert to molec/cm3
out.nh_oh           = F(:,9) * params.n_air/1e9;                    % Convert to molec/cm3
out.sh_oh           = F(:,10) * params.n_air/1e9;                   % Convert to molec/cm3
out.co              = mean([F(:,11),F(:,12)],2);
out.nh_co           = F(:,11);
out.sh_co           = F(:,12);
out.ch4_strat       = mean([F(:,13),F(:,14)],2);
out.nh_ch4_strat    = F(:,9);
out.sh_ch4_strat    = F(:,10);
out.d13c_strat      = mean([F(:,15),F(:,16)],2);
out.nh_d13c_strat   = F(:,15);
out.sh_d13c_strat   = F(:,16);
out.d14c_strat      = mean([F(:,17),F(:,18)],2);
out.nh_d14c_strat   = F(:,17);
out.sh_d14c_strat   = F(:,18);
out.dD_strat        = mean([F(:,19),F(:,20)],2);
out.nh_dD_strat     = F(:,19);
out.sh_dD_strat     = F(:,20);
out.oh_strat        = mean([F(:,21),F(:,22)],2) * params.n_air/1e9; % Convert to molec/cm3
out.nh_oh_strat     = F(:,21) * params.n_air/1e9;                   % Convert to molec/cm3
out.sh_oh_strat     = F(:,22) * params.n_air/1e9;                   % Convert to molec/cm3
out.co_strat        = mean([F(:,23),F(:,24)],2);
out.nh_co_strat     = F(:,23);
out.sh_co_strat     = F(:,24);
out.cl              = cl_conc; % molec/cm3
% Also store the inputs
out.ch4_ems       = ch4_ems.nh+ch4_ems.sh;
out.nh_ch4_ems    = ch4_ems.nh;
out.sh_ch4_ems    = ch4_ems.sh;
out.ch4c13_ems    = mean([ch4c13_ems.nh,ch4c13_ems.sh],2); % we don't save these ems for other isotopologues?
out.nh_ch4c13_ems = ch4c13_ems.nh;
out.sh_ch4c13_ems = ch4c13_ems.sh;
out.oh_ems        = oh_ems.nh+oh_ems.sh;
out.nh_oh_ems     = oh_ems.nh;
out.sh_oh_ems     = oh_ems.sh;
out.co_ems        = co_ems.nh+co_ems.sh;
out.nh_co_ems     = co_ems.nh;
out.sh_co_ems     = co_ems.sh;
% Save out sectoral emissions
out.ch4_ems_wet_tropical = ch4_ems.nh_wet_trop + ch4_ems.sh_wet_trop;
out.ch4_ems_wet_boreal   = ch4_ems.nh_wet_bor  + ch4_ems.sh_wet_bor;
out.ch4_ems_fire         = ch4_ems.nh_fire     + ch4_ems.sh_fire;
out.ch4_ems_fossil       = ch4_ems.nh_fossil   + ch4_ems.sh_fossil;
out.ch4_ems_animal       = ch4_ems.nh_animal   + ch4_ems.sh_animal;
out.co_ems_fire          = co_ems.nh_fire      + co_ems.sh_fire;
out.co_ems_ocean         = co_ems.nh_ocean     + co_ems.sh_ocean;
% Compute OH loss rate
out.oh_loss_12ch4 = params.mm_oh.*params.YrToDay/(params.n_air/1e9) * params.k_12ch4                .* out.oh .* out.ch4;
out.oh_loss_13ch4 = params.mm_oh.*params.YrToDay/(params.n_air/1e9) * params.k_12ch4/params.KIE_oh  .* out.oh .* out.ch4.*(1 + out.d13c/1000);
out.oh_loss_ch4   = out.oh_loss_12ch4 + out.oh_loss_13ch4;
out.oh_loss_co    = params.mm_oh.*params.YrToDay/(params.n_air/1e9) * params.k_co                   .* out.oh .* out.co;
out.oh_loss_other = params.mm_oh.*params.YrToDay/(params.n_air/1e9) * mean([kX_NH,kX_SH],2)         .* out.oh;
out.oh_loss       = out.oh_loss_ch4 + out.oh_loss_co + out.oh_loss_other;
% Get methane lifetime
out.ch4_lifetime_nh = 1 ./ (F(:, 9).*k_ch4 + k_cl + params.k_ch4_soil);
out.ch4_lifetime_sh = 1 ./ (F(:,10).*k_ch4 +k_cl + params.k_ch4_soil);
out.ch4_lifetime    = mean([out.ch4_lifetime_nh,out.ch4_lifetime_sh],2);
% Add the strat-trop exchange rate
out.tau_TS = tau_TS / params.YrToDay;

%%% Set any 14C to NaN for older than 50ky BP
nan14c = St < datenum(-50*1000,1,1);
out.d14c(nan14c)          = NaN;
out.nh_d14c(nan14c)       = NaN;
out.sh_d14c(nan14c)       = NaN;
out.d14c_strat(nan14c)    = NaN;
out.nh_d14c_strat(nan14c) = NaN;
out.sh_d14c_strat(nan14c) = NaN;

%%% Are we plotting?
if params.diagnostics
    close all; plotResults( params, out, params.obs ); pause(0.001);
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================
