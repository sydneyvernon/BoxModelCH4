%%% =======================================================================
%%% = getLifetime.m
%%% = Alex Turner
%%% = 06/03/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Computes the methane lifetime.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St  -- Our time vector.
%%% =  ( 2): ems -- Emission sources (and OH) for the box model.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): tau -- Structure with the methane lifetimes.
%%% =======================================================================

function [ tau ] = getLifetime( St, ems )

%%% Get some parameters for computing the methane lifetime
params  = getParameters(St);
k_12ch4 = 3.395e-15;

%%% Accounting for non-OH reactions (Holmes et al., 2018)
R = -0.27;

%%% Get time-dependent OH reaction rate
R_NH = ems.nh_oh * k_12ch4; % NH
R_SH = ems.sh_oh * k_12ch4; % SH

%%% Get the lifetimes
% Total
tau_NH  = 1 ./ ( R_NH * (1 + R) ) ./ (params.YrToDay*24*60*60);     % yr
tau_SH  = 1 ./ ( R_SH * (1 + R) ) ./ (params.YrToDay*24*60*60);     % yr
tau_GLO = mean([tau_NH,tau_SH],2);                                  % yr
% OH
tau_NHoh  = 1 ./ ( R_NH ) ./ (params.YrToDay*24*60*60);     % yr
tau_SHoh  = 1 ./ ( R_SH ) ./ (params.YrToDay*24*60*60);     % yr
tau_GLOoh = mean([tau_NHoh,tau_SHoh],2);                        % yr

%%% Make an output structure
tau.nh  = tau_NH;
tau.sh  = tau_SH;
tau.glo = tau_GLO;
tau.nhoh  = tau_NHoh;
tau.shoh  = tau_SHoh;
tau.glooh = tau_GLOoh;

end


%%% =======================================================================
%%% = END
%%% =======================================================================
