%%% =======================================================================
%%% = disassembleStateVector.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Break the state vector into a matrix of emission sources and a
%%% =        vector of ICs.
%%% =  ( 2): "assembleStateVector" will do the opposite.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): x  -- State vector.
%%% =  ( 2): nT -- Number of time steps in our time vector.
%%% =  ( 3): nE -- Number of different emission sources.
%%% =  ( 4): nI -- Number of ICs.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): ems -- Matrix with the emission sources (and OH).
%%% =  ( 2): IC  -- Vector with the initial conditions.
%%% =======================================================================

function [ emsParams ] = disassembleStateVector( x, inParams )

% Required parameters
emsParams.soilOffset = inParams.soilOffset;
emsParams.fireFacCO  = inParams.fireFacCO;
emsParams.forcings   = inParams.forcings;
emsParams.base_ocean = inParams.base_ocean;
% Parameters
emsParams.oh_scale         = x(1);
emsParams.Q10_tropical     = x(2);
emsParams.amp_wet_tropical = x(3);
emsParams.Q10_boreal       = x(4);
emsParams.amp_wet_boreal   = x(5);
emsParams.stepChange       = x(6);
% Baseline emissions
emsParams.base_WT    = x(7);
emsParams.base_FF    = x(8);
emsParams.base_BB    = x(9);
emsParams.base_OH    = x(10);
emsParams.base_tauTS = x(11);
emsParams.base_AN    = x(12);
emsParams.base_CL    = x(13);
% Scalings
ii = 14;
nI = length(inParams.WT_scalings)-1;
emsParams.WT_scalings  = x(ii:ii+nI); ii = ii+nI+1;
emsParams.FF_scalings  = x(ii:ii+nI); ii = ii+nI+1;
emsParams.BB_scalings  = x(ii:ii+nI); ii = ii+nI+1;
emsParams.OH_scalings  = x(ii:ii+nI); ii = ii+nI+1;
emsParams.tau_scalings = x(ii:ii+nI); ii = ii+nI+1;
emsParams.AN_scalings  = x(ii:ii+nI); ii = ii+nI+1;
emsParams.CL_scalings  = x(ii:ii+nI); ii = ii+nI+1;
% ICs
nIC = length(inParams.IC)-1;
emsParams.IC = x(ii:ii+nIC);

end


%%% =======================================================================
%%% = END
%%% =======================================================================
