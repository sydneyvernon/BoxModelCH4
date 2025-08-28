%%% =======================================================================
%%% = assembleStateVector.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Transform the emission sources and ICs from matrice/vectors to
%%% =        a state vector for the inversion.
%%% =  ( 2): "disassembleStateVector" will do the opposite.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): ems -- Matrix with the emission sources (and OH).
%%% =  ( 2): IC  -- Vector with the initial conditions.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): x -- State vector.
%%% =======================================================================

function [ x ] = assembleStateVector( emsParams )
% State vector
x = [emsParams.oh_scale; emsParams.Q10_tropical; emsParams.amp_wet_tropical; emsParams.Q10_boreal; emsParams.amp_wet_boreal; emsParams.stepChange;...
     emsParams.base_WT; emsParams.base_FF; emsParams.base_BB; emsParams.base_OH; emsParams.base_tauTS; emsParams.base_AN; emsParams.base_CL;...
     emsParams.WT_scalings; emsParams.FF_scalings; emsParams.BB_scalings; emsParams.OH_scalings; emsParams.tau_scalings; emsParams.AN_scalings; emsParams.CL_scalings;...
     emsParams.IC];
end


%%% =======================================================================
%%% = END
%%% =======================================================================
