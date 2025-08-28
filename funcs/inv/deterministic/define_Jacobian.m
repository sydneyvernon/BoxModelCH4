%%% =======================================================================
%%% = define_Jacobian.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Constructs the Jacobian for a deterministic inversion.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St           -- Our time vector.
%%% =  ( 2): ems          -- Matrix with the emission sources (and OH).
%%% =  ( 3): IC           -- Vector with the initial conditions.
%%% =  ( 4): params       -- Structure with parameters for the box model.
%%% =  ( 5): run_parallel -- Are we running in parallel?  (True/False).
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): jacobian_ems -- Jacobian for the emission sources.
%%% =  ( 2): jacobian_IC  -- Jacobian for the ICs.
%%% =======================================================================

function [ jacobian ] = define_Jacobian( St, params, emsParams, IC, run_parallel)

%%% Get the unperturbed output
out_base = assembleObs(boxModel_wrapper(params,St,emsParams,IC));

%%% Assemble the state vector
state_vec = assembleStateVector(emsParams);

%%% Initialize the jacobians
nX = length(state_vec);     % Number of state vector elements
nY = length(out_base);      % Maximum number of observations
% Jacobian
jacobian = zeros(nY,nX);

%%% Perturbations for the Jacobian
perturb.oh_scale    = 0.01;
perturb.Q10         = 0.05;
perturb.amp_wet     = 1;
perturb.base_WT     = 1;
perturb.base_FF     = 1;
perturb.base_BB     = 1;
perturb.base_OH     = 20;
perturb.WT_scalings = 0.1*ones(size(emsParams.WT_scalings));
perturb.FF_scalings = 0.1*ones(size(emsParams.FF_scalings));
perturb.BB_scalings = 0.1*ones(size(emsParams.BB_scalings));
perturb.OH_scalings = 0.1*ones(size(emsParams.OH_scalings));
% Reshape into state vector format
perturb = assembleStateVector(perturb);

%%% Jacobian
if run_parallel % Parallel
    parfor i = 1:nX
        state_pos     = state_vec;
        state_pos(i)  = state_pos(i) + perturb(i);
        state_pos     = disassembleStateVector(state_pos,emsParams);
        out_plus      = assembleObs(boxModel_wrapper(params,St,state_pos,IC));
        state_neg     = state_vec;
        state_neg(i)  = state_neg(i) - perturb(i);
        state_neg     = disassembleStateVector(state_neg,emsParams);
        out_neg       = assembleObs(boxModel_wrapper(params,St,state_neg,IC));
        jacobian(:,i) = (out_plus - out_neg)./(2*perturb(i));
    end
else
    for i = 1:nX
        state_pos     = state_vec;
        state_pos(i)  = state_pos(i) + perturb(i);
        state_pos     = disassembleStateVector(state_pos,emsParams);
        out_plus      = assembleObs(boxModel_wrapper(params,St,state_pos,IC));
        state_neg     = state_vec;
        state_neg(i)  = state_neg(i) - perturb(i);
        state_neg     = disassembleStateVector(state_neg,emsParams);
        out_neg       = assembleObs(boxModel_wrapper(params,St,state_neg,IC));
        jacobian(:,i) = (out_plus - out_neg)./(2*perturb(i));
    end
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================
