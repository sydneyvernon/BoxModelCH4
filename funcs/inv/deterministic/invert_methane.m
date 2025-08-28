%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Performs a deterministic inversion.
%%% =  ( 2): Depending on the flags, it can do a linear inversion or a
%%% =        non-linear inversion.  The non-linear inversion can use either
%%% =        Gauss-Newton or Levenberg-Marquardt.
%%% =  ( 3): Prior covariance matrix is defined in the subfunction called:
%%% =        "update_solution".cd ,,.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St           -- Our time vector.
%%% =  ( 2): obs          -- Structure with the observations.
%%% =  ( 3): ems_p        -- Matrix with the prior emissions (and OH).
%%% =  ( 4): IC_p         -- Vector with the prior initial conditions.
%%% =  ( 5): params       -- Structure with parameters for the box model.
%%% =  ( 6): linear       -- Are we doing a linear inversion?  (True/False)
%%% =  ( 7): run_parallel -- Are we running in parallel?  (True/False).
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): soln    -- Cell array with the emissions and ICs.
%%% =  ( 2): K_ems   -- Jacobian for the emission sources (and OH).
%%% =  ( 3): K_IC    -- Jacobian for the ICs.
%%% =  ( 4): reldiff -- Vector of relative differences while iterating.
%%% =  ( 5): absdiff -- Vector of absolute differences while iterating.
%%% =======================================================================

function [ soln, K, reldiff, absdiff] = invert_methane( St, params, emsParams_p, IC, obs, linear, run_parallel )

%%% Define tolerances
reltol  = 1e-6;
abstol  = 1e-6;
kmax    = 5;
k       = 1;
iter    = true;
reldiff = nan(kmax,1);
absdiff = nan(kmax,1);

%%% Initialize the Levenberg-Marquardt parameters (gamma = zero reverts to IMAP)
LM_param.chi2           = NaN;
LM_param.gamma          = 1e4;
%LM_param.gamma          = 0;
LM_param.cfn_true_old   = NaN;
LM_param.cfn_linear_old = NaN;

%%% Are we doing a linear or nonlinear inversion?
if linear
    LM_param.gamma = 0;
    iter           = false;
end

%%% Diagnostics
if linear
    fprintf('   * STARTING LINEAR INVERSION\n');
else
    if LM_param.gamma == 0
        fprintf('   * STARTING NONLINEAR IMAP INVERSION\n');
    else
        fprintf('   * STARTING NONLINEAR LEVENBERG-MARQUARDT INVERSION\n');
    end
end

%%% Get the solution for the first step
K               = define_Jacobian(St,params,emsParams_p,IC,run_parallel);
[soln,LM_param] = update_solution(St,params,IC,emsParams_p,emsParams_p,LM_param,K,obs);

% Diagnostic
%out_anal = boxModel_wrapper(params,St,soln{1},IC);
%plotResults( params.pSt, out_anal )

while iter
    tic;
    % Only update the jacobian if we accepted the previous step
    if ~isnan(LM_param.chi2) && LM_param.chi2 ~= chi2o
        K = define_Jacobian(St,params,emsParams_i,IC,run_parallel);
    end
    % Update parameters
    emsParams_i = soln{1};
    chi2o       = LM_param.chi2;
    % Get new solution
    [soln,LM_param] = update_solution(St,params,IC,emsParams_i,emsParams_p,LM_param,K,obs);
    % Get new relative difference
    state_vec_old = assembleStateVector(soln{1});
    state_vec_new = assembleStateVector(emsParams_i);
    absdiff       = abs(state_vec_old - state_vec_new);
    reldiff       = absdiff ./ min(abs(state_vec_old),abs(state_vec_new));
    absdiff(k)    = sum(absdiff);
    reldiff(k)    = rms(reldiff);
    % Diagnostic
    fprintf('ITERATION %3i/%3i: reldiff = %5.2e & absdiff = %5.2e\n',k,kmax,reldiff(k),absdiff(k));
    % Check convergence
    if (kmax <= k)
        iter = false;
    elseif ~isnan(LM_param.chi2) && LM_param.chi2 ~= chi2o % Did we accept the previous solution?
        k = k + 1;
        if (reldiff(k) < reltol) || (absdiff(k) < abstol)
            iter = false;
        end
    end
    toc
end

end

function [ out, LM_param ] = update_solution( St, params, IC, emsParams_i, emsParams_p, LM_param, K, obs )

%%% Get dimensions
nY = size(K,1);
nX = size(K,2);

%%% Assemble the state and observation vectors
F  = assembleObs(boxModel_wrapper(params,St,emsParams_i,IC));
y  = assembleObs(obs);
xp = assembleStateVector(emsParams_p);
xi = assembleStateVector(emsParams_i);

%%% Determine which rows to throw out (no observations!)
indGood = ~isnan(y);
y       = y(indGood);
F       = F(indGood);
K       = K(indGood,:);

%%% Create the prior error covariance matrix
% Components
prior_cov.oh_scale    = 0.1;
prior_cov.Q10         = 0.2;
prior_cov.amp_wet     = 2;
prior_cov.base_WT     = 2;
prior_cov.base_FF     = 2;
prior_cov.base_BB     = 2;
prior_cov.base_OH     = 200;
prior_cov.WT_scalings = 5*ones(size(emsParams_i.WT_scalings));
prior_cov.FF_scalings = 5*ones(size(emsParams_i.FF_scalings));
prior_cov.BB_scalings = 5*ones(size(emsParams_i.BB_scalings));
prior_cov.OH_scalings = 5*ones(size(emsParams_i.OH_scalings));
% Reshape into state vector format
Sa  = assembleStateVector(prior_cov);
Sa  = diag(Sa.^2);
tau = 0;

% Construct matrix
% Sa_ems = [Sa_ch4,Sa_ch4,Sa_ch4c13,  Sa_ch4c13,
%            Sa_oh, Sa_oh,    Sa_co, 0.07*Sa_co, Sa_tau, Sa_kx_NH, Sa_kx_SH];
% Sa     = diag(assembleStateVector(Sa_ems,Sa_IC));
% tau    = [tau_ch4,tau_ch4,tau_ch4c13,tau_ch4c13,
%           tau_oh, tau_oh,     tau_co,    tau_co, tau_tau]*365.25;
% Sa     = fillDiagonalsAnal(Sa,tau,St);

%%% Construct the observational error covariance matrix
So = [obs.ch4_err;obs.ch4_err;obs.ch4c13_err;obs.ch4c13_err].^2;
So = diag(So(indGood));

%%% Get inverse covariance matrices
if any(tau > 0)
    SaI = inv(Sa);
else
    SaI = diag(Sa);
    SaI = diag(1./SaI);
end
SoI = diag(So);
SoI = diag(1./SoI);

%%% Invert with the n-form from Rodgers
gamma = LM_param.gamma;
LHS   = (SaI + K'*SoI*K + gamma*SaI);
RHS   = (K'*SoI * (y - F) - SaI*(xi - xp));
dx    = LHS \ RHS;
xhat  = xi + dx;

%%% Are we doing Levenburg-Marquardt?
if gamma > 0
    % Is this the first step
    if isnan(LM_param.cfn_true_old) || isnan(LM_param.cfn_linear_old)
    	LM_param.cfn_true_old   = (y -    F)'*SoI*(y -    F) + (xi - xp)'*SaI*(xi - xp);
        LM_param.cfn_linear_old = (y - K*xi)'*SoI*(y - K*xi) + (xi - xp)'*SaI*(xi - xp);
    end
end

%%% Get the posterior covariance matrix and the averaging kernel matrix
S_hat  = inv(K'*SoI*K + SaI);
AK_hat = eye(size(Sa)) - S_hat*SaI;
C_hat  = sqrt(diag(S_hat));
C_hat  = diag(1./C_hat);
C_hat  = C_hat * S_hat * C_hat;

%%% Simulate the new concentrations
out_new = boxModel_wrapper(params,St,disassembleStateVector(xhat,emsParams_i),IC); plotResults( params.pSt, out_new );
F_new = assembleObs(out_new);
F_new = F_new(indGood);

%%% Measure performance of inversion step
% Get chi2
chi2_old = LM_param.chi2;
chi2_new = (xhat - xp)' / (SaI + K'*SoI*K) * (xhat - xp);
% Compute the change in the cost function
cfn_true_new   = (y -  F_new)'*SoI*(y -  F_new) + (xhat - xp)'*SaI*(xhat - xp);
cfn_linear_new = (y - K*xhat)'*SoI*(y - K*xhat) + (xhat - xp)'*SaI*(xhat - xp);
% Ratio of the true change in cost function to the linear approximation of the cost function
cfn_ratio = (cfn_true_new - LM_param.cfn_true_old) ./ (cfn_linear_new - LM_param.cfn_linear_old);

%%% How should we proceed?
if cfn_ratio < 0 || cfn_true_new > LM_param.cfn_true_old % This was a bad step
    gamma    = gamma * 10;
    xhat     = xi;
    chi2_new = chi2_old;
else % We accept the step, now change parameters
    LM_param.cfn_true_old   = cfn_true_new;
    LM_param.cfn_linear_old = cfn_linear_old;
    % Now let's evaluate gamma
    if cfn_ratio == 1 % Perfect linear approximation, get rid of gamma
        gamma = 0;
    elseif cfn_ratio < 0.25 % Need to increase gamma
        gamma = gamma * 2;
    elseif cfn_ratio > 0.75 % Can decrease gamma
        gamma = gamma / 3;
    end
end

%%% Store the Levenberg-Marquardt parameters
LM_param.gamma = gamma;
LM_param.chi2  = chi2_new;

%%% Disassemble the state vector
emsParams_hat = disassembleStateVector(xhat,emsParams_i);

%%% Make the other output
out = {emsParams_hat,S_hat,AK_hat,C_hat,Sa};

end

%%% Fill diagonals
function [ sig ] = fillDiagonalsAnal(sigO,tau,St)

%%% Populate the diagonals of the covariance matrices (temporal correlations)
sig = sigO;
for k = 1:length(tau)
    if tau(k) > 0
        kk = (k-1)*length(St) + 1;
        for i = kk:(kk+length(St)-1)
            ii = mod(i-1,length(St)) + 1;
            for j = kk:(kk+length(St)-1)
                jj = mod(j-1,length(St)) + 1;
                sig(i,j) = sqrt(sigO(i,i))*sqrt(sigO(j,j)) * exp(-abs(St(ii) - St(jj))/tau(k));
            end
        end
    end
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================

