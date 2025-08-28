%%% =======================================================================
%%% = define_prior.m
%%% = Alex Turner
%%% = 03/30/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Defines the prior distribution.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): ems         -- Matrix with the emission sources (and OH).
%%% =  ( 2): IC          -- Vector with the ICs.
%%% =  ( 3): input_param -- A structure containing inputs to the inversion.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): p_prior -- Probability of the prior distribution.
%%% =======================================================================

function [ p_prior ] = define_prior( emsParams, input_param )

%%% Extract information from the input params
use_log     = input_param.use_log;
prior_state = input_param.emsParams;

%%% Components of the prior
% OH scaling factor (bounded normal)
%mode               = prior_state.oh_scale;
%sig                = exp(0.1).^2;
%func               = @(x) p_lognormal(x,mode,sig,use_log);
%prior_cov.oh_scale = func(emsParams.oh_scale);
mu                 = prior_state.oh_scale;
sig                = 0.1.^2;
xL                 = 0.01;
xU                 = 2;
func               = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
prior_cov.oh_scale = func(emsParams.Q10_tropical);
% Q10 (bounded normal)
mu                     = prior_state.Q10_tropical;
sig                    = 1.^2;
xL                     = 1;
xU                     = 25;
func                   = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
prior_cov.Q10_tropical = func(emsParams.Q10_tropical);
mu                     = prior_state.Q10_boreal;
sig                    = 1.^2;
xL                     = 1;
xU                     = 25;
func                   = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
prior_cov.Q10_boreal   = func(emsParams.Q10_boreal);
% Amplitude of the wetland emissions (normal)
mu                         = prior_state.amp_wet_tropical;
sig                        = 2.^2;
func                       = @(x) p_normal(x,mu,sig,use_log);
prior_cov.amp_wet_tropical = func(emsParams.amp_wet_boreal);
mu                         = prior_state.amp_wet_boreal;
sig                        = 2.^2;
func                       = @(x) p_normal(x,mu,sig,use_log);
prior_cov.amp_wet_boreal   = func(emsParams.amp_wet_boreal);
% Baseline wetland emissions (normal)
mu                = prior_state.base_WT;
sig               = 10.^2;
func              = @(x) p_normal(x,mu,sig,use_log);
prior_cov.base_WT = func(emsParams.base_WT);
% Baseline fossil emissions (normal)
mu                = prior_state.base_FF;
sig               = 10.^2;
func              = @(x) p_normal(x,mu,sig,use_log);
prior_cov.base_FF = func(emsParams.base_FF);
% Baseline fire emissions (normal)
mu                = prior_state.base_BB;
sig               = 10.^2;
func              = @(x) p_normal(x,mu,sig,use_log);
prior_cov.base_BB = func(emsParams.base_BB);
% Baseline OH emissions (normal)
mu                = prior_state.base_OH;
sig               = 50.^2;
func              = @(x) p_normal(x,mu,sig,use_log);
prior_cov.base_OH = func(emsParams.base_OH);
% Baseline strat-trop exchange (bounded normal)
mu                   = prior_state.base_tauTS;
sig                  = 1.^2;
xL                   = 4;
xU                   = 16;
func                 = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
prior_cov.base_tauTS = func(emsParams.base_tauTS);
% Baseline animal emissions (bounded normal)
mu                   = prior_state.base_AN;
sig                  = 5.^2;
xL                   = 0;
xU                   = 15;
func                 = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
prior_cov.base_AN = func(emsParams.base_AN);
% Baseline chlorine abundance (bounded normal)
mu                   = prior_state.base_CL;
%sig                  = (prior_state.base_CL*.05).^2;
sig                  = 100.^2;
xL                   = 0;
%xU                   = prior_state.base_CL*2;
xU                   = 2050;
func                 = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
prior_cov.base_CL = func(emsParams.base_CL);
% Timing of step function (uniform)
xL                   = min(input_param.St);
xU                   = max(input_param.St);
func                 = @ (x) p_uniform(x,xL,xU,use_log);
prior_cov.stepChange = func(emsParams.stepChange);
% Wetland scalings (normal)
prior_cov.WT_scalings = prior_state.WT_scalings;
for i = 1:length(prior_cov.WT_scalings)
    mu                       = prior_cov.WT_scalings(i);
    sig                      = 5.^2;
    func                     = @(x) p_normal(x,mu,sig,use_log);
    prior_cov.WT_scalings(i) = func(emsParams.WT_scalings(i));
end
% Fossil scalings (normal)
prior_cov.FF_scalings = prior_state.FF_scalings;
for i = 1:length(prior_cov.WT_scalings)
    mu                       = prior_cov.FF_scalings(i);
    sig                      = 5.^2;
    func                     = @(x) p_normal(x,mu,sig,use_log);
    prior_cov.FF_scalings(i) = func(emsParams.FF_scalings(i));
end
% Fire scalings (normal)
prior_cov.BB_scalings = prior_state.BB_scalings;
for i = 1:length(prior_cov.BB_scalings)
    mu                       = prior_cov.BB_scalings(i);
    sig                      = 5.^2;
    func                     = @(x) p_normal(x,mu,sig,use_log);
    prior_cov.BB_scalings(i) = func(emsParams.BB_scalings(i));
end
% OH scalings (normal)
prior_cov.OH_scalings = prior_state.OH_scalings;
for i = 1:length(prior_cov.OH_scalings)
    mu                       = prior_cov.OH_scalings(i);
    sig                      = 5.^2;
    func                     = @(x) p_normal(x,mu,sig,use_log);
    prior_cov.OH_scalings(i) = func(emsParams.OH_scalings(i));
end
% Strat-trop exchange (normal)
prior_cov.tau_scalings = prior_state.tau_scalings;
for i = 1:length(prior_cov.tau_scalings)
    mu                        = prior_state.tau_scalings(i);
    sig                       = 5.^2;
    func                      = @(x) p_normal(x,mu,sig,use_log);
    prior_cov.tau_scalings(i) = func(emsParams.tau_scalings(i));
end
% Animal scalings (normal)
prior_cov.AN_scalings = prior_state.AN_scalings;
for i = 1:length(prior_cov.AN_scalings)
    mu                       = prior_cov.AN_scalings(i);
    sig                      = 5.^2;
    func                     = @(x) p_normal(x,mu,sig,use_log);
    prior_cov.AN_scalings(i) = func(emsParams.AN_scalings(i));
end
% Chlorine scalings (normal)
prior_cov.CL_scalings = prior_state.CL_scalings;
for i = 1:length(prior_cov.CL_scalings)
    mu                       = prior_cov.CL_scalings(i);
    sig                      = 5.^2;
    func                     = @(x) p_normal(x,mu,sig,use_log);
    prior_cov.CL_scalings(i) = func(emsParams.CL_scalings(i));
end

%%% Initial conditions (uniform)
xL           = [ 300,300, -50,-50,    0,   0, -160,-160,  6, 6, 15,15]'; xL = [xL;xL];
xU           = [ 800,800, -44,-44, 1000,1000,  -20, -20, 12,12, 40,40]'; xU = [xU;xU];
func         = @ (x) p_uniform(x,xL,xU,use_log);
prior_cov.IC = func(emsParams.IC);

%%% Reshape into state vector format
priors  = assembleStateVector(prior_cov);

% % Diagnostic
% if any(isnan(priors))
%     fprintf('NH ch4 = %4.0f, ch4c13 = %4.0f, mcf = %4.0f, n2o = %4.0f | SH ch4 = %4.0f, ch4c13 = %4.0f, mcf = %4.0f, n2o = %4.0f | NH oh = %4.0f, SH oh = %4.0f, tau = %4.0f, IC = %4.0f \n',...
%             priors(1),priors(3),priors(5),priors(7),priors(2),priors(4),priors(6),priors(8),priors(9),priors(10),priors(11),priors(12));
% end
if use_log
    p_prior = sum(priors);
else
    p_prior = prod(priors);
end

end

%%% Fill diagonals
function [ sig ] = fillDiagonals(sig,tau,St)

%%% Populate the diagonals of the covariance matrices (temporal correlations)
n = length(sig);
for i = 1:n
for j = 1:n
    if i ~= j
      sig(i,j) = sqrt(sig(i,i))*sqrt(sig(j,j)) * exp(-abs(St(i) - St(j))/tau);
    end
end
end

end

%%% Distributioins
function [ p ] = p_lognormal(x,mode,sig,use_log)
mu = log(mode) + sig.^2;
if use_log
    p = log(logmvnpdf(x,mu,sig));
else
    p = logmvnpdf(x,mu,sig);
end
end
function [ p ] = p_normalB(x,mu,sig,xL,xU,use_log)
if any(x < xL) || any(xU < x)
    if use_log
        p = -1d6; % Make this a large negative number (rather than NaN)
    else
        p = 1d-6; % Make this a very small number (rather than NaN)
    end
else
	if use_log
        p = logmvnpdf(x,mu,sig);
    else
        p = mvnpdf(x,mu,sig);
    end
end
end
function [ p ] = p_normal(x,mu,sig,use_log)
if use_log
    p = logmvnpdf(x,mu,sig);
else
    p = mvnpdf(x,mu,sig);
end
end
function [ p ] = p_uniform(x,xL,xU,use_log)
if any(x < xL) || any(xU < x)
    if use_log
        p = -1d6; % Make this a large negative number (rather than NaN)
    else
        p = 1d-6; % Make this a very small number (rather than NaN)
    end
else
    p = 1;
	if use_log
        p = log(p);
    end
end
end


%%% =======================================================================
%%% = END
%%% =======================================================================
