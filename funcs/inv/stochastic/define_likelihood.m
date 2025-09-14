%%% =======================================================================
%%% = define_likelihood.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Defines the likelihood distribution.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): ems         -- Matrix with the emission sources (and OH).
%%% =  ( 2): IC          -- Vector with the ICs.
%%% =  ( 3): input_param -- A structure containing inputs to the inversion.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): p_like -- Probability of the likelihood distribution.
%%% =======================================================================

function [ p_like ] = define_likelihood( emsParams, input_param )

%%% Read the structure
obs     = input_param.obs;
params  = input_param.params;
use_log = input_param.use_log;

%%% Run the box model
out = boxModel_wrapper(params,params.St,emsParams);

%%% Define different types of distributions to use
if use_log
    p_normal  = @(x,mu,sig) logmvnpdf(x,mu,sig);
else
    p_normal  = @(x,mu,sig) mvnpdf(x,mu,sig);
end

%%% Define the likelihood distributions for each source
% CH4 NH (normal)
y     = out.nh_ch4;
mu    = obs.ch4_NH;
sig   = obs.ch4_NH_err;
ind   = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_ch4_NH = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% dD NH (normal)
y     = out.nh_dD;
mu    = obs.dD_NH;
sig   = obs.dD_NH_err;
ind   = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_dD_NH = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% CH4 (normal)
y     = out.sh_ch4;
mu    = obs.ch4;
sig   = obs.ch4_err;
ind   = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_ch4 = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% d13C (normal)
y      = out.d13c;
mu     = obs.d13c;
sig    = obs.d13c_err;
ind    = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_d13c = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% dD (normal)
y    = out.dD;
mu   = obs.dD;
sig  = obs.dD_err;
ind  = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_dD = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% 14CH4 (normal)
y      = out.d14c;
mu     = obs.d14c;
sig    = obs.d14c_err;
ind    = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_d14c = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% OH (normally distributed about ~10^6)
y    = out.oh;
mu   = params.gmOHplot*ones(size(y));
sig  = mu/8;
ind  = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_oh = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% tau_TS (normally distributed about 8yr)
y     = out.tau_TS;
mu    = 8*ones(size(y));
sig   = mu/3;
ind   = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_tau = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% Chlorine abundance (bounded normal)
y    = out.cl;
mu   = input_param.emsParams.base_CL*ones(size(y));
sig  = mu/6;
ind  = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_cl = p_normal(y(ind),mu(ind),diag(sig(ind).^2));

%%% Construct the full likelihood distribution
likeli = [p_ch4_NH, p_dD_NH, p_ch4, p_d13c, p_dD, p_d14c];
%likeli = [p_ch4, p_d13c, p_dD, p_d14c, p_oh, p_tau, p_cl];

% % Diagnostic
% if any(isnan(likeli))
%     fprintf('NH ch4 = %4.0f, ch4c13 = %4.0f, mcf = %4.0f, n2o = %4.0f | ch4 = %4.0f, ch4c13 = %4.0f, mcf = %4.0f, n2o = %4.0f SH\n',...
%              likeli(1),likeli(3),likeli(5),likeli(7),likeli(2),likeli(4),likeli(6),likeli(8));
% end
if use_log
    p_like = sum(likeli);
else
    p_like = prod(likeli);
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
