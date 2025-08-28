%%% =======================================================================
%%% = getCOems.m
%%% = Alex Turner
%%% = 05/02/2017
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Create the CO emissions and puts them onto our temporal 
%%% =        grid.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St      -- Our time vector.
%%% =  ( 2): tRes    -- String containing the temporal resolution.
%%% =  ( 3): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure containing the methane emissions.
%%% =======================================================================

function [ out ] = getOHems( St, params, emsParams )

%%% Figure out the case
const_ems = true; % This is the fallback case
caseA     = false;
caseB     = false;
caseC     = false;
caseD     = false;
% Sort out the case info
if strcmp(params.caseName,'caseAa')     || strcmp(params.caseName,'caseAb')
    caseA = true;
elseif strcmp(params.caseName,'caseBa') || strcmp(params.caseName,'caseBb')
    caseB = true;
elseif strcmp(params.caseName,'caseCa') || strcmp(params.caseName,'caseCb')
    caseC = true;
elseif strcmp(params.caseName,'caseDa') || strcmp(params.caseName,'caseDb')
    caseD = true;
else
    const_ems = true;
end

%%% Which emissions do we want to use?
const_ems  = true;
if const_ems
    % Total OH source (constant)
    tot_oh  = 2400; % Tg/yr (Murray et al., 2014 estimate is ~2300 Tg/yr)
    frac_nh = 0.50; % Fraction of emissions in the NH
    % OH Emissions
    ems_nh = zeros(size(St)) + tot_oh * frac_nh;
    ems_sh = zeros(size(St)) + tot_oh * (1 - frac_nh);
end
if caseC
    % Fraction of emissions in the NH
    frac_nh = 0.50;
    % OH Emissions
%    ems_oh = emsParams.base_oh + emsParams.amp_oh .* params.veg;
    % Additional perturbations
    e_OH = zeros(length(St),1);
    for i = 1:size(emsParams.forcings,2)
        e_OH = e_OH + emsParams.forcings(:,i) .* emsParams.OH_scalings(i);
    end
    % Combine the baseline emissions and the additional perturbations
    ems_OH = emsParams.base_OH + e_OH; ems_oh = ems_OH;
    % Ensure positivity
    ems_oh(ems_oh < 0) = 0;
    % OH Emissions per hemisphere
    ems_nh = ems_oh  * frac_nh;
    ems_sh = ems_oh  * (1 - frac_nh);
end

% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================