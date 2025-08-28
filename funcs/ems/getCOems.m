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

function [ out ] = getCOems( St, params, emsParams )

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
if const_ems || caseA
    % Total emissions
    anth_ems = 588;  % Tg/yr (Lamarque et al., 2010; Yin et al., 2015)
    nat_ems  = 381;  % Tg/yr (Van der Werf et al., 2010; Yin et al., 2015)
    % Spatial distribution of emissions
    frac_nh_anth = 0.95; % Fraction of anthropogenic emissions in the NH
    frac_nh_nat  = 0.90; % Fraction of natural emissions in the NH
    % Emissions per hemisphere
    ems_nat_nh  = zeros(size(St)) + nat_ems  * frac_nh_nat;
    ems_nat_sh  = zeros(size(St)) + nat_ems  * (1 - frac_nh_nat);
    ems_anth_nh = zeros(size(St)) + anth_ems * frac_nh_anth;
    ems_anth_sh = zeros(size(St)) + anth_ems * (1 - frac_nh_anth);
    % When do we add anthropogenic?
    indAnth = St >= datenum(1980,1,1);
    ems_anth_nh(~indAnth) = 0;
    ems_anth_sh(~indAnth) = 0;
    % CO Emissions
    ems_nh = ems_nat_nh + ems_anth_nh;
    ems_sh = ems_nat_sh + ems_anth_sh;
end
if caseB || caseC
    % Fraction of emissions in the NH
    frac_nh = 0.90;
    % Emissions
%    ems_fire  = emsParams.fireFacCO * (emsParams.base_fire + params.veg  * emsParams.amp_fire);
    ems_ocean = emsParams.base_ocean*ones(size(St));
    % Additional perturbations
    e_BB = zeros(length(St),1);
    for i = 1:size(emsParams.forcings,2)
        e_BB = e_BB + emsParams.forcings(:,i) .* emsParams.BB_scalings(i);
    end
    % Combine the baseline emissions and the additional perturbations
    ems_BB = emsParams.fireFacCO * (emsParams.base_BB + e_BB); ems_fire = ems_BB;
    % Ensure positivity
    ems_fire(ems_fire < 0)   = 0;
    ems_ocean(ems_ocean < 0) = 0;
    % Partition Emissions to the hemispheres
    ems_fire_nh  = ems_fire  * frac_nh;
    ems_ocean_nh = ems_ocean * frac_nh;
    ems_fire_sh  = ems_fire  * (1 - frac_nh);
    ems_ocean_sh = ems_ocean * (1 - frac_nh);
    % CO Emissions
    ems_nh = ems_fire_nh + ems_ocean_nh;
    ems_sh = ems_fire_sh + ems_ocean_sh;
    % Store sectoral emissions
    out.nh_fire  = ems_fire_nh;
    out.nh_ocean = ems_ocean_nh;
    out.sh_fire  = ems_fire_sh;
    out.sh_ocean = ems_ocean_sh;
end
% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================