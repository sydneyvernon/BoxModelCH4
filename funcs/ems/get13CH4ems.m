%%% =======================================================================
%%% = get13CH4ems.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Creates the isotopic composition for the Northern and Southern
%%% =        hemispheric methane emissions.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St      -- Our time vector.
%%% =  ( 2): tRes    -- String containing the temporal resolution.
%%% =  ( 3): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure containing the d13C composition.
%%% =======================================================================

function [ out ] = get13CH4ems( St, params, emsParams, ch4_ems )

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

if const_ems
    %%% Define the delta13C composition
    % delta13C source
    nh_ch4c13 = -51.90; % permil
    sh_ch4c13 = -51.00; % permil
    % delta13C composition
    ems_nh = zeros(size(St)) + nh_ch4c13;
    ems_sh = zeros(size(St)) + sh_ch4c13;
end
if caseB || caseC
    % delta13C source signatures
    sig_wet_trop = -61.4;
    sig_wet_bor  = -62.6;
    sig_fos      = -44;
    sig_ani      = -65.4; % Sherwood et al. (2017)
    sig_fir      = -22.3;
    % Change in source signature during glacials (Moller et al., Nat Geo 2013; Bock et al., PNAS 2017)
    source_change = 2.6 * params.iceVolN;
    sig_wet_trop  = sig_wet_trop + source_change;
    sig_wet_bor   = sig_wet_bor  + source_change;
    sig_fir       = sig_fir      + source_change;
    sig_ani       = sig_ani      + source_change;
    % Fraction of emissions from each sector
    nh_wetTF = ch4_ems.nh_wet_trop ./ ch4_ems.nh;
    nh_wetBF = ch4_ems.nh_wet_bor  ./ ch4_ems.nh;
    nh_firF  = ch4_ems.nh_fire     ./ ch4_ems.nh;
    nh_fosF  = ch4_ems.nh_fossil   ./ ch4_ems.nh;
    nh_aniF  = ch4_ems.nh_animal   ./ ch4_ems.nh;
    sh_wetTF = ch4_ems.sh_wet_trop ./ ch4_ems.sh;
    sh_wetBF = ch4_ems.sh_wet_bor  ./ ch4_ems.sh;
    sh_firF  = ch4_ems.sh_fire     ./ ch4_ems.sh;
    sh_fosF  = ch4_ems.sh_fossil   ./ ch4_ems.sh;
    sh_aniF  = ch4_ems.sh_animal   ./ ch4_ems.sh;
    % Composition of sources
    ems_nh = nh_wetTF.*sig_wet_trop + nh_wetBF.*sig_wet_bor + nh_firF.*sig_fir + nh_fosF.*sig_fos + nh_aniF.*sig_ani;
    ems_sh = sh_wetTF.*sig_wet_trop + sh_wetBF.*sig_wet_bor + sh_firF.*sig_fir + sh_fosF.*sig_fos + sh_aniF.*sig_ani;
end
% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================
