%%% =======================================================================
%%% = get14CH4ems.m
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

function [ out ] = get14CH4ems( St, params, emsParams, ch4_ems )

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

if caseB || caseC
    % Delta14 source signatures
    sig_mod = params.atmos14CO2;
    sig_fos = 0;
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
    ems_nh = nh_wetTF.*sig_mod + nh_wetBF.*sig_mod + nh_firF.*sig_mod + nh_fosF.*sig_fos + nh_aniF.*sig_mod;
    ems_sh = sh_wetTF.*sig_mod + sh_wetBF.*sig_mod + sh_firF.*sig_mod + sh_fosF.*sig_fos + sh_aniF.*sig_mod;
end
% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================
