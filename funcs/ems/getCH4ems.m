%%% =======================================================================
%%% = getCH4ems.m
%%% =
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Create the methane emissions and puts them onto our temporal 
%%% =        grid.  Can either use constant emissions, constant with a 30
%%% =        Tg/yr jump in 2007, or EDGAR v4.2FT2010.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St      -- Our time vector.
%%% =  ( 2): tRes    -- String containing the temporal resolution.
%%% =  ( 3): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure containing the methane emissions.
%%% =======================================================================

function [ out ] = getCH4ems( St, params, emsParams )

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
    anth_ems = 355;
    nat_ems  = 225;
    % Spatial distribution of emissions
    frac_nh_anth = 0.9;     % Fraction of anthropogenic emissions in the NH
    frac_nh_nat  = 0.4;     % Fraction of natural emissions in the NH
    % Emissions per hemisphere
    ems_nat_nh  = zeros(size(St)) + nat_ems  * frac_nh_nat;
    ems_nat_sh  = zeros(size(St)) + nat_ems  * (1 - frac_nh_nat);
    ems_anth_nh = zeros(size(St)) + anth_ems * frac_nh_anth;
    ems_anth_sh = zeros(size(St)) + anth_ems * (1 - frac_nh_anth);
    % When do we add anthropogenic?
    indAnth = St >= datenum(1980,1,1);
    ems_anth_nh(~indAnth) = 0;
    ems_anth_sh(~indAnth) = 0;
    % CH4 Emissions
    ems_nh = ems_nat_nh + ems_anth_nh;
    ems_sh = ems_nat_sh + ems_anth_sh;
end
if caseB || caseC
    % Surface temperature
    sTemp   = emsParams.forcings(:,3) .* params.sTemp ./ params.sTempN;
    iceVolN = emsParams.forcings(:,2);
    veg     = 1 - iceVolN;
    % Fraction of emissions in the NH
    frac_nh = 0.4;      % Fraction of emissions in the NH
    % Emissions
    ems_wet_tropical = emsParams.amp_wet_tropical .* emsParams.Q10_tropical .^ ( (sTemp + emsParams.soilOffset)/10 );
    ems_wet_boreal   = veg .* emsParams.amp_wet_boreal   .* emsParams.Q10_boreal .^ ( (sTemp + emsParams.soilOffset)/10 );
    %ems_fossil = emsParams.base_fossil + emsParams.amp_fossil .*
    %params.iceVolN;  % why commented??
    %ems_fire   = emsParams.base_fire   + emsParams.amp_fire   .* params.veg;
    % Additional perturbations
    e_WT = zeros(length(St),1);
    e_FF = zeros(length(St),1);
    e_BB = zeros(length(St),1);
    e_AN = zeros(length(St),1);
    for i = 1:size(emsParams.forcings,2)
        e_WT = e_WT + emsParams.forcings(:,i) .* emsParams.WT_scalings(i);
        e_FF = e_FF + emsParams.forcings(:,i) .* emsParams.FF_scalings(i);
        e_BB = e_BB + emsParams.forcings(:,i) .* emsParams.BB_scalings(i);
        e_AN = e_AN + emsParams.forcings(:,i) .* emsParams.AN_scalings(i);
    end
    % Combine the baseline emissions and the additional perturbations
    ems_WTT = emsParams.base_WT + e_WT + ems_wet_tropical; ems_tropical = ems_WTT;
    ems_WTB = ems_wet_boreal;                              ems_boreal   = ems_WTB;
    ems_FF  = emsParams.base_FF + e_FF;                    ems_fossil = ems_FF;
    ems_AN  = emsParams.base_AN + e_AN;                    ems_animal = ems_AN;
    ems_BB  = emsParams.base_BB + e_BB;                    ems_fire = ems_BB;
    % Ensure positivity
    ems_tropical(ems_tropical < 0) = 0;
    ems_boreal(ems_boreal < 0)     = 0;
    ems_fossil(ems_fossil < 0)     = 0;
    ems_animal(ems_animal < 0)     = 0;
    ems_fire(ems_fire < 0)         = 0;
    % Partition Emissions to the hemispheres
    ems_wet_trop_nh = ems_tropical * frac_nh;
    ems_wet_bor_nh  = ems_boreal;
    ems_fire_nh     = ems_fire * frac_nh;
    ems_fossil_nh   = ems_fossil * frac_nh;
    ems_animal_nh   = ems_animal * frac_nh;
    ems_wet_trop_sh = ems_tropical * (1 - frac_nh);
    ems_wet_bor_sh  = ems_boreal * 0;
    ems_fire_sh     = ems_fire * (1 - frac_nh);
    ems_fossil_sh   = ems_fossil * (1 - frac_nh);
    ems_animal_sh   = ems_animal * (1 - frac_nh);
    % CH4 Emissions
    ems_nh = ems_wet_trop_nh + ems_wet_bor_nh + ems_fire_nh + ems_fossil_nh + ems_animal_nh;
    ems_sh = ems_wet_trop_sh + ems_wet_bor_sh + ems_fire_sh + ems_fossil_sh + ems_animal_sh;
    % Store sectoral emissions
    out.nh_wet_trop = ems_wet_trop_nh;
    out.nh_wet_bor  = ems_wet_bor_nh;
    out.nh_wet      = ems_wet_trop_nh + ems_wet_bor_nh;
    out.nh_fire     = ems_fire_nh;
    out.nh_fossil   = ems_fossil_nh;
    out.nh_animal   = ems_animal_nh;
    out.sh_wet_trop = ems_wet_trop_sh;
    out.sh_wet_bor  = ems_wet_bor_sh;
    out.sh_wet      = ems_wet_trop_sh + ems_wet_bor_sh;
    out.sh_fire     = ems_fire_sh;
    out.sh_fossil   = ems_fossil_sh;
    out.sh_animal   = ems_animal_sh;
end

% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================