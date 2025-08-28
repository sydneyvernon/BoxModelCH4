%%% =======================================================================
%%% = boxModel.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): 2-box model for methane, delta13C, and methylchloroform.
%%% =  ( 2): Adapted from C. Frankenberg's box model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): t      -- Box model time.
%%% =  ( 2): y      -- Concentrations at time t.
%%% =  ( 3): St     -- Out time vector.
%%% =  ( 4): S      -- Emission sources for the box model.
%%% =  ( 5): tau_ST -- Strat-trop exchange (function of time).
%%% =  ( 6): params -- Structure with the parameters for the box model.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): dy -- Changes in the box model concentrations.
%%% =======================================================================

function [ dy ] = boxModel(t,y,St,S,tau_TS,kX_NH,kX_SH,k_ch4,k_cl,params)

%%% Initialize dy for troposphere & stratosphere (strat is same as trop)
% - Column 1:   12CH4 in the Northern Hemisphere
% - Column 2:   12CH4 in the Southern Hemisphere
% - Column 3:   13CH4 in the Northern Hemisphere
% - Column 4:   13CH4 in the Southern Hemisphere
% - Column 5:   14CH4 in the Northern Hemisphere
% - Column 6:   14CH4 in the Southern Hemisphere
% - Column 7:  12CH3D in the Northern Hemisphere
% - Column 8:  12CH3D in the Southern Hemisphere
% - Column 9:      OH in the Northern Hemisphere
% - Column 10:     OH in the Southern Hemisphere
% - Column 11:     CO in the Northern Hemisphere
% - Column 12:     CO in the Southern Hemisphere
dy = zeros(16,1);

%%% Interpolate sources, strat-trop exchange, and reaction rates in time
ind = find(t == St);        % Find the index for this timestep
if ~isempty(ind)            % Do we already have it?
    S      = S(ind(1),:);
    tau_TS = tau_TS(ind(1));
    kX_NH  = kX_NH(ind(1));
    kX_SH  = kX_SH(ind(1));
    k_ch4  = k_ch4(ind(1));
    k_cl   = k_cl(ind(1));
else                        % If not, we'll interpolate
    [pindex,index,slope] = lininterp1_ind(St,t);
    S      = S(pindex,:)      * (1 - slope) + slope * S(index,:);
    tau_TS = tau_TS(pindex,:) * (1 - slope) + slope * tau_TS(index,:);
    kX_NH  = kX_NH(pindex,:)  * (1 - slope) + slope * kX_NH(index,:);
    kX_SH  = kX_SH(pindex,:)  * (1 - slope) + slope * kX_SH(index,:);
    k_ch4  = k_ch4(pindex,:)  * (1 - slope) + slope * k_ch4(index,:);
    k_cl   = k_cl(pindex,:)   * (1 - slope) + slope * k_cl(index,:);
end
% Define the North-South exchange rates
tau_NS       = params.tau_NS;
tau_NS_strat = params.tau_NS_strat;

%%% Conserve mass between strat & trop, so:
% tau_ST = tau_TS/((p_surface-p_tropopause)/(p_tropopause-p_stratopause))
% where m_fac = (p_surface-p_tropopause)/(p_tropopause-p_stratopause)
%             = 5.7047
tau_ST = tau_TS / 5.7047;

%%% Reaction Rates
% Tropospheric reaction rates (y(9) & y(10) are the OH concentrations)
k_12ch4_NH    = (y( 9) * k_ch4);                    % NH
k_12ch4_SH    = (y(10) * k_ch4);                    % SH
k_13ch4_NH    = (y( 9) * k_ch4 / params.KIE_oh);    % NH
k_13ch4_SH    = (y(10) * k_ch4 / params.KIE_oh);    % SH
k_14ch4_NH    = (y( 9) * k_ch4 / params.KIE_oh^2);	% NH
k_14ch4_SH    = (y(10) * k_ch4 / params.KIE_oh^2);	% SH
k_12ch3d_NH   = (y( 9) * k_ch4 / params.KIE_dD_oh);	% NH
k_12ch3d_SH   = (y(10) * k_ch4 / params.KIE_dD_oh);	% SH
k_co_NH       = (y( 9) * params.k_co   );	% NH
k_co_SH       = (y(10) * params.k_co   );	% SH
k_ch4_cl      = k_cl;                           % Chlorine reaction
k_13ch4_cl    = k_ch4_cl / params.KIE_cl;       % Chlorine reaction
k_14ch4_cl    = k_ch4_cl / params.KIE_cl^2;     % Chlorine reaction
k_12ch3d_cl   = k_ch4_cl / params.KIE_dD_cl;	% Chlorine reaction
k_ch4_soil    = params.k_ch4_soil;                  % Soil uptake
k_13ch4_soil  = k_ch4_soil / params.KIE_soil;       % Soil uptake
k_14ch4_soil  = k_ch4_soil / params.KIE_soil^2;     % Soil uptake
k_12ch3d_soil = k_ch4_soil / params.KIE_dD_soil;	% Soil uptake
k_co_other    = params.k_co_other;	% Other CO loss pathways
k_14ch4_rad   = params.k_14ch4_rad;	% Radioactive decay
% Stratospheric reaction rates
k_ch4_strat_NH    = params.k_ch4_strat_nh;                          % NH
k_ch4_strat_SH    = params.k_ch4_strat_sh;                          % SH
k_13ch4_strat_NH  = params.k_ch4_strat_nh / params.KIE_strat;       % NH
k_13ch4_strat_SH  = params.k_ch4_strat_sh / params.KIE_strat;       % SH
k_14ch4_strat_NH  = params.k_ch4_strat_nh / params.KIE_strat^2;     % NH
k_14ch4_strat_SH  = params.k_ch4_strat_sh / params.KIE_strat^2;     % SH
k_12ch3d_strat_NH = params.k_ch4_strat_nh / params.KIE_dD_strat;    % NH
k_12ch3d_strat_SH = params.k_ch4_strat_sh / params.KIE_dD_strat;    % SH
k_oh_strat        = params.k_oh_strat;	% Single reaction rate
k_co_strat        = params.k_co_strat;	% Single reaction rate


%%% Compute dy in the troposphere
% 12CH4
dy( 1) = S( 1) + (y( 2)-y( 1))/tau_NS + (y(13)-y( 1))/tau_TS - y( 1)*(k_12ch4_NH+k_ch4_cl+k_ch4_soil);  % NH
dy( 2) = S( 2) + (y( 1)-y( 2))/tau_NS + (y(14)-y( 2))/tau_TS - y( 2)*(k_12ch4_SH+k_ch4_cl+k_ch4_soil);  % SH
% 13CH4
dy( 3) = S( 3) + (y( 4)-y( 3))/tau_NS + (y(15)-y( 3))/tau_TS - y( 3)*(k_13ch4_NH+k_13ch4_cl+k_13ch4_soil);  % NH
dy( 4) = S( 4) + (y( 3)-y( 4))/tau_NS + (y(16)-y( 4))/tau_TS - y( 4)*(k_13ch4_SH+k_13ch4_cl+k_13ch4_soil);  % SH
% 14CH4
dy( 5) = S( 5) + (y( 6)-y( 5))/tau_NS + (y(17)-y( 5))/tau_TS - y( 5)*(k_14ch4_NH+k_14ch4_cl+k_14ch4_soil+k_14ch4_rad);  % NH
dy( 6) = S( 6) + (y( 5)-y( 6))/tau_NS + (y(18)-y( 6))/tau_TS - y( 6)*(k_14ch4_SH+k_14ch4_cl+k_14ch4_soil+k_14ch4_rad);  % SH
% 12CH3D
dy( 7) = S( 7) + (y( 8)-y( 7))/tau_NS + (y(19)-y( 7))/tau_TS - y( 7)*(k_12ch3d_NH+k_12ch3d_cl+k_12ch3d_soil);  % NH
dy( 8) = S( 8) + (y( 7)-y( 8))/tau_NS + (y(20)-y( 8))/tau_TS - y( 8)*(k_12ch3d_SH+k_12ch3d_cl+k_12ch3d_soil);  % SH
% OH (allow the feedback?)
if params.interactive_OH
dy( 9) = S( 9) + (y(10)-y( 9))/tau_NS + (y(21)-y( 9))/tau_TS - y( 9)*k_12ch4_NH - y( 3)*k_13ch4_NH - y(11)*k_co_NH - y( 9)*kX_NH; % NH
dy(10) = S(10) + (y( 9)-y(10))/tau_NS + (y(22)-y(10))/tau_TS - y(10)*k_12ch4_SH - y( 4)*k_13ch4_SH - y(12)*k_co_SH - y(10)*kX_SH; % SH
end
% CO
dy(11) = S(11) + (y(12)-y(11))/tau_NS + (y(23)-y(11))/tau_TS - y(11)*(k_co_NH+k_co_other) + y( 1)*k_12ch4_NH + y( 3)*k_13ch4_NH; % NH
dy(12) = S(12) + (y(11)-y(12))/tau_NS + (y(24)-y(12))/tau_TS - y(12)*(k_co_SH+k_co_other) + y( 2)*k_12ch4_SH + y( 4)*k_13ch4_SH; % SH

%%% Compute dy in the stratosphere (assuming that all loss can be representated as first-order)
% Are we using the stratosphere?
if params.use_strat
% 12CH4
dy(13) = (y(14)-y(13))/tau_NS_strat + (y( 1)-y(13))/tau_ST - y(13)*k_ch4_strat_NH;   % NH
dy(14) = (y(13)-y(14))/tau_NS_strat + (y( 2)-y(14))/tau_ST - y(14)*k_ch4_strat_SH;   % SH
% 13CH4
dy(15) = (y(16)-y(15))/tau_NS_strat + (y( 3)-y(15))/tau_ST - y(15)*k_13ch4_strat_NH; % NH
dy(16) = (y(15)-y(16))/tau_NS_strat + (y( 4)-y(16))/tau_ST - y(16)*k_13ch4_strat_SH; % SH
% 14CH4
dy(17) = (y(18)-y(17))/tau_NS_strat + (y( 5)-y(17))/tau_ST - y(17)*(k_14ch4_strat_NH+k_14ch4_rad); % NH
dy(18) = (y(17)-y(18))/tau_NS_strat + (y( 6)-y(18))/tau_ST - y(18)*(k_14ch4_strat_SH+k_14ch4_rad); % SH
% 12CH3D
dy(19) = (y(20)-y(19))/tau_NS_strat + (y( 7)-y(19))/tau_ST - y(19)*k_12ch3d_strat_NH; % NH
dy(20) = (y(19)-y(20))/tau_NS_strat + (y( 8)-y(20))/tau_ST - y(20)*k_12ch3d_strat_SH; % SH
% OH
dy(21) = (y(22)-y(21))/tau_NS_strat + (y( 9)-y(21))/tau_ST - y(21)*k_oh_strat;       % NH
dy(22) = (y(21)-y(22))/tau_NS_strat + (y(10)-y(22))/tau_ST - y(22)*k_oh_strat;       % SH
% CO
dy(23) = (y(24)-y(23))/tau_NS_strat + (y(11)-y(23))/tau_ST - y(23)*k_co_strat;       % NH
dy(24) = (y(23)-y(24))/tau_NS_strat + (y(12)-y(24))/tau_ST - y(24)*k_co_strat;       % SH
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================