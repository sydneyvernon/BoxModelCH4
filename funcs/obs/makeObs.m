%%% =======================================================================
%%% = makeObs.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Aggregates the observations into Northern and Southern
%%% =        hemispheric averages and puts them onto my temporal grid.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St         -- Our time vector.
%%% =  ( 2): tAvg       -- String defining our averaging period.
%%% =  ( 3): ch4_obs    -- Structure containing the methane obs.
%%% =  ( 4): ch4c13_obs -- Structure containing the d13C obs.
%%% =  ( 5): co_obs     -- Structure containing the carbon monoxide obs.
%%% =  ( 6): ch4c13_obs -- Structure containing the dD obs.
%%% =  ( 7): dataDir    -- Directory containing the data.
%%% =  ( 8): reread     -- Structure that says if we'll re-read the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = makeObs( St, tAvg, ch4_NH_obs, ch4_obs, d13c_obs, dD_obs, d14c_obs, reread )

%%% Diagnostic
fprintf('\n *** MAKING THE OBSERVATION STRUCTURE *** \n');

%%% Create the output structure
out = struct;

%%% Get the years
yrs = datevec(St);
yrs = yrs(:,1);


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/InputData_%4i-%4i_%s-%s.mat',...
                  reread.dir,reread.sYear,reread.eYear,reread.tRes,reread.tAvg);

%%% Load pre-existing data file
if ~reread.flag
    % Check if a file exists
    if exist(OutName, 'file') == 2
        fprintf('   * LOADING OLD OBS STRUCTURE\n');
        load(OutName);
        return % Don't need to read the data
    end
end

%%% Get the averaging window
fDays = 365.25; % Number of days in the block average
if strcmp(tAvg,'month') || strcmp(tAvg,'MONTH') || strcmp(tAvg,'monthly')
    fDays = fDays / 12;
elseif strcmp(tAvg,'kyear') || strcmp(tAvg,'KYR') || strcmp(tAvg,'kyr')
    fDays = fDays * 1000;
end


%%% =======================================================================
%%% CH4 (NH)
%%% =======================================================================

    %%% Read the observation structure
    fprintf('   * CH4_NH\n');
    
    %%% Load the data
    tDatI = ch4_NH_obs.tim.icecore;
    yDatI = ch4_NH_obs.obs.icecore;
    eDatI = ch4_NH_obs.sig.icecore;
    
    %%% Block average
    [tDat, yDat, eDat] = BlockAverage_AltError(tDatI,yDatI,ones(size(tDatI)),fDays);
    
    %%% Put the data on our temporal grid
    oDat = interp1(tDat,yDat,St);
    eDat = interp1(tDat,eDat,St);
    % Force the minimum uncertainty to 10 ppb
    eDat(eDat < 10)   = 10;
    eDat(isnan(eDat)) = 10;
    
    %%% Scale error by factor of 2
    eDat = 2*eDat;
    
    %%% Make sure we didn't lose any data
    indMiss = find(isnan(oDat));
    tDat    = tDat(~isnan(yDat));
    yDat    = yDat(~isnan(yDat));
    for i = 1:length(indMiss)
        [tOff,~] = min(abs(St(indMiss(i)) - tDat));
        if tOff < fDays
            oDat(indMiss(i)) = interp1(tDat,yDat,St(indMiss(i)));
        end
    end
    
    %%% Put the data in the output structure
    out.ch4_NH     = oDat;
    out.ch4_NH_err = eDat;


% %%% =======================================================================
% %%% dD (NH)
% %%% =======================================================================
% 
% %%% Read the observation structure
% fprintf('   * dD_NH\n');
% 
% %%% Load the data
% tDatI = dD_NH_obs.tim.icecore;
% yDatI = dD_NH_obs.obs.icecore;
% eDatI = dD_NH_obs.sig.icecore;
% 
% %%% Block average
% [tDat, yDat, eDat] = BlockAverage_AltError(tDatI,yDatI,ones(size(tDatI)),fDays);
% 
% %%% Put the data on our temporal grid
% oDat = interp1(tDat,yDat,St);
% eDat = interp1(tDat,eDat,St);
% % Force the minimum uncertainty
% min_uncert              = min(eDatI);
% eDat(eDat < min_uncert) = min_uncert;
% eDat(isnan(eDat))       = max(eDatI);
% 
% %%% Scale error by factor of 3
% eDat = 3*eDat;
% 
% %%% Make sure we didn't lose any data
% indMiss = find(isnan(oDat));
% tDat    = tDat(~isnan(yDat));
% yDat    = yDat(~isnan(yDat));
% for i = 1:length(indMiss)
%     [tOff,~] = min(abs(St(indMiss(i)) - tDat));
%     if tOff < fDays
%         oDat(indMiss(i)) = interp1(tDat,yDat,St(indMiss(i)));
%     end
% end
% 
% %%% Put the data in the output structure
% out.dD_NH     = oDat;
% out.dD_NH_err = eDat;
% 

%%% =======================================================================
%%% CH4 (SH)
%%% =======================================================================

%%% Read the observation structure
fprintf('   * CH4\n');

%%% Load the data
tDatI = ch4_obs.tim.icecore;
yDatI = ch4_obs.obs.icecore;
eDatI = ch4_obs.sig.icecore;

%%% Block average
[tDat, yDat, eDat] = BlockAverage_AltError(tDatI,yDatI,ones(size(tDatI)),fDays);

%%% Put the data on our temporal grid
oDat = interp1(tDat,yDat,St);
eDat = interp1(tDat,eDat,St);
% Force the minimum uncertainty to 10 ppb
eDat(eDat < 10)   = 10;
eDat(isnan(eDat)) = 10;

%%% Scale error by factor of 2
eDat = 2*eDat;

%%% Make sure we didn't lose any data
indMiss = find(isnan(oDat));
tDat    = tDat(~isnan(yDat));
yDat    = yDat(~isnan(yDat));
for i = 1:length(indMiss)
    [tOff,~] = min(abs(St(indMiss(i)) - tDat));
    if tOff < fDays
        oDat(indMiss(i)) = interp1(tDat,yDat,St(indMiss(i)));
    end
end

%%% Put the data in the output structure
out.ch4     = oDat;
out.ch4_err = eDat;


%%% =======================================================================
%%% d13C
%%% =======================================================================

%%% Read the observation structure
fprintf('   * d13C\n');

%%% Load the data
tDatI = d13c_obs.tim.icecore;
yDatI = d13c_obs.obs.icecore;
eDatI = d13c_obs.sig.icecore;

%%% Block averave
[tDat, yDat, eDat] = BlockAverage_AltError(tDatI,yDatI,ones(size(tDatI)),fDays);

%%% Put the data on our temporal grid
oDat = interp1(tDat,yDat,St);
eDat = interp1(tDat,eDat,St);
% Force the minimum uncertainty
min_uncert              = min(eDatI);
eDat(eDat < min_uncert) = min_uncert;
eDat(isnan(eDat))       = max(eDatI);

%%% Scale error by factor of 4
eDat = 4*eDat;

%%% Make sure we didn't lose any data
indMiss = find(isnan(oDat));
tDat    = tDat(~isnan(yDat));
yDat    = yDat(~isnan(yDat));
for i = 1:length(indMiss)
    [tOff,~] = min(abs(St(indMiss(i)) - tDat));
    if tOff < fDays
        oDat(indMiss(i)) = interp1(tDat,yDat,St(indMiss(i)));
    end
end

%%% Put the data in the output structure
out.d13c     = oDat;
out.d13c_err = eDat;


%%% =======================================================================
%%% dD
%%% =======================================================================

%%% Read the observation structure
fprintf('   * dD\n');

%%% Load the data
tDatI = dD_obs.tim.icecore;
yDatI = dD_obs.obs.icecore;
eDatI = dD_obs.sig.icecore;

%%% Block average
[tDat, yDat, eDat] = BlockAverage_AltError(tDatI,yDatI,ones(size(tDatI)),fDays);

%%% Put the data on our temporal grid
oDat = interp1(tDat,yDat,St);
eDat = interp1(tDat,eDat,St);
% Force the minimum uncertainty
min_uncert              = min(eDatI);
eDat(eDat < min_uncert) = min_uncert;
eDat(isnan(eDat))       = max(eDatI);

%%% Scale error by factor of 3
eDat = 3*eDat;

%%% Make sure we didn't lose any data
indMiss = find(isnan(oDat));
tDat    = tDat(~isnan(yDat));
yDat    = yDat(~isnan(yDat));
for i = 1:length(indMiss)
    [tOff,~] = min(abs(St(indMiss(i)) - tDat));
    if tOff < fDays
        oDat(indMiss(i)) = interp1(tDat,yDat,St(indMiss(i)));
    end
end

%%% Put the data in the output structure
out.dD     = oDat;
out.dD_err = eDat;


%%% =======================================================================
%%% d14C
%%% =======================================================================

%%% Read the observation structure
fprintf('   * d14C\n');

%%% Load the data
tDatI = d14c_obs.tim.icecore;
yDatI = d14c_obs.obs.icecore;
eDatI = d14c_obs.sig.icecore;

%%% Block average
%[tDat, yDat, eDat] = BlockAverage_AltError(tDatI,yDatI,ones(size(tDatI)),1000);

%%% Put the data on our temporal grid
dT = St(2) - St(1);
oDat = nan(size(St));
eDat = nan(size(St));
for i = 1:length(St)
    iUse    = (St(i) - dT/2) < tDatI & tDatI <= (St(i) + dT/2);
    if any(iUse)
        oDat(i) = nansum(yDatI(iUse)./eDatI(iUse))./nansum(1./eDatI(iUse));
        eDat(i) = 1/sqrt(sum(iUse))*sqrt(nansum(eDatI(iUse).^2));
    end
end
% Force the minimum uncertainty
min_uncert              = min(eDatI);
eDat(eDat < min_uncert) = min_uncert;
eDat(isnan(eDat))       = max(eDatI);

%%% Scale error by factor of 2
%eDat = 2*eDat;

%%% Put the data in the output structure
out.d14c     = oDat;
out.d14c_err = eDat;


%%% =======================================================================
%%% SAVE THIS OBSERVATION FILE
%%% =======================================================================

%%% Save the structure
fprintf('   * SAVING OBS STRUCTURE\n');
if exist(OutName, 'file') == 2
    delete(OutName);
end
save(OutName,'out');


end

%%% Read the observation structure
function [ oDat_NH, oDat_SH, eDat_NH, eDat_SH ] = ReadObsStruct( t, tAvg, obs )

%%% Define parameters for the block averaging
fDays = 365.25; % Number of days in the block average
if strcmp(tAvg,'month') || strcmp(tAvg,'MONTH') || strcmp(tAvg,'monthly')
    fDays = fDays / 12;
end

%%% Get the site names
sNames = fieldnames(obs.obs);

%%% Initialize a vector for the observation times and 
tDat_NH = [];
yDat_NH = [];
iDat_NH = [];
tDat_SH = [];
yDat_SH = [];
iDat_SH = [];

%%% Get the data from the structure and sort it
% Get the data
for i = 1:length(sNames);
    % Get data for this site
    lat  = obs.lat.(sNames{i});
    tDat = obs.tim.(sNames{i});
    yDat = obs.obs.(sNames{i});
    % Remove the seasonal cycle
    yDat_noSeas = DeseasonalizeData(tDat,yDat,fDays);
    % Throw out NaNs
    ind = ~isnan(tDat) & ~isnan(yDat_noSeas);
    % Require at least a 5-year record
    tDat        = tDat(ind);
    yDat_noSeas = yDat_noSeas(ind);
    if sum(ind) > 0
        if ( abs(tDat(end) - tDat(1)) > 365.25*5 )
            % Sort by latitude
            if lat > 0 % NH
                tDat_NH = [tDat_NH;tDat];
                yDat_NH = [yDat_NH;yDat_noSeas];
                iDat_NH = [iDat_NH;repmat(i,size(tDat))];
            end
            if lat < 0 % SH
                tDat_SH = [tDat_SH;tDat];
                yDat_SH = [yDat_SH;yDat_noSeas];
                iDat_SH = [iDat_SH;repmat(i,size(tDat))];
            end
        end
    end
end

%%% Allow dD measurements to have a 3 year record
if length(yDat_NH) < 1
    tDat_NH = [];
    yDat_NH = [];
    iDat_NH = [];
    tDat_SH = [];
    yDat_SH = [];
    iDat_SH = [];
    for i = 1:length(sNames);
        % Get data for this site
        lat  = obs.lat.(sNames{i});
        tDat = obs.tim.(sNames{i});
        yDat = obs.obs.(sNames{i});
        % Remove the seasonal cycle
        yDat_noSeas = DeseasonalizeData(tDat,yDat,fDays);
        % Throw out NaNs
        ind = ~isnan(tDat) & ~isnan(yDat_noSeas);
        % Require at least a 5-year record
        tDat        = tDat(ind);
        yDat_noSeas = yDat_noSeas(ind);
        if sum(ind) > 0
            if ( abs(tDat(end) - tDat(1)) > 365.25*3 )
                % Sort by latitude
                if lat > 0 % NH
                    tDat_NH = [tDat_NH;tDat];
                    yDat_NH = [yDat_NH;yDat_noSeas];
                    iDat_NH = [iDat_NH;repmat(i,size(tDat))];
                end
                if lat < 0 % SH
                    tDat_SH = [tDat_SH;tDat];
                    yDat_SH = [yDat_SH;yDat_noSeas];
                    iDat_SH = [iDat_SH;repmat(i,size(tDat))];
                end
            end
        end
    end
end

% Sort
[~,ind] = sort(tDat_NH);
tDat_NH = tDat_NH(ind);
yDat_NH = yDat_NH(ind);
iDat_NH = iDat_NH(ind);
[~,ind] = sort(tDat_SH);
tDat_SH = tDat_SH(ind);
yDat_SH = yDat_SH(ind);
iDat_SH = iDat_SH(ind);

%%% Bootstrap the mean and error
nBoot    = 50;
tDatO_NH = [];
yDatO_NH = [];
tDatO_SH = [];
yDatO_SH = [];
for i = 1:nBoot
    
    %%% NH
    uniqID = unique(iDat_NH);   % Get the unique IDs
    %uniqID
    nDraw  = length(uniqID);    % How many different timeseries do we have?
    % Randomly sample from our sites
    tDat = [];
    yDat = [];
    for j = 1:nDraw
        ind  = randsample(uniqID,1) == iDat_NH;
        tDat = [tDat;tDat_NH(ind)];
        yDat = [yDat;yDat_NH(ind)];
    end
    %size(tDat)
    % Block average
    [tDat, yDat] = BlockAverage(tDat,yDat,ones(size(tDat)),fDays);
    % Store the data
    tDatO_NH = [tDatO_NH;tDat];
    yDatO_NH = [yDatO_NH;yDat];
    
    %%% SH
    uniqID = unique(iDat_SH);   % Get the unique IDs
    nDraw  = length(uniqID);    % How many different timeseries do we have?
    % Randomly sample from our sites
    tDat = [];
    yDat = [];
    for j = 1:nDraw
        ind  = randsample(uniqID,1) == iDat_SH;
        tDat = [tDat;tDat_SH(ind)];
        yDat = [yDat;yDat_SH(ind)];
    end
    % Block average
    [tDat, yDat] = BlockAverage(tDat,yDat,ones(size(tDat)),fDays);
    % Store the data
    tDatO_SH = [tDatO_SH;tDat];
    yDatO_SH = [yDatO_SH;yDat];
    
end

%%% Block average to get the mean and error
% Remove NaNs
ind      = ~isnan(tDatO_NH) & ~isnan(yDatO_NH);
tDatO_NH = tDatO_NH(ind);
yDatO_NH = yDatO_NH(ind);
ind      = ~isnan(tDatO_SH) & ~isnan(yDatO_SH);
tDatO_SH = tDatO_SH(ind);
yDatO_SH = yDatO_SH(ind);
% Sort
[~,ind]  = sort(tDatO_NH);
tDatO_NH = tDatO_NH(ind);
yDatO_NH = yDatO_NH(ind);
[~,ind]  = sort(tDatO_SH);
tDatO_SH = tDatO_SH(ind);
yDatO_SH = yDatO_SH(ind);
% Block average
[tDat_NH, yDat_NH, eDat_NH] = BlockAverage_AltError(tDatO_NH,yDatO_NH,ones(size(tDatO_NH)),fDays);
[tDat_SH, yDat_SH, eDat_SH] = BlockAverage_AltError(tDatO_SH,yDatO_SH,ones(size(tDatO_SH)),fDays);

%%% Put the data on our temporal grid
oDat_NH = interp1(tDat_NH,yDat_NH,t);
eDat_NH = interp1(tDat_NH,eDat_NH,t);
oDat_SH = interp1(tDat_SH,yDat_SH,t);
eDat_SH = interp1(tDat_SH,eDat_SH,t);

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================
