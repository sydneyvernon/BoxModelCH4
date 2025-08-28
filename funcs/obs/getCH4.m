%%% =======================================================================
%%% = getCH4.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the atmospheric methane observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =  ( 2): reread  -- Structure that says if we'll re-read the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getCH4( dataDir, reread )

%%% Diagnostic
fprintf('   * CH4\n');


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/ch4_%4i-%4i_%s-%s.mat',...
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


%%% =======================================================================
%%% READ DATA
%%% =======================================================================

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;


%%% =======================================================================
%%% EDC CH4 data
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/iceCore',dataDir);

%%% Load the data
iceCoreDat = readtable(sprintf('%s/41586_2008_BFnature06950_MOESM33_ESM.xls',dataDirU));
out.tim.icecore = datenum(1950 - table2array(iceCoreDat(:,2)),1,1);
out.obs.icecore = table2array(iceCoreDat(:,3));
out.sig.icecore = 10*ones(size(out.tim.icecore));
% Add the Yan data?
iceCoreDat   = readtable(sprintf('%s/AllanHillsCH4.csv',dataDirU));
yanDat.tim   = datenum(1950 - table2array(iceCoreDat(:,5))*1000,1,1);
yanDat.ch4   = table2array(iceCoreDat(:,3));
yanDat.ch4E  = table2array(iceCoreDat(:,4));
yanDat.sFlag = table2array(iceCoreDat(:,7));
yanDat.flag  = zeros(size(yanDat.tim));
for i = 1:length(yanDat.sFlag)
    yanDat.flag(i) = strcmp('Yes',yanDat.sFlag{i}) || isnan(yanDat.tim(i)) || isnan(yanDat.ch4(i)) || isnan(yanDat.ch4E(i));
end
yanDat.tim   = yanDat.tim(~yanDat.flag);
yanDat.ch4   = yanDat.ch4(~yanDat.flag);
yanDat.ch4E  = yanDat.ch4E(~yanDat.flag);
yanDat.sFlag = yanDat.sFlag(~yanDat.flag);
yanDat.flag  = yanDat.flag(~yanDat.flag);
% Add the Yan data
out.tim.icecore = [out.tim.icecore(:);yanDat.tim];
out.obs.icecore = [out.obs.icecore(:);yanDat.ch4];
out.sig.icecore = [out.sig.icecore(:);yanDat.ch4E];

% Sort all the data
[~,ind]         = sort(out.tim.icecore,'ascend');
out.tim.icecore = out.tim.icecore(ind);
out.obs.icecore = out.obs.icecore(ind);
out.sig.icecore = out.sig.icecore(ind);
out.lat.icecore = -90;


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


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================