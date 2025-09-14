%%% =======================================================================
%%% = getdD.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the N Hemisphere dD observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =  ( 2): reread  -- Structure that says if we'll re-read the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getdD_NH( dataDir, reread )

%%% Diagnostic
fprintf('   * dD\n');


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/dD_NH_%4i-%4i_%s-%s.mat',...
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
%%% Ice Core data from Bock et al
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/iceCore',dataDir);

%%% EDC dD
fdat   = readtable(sprintf('%s/Bock_NGRIP_dD.csv',dataDirU));
dat.tim.NGRIP = datenum(1950-table2array(fdat(:,3)),1,1);
dat.obs.NGRIP = table2array(fdat(:,6));  % "corrected" measurement column
dat.sig.NGRIP = 3 * ones(size(table2array(fdat(:,6)))); % guess - slightly larger than max err on the EDC dD.
%fclose(fileID);

% %%% EDML dD
% fspec           = '%f %f %f %f %f %f %f %f %f %f';
% fileID          = fopen(sprintf('%s/EDML_dDCH4.tab',dataDirU),'r');
% fdat            = textscan(fileID,fspec,'HeaderLines',19);
% dat.tim.EDML = datenum(1950-fdat{2}(:)*1000,1,1);
% dat.obs.EDML = fdat{6}(:);
% dat.sig.EDML = fdat{4}(:);
% fclose(fileID);

%%% Combine the records
out.tim.icecore = [dat.tim.NGRIP];
out.obs.icecore = [dat.obs.NGRIP];
out.sig.icecore = [dat.sig.NGRIP];
[~,ind]         = sort(out.tim.icecore,'ascend'); % should be redundant
out.tim.icecore = out.tim.icecore(ind);
out.obs.icecore = out.obs.icecore(ind);
out.sig.icecore = out.sig.icecore(ind);
out.lat.icecore = 90;

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
