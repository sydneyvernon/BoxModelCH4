%%% =======================================================================
%%% = getdD.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the d13C observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =  ( 2): reread  -- Structure that says if we'll re-read the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getdD( dataDir, reread )

%%% Diagnostic
fprintf('   * dD\n');


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/dD_%4i-%4i_%s-%s.mat',...
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
fspec       = '%f %f %f %f %f %f %f %f %f %f';
fileID      = fopen(sprintf('%s/EDC_dDCH4.tab',dataDirU),'r');
fdat        = textscan(fileID,fspec,'HeaderLines',19);
dat.tim.EDC = datenum(1950-fdat{2}(:)*1000,1,1);
dat.obs.EDC = fdat{6}(:);
dat.sig.EDC = fdat{4}(:);
fclose(fileID);

%%% EDML dD
fspec           = '%f %f %f %f %f %f %f %f %f %f';
fileID          = fopen(sprintf('%s/EDML_dDCH4.tab',dataDirU),'r');
fdat            = textscan(fileID,fspec,'HeaderLines',19);
dat.tim.EDML = datenum(1950-fdat{2}(:)*1000,1,1);
dat.obs.EDML = fdat{6}(:);
dat.sig.EDML = fdat{4}(:);
fclose(fileID);

%%% Combine the records
out.tim.icecore = [dat.tim.EDC;dat.tim.EDML];
out.obs.icecore = [dat.obs.EDC;dat.obs.EDML];
out.sig.icecore = [dat.sig.EDC;dat.sig.EDML];
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
