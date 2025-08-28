%%% =======================================================================
%%% = getd13C.m
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

function [ out ] = getd13C( dataDir, reread )

%%% Diagnostic
fprintf('   * d13C\n');


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/d13c_%4i-%4i_%s-%s.mat',...
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

%%% Conversion from VPDB to NIWA scale (d13C(NIWA) = d13C + 0.169 permil)
NIWAconv = 0.169; % permi


%%% =======================================================================
%%% Ice Core data from Bock et al
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/iceCore',dataDir);

%%% EDC d13C
fspec       = '%f %f %f %f %f %f %f %f %f %f';
fileID      = fopen(sprintf('%s/EDC_d13CH4.tab',dataDirU),'r');
fdat        = textscan(fileID,fspec,'HeaderLines',25);
dat.tim.EDC = datenum(1950-fdat{2}(:)*1000,1,1);
dat.obs.EDC = fdat{10}(:);
dat.sig.EDC = fdat{6}(:);
fclose(fileID);

%%% TALDICE d13C
fspec           = '%f %f %f %f %f %f %f %f %f %f';
fileID          = fopen(sprintf('%s/TALDICE_d13CH4.tab',dataDirU),'r');
fdat            = textscan(fileID,fspec,'HeaderLines',25);
dat.tim.TALDICE = datenum(1950-fdat{2}(:)*1000,1,1);
dat.obs.TALDICE = fdat{10}(:);
dat.sig.TALDICE = fdat{6}(:);
fclose(fileID);

%%% Vostok d13C
fspec          = '%f %f %f %f %f %f %f %f %f %f';
fileID         = fopen(sprintf('%s/Vostok_d13CH4.tab',dataDirU),'r');
fdat           = textscan(fileID,fspec,'HeaderLines',26);
dat.tim.Vostok = datenum(1950-fdat{2}(:)*1000,1,1);
dat.obs.Vostok = fdat{10}(:);
dat.sig.Vostok = fdat{6}(:);
fclose(fileID);

%%% Combine the records
out.tim.icecore = [dat.tim.EDC;dat.tim.TALDICE;dat.tim.Vostok];
out.obs.icecore = [dat.obs.EDC;dat.obs.TALDICE;dat.obs.Vostok];
out.sig.icecore = [dat.sig.EDC;dat.sig.TALDICE;dat.sig.Vostok];
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
