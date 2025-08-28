%%% =======================================================================
%%% = getd14C.m
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

function [ out ] = getd14C( dataDir, reread )

%%% Diagnostic
fprintf('   * d14C\n');


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/d14c_%4i-%4i_%s-%s.mat',...
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

%%% Taylor glacier data from Dynosius et al., Science (2020): TG_deglacial14C_finaldata.xlsx
dat.tim.Taylor = datenum(1950 - [
    14.92, 14.86, 14.58, 14.54, 14.421, 14.42, 13.00, 10.131, 10.13,  9.21,  7.94]'*1000,1,1);
dat.obs.Taylor = [
      317,   327,   288,   287,   216,   204,   246,   126,   139,   141,    92]';
dat.sig.Taylor = [
      166,   151,   128,   112,   109,   111,    90,    58,    57,    63,    74]';
dat.lat.Taylor = -90;

% %%% Pakitsoq Greenland data from Petrenko et al., Science (2009): pakitsoq2009-14ch4.txt
% dat.tim.Greenland = datenum(1950 - [
%     11637, 11631, 11480, 11470, 11364, 11354]',1,1));
% dat.obs.Greenland = [
%     170.4, 161.8, 134.1, 159.8,  32.1, 112.2]';
% dat.sig.Greenland = [
%      50.3,  49.2,  43.0,  42.0,  40.4,  41.6]';
% dat.lat.Greenland = 69.43050;

%%% Taylor glacier data from Petrenko et al., Science (2017): Petrenko_TG_YD-PB_14CH4.xlsx
dat.tim.petrenko = datenum(1950 - [
    11715, 11559, 11515, 11453, 11357]',1,1);
dat.obs.petrenko = [
    192.2, 130.6, 125.8, 150.2, 132.3]';
dat.sig.petrenko = [
     52.4,  38.1,  35.3,  40.5,  35.1]';
dat.lat.petrenko = 69.43050;

%%% Combine the records
out.tim.icecore = [dat.tim.Taylor;dat.tim.petrenko];
out.obs.icecore = [dat.obs.Taylor;dat.obs.petrenko];
out.sig.icecore = [dat.sig.Taylor;dat.sig.petrenko];
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
