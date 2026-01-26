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

function [ out ] = getd14C_NH( dataDir, reread )

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

% three data pts from Summit, Greenland % is this just GRIP?
% we only have top depth and bottom depth, no gas age :(
% Measurements of 14CO and 14CH4 from large volume deep ice samples collecd at Summit, Greenland in May-June 2015.  
% All measurements corrected for extraneous carbon addition during sample processing and AMS measurment and normalized 
% to sample Delta13C.  Reported measuements for 'C14CO' and 'C14CH4' represents a combination of in situ cosmogenic 
% and paleoatmopsheric 14C, 'C14CH4_cosm_corr' removes the contribution from in situ cosmogenic 14CH4 while 'insitu 
% cosm C14CO' reports the in situ only component of C14CO. Units for all quantities are given in square brackets 
% and all reported errors are ± 1sigma
out.tim.icecore = datenum([0,0,0]'*1000,1,1);
out.obs.icecore = [89.65, 98.76, 99.52]'; % this is C14CH4 -- includes in situ cosmogenic.
out.sig.Taylor = [0.86, 0.87, 1]';

%%% sort - should be redundant
[~,ind]         = sort(out.tim.icecore,'ascend');
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
