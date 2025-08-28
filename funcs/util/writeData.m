%%% =======================================================================
%%% = writeData.m
%%% = Alex Turner
%%% = 06/06/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Saves the output as a csv file.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St       -- Our time vector.
%%% =  ( 2): obs      -- Structure containing the observations.
%%% =  ( 3): ems      -- Emission sources (and OH) for the box model.
%%% =  ( 4): mod      -- Structure containing the model results.
%%% =  ( 5): IC       -- Initial conditions.
%%% =  ( 6): baseName -- Prefix for the plots.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = writeData( pSt, out, outDir, caseName )

%%% Write out the results
% Simulated CH4
head  = {'kyr since 1950','ch4 (ppb)'};
fname = sprintf('%s/csv/%s/ch4.csv',outDir,caseName);
dat   = [pSt,out.ch4];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [head;num2cell(dat)];
% Simulated OH
head  = {'kyr since 1950','oh (molec/cm3)'};
fname = sprintf('%s/csv/%s/oh.csv',outDir,caseName);
dat   = [pSt,out.oh];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Simulated CO
head  = {'kyr since 1950','co (ppb)'};
fname = sprintf('%s/csv/%s/co.csv',outDir,caseName);
dat   = [pSt,out.co];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Simulated d13c
head  = {'kyr since 1950','d13c (permill)'};
fname = sprintf('%s/csv/%s/d13c.csv',outDir,caseName);
dat   = [pSt,out.d13c];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Simulated dD
head  = {'kyr since 1950','dD (permill)'};
fname = sprintf('%s/csv/%s/dD.csv',outDir,caseName);
dat   = [pSt,out.dD];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Simulated d14c
head  = {'kyr since 1950','d14c (permill)'};
fname = sprintf('%s/csv/%s/d14c.csv',outDir,caseName);
dat   = [pSt,out.d14c];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Total CH4 emissions
head  = {'kyr since 1950','ch4_ems (Tg/yr)'};
fname = sprintf('%s/csv/%s/ems_ch4.csv',outDir,caseName);
dat   = [pSt,out.ch4_ems_wet_tropical + out.ch4_ems_wet_boreal + out.ch4_ems_fire + out.ch4_ems_fossil];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Wetland CH4 emissions
head  = {'kyr since 1950','ch4_ems_wet (Tg/yr)'};
fname = sprintf('%s/csv/%s/ems_ch4_wet.csv',outDir,caseName);
dat   = [pSt,out.ch4_ems_wet_tropical + out.ch4_ems_wet_boreal];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Tropical wetland CH4 emissions
head  = {'kyr since 1950','ch4_ems_wet_tropical (Tg/yr)'};
fname = sprintf('%s/csv/%s/ems_ch4_wet_tropical.csv',outDir,caseName);
dat   = [pSt,out.ch4_ems_wet_tropical];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Boreal wetland CH4 emissions
head  = {'kyr since 1950','ch4_ems_boreal (Tg/yr)'};
fname = sprintf('%s/csv/%s/ems_ch4_wet_boreal.csv',outDir,caseName);
dat   = [pSt,out.ch4_ems_wet_boreal];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Fire CH4 emissions
head  = {'kyr since 1950','ch4_ems_fire (Tg/yr)'};
fname = sprintf('%s/csv/%s/ems_ch4_fire.csv',outDir,caseName);
dat   = [pSt,out.ch4_ems_fire];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Fossil CH4 emissions
head  = {'kyr since 1950','ch4_ems_fossil (Tg/yr)'};
fname = sprintf('%s/csv/%s/ems_ch4_fossil.csv',outDir,caseName);
dat   = [pSt,out.ch4_ems_fossil];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Total OH production
head  = {'kyr since 1950','oh_prod (Tg/yr)'};
fname = sprintf('%s/csv/%s/prod_oh.csv',outDir,caseName);
dat   = [pSt,out.oh_ems];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Total OH loss
head  = {'kyr since 1950','oh_loss (Tg/yr)'};
fname = sprintf('%s/csv/%s/loss_oh.csv',outDir,caseName);
dat   = [pSt,out.oh_loss];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% OH loss to CH4
head  = {'kyr since 1950','oh_loss_ch4 (Tg/yr)'};
fname = sprintf('%s/csv/%s/loss_oh_ch4.csv',outDir,caseName);
dat   = [pSt,out.oh_loss_ch4];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% OH loss to CO
head  = {'kyr since 1950','oh_loss_co (Tg/yr)'};
fname = sprintf('%s/csv/%s/loss_oh_co.csv',outDir,caseName);
dat   = [pSt,out.oh_loss_co];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% OH loss to other
head  = {'kyr since 1950','oh_loss_other (Tg/yr)'};
fname = sprintf('%s/csv/%s/loss_oh_other.csv',outDir,caseName);
dat   = [pSt,out.oh_loss_other];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Methane lifetime
head  = {'kyr since 1950','ch4_lifetime (yr)'};
fname = sprintf('%s/csv/%s/ch4_lifetime.csv',outDir,caseName);
dat   = [pSt,out.ch4_lifetime];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];
% Stratospheric residence time
head  = {'kyr since 1950','strat_res_time (yr)'};
fname = sprintf('%s/csv/%s/tau_ST.csv',outDir,caseName);
dat   = [pSt,out.tau_TS / 5.7047];
writetable(cell2table([head;num2cell(dat(~isnan(dat(:,2)),:))]),fname,'writevariablenames',0);
full_struct = [full_struct(:,:),[head(2);num2cell(dat(:,2))]];

%%% Save the full output
fname = sprintf('%s/csv/%s/allOutput.csv',outDir,caseName);
writetable(cell2table(full_struct),fname,'writevariablenames',0);

end


%%% =======================================================================
%%% = END
%%% =======================================================================
