%%% =======================================================================
%%% = plotAllObs.m
%%% = Alex Turner
%%% = 06/02/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Plot the raw observations and the box-model output.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): obs  -- Structure containing the observations.
%%% =  ( 2): type -- String indicating the type of observations.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = plotResults_alt( params, out, obs, add_name, OH, CO)

%%% Close everything
close all;

%%% Build our plotting structure
St                = params.pSt;
xLims             = [St(1) + double(params.spinup)/1000,St(end)];
dat.t             = St;
dat.ch4_NH        = out.nh_ch4;
dat.ch4           = out.sh_ch4;
dat.dch4          = out.ch4(2:end) - out.ch4(1:end-1);
dat.oh            = out.oh;
dat.co            = out.co;
dat.d13c          = out.d13c;
dat.dD            = out.dD;
dat.d14c          = out.d14c;
dat.cl            = out.cl;
dat.ems_wet       = out.ch4_ems_wet_tropical + out.ch4_ems_wet_boreal;
dat.ems_fire      = out.ch4_ems_fire;
dat.ems_fossil    = out.ch4_ems_fossil;
dat.ems_animal    = out.ch4_ems_animal;
dat.ems_oh        = out.oh_ems;
dat.loss_oh       = out.oh_loss;
dat.loss_oh_ch4   = out.oh_loss_ch4;
dat.loss_oh_co    = out.oh_loss_co;
dat.loss_oh_other = out.oh_loss_other;
dat.lifetime      = out.ch4_lifetime;
dat.tau_TS        = out.tau_TS;
dat.tau_ST        = out.tau_TS / 5.7047;
dat.temp          = params.sTemp;

%%% Get the ice core data
% CH4
iceCoreDat = readtable('./data/obs/iceCore/41586_2008_BFnature06950_MOESM33_ESM.xls');
iceCH4.t   = (1950 - table2array(iceCoreDat(:,2)))/1000;
iceCH4.ch4 = table2array(iceCoreDat(:,3));
% Add the Yan data
if St(1) < -2000
    iceCoreDat   = readtable('./data/obs/iceCore/AllanHillsCH4.csv');
    yanDat.tim   = (1950/1000 - table2array(iceCoreDat(:,5)));
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
    % Combine the records
    iceCH4.t   = [iceCH4.t;yanDat.tim];
    iceCH4.ch4 = [iceCH4.ch4;yanDat.ch4];
end

%%% Temperature
fspec     = '%f %f %f %f %f %f';
fileID    = fopen('./data/obs/iceCore/edc3deuttemp2007.txt','r');
fdat      = textscan(fileID,fspec,'HeaderLines',92);
iceT.t    = (1950 - fdat{3}(:))/1000;
iceT.temp = fdat{5}(:);
fclose(fileID);
% % Sync
% iceCore.t    = dat.t;
% iceCore.temp = interp1(iceT.t,iceT.temp,iceCore.t,'linear','extrap');
% iceCore.ch4  = interp1(iceCH4.t,iceCH4.ch4,iceCore.t,'linear','extrap');

%%% Read in the Bock et al isotope data (d13C)
% EDC d13C
fspec    = '%f %f %f %f %f %f %f %f %f %f';
fileID   = fopen('./data/obs/iceCore/EDC_d13CH4.tab','r');
fdat     = textscan(fileID,fspec,'HeaderLines',25);
edc.t    = (1950/1000-fdat{2}(:));
edc.ch4  = fdat{3}(:);
edc.d13c = fdat{10}(:);
edc.sig  = fdat{6}(:);
fclose(fileID);
% TALDICE d13C
fspec        = '%f %f %f %f %f %f %f %f %f %f';
fileID       = fopen('./data/obs/iceCore/TALDICE_d13CH4.tab','r');
fdat         = textscan(fileID,fspec,'HeaderLines',25);
taldice.t    = (1950/1000-fdat{2}(:));
taldice.ch4  = fdat{3}(:);
taldice.d13c = fdat{10}(:);
taldice.sig  = fdat{6}(:);
fclose(fileID);
% VOSTOK d13C
fspec       = '%f %f %f %f %f %f %f %f %f %f';
fileID      = fopen('./data/obs/iceCore/Vostok_d13CH4.tab','r');
fdat        = textscan(fileID,fspec,'HeaderLines',26);
vostok.t    = (1950/1000-fdat{2}(:));
vostok.ch4  = fdat{3}(:);
vostok.d13c = fdat{10}(:);
vostok.sig  = fdat{6}(:);
fclose(fileID);
% Combine the records
d13c.t    = [edc.t;taldice.t;vostok.t];
d13c.ch4  = [edc.ch4;taldice.ch4;vostok.ch4];
d13c.d13c = [edc.d13c;taldice.d13c;vostok.d13c];
d13c.sig  = [edc.sig;taldice.sig;vostok.sig];
[~,ind]   = sort(d13c.t,'descend');
d13c.t    = d13c.t(ind);
d13c.ch4  = d13c.ch4(ind);
d13c.d13c = d13c.d13c(ind);
d13c.sig  = d13c.sig(ind);

%%% Read in the Bock et al isotope data (dD)
% EDC dD
fspec   = '%f %f %f %f %f %f %f %f %f %f';
fileID  = fopen('./data/obs/iceCore/EDC_dDCH4.tab','r');
fdat    = textscan(fileID,fspec,'HeaderLines',19);
edc.t   = (1950/1000-fdat{2}(:));
edc.dD  = fdat{6}(:);
edc.sig = fdat{4}(:);
fclose(fileID);
% EDML dD
fspec    = '%f %f %f %f %f %f %f %f %f %f';
fileID   = fopen('./data/obs/iceCore/EDML_dDCH4.tab','r');
fdat     = textscan(fileID,fspec,'HeaderLines',19);
edml.t   = (1950/1000-fdat{2}(:));
edml.dD  = fdat{6}(:);
edml.sig = fdat{4}(:);
fclose(fileID);
% Combine the records
dD.t    = [edc.t;edml.t];
dD.dD   = [edc.dD;edml.dD];
dD.sig  = [edc.sig;edml.sig];
[~,ind] = sort(dD.t,'descend');
dD.t    = dD.t(ind);
dD.dD   = dD.dD(ind);
dD.sig  = dD.sig(ind);

%%% 14CH4 data
% Taylor glacier data from Dynosius et al., Science (2020): TG_deglacial14C_finaldata.xlsx
taylor.t = (1950/1000 - [
    14.92, 14.86, 14.58, 14.54, 14.42, 14.42, 13.00, 10.13, 10.13,  9.21,  7.94]');
taylor.d14c = [
      317,   327,   288,   287,   216,   204,   246,   126,   139,   141,    92]';
taylor.sig = [
      166,   151,   128,   112,   109,   111,    90,    58,    57,    63,    74]';
% Taylor glacier data from Petrenko et al., Science (2017): Petrenko_TG_YD-PB_14CH4.xlsx
petrenko.t = (1950 - [
    11715, 11559, 11515, 11453, 11357]')/1000;
petrenko.d14c = [
    192.2, 130.6, 125.8, 150.2, 132.3]';
petrenko.sig = [
     52.4,  38.1,  35.3,  40.5,  35.1]';
% Combine the records
d14c.t    = [taylor.t;petrenko.t];
d14c.d14c = [taylor.d14c;petrenko.d14c];
d14c.sig  = [taylor.sig;petrenko.sig];
[~,ind]   = sort(d14c.t,'descend');
d14c.t    = d14c.t(ind);
d14c.d14c = d14c.d14c(ind);
d14c.sig  = d14c.sig(ind);


%%
%%% =======================================================================
%%% 3. Plotting Methane
%%% =======================================================================

%%% Colors to use
cols   = lines(5);
modCol = cols(1,:);
emsCol = cols([4,2,5,3],:);
obsCol = [0,0,0];

%%% Make the methane uncertainty patch
% Methane
unq_ch4 = zeros(size(obs.ch4));
counter = 0;
for i = 1:length(unq_ch4)
    if i == 1 && ~isnan(obs.ch4(i))
        counter    = counter + 1;
        unq_ch4(i) = counter;
    else
        if ~isnan(obs.ch4(i))
            if unq_ch4(i-1) <= 0
                 counter = counter + 1;
            end
            unq_ch4(i) = counter;
        end
    end
end
unq_ch4(unq_ch4 == 0) = NaN;
% d13C
unq_d13C = zeros(size(obs.d13c));
counter  = 0;
for i = 1:length(unq_d13C)
    if i == 1 && ~isnan(obs.d13c(i))
        counter     = counter + 1;
        unq_d13C(i) = counter;
    else
        if ~isnan(obs.d13c(i))
            if unq_d13C(i-1) <= 0
                 counter = counter + 1;
            end
            unq_d13C(i) = counter;
        end
    end
end
unq_d13C(unq_d13C == 0) = NaN;
% dD
unq_dD = zeros(size(obs.dD));
counter  = 0;
for i = 1:length(unq_dD)
    if i == 1 && ~isnan(obs.dD(i))
        counter     = counter + 1;
        unq_d13C(i) = counter;
    else
        if ~isnan(obs.dD(i))
            if unq_dD(i-1) <= 0
                 counter = counter + 1;
            end
            unq_dD(i) = counter;
        end
    end
end
unq_dD(unq_dD == 0) = NaN;
% d14C
unq_d14C = zeros(size(obs.d14c));
counter  = 0;
for i = 1:length(unq_d14C)
    if i == 1 && ~isnan(obs.d14c(i))
        counter     = counter + 1;
        unq_d14C(i) = counter;
    else
        if ~isnan(obs.d14c(i))
            if unq_d14C(i-1) <= 0
                 counter = counter + 1;
            end
            unq_d14C(i) = counter;
        end
    end
end
unq_d14C(unq_d14C == 0) = NaN;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% TEST
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dat.d14c_N = out.nh_d14c;
% dat.d14c_S = out.sh_d14c;
% dat.d14c_D = out.nh_d14c - out.sh_d14c;
% 
% % Plot
% f0 = figure(1);pos = get(f0,'position');set(gcf,'color','w');
% set(f0,'position',[0 pos(2) pos(3)*1.7 pos(4)*2.7])
% h1 = subplot(2,1,1);
% hold on; box on;
% set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
% p1 = plot(dat.t,dat.d14c_N,dat.t,dat.d14c_S);
% set(p1,'LineWidth',3);
% ylabel('\Delta^{14}C of CH_4 (‰)');
% xlim([-50,0])
% xlabel('Time (kyr BP)')
% h2 = subplot(2,1,2);
% hold on; box on;
% set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
% plot(dat.t,dat.d14c_D,'-','Color',obsCol,'LineWidth',3);
% ylabel('Difference (NH - SH; ‰)');
% xlim([-50,0])
% xlabel('Time (kyr BP)')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot the results
% Initialize
f1 = figure(1);pos = get(f1,'position');set(gcf,'color','w');  % changed set gcf to set f1
set(f1,'position',[0 pos(2) pos(3)*1.7 pos(4)*2.7]);
% Methane

h1 = subplot(5,1,1);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
for i = 1:nanmax(unq_ch4)
    iUse   = unq_ch4 == i;
    tPatch = [St(iUse);flipud(St(iUse))];
    uPatch = [obs.ch4(iUse) - obs.ch4_err(iUse);flipud(obs.ch4(iUse) + obs.ch4_err(iUse))];
    patch(tPatch,uPatch,(1-obsCol)*.6+obsCol,'EdgeColor','none');
end
plot(dat.t,dat.ch4,'-','Color',modCol,'LineWidth',3);
%plot(dat.t,dat.ch4_NH,'-','Color','green','LineWidth',3);  % greenland -- should plot simulated IPD
plot(St,obs.ch4,'-','Color',obsCol);
%plot(St,obs.ch4_NH,'-','Color','cyan');  % greenland

plot(iceCH4.t,iceCH4.ch4,'.','Color',obsCol);
set(gca,'XTickLabel',{},'YGrid','on');
ylabel('CH_4 (ppb)');YYCol = get(gca,'YColor');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol);

% d13C
h2 = subplot(5,1,2);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on');
for i = 1:nanmax(unq_d13C)
    iUse   = unq_d13C == i;
    tPatch = [St(iUse);flipud(St(iUse))];
    uPatch = [obs.d13c(iUse) - obs.d13c_err(iUse);flipud(obs.d13c(iUse) + obs.d13c_err(iUse))];
    patch(tPatch,uPatch,(1-obsCol)*.6+obsCol,'EdgeColor','none');
end
plot(dat.t,dat.d13c,'-','Color',modCol,'LineWidth',3);
plot(d13c.t,d13c.d13c,'.','Color',obsCol);
for i = 1:length(d13c.t)
    plot(d13c.t(i)*[1,1],d13c.d13c(i) + d13c.sig(i)*[-1,1],'k-');
end
set(gca,'XTickLabel',{},'YGrid','on','YAxisLocation','right')
%axis ij;
ylabel('\delta^{13}C of CH_4 (‰)');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol);

% dD
h3 = subplot(5,1,3);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
for i = 1:nanmax(unq_dD)
    iUse   = unq_dD == i;
    tPatch = [St(iUse);flipud(St(iUse))];
    uPatch = [obs.dD(iUse) - obs.dD_err(iUse);flipud(obs.dD(iUse) + obs.dD_err(iUse))];
    patch(tPatch,uPatch,(1-obsCol)*.6+obsCol,'EdgeColor','none');
end
plot(dat.t,dat.dD,'-','Color',modCol,'LineWidth',3);
plot(dD.t,dD.dD,'.','Color',obsCol);
for i = 1:length(dD.t)
    plot(dD.t(i)*[1,1],dD.dD(i) + dD.sig(i)*[-1,1],'k-');
end
set(gca,'XTickLabel',{},'YGrid','on')
ylabel('\delta{}D of CH_4 (‰)');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol);

% d14C
h4 = subplot(5,1,4);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
for i = 1:nanmax(unq_d14C)
    iUse   = unq_d14C == i;
    tPatch = [St(iUse);flipud(St(iUse))];
    uPatch = [obs.d14c(iUse) - obs.d14c_err(iUse);flipud(obs.d14c(iUse) + obs.d14c_err(iUse))];
    patch(tPatch,uPatch,(1-obsCol)*.6+obsCol,'EdgeColor','none');
end
plot(dat.t,dat.d14c,'-','Color',modCol,'LineWidth',3);
plot(d14c.t,d14c.d14c,'.','Color',obsCol);
for i = 1:length(d14c.t)
    plot(d14c.t(i)*[1,1],d14c.d14c(i) + d14c.sig(i)*[-1,1],'k-');
end
set(gca,'XTickLabel',{},'YGrid','on','YAxisLocation','right')
ylabel('\Delta^{14}C of CH_4 (‰)');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol);

% Emissions
h5 = subplot(5,1,5);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
plot(dat.t,dat.ems_animal,'-','Color',emsCol(4,:),'LineWidth',3);
plot(dat.t,dat.ems_fossil,'-','Color',emsCol(1,:),'LineWidth',3);
plot(dat.t,dat.ems_fire,'-','Color',emsCol(2,:),'LineWidth',3);
plot(dat.t,dat.ems_wet,'-','Color',emsCol(3,:),'LineWidth',3);
set(gca,'YGrid','on','YAxisLocation','right');
ylabel('CH_4 emissions (Tg/yr)');
xlabel('Time (kyr BP)');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol);
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol);

% Reposition
set(h1,'Position',[.13,.82,.7750,.16]);
set(h2,'Position',[.13,.63,.7750,.16]);
set(h3,'Position',[.13,.44,.7750,.16]);
set(h4,'Position',[.13,.25,.7750,.16]);
set(h5,'Position',[.13,.06,.7750,.16]);

% Save
export_fig(f1,['./output/' add_name 'ch4_simulated.png'],'-png','-m2','-painters','-cmyk');
hold off;

% quick check
% f_NHSH = figure(1);
% clf(f_NHSH);
% hold on;
% plot(out.ch4, 'DisplayName', 'ch4');
% plot(out.ch4_ems, 'DisplayName', 'ch4\_ems');
% plot(out.ch4_lifetime, 'DisplayName', 'ch4\_lifetime');
% plot(out.ch4_lifetime_nh, 'DisplayName', '...\_nh');
% plot(out.ch4_lifetime_sh, 'DisplayName', '...\_sh');
% hold off;
% legend('show');
% title('comparing simulated CH4 metrics');
% export_fig(f_NHSH, ['./output/' add_name 'NHSH.png'], '-png', '-r300');

% plot IPD
% Create figure with exact size
f_IPD = figure;
set(f_IPD, 'Units', 'inches', 'Position', [1, 1, 4, 6]);
clf(f_IPD);
inds = ~isnan(obs.ch4_NH); % where we have Greenland data
IPD_sim = (dat.ch4_NH(inds) - dat.ch4(inds))./(0.5*(dat.ch4(inds)+dat.ch4_NH(inds)));
IPD_obs = (obs.ch4_NH(inds) - obs.ch4(inds))./(0.5*(obs.ch4(inds)+obs.ch4_NH(inds)));
ax1 = subplot(2,1,1);
plot(dat.t(inds), IPD_sim, '-', 'Color', modCol, 'LineWidth', 3); 
hold on; box on;
% plot(St(inds), IPD_obs, '-', 'Color', obsCol);
% plot(St(inds), IPD_obs, '.', 'Color', obsCol);
ylabel('rIPD (following Brook)');
hold off;
ax2 = subplot(2,1,2);
plot(dat.t(inds), dat.ch4_NH(inds), '-', 'Color', 'green', 'LineWidth', 1); hold on; box on;
plot(St(inds), obs.ch4_NH(inds), '.', 'Color', 'green', 'LineWidth', 1);
plot(dat.t(inds), dat.ch4(inds), '-', 'Color', 'blue', 'LineWidth', 1);
plot(St(inds), obs.ch4(inds), '.', 'Color', 'blue', 'LineWidth', 1);

xlabel('ka');
ylabel('CH4 (ppb)');
hold off;
%set([ax1, ax2], 'FontSize', 10);  % Optional: set consistent font size
export_fig(f_IPD, ['./output/' add_name 'rIPD.png'], '-png', '-r300');


if OH     % plot OH concentration
    f_oh = figure(1);
    set(f_oh, 'Units', 'inches', 'Position', [1, 1, 8, 6]);
    % First subplot: OH concentration
    subplot(3,1,1);
    plot(out.oh);
    title('OH Concentration');
    ylabel('OH');
    grid on;
    % Second subplot: OH emissions
    subplot(3,1,2);
    plot(out.oh_ems);
    title('OH Emissions');
    ylabel('OH Emissions');
    grid on;
    % Third subplot: OH loss
    subplot(3,1,3);
    plot(out.oh_loss);
    title('OH Loss');
    ylabel('OH Loss');
    xlabel('Time');
    grid on;
    export_fig(f_oh, ['./output/' add_name 'oh_simulated.png'], '-png', '-r300');
end

if CO     % plot OH concentration
    f_co = figure(1);
    set(f_co, 'Units', 'inches', 'Position', [1, 1, 8, 4]);
    % First subplot: OH concentration
    subplot(2,1,1);
    plot(out.co);
    title('CO Concentration');
    ylabel('CO');
    grid on;

    subplot(2,1,2);
    plot(out.co_ems);
    title('CO Emissions');
    ylabel('CO Emissions');
    grid on;

    export_fig(f_co, ['./output/' add_name 'co_simulated.png'], '-png');
end




%%
%%% =======================================================================
%%% 3. Plotting Drivers
%%% =======================================================================

iceCol = [119,170,221]/256;

%%% Plot the results
% Initialize
f2 = figure(2);pos = get(f2,'position');set(gcf,'color','w');
set(f2,'position',[pos(3)*1.7 pos(2) pos(3)*1.7 pos(4)*2.7])
% Emissions
h1 = subplot(5,1,1);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
plot(dat.t,dat.ems_animal,'-','Color',emsCol(4,:),'LineWidth',3);
plot(dat.t,dat.ems_fossil,'-','Color',emsCol(1,:),'LineWidth',3);
plot(dat.t,dat.ems_fire,'-','Color',emsCol(2,:),'LineWidth',3);
plot(dat.t,dat.ems_wet,'-','Color',emsCol(3,:),'LineWidth',3);
set(gca,'XTickLabel',{},'YGrid','on')
ylabel('CH_4 emissions (Tg/yr)');
xlim(xLims);xLL = xlim;yLL = ylim;ylim(yLL);
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% OH
h2 = subplot(5,1,2);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
plot(dat.t,dat.oh*1e-5,'-','Color',modCol,'LineWidth',3);
set(gca,'XTickLabel',{},'YGrid','on','YAxisLocation','right')
ylabel('OH (10^5 molec/cm^3)');
xlim(xLims);xLL = xlim;yLL = ylim;ylim(yLL);
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% Production/Loss
h3 = subplot(5,1,3);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
tPatch = [dat.t;flipud(dat.t)];
patch(tPatch,[dat.ems_oh;flipud(zeros(size(dat.t)))],modCol,'EdgeColor','none')
patch(tPatch,[-1*(dat.loss_oh_ch4 + dat.loss_oh_co + dat.loss_oh_other);flipud(zeros(size(dat.t)))],emsCol(3,:),'EdgeColor','none')
patch(tPatch,[-1*(dat.loss_oh_ch4 + dat.loss_oh_co);flipud(zeros(size(dat.t)))],emsCol(2,:),'EdgeColor','none')
patch(tPatch,[-1*(dat.loss_oh_ch4);flipud(zeros(size(dat.t)))],emsCol(1,:),'EdgeColor','none')
plot(dat.t,dat.ems_oh,'-','Color',modCol,'LineWidth',3);
plot(dat.t,-1*(dat.loss_oh_ch4),'-','Color',emsCol(1,:),'LineWidth',3);
plot(dat.t,-1*(dat.loss_oh_ch4 + dat.loss_oh_co),'-','Color',emsCol(2,:),'LineWidth',3);
plot(dat.t,-1*(dat.loss_oh_ch4 + dat.loss_oh_co + dat.loss_oh_other),'-','Color',emsCol(3,:),'LineWidth',3);
plot(dat.t,zeros(size(dat.t)),'k-','LineWidth',2,'LineWidth',2,'Color',YYCol);
set(gca,'XTickLabel',{},'YGrid','on')
ylabel('OH Source (Tg/yr)');
%ylim([-1500,1500]);
xlim(xLims);xLL = xlim;yLL = ylim;yLL = [-max(abs(yLL)),max(abs(yLL))];ylim(yLL);
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% Chlorine
h4 = subplot(5,1,4);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
plot(dat.t,dat.cl,'-','Color',modCol,'LineWidth',3);
set(gca,'XTickLabel',{},'YGrid','on','YAxisLocation','right')
ylabel('Cl (molec/cm^3)');
xlim(xLims);xLL = xlim;yLL = ylim;ylim(yLL);
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% Strat-trop exchange
h5 = subplot(5,1,5);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
yyaxis right
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on','YColor',iceCol)
%plot(dat.t,params.iceVol,'-','LineWidth',3,'Color',iceCol);
patch([dat.t;flipud(dat.t)],[params.iceVol;flipud(zeros(size(params.iceVol)))],iceCol,'EdgeColor','none');
xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',iceCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
ylabel('Ice volume (m)');
yyaxis left
plot(dat.t,dat.tau_ST,'-','Color',modCol,'LineWidth',3);
set(gca,'YGrid','on','YColor',YYCol)
ylabel('Strat-trop exchange (yr)');
xlabel('Time (kyr BP)')
%ylim([1.35,1.45])
xlim(xLims);
xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',iceCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% Reposition
set(h1,'Position',[.13,.82,.7750,.16])
set(h2,'Position',[.13,.63,.7750,.16])
set(h3,'Position',[.13,.44,.7750,.16])
set(h4,'Position',[.13,.25,.7750,.16])
set(h5,'Position',[.13,.06,.7750,.16])
% Save
exportgraphics(f2,['./output/' add_name 'ch4_drivers.png']) %,'-png','-m2','-painters','-cmyk')


%%
%%% =======================================================================
%%% 4. Last Glacial Maximum
%%% =======================================================================

%%% Compute the residuals in methane
iceCH4.ch4_resid = iceCH4.ch4;
for i = 1:length(iceCH4.ch4_resid)
    iceCH4.ch4_resid(i) = iceCH4.ch4(i) - interp1(dat.t,dat.ch4,iceCH4.t(i));
end
xLims = [-150,0];Fsize = 5;modCol = [0,0,0];

%%% Plot the results
% Initialize
f3 = figure(3);pos = get(f3,'position');set(gcf,'color','w');
set(f3,'position',[pos(3)*3.4 pos(2) pos(3)*1.3 pos(4)*2.7])
% Methane
h1 = subplot(5,1,1);FaceCol = [51,187,238]/256;
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
for i = 1:nanmax(unq_ch4)
    iUse   = unq_ch4 == i;
    tPatch = [St(iUse);flipud(St(iUse))];
    uPatch = [obs.ch4(iUse) - obs.ch4_err(iUse);flipud(obs.ch4(iUse) + obs.ch4_err(iUse))];
    patch(tPatch,uPatch,(1-FaceCol)*.6+FaceCol,'EdgeColor','none')
end
for i = 1:length(iceCH4.t)
    plot(iceCH4.t(i)*[1,1],iceCH4.ch4(i) + mean(obs.ch4_err)*[-1,1],'k-','LineWidth',1,'Color',FaceCol)
end
plot(dat.t,dat.ch4,'-','Color',modCol,'LineWidth',3);
plot(St,obs.ch4,'-','Color',obsCol);
plot(iceCH4.t,iceCH4.ch4,'o','Color',obsCol,'MarkerFaceColor',FaceCol,'MarkerEdgeColor',obsCol,'MarkerSize',Fsize);
set(gca,'XTickLabel',{},'YGrid','on')
ylabel('CH_4 (ppb)');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% residuals
h2 = subplot(5,1,2);
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
for i = 1:length(iceCH4.t)
    plot(iceCH4.t(i)*[1,1],iceCH4.ch4_resid(i) + mean(obs.ch4_err)*[-1,1],'k-','LineWidth',1,'Color',FaceCol)
end
plot(iceCH4.t,iceCH4.ch4_resid,'o','Color',obsCol,'MarkerFaceColor',FaceCol,'MarkerEdgeColor',obsCol,'MarkerSize',Fsize);
set(gca,'XTickLabel',{},'YGrid','on','YAxisLocation','right')
%axis ij;
ylabel('\delta^{13}C of CH_4 (‰)');
ylabel('CH_4 residuals (ppb)');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% d13C
h3 = subplot(5,1,3);FaceCol = [238,119,51]/256;
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
for i = 1:nanmax(unq_d13C)
    iUse   = unq_d13C == i;
    tPatch = [St(iUse);flipud(St(iUse))];
    uPatch = [obs.d13c(iUse) - obs.d13c_err(iUse);flipud(obs.d13c(iUse) + obs.d13c_err(iUse))];
    patch(tPatch,uPatch,(1-FaceCol)*.6+FaceCol,'EdgeColor','none')
end
for i = 1:length(d13c.t)
    plot(d13c.t(i)*[1,1],d13c.d13c(i) + d13c.sig(i)*[-1,1],'k-','LineWidth',1,'Color',FaceCol)
end
plot(dat.t,dat.d13c,'-','Color',modCol,'LineWidth',3);
plot(d13c.t,d13c.d13c,'o','Color',obsCol,'MarkerFaceColor',FaceCol,'MarkerEdgeColor',obsCol,'MarkerSize',Fsize);
set(gca,'XTickLabel',{},'YGrid','on')
ylabel('\delta^{13}C of CH_4 (‰)');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% dD
h4 = subplot(5,1,4);FaceCol = [238,51,119]/256;
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
for i = 1:nanmax(unq_dD)
    iUse   = unq_dD == i;
    tPatch = [St(iUse);flipud(St(iUse))];
    uPatch = [obs.dD(iUse) - obs.dD_err(iUse);flipud(obs.dD(iUse) + obs.dD_err(iUse))];
    patch(tPatch,uPatch,(1-FaceCol)*.6+FaceCol,'EdgeColor','none')
end
for i = 1:length(dD.t)
    plot(dD.t(i)*[1,1],dD.dD(i) + dD.sig(i)*[-1,1],'k-','LineWidth',1,'Color',FaceCol)
end
plot(dat.t,dat.dD,'-','Color',modCol,'LineWidth',3);
plot(dD.t,dD.dD,'o','Color',obsCol,'MarkerFaceColor',FaceCol,'MarkerEdgeColor',obsCol,'MarkerSize',Fsize);
set(gca,'XTickLabel',{},'YGrid','on','YAxisLocation','right')
ylabel('\delta{}D of CH_4 (‰)');
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% d14C
h5 = subplot(5,1,5);FaceCol = [34,136,51]/256;
hold on; box on;
set(gca,'YDir','normal','FontName','Helvetica','FontSize',16,'TickDir','out','LineWidth',2,'XMinorTick','on','YMinorTick','on')
for i = 1:nanmax(unq_d14C)
    iUse   = unq_d14C == i;
    tPatch = [St(iUse);flipud(St(iUse))];
    uPatch = [obs.d14c(iUse) - obs.d14c_err(iUse);flipud(obs.d14c(iUse) + obs.d14c_err(iUse))];
    patch(tPatch,uPatch,(1-FaceCol)*.6+FaceCol,'EdgeColor','none')
end
for i = 1:length(d14c.t)
    plot(d14c.t(i)*[1,1],d14c.d14c(i) + d14c.sig(i)*[-1,1],'k-','LineWidth',1,'Color',FaceCol)
end
plot(dat.t,dat.d14c,'-','Color',modCol,'LineWidth',3);
plot(d14c.t,d14c.d14c,'o','Color',obsCol,'MarkerFaceColor',FaceCol,'MarkerEdgeColor',obsCol,'MarkerSize',Fsize);
set(gca,'YGrid','on','YAxisLocation','right')
ylabel('\Delta^{14}C of CH_4 (‰)');
xlabel('Time (kyr BP)')
xlim(xLims);xLL = xlim;yLL = ylim;
plot([xLL(1),xLL(2)],[yLL(1),yLL(1)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(2)],[yLL(1),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(2),xLL(1)],[yLL(2),yLL(2)],'-','LineWidth',2,'Color',YYCol)
plot([xLL(1),xLL(1)],[yLL(2),yLL(1)],'-','LineWidth',2,'Color',YYCol)
% Reposition
set(h1,'Position',[.11,.82,.7750,.16])
set(h2,'Position',[.11,.63,.7750,.16])
set(h3,'Position',[.11,.44,.7750,.16])
set(h4,'Position',[.11,.25,.7750,.16])
set(h5,'Position',[.11,.06,.7750,.16])
% Save
exportgraphics(f3,['./output/' add_name 'ch4_lgm.png']) %,'-png','-m2','-painters','-cmyk')


end


%%% =======================================================================
%%% = END
%%% =======================================================================
