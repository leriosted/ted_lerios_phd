close all; clear all; clc
atWork = 0;
print_on = 1;

colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
colorGreen = '#009A17';
colorOrange = '#FF6600';
pt = 1.0;
FS = 10;

%% Plot Setup:
if atWork == 1
    figPos1 = [-1900, 100, 600, 700]; % Work
    figPos2 = [-1250, 400, 600, 600];
    figPos3 = [-600, 400, 600, 600];
    figPos4 = [-1900, 100, 600, 600];
    figPos5 = [-1250, 100, 600, 600];
    figPos6 = [-600, 100, 600, 300];
else
    figPos1 = [2600, 150, 600, 700]; % Home
    figPos2 = [-1250, 500, 600, 600];
    figPos3 = [-600, 500, 600, 600];
    figPos4 = [-1900, 150, 600, 600];
    figPos5 = [-1250, 150, 600, 600];
    figPos6 = [-600, 150, 600, 300];
end

%%
load young_results.mat
load old_results.mat
load FL_results.mat
load NFL_results.mat

youngIndex = 3;
oldIndex = 4;
NFLIndex = 3;
FLIndex = 2;

tlim = [0 1.5];
qlim = [-1.0 1.0];
plim = [-17 15];
young.index = 2;
old.index = 2;
NFL.index = 2;
FL.index = 2;
young.work = young.work_resistance_dV_Array;
old.work = old.work_resistance_dV_Array;
NFL.work = NFL.work_resistance_dV_Array;
FL.work = FL.work_resistance_dV_Array;

%%
figure('Position',figPos1)
% young
subplot(4,4,1)
plot(young.time(1:young.exp_point),young.pressure(young.index,1:young.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
plot(young.time(1:young.exp_point(young.index)-1),young.Peff_Array{1, young.index},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
xlabel('t (s)')
ylabel('P (cmH_2O)')
leg1 = legend('P_{alv}','P_{eff}','FontSize',12,...
    'Orientation','horizontal');
legend boxoff
leg1.Position(1:2) = [.65 .24];
xlim(tlim)
ylim(plim)
grid on
set(gca,'FontSize',FS);
subplot(4,4,5)
plot(young.time(1:young.exp_point),young.flow(young.index,1:young.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
xlabel('t (s)')
ylabel('Q (L/s)')
xlim(tlim)
ylim(qlim)
grid on
set(gca,'FontSize',FS);
subplot(4,4,9)
plot(young.time(1:young.exp_point),young.volume(young.index,1:young.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
xlabel('t (s)')
ylabel('V (L/s)')
xlim(tlim)
ylim([0 0.8])
grid on
set(gca,'FontSize',FS);

% old
subplot(4,4,2)
plot(old.time(1:young.exp_point),old.pressure(old.index,1:old.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
plot(old.time(1:old.exp_point(old.index)-1),old.Peff_Array{1, old.index},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
% plot(old.time(1:young.exp_point(young.index)-1),-old.R1insp(old.index)*old.flow(1:young.exp_point(young.index)-1),'Color',colorBlack,'Linewidth',pt,'Linestyle','--')
% plot(old.time(old.exp_point(old.index):end),-old.resistance(old.index,old.exp_point(old.index):end).*old.flow(old.index,old.exp_point(old.index):end),'Color',colorRed,'Linewidth',pt,'Linestyle','--')
xlabel('t (s)')
ylabel('P (cmH_2O)')
xlim(tlim)
ylim(plim)
grid on
set(gca,'FontSize',FS);
subplot(4,4,6)
plot(old.time(1:young.exp_point),old.flow(old.index,1:old.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
xlabel('t (s)')
ylabel('Q (L/s)')
xlim(tlim)
ylim(qlim)
grid on
set(gca,'FontSize',FS);
subplot(4,4,10)
plot(old.time(1:young.exp_point),old.volume(old.index,1:old.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
xlabel('t (s)')
ylabel('V (L/s)')
xlim(tlim)
ylim([0 0.8])
grid on
set(gca,'FontSize',FS);


% NFL
subplot(4,4,3)
plot(NFL.time(1:young.exp_point),NFL.pressure(NFL.index,1:NFL.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
plot(NFL.time(1:NFL.exp_point(NFL.index)-1),NFL.Peff_Array{1, NFL.index},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
% plot(NFL.time(1:young.exp_point(young.index)-1),-NFL.R1insp(NFL.index)*NFL.flow(1:young.exp_point(young.index)-1),'Color',colorBlack,'Linewidth',pt,'Linestyle','--')
% plot(NFL.time(NFL.exp_point(NFL.index):end),-NFL.resistance(NFL.index,NFL.exp_point(NFL.index):end).*NFL.flow(NFL.index,NFL.exp_point(NFL.index):end),'Color',colorRed,'Linewidth',pt,'Linestyle','--')
xlabel('t (s)')
ylabel('P_{alv} (cmH_2O)')
xlim(tlim)
ylim(plim)
grid on
set(gca,'FontSize',FS);
subplot(4,4,7)
plot(NFL.time(1:young.exp_point),NFL.flow(NFL.index,1:NFL.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
xlabel('t (s)')
ylabel('Q (L/s)')
xlim(tlim)
ylim(qlim)
grid on
set(gca,'FontSize',FS);
subplot(4,4,11)
plot(NFL.time(1:young.exp_point),NFL.volume(NFL.index,1:NFL.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
xlabel('t (s)')
ylabel('V (L/s)')
xlim(tlim)
ylim([0 0.8])
grid on
set(gca,'FontSize',FS);


% FL
subplot(4,4,4)
plot(FL.time(1:young.exp_point),FL.pressure(FL.index,1:FL.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
plot(FL.time(1:FL.exp_point(FL.index)-1),FL.Peff_Array{1, FL.index},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
% plot(FL.time(1:young.exp_point(young.index)-1),-FL.R1insp(FL.index)*FL.flow(1:young.exp_point(young.index)-1),'Color',colorBlack,'Linewidth',pt,'Linestyle','--')
% plot(FL.time(FL.exp_point(FL.index):end),-FL.resistance(FL.index,FL.exp_point(FL.index):end).*FL.flow(FL.index,FL.exp_point(FL.index):end),'Color',colorRed,'Linewidth',pt,'Linestyle','--')
xlabel('t (s)')
ylabel('P_{alv} (cmH_2O)')
xlim(tlim)
ylim(plim)
grid on
set(gca,'FontSize',FS);
subplot(4,4,8)
plot(FL.time(1:young.exp_point),FL.flow(FLIndex,1:FL.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
xlabel('t (s)')
ylabel('Q (L/s)')
xlim(tlim)
ylim(qlim)
grid on
set(gca,'FontSize',FS);
subplot(4,4,12)
plot(FL.time(1:young.exp_point),FL.volume(FL.index,1:FL.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
xlabel('t (s)')
ylabel('V (L/s)')
xlim(tlim)
ylim([0 0.8])
grid on
set(gca,'FontSize',FS);




titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',12};
annotation('textbox','Position',[0.0 0.1 0.440 0.880],'String','young',titleSettings{:})

titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',12};
annotation('textbox','Position',[0.0 0.1 0.86 0.880],'String','old',titleSettings{:})

titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',12};
annotation('textbox','Position',[0.0 0.1 1.28 0.880],'String','NFL',titleSettings{:})

titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',12};
annotation('textbox','Position',[0.0 0.1 1.700 0.880],'String','FL',titleSettings{:})

% titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',12};
% annotation('textbox','Position',[0.0 0.1 1.230 0.880],'String','NFL',titleSettings{:})
% 
% titleSettings = {'HorizontalAlignment','center','EdgeColor','none','FontSize',12};
% annotation('textbox','Position',[0.0 0.1 1.630 0.880],'String','FL',titleSettings{:})
%        
if print_on ==1
    print('-r400','-dpng','pressuresALL.png');
else 
end


       