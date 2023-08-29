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
    figPos1 = [2580, -190, 1000, 1270]; % Home
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

tlim = [0 4];
qlim = [-1.2 0];
plim = [-15 15];
rlim = [0 60];

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
% FL
subplot(8,5,1)
plot(NFL.flow(1,1:end),NFL.resistance(1,1:end),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
% plot(FL.time(1:FL.exp_point(1)-1),FL.Peff_Array{1, 1},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
% plot(FL.time,FL.Peff_Array{1, 1},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
xlabel('Q (L/s)')
ylabel('R (cmH_2O/L/s)')
% leg1 = legend('P_{alv}','P_{eff}','FontSize',12,...
%     'Orientation','horizontal');
% legend boxoff
% leg1.Position(1:2) = [.75 .16];
ylim(rlim)
xlim(qlim)
grid on
set(gca,'FontSize',FS);

for ii = 2:25
    subplot(8,5,ii)
    plot(NFL.flow(ii,1:end),NFL.resistance(ii,1:end),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
%     plot(FL.time(1:FL.exp_point(ii)-1),FL.Peff_Array{1, ii},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
%     plot(FL.time,FL.Peff_Array{1, ii},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
    xlabel('Q (L/s)')
    ylabel('R (cmH_2O/L/s)')
    ylim(rlim)
    xlim(qlim)
    grid on
    set(gca,'FontSize',FS);
end
  
if print_on ==1
    print('-r400','-dpng','R-Q_NFL.png');
else 
end


       