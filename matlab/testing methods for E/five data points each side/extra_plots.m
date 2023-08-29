close all; clear all; clc
atWork = 0;
print_on = 1;

colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
colorGreen = '#009A17';
colorOrange = '#FF6600';
pt = 1.0;
FS = 14;

%% Plot Setup:
if atWork == 1
    figPos1 = [-1900, 100, 500, 400]; % Work
    figPos2 = [-1250, 400, 600, 600];
    figPos3 = [-600, 400, 600, 600];
    figPos4 = [-1900, 100, 600, 600];
    figPos5 = [-1250, 100, 600, 600];
    figPos6 = [-600, 100, 600, 300];
else
    figPos1 = [-1900, 150, 500, 400]; % Home
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



%% R EXP
figure('Position',figPos1)
% young

sz = 40;
scatter(young.R1exp,young.flow_min,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1exp,old.flow_min,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1exp,NFL.flow_min,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.R1exp,FL.flow_min,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('R_{1,exp} [cmH_2O/L/s]')
ylabel('peak Q insp [L/s]')
legend('young','old','NFL','FL','Location','southwest')
xlim([0 25])
ylim([-3 0])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','Rexp_peak_Q_insp.png');
else 
end

figure('Position',figPos1)
% young

sz = 40;
scatter(young.R1exp,young.flow_max,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1exp,old.flow_max,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1exp,NFL.flow_max,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.R1exp,FL.flow_max,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('R_{1,exp} [cmH_2O/L/s]')
ylabel('peak Q exp [L/s]')
legend('young','old','NFL','FL','Location','northwest')
xlim([0 25])
ylim([0 3])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','Rexp_peak_Q_exp.png');
else 
end

figure('Position',figPos1)
% young

sz = 40;
scatter(young.R1exp,young.volume_max,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1exp,old.volume_max,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1exp,NFL.volume_max,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.R1exp,FL.volume_max,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('R_{1,exp} [cmH_2O/L/s]')
ylabel('peak V [L]')
legend('young','old','NFL','FL','Location','northwest')
xlim([0 25])
ylim([0 3])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','Rexp_peak_V.png');
else 
end

%% E Plots

figure('Position',figPos1)
% young

sz = 40;
scatter(young.E,young.flow_min,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.E,old.flow_min,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.E,NFL.flow_min,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.E,FL.flow_min,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('E [cmH_2O/L]')
ylabel('peak Q insp [L/s]')
legend('young','old','NFL','FL','Location','southeast')
xlim([0 15])
ylim([-3 0])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','E_peak_Q_insp.png');
else 
end

figure('Position',figPos1)
% young

sz = 40;
scatter(young.E,young.flow_max,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.E,old.flow_max,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.E,NFL.flow_max,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.E,FL.flow_max,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('E [cmH_2O/L]')
ylabel('peak Q exp [L/s]')
legend('young','old','NFL','FL','Location','northeast')
xlim([0 15])
ylim([0 3])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','E_peak_Q_exp.png');
else 
end

figure('Position',figPos1)
% young

sz = 40;
scatter(young.E,young.volume_max,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.E,old.volume_max,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.E,NFL.volume_max,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.E,FL.volume_max,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('E [cmH_2O/L]')
ylabel('peak V [L]')
legend('young','old','NFL','FL','Location','northeast')
xlim([0 15])
ylim([0 3])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','E_peak_V.png');
else 
end

%% R INSP
figure('Position',figPos1)
% young

sz = 40;
scatter(young.R1insp,young.flow_min,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1insp,old.flow_min,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1insp,NFL.flow_min,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.R1insp,FL.flow_min,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('R_{1,insp} [cmH_2O/L/s]')
ylabel('peak Q insp [L/s]')
legend('young','old','NFL','FL','Location','southeast')
xlim([0 15])
ylim([-3 0])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','Rinsp_peak_Q_insp.png');
else 
end

figure('Position',figPos1)
% young

sz = 40;
scatter(young.R1insp,young.flow_max,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1insp,old.flow_max,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1insp,NFL.flow_max,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.R1insp,FL.flow_max,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('R_{1,insp} [cmH_2O/L/s]')
ylabel('peak Q exp [L/s]')
legend('young','old','NFL','FL','Location','northeast')
xlim([0 15])
ylim([0 3])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','Rinsp_peak_Q_exp.png');
else 
end

figure('Position',figPos1)
% young

sz = 40;
scatter(young.R1insp,young.volume_max,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1insp,old.volume_max,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1insp,NFL.volume_max,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.R1insp,FL.volume_max,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('R_{1,insp} [cmH_2O/L/s]')
ylabel('peak V [L]')
legend('young','old','NFL','FL','Location','northeast')
xlim([0 15])
ylim([0 3])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','Rinsp_peak_V.png');
else 
end
       