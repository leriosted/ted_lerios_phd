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
    figPos1 = [-1900, 100, 1100, 400]; % Work
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


%%
figure('Position',figPos1)
% young
% subplot(2,1,1)
% [h1.young,stats.young] = cdfplot(young.int_Q_dR_Array); hold on
% [h1.old,stats.old] = cdfplot(old.int_Q_dR_Array); hold on
[h1.NFL,stats.NFL] = cdfplot(NFL.int_Q_dR_Array); hold on
[h1.FL,stats.FL] = cdfplot(FL.int_Q_dR_Array); hold on

% set(h1.young, 'LineStyle', '-', 'Color', colorBlue,'LineWidth',pt);
% set(h1.old, 'LineStyle', '-', 'Color', colorGreen,'LineWidth',pt);
set(h1.NFL, 'LineStyle', '--', 'Color', colorBlue,'LineWidth',pt);
set(h1.FL, 'LineStyle', '-', 'Color', colorRed,'LineWidth',pt);

xlabel('$ \int_{}^{} Q \,dR $','Interpreter','latex')
ylabel('Proportion')
xlim([0 20])
title('')
legend('NFL','FL','Location','southeast','Fontsize',FS)
grid on
set(gca,'FontSize',FS);

if print_on ==1
    print('-r400','-dpng','CDF_QdR.png');
else 
end

figure('Position',figPos1)
sz = 40;
% scatter(young.R1insp,young.int_Q_dR_Array,sz,...
%                 'MarkerEdgeColor',colorBlue,...
%                 'MarkerFaceColor',colorBlue,...
%                 'LineWidth',1.5); hold on
% scatter(old.R1insp,old.int_Q_dR_Array,sz,...
%                 'MarkerEdgeColor',colorGreen,...
%                 'MarkerFaceColor',colorGreen,...
%                 'LineWidth',1.5); hold on
scatter(NFL.R1insp,NFL.int_Q_dR_Array,sz,'filled','^',...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(FL.R1insp,FL.int_Q_dR_Array,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
xlabel('R_{1,insp}')
ylabel('$ \int_{}^{} Q \,dR $','Interpreter','latex')
legend('NFL','FL','Location','northwest')
xlim([0 20])
ylim([0 30])
grid on
set(gca,'FontSize',14);

       
if print_on ==1
    print('-r400','-dpng','R1_insp_QdR.png');
else 
end


       