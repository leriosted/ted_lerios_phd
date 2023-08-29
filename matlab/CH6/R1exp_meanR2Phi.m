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
    figPos1 = [-1900, 150, 1100, 400]; % Home
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
subplot(1,2,1)
sz = 40;
scatter(young.R1exp,young.meanR2Phi,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1exp,old.meanR2Phi,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1exp,NFL.meanR2Phi,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.R1exp,FL.meanR2Phi,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('R_{1,exp}')
ylabel('mean $R_{2} \Phi$','Interpreter','latex')
legend('young','old','NFL','FL','Location','northwest')
xlim([0 25])
ylim([0 35])
grid on
set(gca,'FontSize',14);

subplot(1,2,2)


sz = 40;
scatter(young.R1exp,young.meanR2Phi./young.R1exp,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1exp,old.meanR2Phi./old.R1exp,sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1exp,NFL.meanR2Phi./NFL.R1exp,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(FL.R1exp,FL.meanR2Phi./FL.R1exp,sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
xlabel('R_{1,exp}')
ylabel('$ \frac{mean R_{2} \Phi}{R_{1,exp}} $','Interpreter','latex')
legend('young','old','NFL','FL','Location','northwest')
xlim([0 25])
ylim([0 4])
grid on
set(gca,'FontSize',14);

if print_on ==1
    print('-r400','-dpng','R1exp_meanR2Phi.png');
else 
end


       