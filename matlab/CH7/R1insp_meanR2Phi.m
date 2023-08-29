close all; clear all; clc
atWork = 0;
print_on = 0;

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

%% Linear regression line
% dlm = fitlm(NFL.R1insp,NFL.meanR2Phi,'Intercept', false)

%%
figure('Position',figPos1);

sz = 40;
scatter(NFL.R1insp,NFL.meanR2Phi,sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
h1 = lsline();
h1.Color = 'k';
h1.LineWidth = 1.5;
scatter(FL.R1insp,FL.meanR2Phi,sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
xlabel('R_{1,insp}')
ylabel('mean $R_{2} \Phi$','Interpreter','latex')
legend('NFL','','FL','Location','northwest')
xlim([0 20])
ylim([0 35])
grid on
set(gca,'FontSize',14);


if print_on ==1
    print('-r400','-dpng','R1insp_meanR2Phi.png');
else 
end


       