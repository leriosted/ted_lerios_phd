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
for ii = 1:20
    young.WOB(ii,1) = abs(young.WOB_dV{1, ii}(end,1));
end
for ii = 1:20
    old.WOB(ii,1) = abs(old.WOB_dV{1, ii}(end,1));
end
for ii = 1:25
    NFL.WOB(ii,1) = abs(NFL.WOB_dV{1, ii}(end,1));
end
for ii = 1:35
    FL.WOB(ii,1) = abs(FL.WOB_dV{1, ii}(end,1));
end
%% Linear regression line
% dlm = fitlm(NFL.R1insp,NFL.meanR2Phi,'Intercept', false)

%%
figure('Position',figPos1);

sz = 40;
scatter(young.R1insp,young.WOB,sz,'^',...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1insp,old.WOB,sz,'v',...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1insp,NFL.WOB,sz,...
                'MarkerEdgeColor',colorOrange,...
                'MarkerFaceColor',colorOrange,...
                'LineWidth',1.5); hold on
scatter(FL.R1insp,FL.WOB,sz,'d',...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on

totalDataX = [young.R1insp old.R1insp NFL.R1insp FL.R1insp]';
totalDataY = [young.WOB; old.WOB; NFL.WOB; FL.WOB];
[p,S]=polyfit(totalDataX,totalDataY,1);
[yfit,delta]=polyval(p,totalDataX,S);
plot(totalDataX,yfit,'k');
RSS = sum((totalDataY-yfit).^2);
TSS = sum((totalDataY-mean(totalDataY)).^2);
Rsquared = 1-RSS/TSS;
anno = annotation('textbox', [0.70, 0.6, 0.1, 0.1], 'String', "R^2 = " + round(Rsquared,2));
anno.EdgeColor = 'w';
anno.BackgroundColor = 'w';
anno.FontSize = 12;
% h1.Color = 'k';
% h1.LineWidth = 1.5;
xlabel('R_{1,insp}')
ylabel('$ \int_{}^{} P_{eff} \,dV $','Interpreter','latex')
legend('young','old','NFL','FL','Location','northwest')
% xlim([0 20])
% ylim([0 35])
grid on
set(gca,'FontSize',14);

if print_on ==1
    print('-r400','-dpng','WOB.png');
else 
end

figure('Position',figPos1);
sz = 40;
scatter(young.R1insp,young.E_Array,sz,'^',...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.R1insp,old.E_Array,sz,'v',...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.R1insp,NFL.E_Array,sz,...
                'MarkerEdgeColor',colorOrange,...
                'MarkerFaceColor',colorOrange,...
                'LineWidth',1.5); hold on
scatter(FL.R1insp,FL.E_Array,sz,'d',...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on

totalDataX = [young.R1insp old.R1insp NFL.R1insp FL.R1insp]';
totalDataY = [young.E_Array'; old.E_Array'; NFL.E_Array'; FL.E_Array'];
[p,S]=polyfit(totalDataX,totalDataY,1);
[yfit,delta]=polyval(p,totalDataX,S);
plot(totalDataX,yfit,'k');
RSS = sum((totalDataY-yfit).^2);
TSS = sum((totalDataY-mean(totalDataY)).^2);
Rsquared = 1-RSS/TSS;
anno = annotation('textbox', [0.70, 0.2, 0.1, 0.1], 'String', "R^2 = " + round(Rsquared,2));
anno.EdgeColor = 'w';
anno.BackgroundColor = 'w';
anno.FontSize = 12;
% h1.Color = 'k';
% h1.LineWidth = 1.5;
xlabel('R_{1,insp}')
ylabel('E (cmH_2O/L)')
% legend('young','old','NFL','FL','Location','southwest')
% xlim([0 20])
% ylim([0 35])
grid on
set(gca,'FontSize',14);


if print_on ==1
    print('-r400','-dpng','E.png');
else 
end


       