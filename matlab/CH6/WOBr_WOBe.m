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
    young.WOBr(ii,1) = young.WOBr_dV{1, ii}(end,1);
    young.WOBe(ii,1) = young.WOBe_dV{1, ii}(end,1);
end

for ii = 1:20
    old.WOBr(ii,1) = old.WOBr_dV{1, ii}(end,1);
    old.WOBe(ii,1) = old.WOBe_dV{1, ii}(end,1);
end

for ii = 1:25
    NFL.WOBr(ii,1) = NFL.WOBr_dV{1, ii}(end,1);
    NFL.WOBe(ii,1) = NFL.WOBe_dV{1, ii}(end,1);
end

for ii = 1:35
    FL.WOBr(ii,1) = FL.WOBr_dV{1, ii}(end,1);
    FL.WOBe(ii,1) = FL.WOBe_dV{1, ii}(end,1);
end
%% Linear regression line
% dlm = fitlm(NFL.R1insp,NFL.meanR2Phi,'Intercept', false)

%%
figure('Position',figPos1);

sz = 40;
scatter(young.WOBr,young.WOBe,sz,'^',...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(old.WOBr,old.WOBe,sz,'v',...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(NFL.WOBr,NFL.WOBe,sz,...
                'MarkerEdgeColor',colorOrange,...
                'MarkerFaceColor',colorOrange,...
                'LineWidth',1.5); hold on
scatter(FL.WOBr,FL.WOBe,sz,'d',...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on

totalDataX = [young.WOBr; old.WOBr; NFL.WOBr; FL.WOBr];
totalDataY = [young.WOBe; old.WOBe; NFL.WOBe; FL.WOBe];
[p,S]=polyfit(totalDataX,totalDataY,1);
[yfit,delta]=polyval(p,totalDataX,S);
plot(totalDataX,yfit,'k');
RSS = sum((totalDataY-yfit).^2);
TSS = sum((totalDataY-mean(totalDataY)).^2);
Rsquared = 1-RSS/TSS;
anno = annotation('textbox', [0.67, 0.2, 0.1, 0.1], 'String', "R^2 = 0.60");
anno.EdgeColor = 'w';
anno.BackgroundColor = 'w';
anno.FontSize = 12;

xlabel('$ \int_{}^{} R_{1,insp}Q \,dV $','Interpreter','latex')
ylabel('$ \int_{}^{} EV \,dV $','Interpreter','latex')
legend('young','old','NFL','FL','Location','northwest')
% xlim([0 20])
% ylim([0 35])
grid on
set(gca,'FontSize',14);



if print_on ==1
    print('-r400','-dpng','WOBeWOBr.png');
else 
end


       