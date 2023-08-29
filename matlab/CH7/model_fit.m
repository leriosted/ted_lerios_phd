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

young.linear_model = young.linear_model';
FL.linear_model = FL.linear_model';

tlim = [0 1.4];
qlim = [-1.0 1.0];
plim = [-17 15];
young.index = 1;
old.index = 2;
NFL.index = 2;
FL.index = 2;
% for ii = 1:10
%     young.WOBr(ii,1) = young.WOBr_dV{1, ii}(end,1);
%     old.WOBr(ii,1) = old.WOBr_dV{1, ii}(end,1);
%     NFL.WOBr(ii,1) = NFL.WOBr_dV{1, ii}(end,1);
%     FL.WOBr(ii,1) = FL.WOBr_dV{1, ii}(end,1);
%     young.WOBe(ii,1) = young.WOBe_dV{1, ii}(end,1);
%     old.WOBe(ii,1) = old.WOBe_dV{1, ii}(end,1);
%     NFL.WOBe(ii,1) = NFL.WOBe_dV{1, ii}(end,1);
%     FL.WOBe(ii,1) = FL.WOBe_dV{1, ii}(end,1);
% end
%% Linear regression line
% dlm = fitlm(NFL.R1insp,NFL.meanR2Phi,'Intercept', false)

%%
figure('Position',figPos1);
plot(young.time(1:young.exp_point),young.pressure(young.index,1:young.exp_point),'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
plot(young.time(1:young.exp_point(young.index)-1),young.linear_model(1:young.exp_point(young.index)-1),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.')

plot(FL.time(1:FL.exp_point),FL.pressure(FL.index,1:FL.exp_point),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
plot(FL.time(1:FL.exp_point(FL.index)-1),FL.linear_model(1:FL.exp_point(FL.index)-1),'Color',colorBlack,'Linewidth',pt,'Linestyle','-.')
xlabel('t (s)')
ylabel('P (cmH_2O)')
legend('young P_{alv}','young P_{model}','FL P_{alv}','FL P_{model}','FontSize',12);
% leg1 = legend('P_{alv}','P_{eff}','FontSize',12,...
%     'Orientation','horizontal');
% legend boxoff
% leg1.Position(1:2) = [.65 .24];
xlim(tlim)
ylim(plim)
grid on
set(gca,'FontSize',FS);

% totalDataX = [young.WOBr; old.WOBr; NFL.WOBr; FL.WOBr];
% totalDataY = [young.WOBe; old.WOBe; NFL.WOBe; FL.WOBe];
% [p,S]=polyfit(totalDataX,totalDataY,1);
% [yfit,delta]=polyval(p,totalDataX,S);
% plot(totalDataX,yfit,'k');
% RSS = sum((totalDataY-yfit).^2);
% TSS = sum((totalDataY-mean(totalDataY)).^2);
% Rsquared = 1-RSS/TSS;
% anno = annotation('textbox', [0.67, 0.2, 0.1, 0.1], 'String', "R^2 = " + round(Rsquared,2));
% anno.EdgeColor = 'w';
% anno.BackgroundColor = 'w';
% anno.FontSize = 12;





if print_on ==1
    print('-r400','-dpng','model_fit.png');
else 
end


       