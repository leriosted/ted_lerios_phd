close all; clear all; clc

figPos1 = [2580, -190, 1000, 1270];
colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
colorGreen = '#009A17';
colorOrange = '#FF6600';
pt = 1.0;
FS = 10;

tlim = [0 4];
qlim = [0 1.5];
vlim = [0 2.0];
plim = [-40 15];

load young_results.mat
load old_results.mat
load FL_results.mat
load NFL_results.mat
%% Analysis

for ii = 1:20
    young.passiveQ_RMSE(ii,1) = rmse(young.PassiveQ(:,ii),young.ExpiredQ(:,ii));
    R = corrcoef(young.PassiveQ(:,ii),young.ExpiredQ(:,ii));
    young.CorrCoeff(ii,1) = R(2,1);
    young.int_V_diff(ii,1) = young.PassiveV(end,ii) - young.ExpiredV(end,ii);
end

for ii = 1:20
    old.passiveQ_RMSE(ii,1) = rmse(old.PassiveQ(:,ii),old.ExpiredQ(:,ii));
    R = corrcoef(old.PassiveQ(:,ii),old.ExpiredQ(:,ii));
    old.CorrCoeff(ii,1) = R(2,1);
    old.int_V_diff(ii,1) = old.PassiveV(end,ii) - old.ExpiredV(end,ii);
end

for ii = 1:25
    NFL.passiveQ_RMSE(ii,1) = rmse(NFL.PassiveQ(:,ii),NFL.ExpiredQ(:,ii));
    R = corrcoef(NFL.PassiveQ(:,ii),NFL.ExpiredQ(:,ii));
    NFL.CorrCoeff(ii,1) = R(2,1);
    NFL.int_V_diff(ii,1) = NFL.PassiveV(end,ii) - NFL.ExpiredV(end,ii);
end

for ii = 1:35
    FL.passiveQ_RMSE(ii,1) = rmse(FL.PassiveQ(:,ii),FL.ExpiredQ(:,ii));
    R = corrcoef(FL.PassiveQ(:,ii),FL.ExpiredQ(:,ii));
    FL.CorrCoeff(ii,1) = R(2,1);
    FL.int_V_diff(ii,1) = FL.PassiveV(end,ii) - FL.ExpiredV(end,ii);
end

% young.RMSE_quartiles = quantile(young.RMSE,3);
% old.RMSE_quartiles = quantile(old.RMSE,3);
% NFL.RMSE_quartiles = quantile(NFL.RMSE,3);
% FL.RMSE_quartiles = quantile(FL.RMSE,3);

VarNames = ["Patient Catagory","Passive Q RMSE",'Passive Q CorrCoeff','Vp-Ve'];

catagory = cell(4,1);
catagory(1) = {'Y'};
catagory(2) = {'O'};
catagory(3) = {'N'};
catagory(4) = {'F'};

RMSE_PassiveQ = cell(4,1);
RMSE_PassiveQ(1) = {quantile(young.passiveQ_RMSE,3)};
RMSE_PassiveQ(2) = {quantile(old.passiveQ_RMSE,3)};
RMSE_PassiveQ(3) = {quantile(NFL.passiveQ_RMSE,3)};
RMSE_PassiveQ(4) = {quantile(FL.passiveQ_RMSE,3)};

cc_PassiveQ = cell(4,1);
cc_PassiveQ(1) = {quantile(young.CorrCoeff,3)};
cc_PassiveQ(2) = {quantile(old.CorrCoeff,3)};
cc_PassiveQ(3) = {quantile(NFL.CorrCoeff,3)};
cc_PassiveQ(4) = {quantile(FL.CorrCoeff,3)};

Int_V_diff = cell(4,1);
Int_V_diff(1) = {quantile(young.int_V_diff*1000,3)};
Int_V_diff(2) = {quantile(old.int_V_diff*1000,3)};
Int_V_diff(3) = {quantile(NFL.int_V_diff*1000,3)};
Int_V_diff(4) = {quantile(FL.int_V_diff*1000,3)};

T = table(catagory,RMSE_PassiveQ,cc_PassiveQ,Int_V_diff,'VariableNames',VarNames);
disp(T)

T2 = rows2vars(T);
disp(T2)

writetable(T,'table_PassiveQ_IQR_MED.csv')
writetable(T2,'table_PassiveQ_IQR_MED_2.csv')

% save('young_results.mat','young','-append')
% save('old_results.mat','old','-append')
% save('NFL_results.mat','NFL','-append')
% save('FL_results.mat','FL','-append')

%% Plots
figure('Position',figPos1);
for ii = 1:20
    subplot(8,5,ii)
    plot(young.PassiveV(:,ii),young.PassiveQ(:,ii),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
    plot(young.ExpiredV(:,ii),young.ExpiredQ(:,ii),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.'); hold on
    xlabel('V (L)')
    ylabel('Q (L/s)')
    
    xlim(vlim)
    ylim(qlim)
    grid on
    set(gca,'FontSize',FS);
    title(['Patient ',num2str(ii)])
%     title(['RMSE: ',num2str(round(young.passiveQ_RMSE(ii),2))]);
end
leg1 = legend('Passive V - Q','Expired V - Q','FontSize',12,...
    'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.62 .48];
%     print('-r400','-dpng','V-Q_young.png');

figure('Position',figPos1);
for ii = 1:20
    subplot(8,5,ii)
    plot(old.PassiveV(:,ii),old.PassiveQ(:,ii),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
    plot(old.ExpiredV(:,ii),old.ExpiredQ(:,ii),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.'); hold on
    xlabel('V (L)')
    ylabel('Q (L/s)')
    
    xlim(vlim)
    ylim(qlim)
    grid on
    set(gca,'FontSize',FS);
    title(['Patient ',num2str(ii)])
%     title(['RMSE: ',num2str(round(old.passiveQ_RMSE(ii),2))]);
end
leg1 = legend('Passive V - Q','Expired V - Q','FontSize',12,...
    'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.62 .48];
%     print('-r400','-dpng','V-Q_old.png');

figure('Position',figPos1);
for ii = 1:25
    subplot(8,5,ii)
    plot(NFL.PassiveV(:,ii),NFL.PassiveQ(:,ii),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
    plot(NFL.ExpiredV(:,ii),NFL.ExpiredQ(:,ii),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.'); hold on
    xlabel('V (L)')
    ylabel('Q (L/s)')

    xlim(vlim)
    ylim(qlim)
    grid on
    set(gca,'FontSize',FS);
    title(['Patient ',num2str(ii)])
%     title(['RMSE: ',num2str(round(NFL.passiveQ_RMSE(ii),2))]);
end

    leg1 = legend('Passive V - Q','Expired V - Q','FontSize',12,...
    'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.62 .38];
%     print('-r400','-dpng','V-Q_NFL.png');

figure('Position',figPos1);
for ii = 1:35
    subplot(8,5,ii)
    plot(FL.PassiveV(:,ii),FL.PassiveQ(:,ii),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
    plot(FL.ExpiredV(:,ii),FL.ExpiredQ(:,ii),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.'); hold on
    xlabel('V (L)')
    ylabel('Q (L/s)')

    xlim(vlim)
    ylim(qlim)
    grid on
    set(gca,'FontSize',FS);
    title(['Patient ',num2str(ii)])
%     title(['RMSE: ',num2str(round(FL.passiveQ_RMSE(ii),2))]);
end

    leg1 = legend('Passive V - Q','Expired V - Q','FontSize',12,...
    'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.62 .16];
%     print('-r400','-dpng','V-Q_FL.png');
%%

figure('Position',[2580, -190, 1000, 600])
subplot(3,4,1)    
plot(young.PassiveV(:,6),young.PassiveQ(:,6),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(young.ExpiredV(:,6),young.ExpiredQ(:,6),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.');
xlabel('V (L)')
ylabel('Q (L/s)')
leg1 = legend('Passive V - Q','Expired V - Q','FontSize',12,...
    'Orientation','horizontal');
legend boxoff
leg1.Position(1:2) = [.62 .30];
xlim(vlim)
ylim(qlim)
title('a)')
grid on
set(gca,'FontSize',FS);

subplot(3,4,2)    
plot(old.PassiveV(:,2),old.PassiveQ(:,2),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(old.ExpiredV(:,2),old.ExpiredQ(:,2),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('b)')
grid on
set(gca,'FontSize',FS);

subplot(3,4,3)    
plot(NFL.PassiveV(:,1),NFL.PassiveQ(:,1),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(NFL.ExpiredV(:,1),NFL.ExpiredQ(:,1),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('c)')
grid on
set(gca,'FontSize',FS);

subplot(3,4,4)    
plot(FL.PassiveV(:,1),FL.PassiveQ(:,1),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(FL.ExpiredV(:,1),FL.ExpiredQ(:,1),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('d)')
grid on
set(gca,'FontSize',FS);

subplot(3,4,5)    
plot(young.PassiveV(:,15),young.PassiveQ(:,15),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(young.ExpiredV(:,15),young.ExpiredQ(:,15),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('e)')
grid on
set(gca,'FontSize',FS);

subplot(3,4,6)    
plot(old.PassiveV(:,17),old.PassiveQ(:,17),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(old.ExpiredV(:,17),old.ExpiredQ(:,17),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('f)')
grid on
set(gca,'FontSize',FS);

subplot(3,4,7)    
plot(NFL.PassiveV(:,17),NFL.PassiveQ(:,17),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(NFL.ExpiredV(:,17),NFL.ExpiredQ(:,17),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('g)')
grid on
set(gca,'FontSize',FS);

subplot(3,4,8)    
plot(FL.PassiveV(:,10),FL.PassiveQ(:,10),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(FL.ExpiredV(:,10),FL.ExpiredQ(:,10),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('h)')
grid on
set(gca,'FontSize',FS);

print('-r400','-dpng','forced_exp_typical.png');

%% Box Plot

cohorts = ["Young","Young","Young","Young","Young",...
    "Young","Young","Young","Young","Young",...
    "Young","Young","Young","Young","Young",...
    "Young","Young","Young","Young","Young",...
    "Elderly","Elderly","Elderly","Elderly","Elderly",...
    "Elderly","Elderly","Elderly","Elderly","Elderly",...
    "Elderly","Elderly","Elderly","Elderly","Elderly",...
    "Elderly","Elderly","Elderly","Elderly","Elderly",...
    "NFL","NFL","NFL","NFL","NFL",...
    "NFL","NFL","NFL","NFL","NFL",...
    "NFL","NFL","NFL","NFL","NFL",...
    "NFL","NFL","NFL","NFL","NFL",...
    "NFL","NFL","NFL","NFL","NFL",...
    "FL","FL","FL","FL","FL",...
    "FL","FL","FL","FL","FL",...
    "FL","FL","FL","FL","FL",...
    "FL","FL","FL","FL","FL",...
    "FL","FL","FL","FL","FL",...
    "FL","FL","FL","FL","FL",...
    "FL","FL","FL","FL","FL"];

RMSE_ALL = [young.passiveQ_RMSE;old.passiveQ_RMSE;NFL.passiveQ_RMSE;FL.passiveQ_RMSE]';

figure('Position',[2580, -190, 400, 300])
boxplot(RMSE_ALL,cohorts)
ylabel('Q_P RMSE')
grid on
set(gca,'FontSize',FS);
print('-r400','-dpng','boxplot_RMSE.png');

CorrCoeffALL = [young.CorrCoeff;old.CorrCoeff;NFL.CorrCoeff;FL.CorrCoeff]';
figure('Position',[2580, -190, 400, 300])
boxplot(CorrCoeffALL,cohorts)
ylabel('Q_P r')
grid on
set(gca,'FontSize',FS);
print('-r400','-dpng','boxplot_CC.png');
