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

young.index =   13;
old.index =     13;
NFL.index =     24;
FL.index =      30;

tlim = [0 4];
qlim = [0 1.5];
vlim = [0 2.0];
plim = [-40 15];

figure('Position',figPos1)
subplot(3,2,1)    
plot(young.PassiveV(:,young.index),young.PassiveQ(:,young.index),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(young.ExpiredV(:,young.index),young.ExpiredQ(:,young.index),'Color',colorBlue,'Linewidth',pt,'Linestyle','--');
xlabel('V (L)')
ylabel('Q (L/s)')
leg1 = legend('Passive V - Q','Expired V - Q','FontSize',12,...
    'Orientation','horizontal');
legend boxoff
leg1.Position(1:2) = [.62 .34];
xlim(vlim)
ylim(qlim)
title('Young')
grid on
set(gca,'FontSize',FS);

subplot(3,2,2)    
plot(old.PassiveV(:,old.index),old.PassiveQ(:,old.index),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(old.ExpiredV(:,old.index),old.ExpiredQ(:,old.index),'Color',colorBlue,'Linewidth',pt,'Linestyle','--');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('Elderly')
grid on
set(gca,'FontSize',FS);

subplot(3,2,3)    
plot(NFL.PassiveV(:,NFL.index),NFL.PassiveQ(:,NFL.index),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(NFL.ExpiredV(:,NFL.index),NFL.ExpiredQ(:,NFL.index),'Color',colorBlue,'Linewidth',pt,'Linestyle','--');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('NFL')
grid on
set(gca,'FontSize',FS);

subplot(3,2,4)    
plot(FL.PassiveV(:,FL.index),FL.PassiveQ(:,FL.index),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(FL.ExpiredV(:,FL.index),FL.ExpiredQ(:,FL.index),'Color',colorBlue,'Linewidth',pt,'Linestyle','--');
xlabel('V (L)')
ylabel('Q (L/s)')
xlim(vlim)
ylim(qlim)
title('FL')
grid on
set(gca,'FontSize',FS);

print('-r400','-dpng','forced_exp_typical.png');


    %%
% figure('Position',figPos1)
% % FL
% 
% subplot(8,5,1)
% plot(FL.volume(1,FL.exp_point(1):end),-FL.flow(1,FL.exp_point(1):end),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
% % plot(FL.time(1:FL.exp_point(1)-1),FL.Peff_Array{1, 1},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
% % plot(FL.time,FL.Peff_Array{1, 1},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
% xlabel('V (L)')
% ylabel('Q (L/s)')
% % leg1 = legend('P_{alv}','P_{eff}','FontSize',12,...
% %     'Orientation','horizontal');
% % legend boxoff
% % leg1.Position(1:2) = [.75 .16];
% xlim(vlim)
% ylim(qlim)
% grid on
% set(gca,'FontSize',FS);
% 
% for ii = 2:35
%     subplot(8,5,ii)
%     plot(FL.volume(ii,FL.exp_point(ii):end),-FL.flow(ii,FL.exp_point(ii):end),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
% %     plot(FL.time(1:FL.exp_point(ii)-1),FL.Peff_Array{1, ii},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
% %     plot(FL.time,FL.Peff_Array{1, ii},'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
%     xlabel('V (L)')
%     ylabel('Q (L/s)')
%     xlim(vlim)
%     ylim(qlim)
%     grid on
%     set(gca,'FontSize',FS);
% end
%   
% if print_on ==1
%     print('-r400','-dpng','V-Q_FL.png');
% else 
% end