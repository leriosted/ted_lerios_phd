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

tlim = [0 4];
qlim = [0 1.5];
vlim = [0 2.0];
plim = [-40 15];


%% Passive expiration

for ii = 1:20

    Pa = young.pressure(ii,:)';
    V = young.volume(ii,:)';
    Q = young.flow(ii,:)';
    T = young.time(ii,:)';
   
    aa = find(Q<=0,1,"first"):length(V);
    V2 = V(aa) - V(end);
    Vexp = max(V)-V(end);
    Texp = T(end)-T(find(Q<=0,1,"first"));
    A = pi*(Vexp)/(2*Texp);
    t(:,ii) = T(aa) - T(find(Q<=0,1,"first"));
    PassiveQ(:,ii) = A*sin(pi*t(:,ii)/Texp);
    PassiveV(:,ii) = cumtrapz(t(:,ii),PassiveQ(:,ii));
    ExpiredQ(:,ii) = -Q(aa);
    ExpiredV(:,ii) = -(V2-Vexp);
end

%%
figure('Position',figPos1)
subplot(8,5,1)
plot(PassiveV(:,1),PassiveQ(:,1),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
plot(ExpiredV(:,1),ExpiredQ(:,1),'Color',colorBlue,'Linewidth',pt,'Linestyle','--'); hold on
xlabel('V (L)')
ylabel('Q (L/s)')
leg1 = legend('Passive V - Q','Expired V - Q','FontSize',12,...
    'Orientation','horizontal');
legend boxoff
leg1.Position(1:2) = [.62 .48];
xlim(vlim)
ylim(qlim)
grid on
set(gca,'FontSize',FS);

for ii = 2:20
    subplot(8,5,ii)
    plot(PassiveV(:,ii),PassiveQ(:,ii),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
    plot(ExpiredV(:,ii),ExpiredQ(:,ii),'Color',colorBlue,'Linewidth',pt,'Linestyle','--'); hold on
    xlabel('V (L)')
    ylabel('Q (L/s)')
    xlim(vlim)
    ylim(qlim)
    grid on
    set(gca,'FontSize',FS);
end
  
if print_on ==1
    print('-r400','-dpng','V-Q_young.png');
else 
end
% figure('Position',figPos1)
% for ii = 1:20
%  
%     subplot(8,5,ii)
%     plot(t(:,ii),PassiveQ(:,ii),'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
%     plot(t(:,ii),ExpiredQ(:,ii),'Color',colorBlue,'Linewidth',pt,'Linestyle','--'); hold on
%     xlabel('t (s)')
%     ylabel('Q (L/s)')
% %     xlim(vlim)
% %     ylim(qlim)
%     grid on
%     set(gca,'FontSize',FS);
% end


% young.PassiveQ = PassiveQ;
% young.PassiveV = PassiveV;
% young.ExpiredQ = ExpiredQ;
% young.ExpiredV = ExpiredV;
% 
% save('young_results.mat','young','-append')

%% Analysis RMSE and CORRCOEF

for ii = 1:20
young.RMSE(ii,:) = rmse(young.PassiveQ(:,ii),young.ExpiredQ(:,ii));
end

       