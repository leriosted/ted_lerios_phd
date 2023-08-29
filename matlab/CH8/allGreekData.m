clear all; close all; clc

colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
pt = 0.5;
FS = 10;

load('Patient_1.mat')
figure('Position',[-1900 150 1800 600])
for ii = 1:length(breaths)
    Paw = breaths(ii).Paw;
    Treal = breaths(ii).Treal;
    V = breaths(ii).V - leak*breaths(ii).T; P = breaths(ii).Paw; Q = breaths(ii).Q-leak;
    T = breaths(ii).T;
    Pes = breaths(ii).Pes;
    Paw = breaths(ii).Paw;
    PEEP = P(1);
    [val,ind] = max(V);
    insp_end = ind;
    [val2,ind2] = max(Q);
    Q_peak = ind2;
    E = (Paw(insp_end) - PEEP)/V(insp_end);
    R = Q(Q_peak)\(Paw(Q_peak) - PEEP - E*V(Q_peak));
    R_Array(ii,1) = R;
    R_med = median(R); 
    M = 50;
    [t_spline,y_spline_Peff] = b_spline_basis_functions(M,2,length((T))/100);
    Psi = zeros(length(T),M); % \Psi = Splines_Overall_Peff
    Psi = y_spline_Peff(1:end,:);
    if length(Psi(:,1)) == length(Pes)+1
        Psi = Psi(1:end-1,:);
    end
%     Peff_coeff = lsqlin(Psi, (Paw-PEEP)-E*(V)-R_med*Q);
    Peff_coeff = lsqlin(Psi, (Paw-PEEP)-E*(V)-R_med*Q);
    Peff = Psi(:,:)*Peff_coeff;
    Peff_mean_Array(ii,1) = mean(Peff);
    E_Array(ii,1) = E;
    PEEP_Array(ii,1) = PEEP;
    RMSE(ii,1) = rmse(Peff-PEEP,Pes-PEEP);
    r = corrcoef(Peff-PEEP,Pes-PEEP);
    cc(ii,1) = r(2,1);
    
    minPeff(ii,1) = min(Peff(1:insp_end));
    minPes(ii,1) = min(Pes(1:insp_end));

    %% PLOTS
    subplot(2,1,1)
    plot(Treal,Pes-mode(data.Pes),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
    plot(Treal,Peff,'Color',colorRed,'Linewidth',pt,'Linestyle','-');
    plot(Treal,Paw,'Color',colorBlue,'Linewidth',pt,'Linestyle','-');
    grid on
    ylabel('P')
    legend('Pes','Peff','Paw')
    subplot(2,1,2)
    plot(Treal,Q,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
    ylabel('Q')
%     xlim([0 10])
%     ylim([-10 20])
    grid on
%     title(['E: ',num2str(E),' | R: ',num2str(R),' | Breath: ',num2str(ii)])
%     set(gca,'FontSize',FS);
end
    set(gca,'FontSize',FS);

    print('-r400','-dpng','GreekData_ALL.png');

    figure('Position',[-1900 700 300 300])
coefficients = polyfit(minPes, minPeff, 1);

xFit = linspace(min(minPes), max(minPes), 1000);

yFit = polyval(coefficients , xFit);

plot(minPes, minPeff, 'b.', 'MarkerSize', 15); 
hold on; 
plot(xFit, yFit, 'r-', 'LineWidth', 2); 
plot(xlim,ylim,'k-')
grid on;
xlabel('min P_{es} [cmH_2O]')
ylabel('min P_{eff} [cmH_2O]')
title('Patient 1')
tmp=corrcoef(minPes,minPeff);
str=sprintf('r= %1.2f',tmp(1,2));
T = text(8,-3.25, str); 
set(T, 'fontsize', 12);

print('-r400','-dpng','Patient_1_cc.png');