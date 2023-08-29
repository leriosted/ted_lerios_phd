clear all; close all; clc
colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
pt = 1.0;
FS = 10;
load('Patient_2.mat')
counter = 1;
% load('Patient_2.mat')
% load('Patient_3.mat')
% load('Patient_4.mat')
% load('Patient_5.mat')
% load('Patient_6.mat')
    figure('Position',[3630 -180 1900 950])
%     [val,ind] = max(cumtrapz(data.Q));
for ii = 1:length(breaths)
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

    Peff_coeff = lsqlin(Psi, Paw-E*(V)-R_med*Q);
    Peff = Psi(:,:)*Peff_coeff;

    Peff_mean_Array(ii,1) = mean(Peff);
    
    E_Array(ii,1) = E;
    PEEP_Array(ii,1) = PEEP;

    RMSE(ii,1) = rmse(Peff-PEEP,Pes-PEEP);

    r = corrcoef(Peff-PEEP,Pes-PEEP);
    cc(ii,1) = r(2,1);

        minPeff(ii,1) = min(Peff-PEEP);
    minPes(ii,1) = min(Pes-PEEP);

    P2(ii).Paw = Paw;
    P2(ii).Q = Q;
    P2(ii).V = V;
    P2(ii).E = E;
    P2(ii).R = R;
    P2(ii).Peff = Peff;
    P2(ii).Pes = Pes;
    P2(ii).PEEP = PEEP;
    P2(ii).insp_end = insp_end;
    P2(ii).r = r(2,1);
    P2(ii).RMSE = RMSE(ii,1);
    P2(ii).cc = cc(ii,1);

    subplot(8,5,ii)
%     plot(T,V,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
    plot(T,Pes-PEEP,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
    plot(T,Peff-PEEP,'Color',colorRed,'Linewidth',pt,'Linestyle','-');
    plot(T,Paw,'Color',colorBlue,'Linewidth',0.6,'Linestyle','-');
%     plot(T,Q,'Color',colorBlack,'Linewidth',0.4,'Linestyle',':')
    
    xlim([0 10])
    ylim([-10 20])
    grid on
    title(['E: ',num2str(E),' | R: ',num2str(R),' | Breath: ',num2str(ii)])
    set(gca,'FontSize',FS);
end
    set(gca,'FontSize',FS);
    leg1 = legend('P_{es}','P_{eff}','P_{aw}','FontSize',12,...
        'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.68 .185];
    sgtitle('Patient 2')

    print('-r400','-dpng','Patient_2_ALL.png');


meanE = mean(E_Array);
meanR = mean(R_Array);
meanPeff = mean(Peff_mean_Array);
meanRMSE = mean(RMSE);
meancc = mean(cc);

    figure('Position',[-1900 150 300 300])
    for ii = 1:length(breaths)
        V = breaths(ii).V - leak*breaths(ii).T; P = breaths(ii).Paw; Q = breaths(ii).Q-leak;
        T = breaths(ii).T;
        Pes = breaths(ii).Pes;
        Paw = breaths(ii).Paw;
        PEEP = P(1);
    
        [val,ind] = max(V);
        insp_end = ind;
    
        [val2,ind2] = max(Q);
        Q_peak = ind2;

        if ii == 4
            T = T-0.95;
       end
    
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
    
        Peff_coeff = lsqlin(Psi, Paw-E*(V)-R_med*Q);
        Peff = Psi(:,:)*Peff_coeff;
    
        Peff_mean_Array(ii,1) = mean(Peff);
        
        E_Array(ii,1) = E;
        PEEP_Array(ii,1) = PEEP;
    
        

        if ii == 3 | ii == 4 
            
            
            subplot(3,1,counter)
            counter = counter + 1;
            
            plot(T,Pes-PEEP,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
            plot(T,Peff-PEEP,'Color',colorRed,'Linewidth',pt,'Linestyle','-');
            plot(T,Paw,'Color',colorBlue,'Linewidth',0.6,'Linestyle','-');
            xlim([0 5])
            ylim([-15 22])
            grid on
%             title(['E: ',num2str(E),' | R: ',num2str(R)])
            set(gca,'FontSize',FS);
        end
    end
        set(gca,'FontSize',FS);
        leg1 = legend('P_{es}','P_{eff}','P_{aw}','FontSize',10,...
            'Orientation','horizontal');
        legend boxoff
        leg1.Position(1:2) = [.20 .185];
        sgtitle('Patient 2')
    
        print('-r400','-dpng','Patient_2.png');


                disp('Patient 2')
        disp(['mean E: ',num2str(meanE)])
        disp(['mean R: ',num2str(meanR)])
        disp(['mean Peff: ',num2str(meanPeff)])

        disp(['mean RMSE: ',num2str(meanRMSE)])
        disp(['mean cc: ',num2str(meancc)])


        figure('Position',[-1900 500 300 300])
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
title('Patient 2')
tmp=corrcoef(minPes,minPeff);
str=sprintf('r= %1.2f',tmp(1,2));
T = text(-5,-12, str); 
set(T, 'fontsize', 12);

print('-r400','-dpng','Patient_2_cc.png');

save('Greek_data.mat','P2','-append')
