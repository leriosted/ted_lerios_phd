clear all; close all; clc

colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
pt = 0.5;
FS = 10;

for jj = 1:6
    file = append('Patient_',num2str(jj));

load(file)

for ii = 1:length(breaths)
    Paw = breaths(ii).Paw;
    Treal = breaths(ii).Treal;
    V = breaths(ii).V - leak*breaths(ii).T; P = breaths(ii).Paw; Q = breaths(ii).Q-leak;
    T = breaths(ii).T;
    Pes = breaths(ii).Pes-mode(data.Pes);
    Paw = breaths(ii).Paw;
    PEEP = P(1);
    [val,ind] = max(V);
    insp_end = ind;
    [val2,ind2] = max(Q);
    Q_peak = ind2;
    E = (Paw(insp_end) - PEEP)/V(insp_end);
    R = Q(Q_peak:end)\(Paw(Q_peak:end) - PEEP - E*V(Q_peak:end));
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
    Peff_coeff = lsqlin(Psi, (Paw)-E*(V)-R_med*Q);
    Peff = Psi(:,:)*Peff_coeff;
%     Peff = Peff - mode(data.Pes);
    Peff_mean_Array(ii,1) = mean(Peff);
    E_Array(ii,1) = E;
    PEEP_Array(ii,1) = PEEP;
    RMSE(ii,1) = rmse(Peff,Pes);
    r = corrcoef(Peff,Pes);
    cc(ii,1) = r(2,1);
    
    minPeff(ii,1) = min(Peff(1:end));
    minPes(ii,1) = min(Pes(1:end));

    %% PLOTS
    figure('Position',[-1900 150 300 300])
    plot(T,Pes,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
    plot(T,Peff,'Color',colorRed,'Linewidth',pt,'Linestyle','-');
    plot(T,Paw,'Color',colorBlue,'Linewidth',pt,'Linestyle','-');
    grid on
    xlabel('t (s)')
    ylabel('P (cmH_2O)')
    ylim([-20 30])
    xlim([0 8])
%     legend('Pes','Peff','Paw')
    set(gca,'FontSize',FS);
    filename1 = append('patient_',num2str(jj),'_breath',num2str(ii));
    print('-r400','-dpng',filename1);
end
end