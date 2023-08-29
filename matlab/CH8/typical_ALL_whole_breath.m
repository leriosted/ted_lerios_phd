clear all; close all; clc
colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
pt = 0.5;
FS = 14;
figure('Position',[-1900 100 700 900])
for kk = 1:6
    file = append('Patient_',num2str(kk));
    load(file)

for ii = 1:length(breaths)
    Paw = breaths(ii).Paw;
    Treal = breaths(ii).Treal;
    V = breaths(ii).V - leak*breaths(ii).T; P = breaths(ii).Paw; Q = breaths(ii).Q-leak;
    T = breaths(ii).T;
    Pes = breaths(ii).Pes-mode(data.Pes);
    if kk == 4
        Pes = breaths(ii).Pes - breaths(ii).Pes(end);
    end
%     Pes = breaths(ii).Pes - breaths(ii).Pes(end);
    
    Paw = breaths(ii).Paw;
    PEEP = P(1);
    PEEPs(ii,1) = PEEP;
    aa = min(find(P==max(P),1,'last'),find(V == max(V),1,"first"));
    aa1 = find(P==max(P),1,'last');
    aa2 = find(V == max(V),1,"first");
    [val,ind] = max(V);
    insp_end = ind;
    [val2,ind2] = min(Q);
    Q_peak = ind2;

%     E(ii,1) = (P(aa2)-PEEP)/V(aa2);                 % Analysis 1

%     E(ii,1) = (P(aa)-PEEP)/V(aa);                   % Analysis 2

    E(ii,1) = (P(aa2) - PEEP - Pes(aa2))./V(aa2);   % Analysis 3

    A = Q(Q_peak:end);
    b = P(Q_peak:end) - PEEP - E(ii)*V(Q_peak:end);
    R(ii,1) = A\b;
    
end

medE = median(E);
medR = median(R);

for jj = 1:length(breaths)
    M = 50;
    Paw = breaths(jj).Paw;
    Treal = breaths(jj).Treal;
    V = breaths(jj).V - leak*breaths(jj).T; P = breaths(jj).Paw; Q = breaths(jj).Q-leak;
    T = breaths(jj).T;
    Pes = breaths(jj).Pes-mode(data.Pes);

        if kk == 4
            Pes = breaths(jj).Pes - breaths(jj).Pes(end);
        end

    Paw = breaths(jj).Paw;
    PEEP = P(1);
    [t_spline,y_spline_Peff] = b_spline_basis_functions(M,2,length((T))/100);
    Psi = zeros(length(T),M); % \Psi = Splines_Overall_Peff
    Psi = y_spline_Peff(1:end,:);
    if length(Psi(:,1)) == length(Pes)+1
        Psi = Psi(1:end-1,:);
    end
    Peff_coeff = lsqlin(Psi, (Paw-PEEP)-medE*(V)-medR*Q);
    Peff = Psi(:,:)*Peff_coeff;

    RMSE(jj,1) = rmse(Peff,Pes);
    r = corrcoef(Peff,Pes);
    cc(jj,1) = r(2,1);
    
    minPeff(jj,1) = min(Peff(1:end));
    minPes(jj,1) = min(Pes(1:end));

%% PLOTS
    subplot(6,1,kk)
    plot(Treal,Paw,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
    plot(Treal,Pes,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
    plot(Treal,Peff,'Color',colorRed,'Linewidth',pt,'Linestyle','-'); hold on
    grid on
    xlim([0 150])
    ylim([-20 30])
    xlabel('t (s)')
    ylabel('P (cmH_2O)')
    set(gca,'FontSize',FS);


%     subplot(2,1,2)
%     plot(Treal,Q,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
%     ylabel('Q')
%     grid on

end

% % Spline Check
% figure()
% for iii = 1:50
%     plot(Psi(:,iii)*Peff_coeff(iii)); hold on
% end
% plot(Peff,'r','Linewidth',2)
% 
% %% Correlation Coefficient Plot
% FigH = figure('Position',[-1900 700 300 300]);
% coefficients = polyfit(minPes, minPeff, 1);
% 
% xFit = linspace(min(minPes), max(minPes), 1000);
% 
% yFit = polyval(coefficients , xFit);
% 
% plot(minPes, minPeff, 'b.', 'MarkerSize', 25); 
% hold on; 
% plot(xFit, yFit, 'r-', 'LineWidth', 3); 
% plot(xlim,ylim,'k-', 'LineWidth',2);
% grid on;
% xlabel('min P_{es} [cmH_2O]')
% ylabel('min P_{eff,3} [cmH_2O]')
% str1 = append('Patient ',num2str(kk));
% title(str1)
% tmp=corrcoef(minPes,minPeff);
% str=sprintf('r= %1.2f',tmp(1,2));
% % T = text(-4,-8, str); 
% % set(T, 'fontsize', 12);
% 
% AxesH = axes('Parent', FigH, ...
%   'Units', 'normalized', ...
%   'Position', [0, 0, 1, 1], ...
%   'Visible', 'off', ...
%   'XLim', [0, 1], ...
%   'YLim', [0, 1], ...
%   'NextPlot', 'add');
% TextH = text(0,1,str, ...
%   'HorizontalAlignment', 'left', ...
%   'VerticalAlignment', 'top');
% fontsize(FigH,20,"pixels")
% filename2 = append('Patient_cc_',num2str(kk),'A3');
% print('-r400','-dpng',filename2);

pp_E(kk,1) = medE;
pp_R(kk,1) = medR;
% pp_r(kk,1) = tmp(1,2);
pp_PEEP(kk,1) = median(PEEPs);
end
leg1 = legend('P_{aw}','P_{es}','P_{eff,3}','FontSize',12,...
        'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.56 .05];
set(gca,'FontSize',FS);
    filename1 = append('typical_ALL_whole_breath_A1');
    print('-r400','-dpng',filename1);

