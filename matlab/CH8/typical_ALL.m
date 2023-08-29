clear all; close all; clc
colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
pt = 0.5;
FS = 10;
figure('Position',[-1900 100 700 900])
for kk = 1:6
    file = append('Patient_',num2str(kk));
    load(file)

for ii = 1:length(breaths)
    Paw = breaths(ii).Paw;
    Treal = breaths(ii).Treal;
    V = breaths(ii).V - leak*breaths(ii).T; P = breaths(ii).Paw; Q = breaths(ii).Q-leak;
    T = breaths(ii).T;
    
    Pes = breaths(ii).Pes;
%     if kk == 4
%         Pes = breaths(ii).Pes - breaths(ii).Pes(end);
%     end
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

    E(ii,1) = (P(aa2)-PEEP)/V(aa2);                 % Analysis 1

%     E(ii,1) = (P(aa)-PEEP)/V(aa);                   % Analysis 2

%     E(ii,1) = (P(aa2) - PEEP - Pes(aa2))./V(aa2);   % Analysis 3

    A_insp = Q(Q_peak:end);
    b_insp = -P(Q_peak:end) + PEEP;
    R_insp(ii,1) = A_insp\b_insp;
    
end

medE = median(E);
medR = median(R_insp);

for jj = 1:length(breaths)
    M = 15;
    Paw = breaths(jj).Paw;
    Treal = breaths(jj).Treal;
    V = breaths(jj).V - leak*breaths(jj).T; P = breaths(jj).Paw; Q = breaths(jj).Q-leak;
    T = breaths(jj).T;
    Pes = breaths(jj).Pes;
    Paw = breaths(jj).Paw;
    PEEP = P(1);

    insp_end = insp_end-1; % MOVE EXP POINT

    [t_spline,y_spline_Peff] = b_spline_basis_functions(M,2,length((T(insp_end-1:end)))/100);
    Psi = zeros(length(T(insp_end-1:end)),M); % \Psi = Splines_Overall_Peff
    Psi = y_spline_Peff(1:end,:);
    if length(Psi(:,1)) == length(Paw(insp_end-1:end))+1
        Psi = Psi(1:end-1,:);
    end

    A_exp = [Q(insp_end-1:end).*Psi];
    b_exp = [(-Paw(insp_end-1:end) + PEEP - E(jj)*V(insp_end-1:end))];
%     x_exp = A_exp\b_exp;
    x_exp = lsqnonneg(A_exp,b_exp);


    R_nlin_coeff = x_exp;

    P_exp_model_nlin = Psi(:,:)*R_nlin_coeff.*(Q(insp_end-1:end)) + E(jj)*V(insp_end-1:end);

    P_insp_model = R_insp(jj)*Q(1:insp_end-2) + E(ii)*V(1:insp_end-2);
    P_exp_model_lin = R_insp(jj)*Q(insp_end-1:end) + E(ii)*V(insp_end-1:end);
     
    P_model_lin = [P_insp_model; P_exp_model_lin];
    P_model_nlin = [P_insp_model; P_exp_model_nlin];

    RMSE_lin(jj,1) = rmse(P_model_lin,Paw-PEEP);
    RMSE_nlin(jj,1) = rmse(P_model_nlin,Paw-PEEP);
    r_lin = corrcoef(P_model_lin,Paw);
    r_nlin = corrcoef(P_model_nlin,Paw);
    cc_lin(jj,1) = r_lin(2,1);
    cc_nlin(jj,1) = r_nlin(2,1);
    
%     minPeff(jj,1) = min(RQ_nlin(1:end));
%     minPes(jj,1) = min(Pes(1:end));

%% PLOTS
    subplot(7,1,kk)
    plot(Treal,Paw-PEEP,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
%     plot(Treal,Pes,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
    plot(Treal,P_model_lin,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); 
    plot(Treal,P_model_nlin,'Color',colorRed,'Linewidth',pt,'Linestyle','-'); 
    grid on
    xlim([0 150])
    ylim([-20 30])
    xlabel('t (s)')
    ylabel('P (cmH_2O)')
    set(gca,'FontSize',FS);

%% Save A1 results
% A1.Treal{:,jj,kk} = Treal;
% A1.T{:,jj,kk} = T;
% A1.Paw{:,jj,kk} = Paw;
% A1.Pes{:,jj,kk} = Pes;
% A1.Peff{:,jj,kk} = Peff;
%% Save A2 results
% A2.Treal{:,jj,kk} = Treal;
% A2.T{:,jj,kk} = T;
% A2.Paw{:,jj,kk} = Paw;
% A2.Pes{:,jj,kk} = Pes;
% A2.Peff{:,jj,kk} = Peff;
%% Save A3 results
% A3.Treal{:,jj,kk} = Treal;
% A3.T{:,jj,kk} = T;
% A3.Paw{:,jj,kk} = Paw;
% A3.Pes{:,jj,kk} = Pes;
% A3.Peff{:,jj,kk} = Peff;

end

%% Spline Check
figure('Position',[-1100 50 1080 900])
for iii = 1:M
    plot(Psi(:,iii)*R_nlin_coeff(iii)); hold on
end
plot(Psi*R_nlin_coeff,'r','Linewidth',2)

%% Correlation Coefficient Plot
% FigH = figure('Position',[-1900 700 300 300]);
% coefficients = polyfit(minPes, minPeff, 1);
% 
% xFit = linspace(min(minPes), max(minPes), 1000);
% 
% yFit = polyval(coefficients , xFit);
% 
% plot(minPes, minPeff, 'b.', 'MarkerSize', 15); 
% hold on; 
% plot(xFit, yFit, 'r-', 'LineWidth', 2); 
% plot(xlim,ylim,'k-')
% grid on;
% xlabel('min P_{es} [cmH_2O]')
% ylabel('min P_{eff} [cmH_2O]')
% % title('Patient 1')
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
% 
% filename2 = append('Patient_cc_',num2str(kk));
% print('-r400','-dpng',filename2);

pp_E(kk,1) = medE;
pp_R(kk,1) = medR;
pp_RMSE_lin(kk,1) = median(RMSE_lin);
pp_RMSE_nlin(kk,1) = median(RMSE_nlin);
% pp_r(kk,1) = tmp(1,2);
pp_PEEP(kk,1) = median(PEEPs);
end
leg1 = legend('P_{aw}','P_{lin}','P_{nlin}','FontSize',12,...
        'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.58 .185];
set(gca,'FontSize',FS);
    filename1 = append('typical_ALL');
    print('-r400','-dpng',filename1);

figure('position',[-1900 50 1800 900])
plot(T,Paw-PEEP,'k'); hold on
plot(T,P_model_lin,'b');
plot(T,P_model_nlin,'r');
plot(T,Q,'k--')
grid on
legend('P','lin model','nlin model','Q (L/s)')
xlabel('t (s)')
ylabel('P (cmH_2O)')
print('-r400','-dpng','test_model_plot');

figure('position',[2800 -200 700 700])
plot(Q(insp_end-1:end),Psi*R_nlin_coeff,'k'); hold on
grid on
xlabel('Q (L/s)')
ylabel('R (cmH_2O/L/s)')


% save('results.mat','A1')
% save('results.mat','A3','-append')

% save('Greek_data.mat','P1')
% save('Greek_data.mat','P1','-append')