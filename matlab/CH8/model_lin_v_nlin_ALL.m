clear all; close all; clc
colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
colorGreen = '#009A17';
colorOrange = '#FF6600';
colorYellow = "#EDB120";
pt = 0.5;
FS = 10;
% figure('Position',[-1900 100 700 900])
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
    M = 35;
    Paw = breaths(jj).Paw;
    Treal = breaths(jj).Treal;
    V = breaths(jj).V - leak*breaths(jj).T; P = breaths(jj).Paw; Q = breaths(jj).Q-leak;
    T = breaths(jj).T;
    Pes = breaths(jj).Pes;
    Paw = breaths(jj).Paw;
    PEEP = P(1);

    insp_end = insp_end +-0;

    insp = 1:insp_end;
    exp = insp_end+1:length(T);


    [t_spline,y_spline] = b_spline_basis_functions(M,2,length((T(exp)))/100);
    Psi = zeros(length(T),M); % \Psi = Splines_Overall_Peff

    if length(y_spline(:,1)) == length(exp) +1
        y_spline = y_spline(1:end-1,:);
    end

    Psi(exp,:) = y_spline;

    A_exp = [Q.*Psi];
    b_exp = [(Paw - PEEP - E(jj)*V)];
    x_exp = lsqnonneg(A_exp,b_exp);

    R_nlin_coeff = x_exp(1:M);
    
    P_model_lin = R_insp(jj)*Q + E(jj)*V + PEEP;

    P_model_nlin = [R_insp(jj)*Q(insp) + E(jj)*V(insp) + PEEP; Psi(exp,:)*R_nlin_coeff.*Q(exp) + E(jj)*V(exp) + PEEP];
% P_model_nlin = [R_insp(jj)*Q(insp) + E(jj)*V(insp) + PEEP; Psi(exp,:)*R_nlin_coeff.*Q(exp) + PEEP ];

    R2Phi(jj,1) = mean(Psi(exp,:)*R_nlin_coeff);
    RMSE_lin(jj,1) = rmse(P_model_lin,Paw);
    RMSE_nlin(jj,1) = rmse(P_model_nlin,Paw);

%% Pick a breath to plot
% if jj == 2 && kk==1
%     figure('position',[-900 50 900 500])
%     for iii = 1:M
%         plot(Psi(:,iii)*R_nlin_coeff(iii)); hold on
%     end
%     plot(Psi*R_nlin_coeff,'r','Linewidth',2)
% 
%     figure('position',[-1900 50 1000 500])
%     plot(T,Paw,'k'); hold on
%     plot(T,P_model_lin,'b');
%     plot(T,P_model_nlin,'r');
%     plot(T,Q,'k--')
%     plot(T,V,'b--')
%     grid on
%     legend('P','lin model','nlin model','Q (L/s)','V (L)')
% %     rrmmssee = rmse(Paw,P_model_nlin);
%     str = append('lin RMSE = ',num2str(rmse(Paw,P_model_lin)),' | nlin RMSE = ',num2str(rmse(Paw,P_model_nlin)));
%     title(str)
%     xlabel('t (s)')
%     ylabel('P (cmH_2O)')
%     print('-r400','-dpng','test_model_plot');
% end
%%

    r_lin = corrcoef(P_model_lin,Paw);
    r_nlin = corrcoef(P_model_nlin,Paw);
    cc_lin(jj,1) = r_lin(2,1);
    cc_nlin(jj,1) = r_nlin(2,1);
    


    if kk == 6 && jj == 1
        figure('position',[-1900 50 400 300])

        plot_ind = find(Q==min(Q),1,'last') + 0;
%         plot_ind_end = find(Q(plot_ind:end)==0,1,'first');
        plot_ind_end = find(T==T(end))-10;

        plot(T(plot_ind:plot_ind_end),Paw(plot_ind:plot_ind_end),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
%         plot(T,P_model_lin,'Color',colorBlue,'Linewidth',pt,'Linestyle','-');
        plot(T(plot_ind:plot_ind_end),P_model_nlin(plot_ind:plot_ind_end),'Color',colorRed,'Linewidth',pt,'Linestyle','-.');
%         plot(T,10*Q,'k')
        % plot(T,Q,'k--')
        grid on
        legend('P_{aw}','nonlinear R model','FontSize',FS)
        xlabel('t (s)')
        ylabel('P(t) (cmH_2O)')
%         xlim([0 2.5])
        ylim([0 15])
        set(gca,'FontSize',FS);
        print('-r400','-dpng','CH8_typical');
        
        figure('position',[2800 -200 400 300])
        plot(Q(plot_ind:plot_ind_end),Psi(plot_ind:plot_ind_end,:)*R_nlin_coeff,'Color',colorBlack,'Linewidth',2,'Linestyle','-'); hold on
        grid on
        xlabel('Q (L/s)')
        ylabel('R (cmH_2O/L/s)')
%         xlim([-1.1 0])
        ylim([0 30])
                grid on
        ylabel('R_{exp}\Phi(t) (cmH_2O/L/s)')
        str = append('Patient ',num2str(kk));
        title(str)
        set(gca,'FontSize',FS);
        print('-r400','-dpng','CH8_RQ');
        
        figure('Position',[-1100 50 400 300])
        for iii = 1:M
            plot(T(plot_ind:plot_ind_end),Psi(plot_ind:plot_ind_end,iii)*R_nlin_coeff(iii)); hold on
        end
        plot(T(plot_ind:plot_ind_end),Psi(plot_ind:plot_ind_end,:)*R_nlin_coeff,'r','Linewidth',2)
        grid on
        str = append('Patient ',num2str(kk));
        title(str)
        % legend('P_{aw}','linear model','nonlinear model','FontSize',FS)
        xlabel('t (s)')
        ylabel('R_{exp}\Phi(t) (cmH_2O/L/s)')
%         xlim([0 2.5])
        ylim([0 30])
        set(gca,'FontSize',FS);
        print('-r400','-dpng','CH8_splines');
    end
%     %% Spline Check


%         for iii = 1:M
%             splines(jj,kk) = Psi(:,iii)*R_nlin_coeff(iii);
%         end
%         sum_splines(jj,kk) = Psi*R_nlin_coeff;


%     minPeff(jj,1) = min(RQ_nlin(1:end));
%     minPes(jj,1) = min(Pes(1:end));

%% PLOTS
%     subplot(7,1,kk)
%     plot(Treal,Paw,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
% %     plot(Treal,Pes,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
%     plot(Treal,P_model_lin,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); 
%     plot(Treal,P_model_nlin,'Color',colorRed,'Linewidth',pt,'Linestyle','-'); 
%     grid on
%     xlim([0 150])
%     ylim([-20 30])
%     xlabel('t (s)')
%     ylabel('P (cmH_2O)')
%     set(gca,'FontSize',FS);

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

meanR2Phi{:,kk} = R2Phi;

%     %% Spline Check
%     if kk ==1
%         figure('Position',[-1100 50 1080 900])
%         for iii = 1:M
%             plot(Psi(:,iii)*R_nlin_coeff(iii)); hold on
%         end
%         plot(Psi*R_nlin_coeff,'r','Linewidth',2)
%     end

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

% pp_Paw(:,kk) = Paw;
% pp_Q(:,kk) = Q;
% pp_V(:,kk) = V;
% pp_P_model_lin(:,kk) = P_model_lin;
% pp_P_model_nlin(:,kk) = P_model_nlin;
% pp_insp(:,kk) = insp;
% pp_exp(:,kk) = exp;
% pp_T(:,kk) = T;
pp_E(kk,1) = medE;
pp_R(kk,1) = medR;
pp_R2Phi(kk,1) = median(R2Phi);
pp_RMSE_lin(kk,1) = median(RMSE_lin);
pp_RMSE_nlin(kk,1) = median(RMSE_nlin);
% pp_r(kk,1) = tmp(1,2);
pp_PEEP(kk,1) = median(PEEPs);
end
% leg1 = legend('P_{aw}','P_{lin}','P_{nlin}','FontSize',12,...
%         'Orientation','horizontal');
%     legend boxoff
%     leg1.Position(1:2) = [.58 .185];
%     set(gca,'FontSize',FS);
%     filename1 = append('typical_ALL');
%     print('-r400','-dpng',filename1);




% save('results.mat','A1')
% save('results.mat','A3','-append')

% save('Greek_data.mat','P1')
% save('Greek_data.mat','P1','-append')



figure('Position',[-1600 50 400 300])

[h2.P1,stats.P1] = cdfplot(meanR2Phi{:,1}); hold on
[h2.P2,stats.P2] = cdfplot(meanR2Phi{:,2}); hold on
[h2.P3,stats.P3] = cdfplot(meanR2Phi{:,3}); hold on
[h2.P4,stats.P4] = cdfplot(meanR2Phi{:,4}); hold on
[h2.P5,stats.P5] = cdfplot(meanR2Phi{:,5}); hold on
[h2.P6,stats.P6] = cdfplot(meanR2Phi{:,6}); hold on

set(h2.P1, 'LineStyle', '-', 'Color', colorBlue,'LineWidth',2);
set(h2.P2, 'LineStyle', '-', 'Color', colorGreen,'LineWidth',2);
set(h2.P3, 'LineStyle', '-', 'Color', colorRed,'LineWidth',2);
set(h2.P4, 'LineStyle', '-', 'Color', colorBlack,'LineWidth',2);
set(h2.P5, 'LineStyle', '-', 'Color', colorOrange,'LineWidth',2);
set(h2.P6, 'LineStyle', '-', 'Color', colorYellow,'LineWidth',2);

xlabel('mean R_{exp}(t) (cmH_2O/L/s)')
ylabel('Proportion')
xlim([0 35])
title(' ')
legend('P1','P2','P3','P4', 'P5','P6','Location','southeast','Fontsize',FS)
grid on
title(' ')
set(gca,'FontSize',FS);
       
print('-r400','-dpng','CH8_CDF_meanR2Phi.png');
