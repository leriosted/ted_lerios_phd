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
    V = V - V(end);
    [val,ind] = max(V);
    insp_end = ind;
    [val2,ind2] = find(Q==min(Q));
    Q_peak = ind2;

    plot_ind = find(Q==min(Q),1,'last') ;
        plot_ind_end = find(T==T(end));
%         plot_ind_end = find(Q(plot_ind:end)>0);
    IND = plot_ind:plot_ind_end;
        IND(P(IND)<PEEP-2)=[];

    E(ii,1) = (P(aa2)-PEEP)/V(aa2);                 % Analysis 1

%     E(ii,1) = (P(aa)-PEEP)/V(aa);                   % Analysis 2

%     E(ii,1) = (P(aa2) - PEEP - Pes(aa2))./V(aa2);   % Analysis 3

    A_roh = [Q(IND), Q(IND).*abs(Q(IND))];
    b_roh = Paw(IND) - PEEP - E(ii)*V(IND);
    x_roh(ii,:) = A_roh\b_roh; 
    R_1(ii,1) = x_roh(ii,1);
    R_2(ii,1) = x_roh(ii,2);
    

    A_lin = Q(IND);
    b_lin = Paw(IND) - PEEP - E(ii)*V(IND);
    x_lin = A_lin\b_lin;
    R_lin = x_lin;
    R_lin_Array(ii,1) = R_lin;


end

medE = median(E);
medR = median(R_1);


for jj = 1:length(breaths)
    M = 35;
    Paw = breaths(jj).Paw;
    Treal = breaths(jj).Treal;
    V = breaths(jj).V - leak*breaths(jj).T; P = breaths(jj).Paw; Q = breaths(jj).Q-leak;
    T = breaths(jj).T;
    Pes = breaths(jj).Pes;
    Paw = breaths(jj).Paw;
    PEEP = P(1);
    aa2 = find(V == max(V),1,"first");
    V = V - V(end);

    insp_end = insp_end +-0;

    insp = 1:insp_end;
    exp = insp_end+1:length(T);

        plot_ind = find(Q==min(Q),1,'last') ;
        plot_ind_end = find(T==T(end));
%         plot_ind_end = find(Q(plot_ind:end)>0);
        IND = plot_ind:plot_ind_end;
        IND(P(IND)<PEEP-2)=[];


    [t_spline,y_spline] = b_spline_basis_functions(M,2,length((T(IND)))/100);
    Psi = zeros(length(T),M); % \Psi = Splines_Overall_Peff

    if length(y_spline(:,1)) == length(IND) +1
        y_spline = y_spline(1:end-1,:);
    end

    Psi(IND,:) = y_spline;

    A_exp = [Q.*Psi];
    b_exp = [(Paw - PEEP - E(jj)*V)];
    x_exp = lsqnonneg(A_exp,b_exp);

    R_nlin_coeff = x_exp(1:M);
    
    P_model_lin = R_lin*Q + E(jj)*V + PEEP;

    R_roh = R_1(jj) + R_2(jj).*abs(Q);
    P_model_roh = R_roh.*Q + E(jj)*V + PEEP;

    P_model_nlin = [Psi(IND,:)*R_nlin_coeff.*Q(IND) + E(jj)*V(IND) + PEEP];


    R_roh_mean(jj,1) = mean(R_roh);
    R2Phi(jj,1) = mean(Psi(IND,:)*R_nlin_coeff);
    RMSE_lin(jj,1) = rmse(P_model_lin(IND),Paw(IND))./(median(Paw(IND)));
    RMSE_roh(jj,1) = rmse(P_model_roh(IND),Paw(IND))./(median(Paw(IND)));
    RMSE_nlin(jj,1) = rmse(P_model_nlin,Paw(IND))./(median(Paw(IND)));

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

%     r_lin = corrcoef(P_model_roh,Paw);
%     r_nlin = corrcoef(P_model_nlin,Paw);
%     cc_lin(jj,1) = r_lin(2,1);
%     cc_nlin(jj,1) = r_nlin(2,1);
    

%%
    if kk == 6 && jj == 5

        figure('position',[-1900 50 400 300])


%         ind(P(IND)<(2))=[];

        plot(T(IND),Paw(IND),'Color',colorBlack,'Linewidth',2,'Linestyle','-'); hold on
        plot(T(IND),P_model_lin(IND),'Color',colorOrange,'Linewidth',2,'Linestyle','-');
        plot(T(IND),P_model_roh(IND),'Color',colorBlue,'Linewidth',2,'Linestyle','-');
        plot(T(IND),P_model_nlin,'Color',colorRed,'Linewidth',2,'Linestyle','-');
%         plot(T(IND),1*Q(IND),'m')
        % plot(T,Q,'k--')
        grid on
        legend('P_{aw} (measured)','P_{aw} (R_{lin})','P_{aw} (R_{Rohrer})','P_{aw} (R_{nlin})','FontSize',FS,'Location','southeast')
        xlabel('t (s)')
        ylabel('P(t) (cmH_2O)')
%         title(append('Patient ',num2str(kk)))
%         xlim([0 2.5])
        ylim([0 18])
        set(gca,'FontSize',14);
        print('-r400','-dpng',append('CH8_models_P',num2str(kk),'_B',num2str(jj)));
        
        figure('position',[2800 -200 400 300])
        plot(Q(IND),R_lin_Array(jj)*ones(length(IND),1),'Color',colorOrange,'Linewidth',2,'Linestyle','-'); hold on
        plot(Q(IND),R_roh((IND)),'Color',colorBlue,'Linewidth',2,'Linestyle','-'); 
        plot(Q(IND),Psi(IND,:)*R_nlin_coeff,'Color',colorRed,'Linewidth',2,'Linestyle','-');  

        grid on
        xlabel('Q (L/s)')
        ylabel('R (cmH_2O/L/s)')
        legend('nonlinear R','Rohrer R', 'linear R')
%         xlim([-1.1 0])
        ylim([0 30])
                grid on
        ylabel('R_{exp}\Phi(t) (cmH_2O/L/s)')
        str = append('Patient ',num2str(kk));
        title(str)
        set(gca,'FontSize',14);
        print('-r400','-dpng',append('CH8_RQ_P',num2str(kk),'_B',num2str(jj)));
        
        R_const = R_lin_Array(jj)*ones(length(IND));

        figure('Position',[-1100 50 400 300])
        for iii = 1:M
            h(iii) = plot(T(IND),Psi(IND,iii)*R_nlin_coeff(iii)); hold on
        end
        h(iii+1) = plot(T(IND),R_roh(IND),'Color',colorBlue,'Linewidth',2,'Linestyle','-');
        h(iii+2) = plot(T(IND),Psi(IND,:)*R_nlin_coeff,'Color',colorRed,'Linewidth',2,'Linestyle','-'); hold on
        hh = plot(T(IND),R_const,'Color',colorOrange,'Linewidth',2,'Linestyle','-');
        
        plot(T(IND),Psi(IND,:)*R_nlin_coeff,'Color',colorRed,'Linewidth',2,'Linestyle','-'); hold on
        plot(T(IND),R_roh(IND),'Color',colorBlue,'Linewidth',2,'Linestyle','-');
        plot(T(IND),R_const,'Color',colorOrange,'Linewidth',2,'Linestyle','-');
%         legend('nonlinear R','Rohrer R','linear R')
        grid on
        str = append('Patient ',num2str(kk));
        
%         title(str)
        
        %       
        %         legend(([hh(1),h(36:37)]),'nonlinear R','Rohrer R','linear R','FontSize',FS)
        %         legend('A','B','C')
        xlabel('t (s)')
        ylabel('R_{exp}\Phi(t) (cmH_2O/L/s)')
        
%         xlim([0 2.5])
        ylim([0 30])
        set(gca,'FontSize',14);
        print('-r400','-dpng',append('CH8_splines_P',num2str(kk),'_B',num2str(jj)));
    end

%%
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

meanR{:,kk} = R_lin_Array;
meanR1{:,kk} = R_1;
meanR2{:,kk} = R_2;
meanR2Phi{:,kk} = R2Phi;
meanRroh{:,kk} = R_roh_mean;

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
pp_RMSE_roh(kk,1) = median(RMSE_roh);
pp_RMSE_lin(kk,1) = median(RMSE_lin);
pp_RMSE_nlin(kk,1) = median(RMSE_nlin);

pp_q_E(kk,:) = quantile(E,3);
pp_q_R(kk,:) = quantile(R_lin,3);
pp_q_R1(kk,:) = quantile(R_1,3);
pp_q_R2(kk,:) = quantile(R_2,3);
pp_q_R2Phi(kk,:) = quantile(R2Phi,3);

pp_RMSE_roh_q(kk,1:3) = quantile(RMSE_roh,3);
pp_RMSE_lin_q(kk,1:3) = quantile(RMSE_lin,3);
pp_RMSE_nlin_q(kk,1:3) = quantile(RMSE_nlin,3);

% quantile(FL.WOBe,3);

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

[h2.P1,stats.P1] = cdfplot(meanR{:,1}); hold on
[h2.P2,stats.P2] = cdfplot(meanR{:,2}); hold on
[h2.P3,stats.P3] = cdfplot(meanR{:,3}); hold on
[h2.P4,stats.P4] = cdfplot(meanR{:,4}); hold on
[h2.P5,stats.P5] = cdfplot(meanR{:,5}); hold on
[h2.P6,stats.P6] = cdfplot(meanR{:,6}); hold on

set(h2.P1, 'LineStyle', '-', 'Color', colorBlue,'LineWidth',2);
set(h2.P2, 'LineStyle', '-', 'Color', colorGreen,'LineWidth',2);
set(h2.P3, 'LineStyle', '-', 'Color', colorRed,'LineWidth',2);
set(h2.P4, 'LineStyle', '-', 'Color', colorBlack,'LineWidth',2);
set(h2.P5, 'LineStyle', '-', 'Color', colorOrange,'LineWidth',2);
set(h2.P6, 'LineStyle', '-', 'Color', colorYellow,'LineWidth',2);

xlabel('R (cmH_2O/L/s)')
ylabel('Proportion')
xlim([0 35])
title(' ')
legend('P1','P2','P3','P4', 'P5','P6','Location','southeast','Fontsize',FS)
grid on
title(' ')
set(gca,'FontSize',FS);
       
print('-r400','-dpng','CH8_CDF_R.png');

figure('Position',[-1600 50 400 300])

[h2.P1,stats.P1] = cdfplot(meanR1{:,1}); hold on
[h2.P2,stats.P2] = cdfplot(meanR1{:,2}); hold on
[h2.P3,stats.P3] = cdfplot(meanR1{:,3}); hold on
[h2.P4,stats.P4] = cdfplot(meanR1{:,4}); hold on
[h2.P5,stats.P5] = cdfplot(meanR1{:,5}); hold on
[h2.P6,stats.P6] = cdfplot(meanR1{:,6}); hold on

set(h2.P1, 'LineStyle', '-', 'Color', colorBlue,'LineWidth',2);
set(h2.P2, 'LineStyle', '-', 'Color', colorGreen,'LineWidth',2);
set(h2.P3, 'LineStyle', '-', 'Color', colorRed,'LineWidth',2);
set(h2.P4, 'LineStyle', '-', 'Color', colorBlack,'LineWidth',2);
set(h2.P5, 'LineStyle', '-', 'Color', colorOrange,'LineWidth',2);
set(h2.P6, 'LineStyle', '-', 'Color', colorYellow,'LineWidth',2);

xlabel('R_1 (cmH_2O/L/s)')
ylabel('Proportion')
xlim([0 35])
title(' ')
legend('P1','P2','P3','P4', 'P5','P6','Location','southeast','Fontsize',FS)
grid on
title(' ')
set(gca,'FontSize',FS);
       
print('-r400','-dpng','CH8_CDF_R1.png');

figure('Position',[-1600 50 400 300])

[h2.P1,stats.P1] = cdfplot(meanR2{:,1}); hold on
[h2.P2,stats.P2] = cdfplot(meanR2{:,2}); hold on
[h2.P3,stats.P3] = cdfplot(meanR2{:,3}); hold on
[h2.P4,stats.P4] = cdfplot(meanR2{:,4}); hold on
[h2.P5,stats.P5] = cdfplot(meanR2{:,5}); hold on
[h2.P6,stats.P6] = cdfplot(meanR2{:,6}); hold on

set(h2.P1, 'LineStyle', '-', 'Color', colorBlue,'LineWidth',2);
set(h2.P2, 'LineStyle', '-', 'Color', colorGreen,'LineWidth',2);
set(h2.P3, 'LineStyle', '-', 'Color', colorRed,'LineWidth',2);
set(h2.P4, 'LineStyle', '-', 'Color', colorBlack,'LineWidth',2);
set(h2.P5, 'LineStyle', '-', 'Color', colorOrange,'LineWidth',2);
set(h2.P6, 'LineStyle', '-', 'Color', colorYellow,'LineWidth',2);

xlabel('R_2 (cmH_2O/L/s)')
ylabel('Proportion')
xlim([0 35])
title(' ')
legend('P1','P2','P3','P4', 'P5','P6','Location','southeast','Fontsize',FS)
grid on
title(' ')
set(gca,'FontSize',FS);
       
print('-r400','-dpng','CH8_CDF_R2.png');

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

figure('Position',[-500 50 400 300])

sz = 35;
scatter(meanRroh{:,1},meanR2Phi{:,1},sz,...
                'MarkerEdgeColor',colorBlue,...
                'MarkerFaceColor',colorBlue,...
                'LineWidth',1.5); hold on
scatter(meanRroh{:,2},meanR2Phi{:,2},sz,...
                'MarkerEdgeColor',colorGreen,...
                'MarkerFaceColor',colorGreen,...
                'LineWidth',1.5); hold on
scatter(meanRroh{:,3},meanR2Phi{:,3},sz,...
                'MarkerEdgeColor',colorRed,...
                'MarkerFaceColor',colorRed,...
                'LineWidth',1.5); hold on
scatter(meanRroh{:,4},meanR2Phi{:,4},sz,...
                'MarkerEdgeColor',colorBlack,...
                'MarkerFaceColor',colorBlack,...
                'LineWidth',1.5); hold on
scatter(meanRroh{:,5},meanR2Phi{:,5},sz,...
                'MarkerEdgeColor',colorOrange,...
                'MarkerFaceColor',colorOrange,...
                'LineWidth',1.5); hold on
scatter(meanRroh{:,6},meanR2Phi{:,6},sz,...
                'MarkerEdgeColor',colorYellow,...
                'MarkerFaceColor',colorYellow,...
                'LineWidth',1.5); hold on
xlabel('mean R_{Rohrer}')
ylabel('mean R_{exp}\Phi')
legend('P1','P2','P3','P4','P5','P6','Location','southeast')
xlim([0 35])
ylim([0 40])
grid on
set(gca,'FontSize',14);

    print('-r400','-dpng','CH8_scatter.png');

