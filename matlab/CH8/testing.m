clear all; close all; clc
colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
pt = 1.0;
FS = 8;
load('Patient_1.mat')
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

    subplot(6,5,ii)
%     plot(T,V,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
    plot(T,Pes-PEEP,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
    plot(T,Peff-PEEP,'Color',colorRed,'Linewidth',pt,'Linestyle','-');
    plot(T,Paw,'Color',colorBlue,'Linewidth',0.6,'Linestyle','-');
%     plot(T,Q,'Color',colorBlack,'Linewidth',0.4,'Linestyle',':')
    
    xlim([0 10])
    ylim([-10 20])
    grid on
    title(['E: ',num2str(E),' | R: ',num2str(R)])
    set(gca,'FontSize',FS);
end
    set(gca,'FontSize',FS);
    leg1 = legend('P_{es}','P_{eff}','P_{aw}','FontSize',12,...
        'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.68 .185];

    print('-r400','-dpng','Greek_all.png');


    figure('Position',[-1900 150 600 600])
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
    
        

        if ii == 5 | ii == 8 | ii == 11 | ii == 14 | ii == 18 | ii == 27
            
            
            subplot(4,2,counter)
            counter = counter + 1;
            
            plot(T,Pes-PEEP,'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
            plot(T,Peff-PEEP,'Color',colorRed,'Linewidth',pt,'Linestyle','-');
            plot(T,Paw,'Color',colorBlue,'Linewidth',0.6,'Linestyle','-');
            xlim([0 8])
            ylim([-10 22])
            grid on
            title(['E: ',num2str(E),' | R: ',num2str(R)])
            set(gca,'FontSize',FS);
        end
    end
        set(gca,'FontSize',FS);
        leg1 = legend('P_{es}','P_{eff}','P_{aw}','FontSize',12,...
            'Orientation','horizontal');
        legend boxoff
        leg1.Position(1:2) = [.52 .185];
    
        print('-r400','-dpng','Greek_typical.png');
