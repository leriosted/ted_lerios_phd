clear all; close all; clc

atWork = 0;
print_on = 0;

save_young = 0;
save_old   = 0;
save_NFL   = 0;
save_FL    = 0;

colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
pt = 1.0;
FS = 8;

%% Plot Setup:
if atWork == 1
    figPos1 = [-1900, 100, 1100, 800]; % Work
    figPos2 = [-1250, 400, 600, 600];
    figPos3 = [-600, 400, 600, 600];
    figPos4 = [-1900, 100, 600, 600];
    figPos5 = [-1250, 100, 600, 600];
    figPos6 = [-600, 100, 600, 300];
else
    figPos1 = [-1900, 150, 1100, 800]; % Home
    figPos2 = [-1250, 500, 600, 600];
    figPos3 = [-600, 500, 600, 600];
    figPos4 = [-1900, 150, 600, 600];
    figPos5 = [-1250, 150, 600, 600];
    figPos6 = [-600, 150, 600, 300];
end

%% Select a folder of data
dataCatagory = 'FL';
files = dir(fullfile(dataCatagory,'*.mat'));

%% Preallocation of Results Arrays


%% Process data (Calculate volume and run model)
for ii = 1:length(files)
    
    thisFileName = fullfile(dataCatagory,files(ii).name);
    fprintf('Processing "%s".\n',thisFileName);
    load(thisFileName)

    data.V_int = cumtrapz(data.T,data.Q);

    [val,ind] = max(cumtrapz(data.Q));
    exp_point = ind;
    [val,ind2] = min(data.Q);
    exp_max = ind2;

%     [P_model_lin,P_model_R2,P_model_Peff,E,R1_insp,R1_exp,R2,phi,Peff,Peff_coeff,psi,...
%         RMSE_lin,RMSE_R2,RMSE_Peff,Ave_R2_Phi,Peak_R2_Phi,Ave_Peff,Max_Peff,AUC_Peff...
%         ,R1_insp_Array,R1_exp_Array,R2Array]...
%         = model_function(data.T,data.Pa,data.T,data.Q,data.V,0,exp_point,exp_max);
   
%     [P_model_lin,P_model_R2,P_model_Peff,E,R1_insp,R1_exp,R2,phi,Peff,Peff_coeff,psi,...
%     RMSE_lin,RMSE_R2,RMSE_Peff,Ave_R2_Phi,Peak_R2_Phi,Ave_Peff,Max_Peff,AUC_Peff...
%     ,R1_insp_Array,R1_exp_Array,R2Array]...
%     = model_function_BACKUP(data.T,data.Pa,data.T,data.Q,data.V,0,exp_point,exp_max);

        [P_model_lin,P_model_R2,P_model_Peff,E,R1_insp,R1_exp,R2,phi,Peff,Peff_coeff,psi,...
    RMSE_lin,RMSE_R2,RMSE_Peff,Ave_R2_Phi,Peak_R2_Phi,Ave_Peff,Max_Peff,AUC_Peff...
    ,R1_insp_Array,R1_exp_Array,R2Array]...
    = model_function_2(data.T,data.Pa,data.T,data.Q,data.V,0,exp_point,exp_max);
    
    %% Store results
    % Linear Data and Model
    pressure(ii,:) = data.Pa;
    flow(ii,:) = data.Q;
    volume(ii,:) = data.V_int;
    linear_model(ii,:) = P_model_lin;
    lin_R1_insp_Array(:,ii) = R1_insp;
    lin_R1_exp_Array(:,ii) = R1_exp;
    lin_RMSE_Array(:,ii) = RMSE_lin;

    nonlinear_model(:,ii) = P_model_R2;
    nlin_R2_phi_Array(:,ii) = phi*R2;
    nlin_RMSE_Array(:,ii) = RMSE_R2;

%     Peff_model = zeros(length(data.T),length(files));
    Peff_model{ii,:} = P_model_Peff;
    Peff_Array{:,ii} = Peff;
    mean_Peff_Array(:,ii) = Ave_Peff;
    Peff_min(:,ii) = min(Peff);

    WOB_dV_Array{:,ii} = cumtrapz(data.V_int(ind2),Peff);
    WOB_dV{:,ii} = WOB_dV_Array{end,ii};

    WOBr_dV_Array{:,ii} = cumtrapz(data.V_int(1:exp_point),R1_insp*data.Q(1:exp_point));
    WOBr_dV{:,ii} = WOBr_dV_Array{end,ii};

    WOBe_dV_Array{:,ii} = cumtrapz(data.V_int(1:exp_point),E*data.V_int(1:exp_point));
    WOBe_dV{:,ii} = WOBe_dV_Array{end,ii};

    E_Array(:,ii) = E;

    resistance(ii,:) = R1_insp_Array;
    resistance(ii,exp_point:end) = R2Array(exp_point:end);
    exp_point_Array(ii,:) = exp_point;

    RMSE_linear(ii,:) = RMSE_lin;
    RMSE_nonlinear(ii,:) = RMSE_R2;
    RMSE_Peff(ii,:) = RMSE_Peff;

    Ave_R2_Phi_Array(:,ii) = Ave_R2_Phi;
    Int_R2_Phi_dT = cumtrapz(data.T,phi*R2);
    Int_R2_Phi_dT_Array(:,ii) = Int_R2_Phi_dT(end);
    Peak_R2_Phi_Array(:,ii) = Peak_R2_Phi;
    int_Q_dR = cumtrapz(data.Q,R2Array);
    int_Q_dR_Array(:,ii) = int_Q_dR(end);

    file_name_list{:,ii} = thisFileName;

    % Slope of RQ
    aa = find(data.Q==min(data.Q));
    ind = aa:length(data.Q);
    R_sub = R2Array(ind);
    Q_sub = data.Q(ind);
    slope = lsqr(Q_sub - min(Q_sub),R_sub-min(R_sub));
    slope_Array(:,ii) = slope;

    work_resistance_dV_Array(:,ii) = cumtrapz(data.V_int,(R2Array.*data.Q));
    work_resistance_dV(:,ii) = work_resistance_dV_Array(end,ii);

    % Export Arrays
    


    %% Plots
%     figure('Position',figPos1)
% 
%     subplot(5,3,[1 2])
%     plot(data.T,data.Pa,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
%     plot(data.T,P_model_lin,'Color',colorBlack,'Linewidth',pt,'Linestyle','--'); hold on
%     plot(data.T,P_model_R2,'Color',colorRed,'Linewidth',pt,'Linestyle','-.'); hold on
%     xlabel('t (s)')
%     ylabel('P_{alv} (cmH_2O)')
%     xlim([0 6])
%     ylim([-15 15])
%     leg1 = legend('data','P_{linear model}','P_{non-linear model}','FontSize',12,...
%         'Orientation','horizontal');
%     legend boxoff
%     leg1.Position(1:2) = [.28 .180];
%     grid on
%     set(gca,'FontSize',FS);
% 
%     subplot(5,3,3)
%     plot(data.Pa,data.Q,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
%     plot(P_model_lin,data.Q,'Color',colorBlack,'Linewidth',pt,'Linestyle','--'); hold on
%     plot(P_model_R2,data.Q,'Color',colorRed,'Linewidth',pt,'Linestyle','-.'); hold on
%     xlabel('P (cmH_2O)')
%     ylabel('Q (L/s)')
%     xlim([-15 15])
%     ylim([-1.5 1.5])
% %     title([thisFileName,'_',num2str(ii)])
%     grid on
%     set(gca,'FontSize',FS);
% 
%     subplot(5,3,[4 5])
%     plot(data.T,data.Q,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
%     xlabel('t (s)')
%     ylabel('Flow (Q) [L/s]')
%     xlim([0 6])
%     ylim([-1.5 1.5])
%     grid on
%     set(gca,'FontSize',FS);
% 
%     subplot(5,3,6)
%     plot(data.Pa,data.V_int,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
%     plot(P_model_lin,data.V_int,'Color',colorBlack,'Linewidth',pt,'Linestyle','--'); hold on
%     plot(P_model_R2,data.V_int,'Color',colorRed,'Linewidth',pt,'Linestyle','-.'); hold on
%     xlabel('P (cmH_2O)')
%     ylabel('V (L)')
%     xlim([-15 15])
%     ylim([0 1.5])
%     grid on
%     set(gca,'FontSize',FS);
% 
%     subplot(5,3,[7 8])
%     plot(data.T,data.V_int,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
% %     plot(data.T,data.V-data.V(1),'Color',colorBlack,'Linewidth',pt,'Linestyle','--'); hold on
%     xlabel('t (s)')
%     ylabel('Volume (V) [L]')
%     xlim([0 6])
%     ylim([0 1.5])
%     grid on
%     set(gca,'FontSize',FS);
% 
%     subplot(5,3,9)
%     plot(data.Q,R2Array,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
%     xlabel('Q (L/s)')
%     ylabel('R (cmH_2O/L/s)')
%     xlim([-1.5 0])
%     ylim([0 50])
%     grid on
%     set(gca,'FontSize',FS);
% 
%     subplot(5,3,[10 11])
%     plot(data.T,R1_insp_Array,'Color',colorBlue,'Linewidth',pt,'Linestyle','--'); hold on
%     plot(data.T,R2Array,'Color',colorRed,'Linewidth',pt,'Linestyle','-'); hold on
%     xlabel('t (s)')
%     ylabel('Resistance (cmH_2O/L/s) [L]')
%     xlim([0 6])
%     ylim([0 50])
%     grid on
%     set(gca,'FontSize',FS);
% 
%     subplot(5,3,12)
%     plot(data.V_int,R2Array,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
%     xlabel('V (L)')
%     ylabel('R (cmH_2O/L/s)')
%     xlim([0 1])
%     ylim([0 50])
%     grid on
%     set(gca,'FontSize',FS);
%     
%     sg = sgtitle([dataCatagory,' ',num2str(ii),...
%         '$ \; \;  R_{1,insp}$ = ',num2str(round(R1_insp,1)),...
%         '$ \; \; R_{1,exp} $ = ',num2str(round(R1_exp,1)),...
%         '$\; \;$ mean $R_{2}\Phi$ = ',num2str(round(Ave_R2_Phi,1)),...
%         '$\; \;$ mean $P_{eff}\Psi$ = ',num2str(round(Ave_Peff,1)),...
%         '$ \; \; \frac{dQ}{dR}$ = ',num2str(round(slope,1))]);
%     sg.Interpreter = 'latex';   
% 
%     fig_savename = ([dataCatagory,'_',num2str(ii),'.png']);
%     if print_on == 1
%         print('-r400','-dpng',fig_savename);
%     else
%     end
%     figure('Position',figPos1)
% 
%     plot(data.T,data.Pa,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
%     plot(data.T(1:exp_point-1),P_model_Peff,'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
% %     plot(data.T,P_model_lin,'Color',colorBlack,'Linewidth',pt,'Linestyle','--'); hold on
% %     plot(data.T,P_model_R2,'Color',colorRed,'Linewidth',pt,'Linestyle','-.'); hold on
%     xlabel('t (s)')
%     ylabel('P_{alv} (cmH_2O)')
%     xlim([0 6])
%     ylim([-15 15])
%     sg = sgtitle([dataCatagory,' ',num2str(ii),...
%         '$ \; \;  R_{1,insp}$ = ',num2str(round(R1_insp,1)),...
%         '$ \; \; R_{1,exp} $ = ',num2str(round(R1_exp,1)),...
%         '$\; \;$ mean $R_{2}\Phi$ = ',num2str(round(Ave_R2_Phi,1)),...
%         '$\; \;$ mean $P_{eff}\Psi$ = ',num2str(round(Ave_Peff,1)),...
%         '$ \; \; \frac{dQ}{dR}$ = ',num2str(round(slope,1))]);
%     sg.Interpreter = 'latex'; 
% %     leg1 = legend('data','P_{linear model}','P_{non-linear model}','FontSize',12,...
% %         'Orientation','horizontal');
% %     legend boxoff
% %     leg1.Position(1:2) = [.28 .180];
%     legend('measured','simulated')
%     grid on
%     set(gca,'FontSize',14);
%     fig_savename = ([dataCatagory,'_Peff_',num2str(ii),'.png']);
%     if print_on == 1
%         print('-r400','-dpng',fig_savename);
%     else
%     end

    figure('Position',[-1900, 150, 400, 400])

    plot(data.T,data.Pa,'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
    plot(data.T,Peff,'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
%     plot(data.T(1:exp_point-1),Peff,'Color',colorRed,'Linewidth',pt,'Linestyle','-.')
    plot(data.T(1:exp_point-1),-R1_insp*data.Q(1:exp_point-1),'Color',colorBlack,'Linewidth',pt,'Linestyle','--')
    plot(data.T(exp_point:end),-R2Array(exp_point:end).*data.Q(exp_point:end),'Color',colorBlue,'Linewidth',pt,'Linestyle','-.')
    xlabel('t (s)')
    ylabel('P_{alv} (cmH_2O)')
    xlim([0 6])
    ylim([-15 15])
%     sg = sgtitle([dataCatagory,' ',num2str(ii),...
%         '$ \; \;  R_{1,insp}$ = ',num2str(round(R1_insp,1)),...
%         '$ \; \; R_{1,exp} $ = ',num2str(round(R1_exp,1)),...
%         '$\; \;$ mean $R_{2}\Phi$ = ',num2str(round(Ave_R2_Phi,1)),...
%         '$\; \;$ mean $P_{eff}\Psi$ = ',num2str(round(Ave_Peff,1)),...
%         '$ \; \; \frac{dQ}{dR}$ = ',num2str(round(slope,1))]);
    sg.Interpreter = 'latex'; 
    legend('P_{alv}','P_{eff}','R_{1,insp} Q','R_{2} \Phi Q')
    grid on
    set(gca,'FontSize',14);
    fig_savename = ([dataCatagory,'_P_',num2str(ii),'.png']);
    if print_on == 1
        print('-r400','-dpng',fig_savename);
    else
    end
end


%% Save results in struct form

if save_young == 1
    young.time = data.T';
    young.pressure = pressure;
    young.flow = flow;
    young.volume = volume;
    young.linear_model = linear_model;
    young.nonlinear_model = nonlinear_model;
    young.RMSE_linear = RMSE_linear;
    young.RMSE_nonlinear = RMSE_nonlinear;
    young.resistance = resistance;
    young.exp_point = exp_point_Array;
    young.meanR2Phi = Ave_R2_Phi_Array;
    young.Int_R2_Phi_dT_Array = Int_R2_Phi_dT_Array;
    young.Peak_R2_Phi_Array = Peak_R2_Phi_Array;
    young.int_Q_dR_Array = int_Q_dR_Array;
    young.R1insp = lin_R1_insp_Array;
    young.R1exp = lin_R1_exp_Array;
    young.slope_Array = slope_Array;
    young.work_resistance_dV = work_resistance_dV;
    young.work_resistance_dV_Array = work_resistance_dV_Array;
    young.Peff_Array = Peff_Array;
    young.WOB_dV = WOB_dV;
    young.E_Array = E_Array;
    young.WOBe_dV = WOBe_dV;
    young.WOBr_dV = WOBr_dV;
    young.Peff_min = Peff_min;
    save('young_results.mat','young')
end

if save_old == 1
    old.time = data.T';
    old.pressure = pressure;
    old.flow = flow;
    old.volume = volume;
    old.linear_model = linear_model;
    old.nonlinear_model = nonlinear_model;
    old.RMSE_linear = RMSE_linear;
    old.RMSE_nonlinear = RMSE_nonlinear;
    old.resistance = resistance;
    old.exp_point = exp_point_Array;
    old.meanR2Phi = Ave_R2_Phi_Array;
    old.Int_R2_Phi_dT_Array = Int_R2_Phi_dT_Array;
    old.int_Q_dR_Array = int_Q_dR_Array;
    old.Peak_R2_Phi_Array = Peak_R2_Phi_Array;
    old.R1insp = lin_R1_insp_Array;
    old.R1exp = lin_R1_exp_Array;
    old.slope_Array = slope_Array;
    old.work_resistance_dV = work_resistance_dV;
    old.work_resistance_dV_Array = work_resistance_dV_Array;
    old.Peff_Array = Peff_Array;
    old.WOB_dV = WOB_dV;
    old.E_Array = E_Array;
    old.WOBe_dV = WOBe_dV;
    old.WOBr_dV = WOBr_dV;
    old.Peff_min = Peff_min;
    save('old_results.mat','old')
end

if save_NFL == 1
    NFL.time = data.T';
    NFL.time = data.T;
    NFL.pressure = pressure;
    NFL.flow = flow;
    NFL.volume = volume;
    NFL.linear_model = linear_model;
    NFL.nonlinear_model = nonlinear_model;
    NFL.RMSE_linear = RMSE_linear;
    NFL.RMSE_nonlinear = RMSE_nonlinear;
    NFL.resistance = resistance;
    NFL.exp_point = exp_point_Array;
    NFL.meanR2Phi = Ave_R2_Phi_Array;
    NFL.Int_R2_Phi_dT_Array = Int_R2_Phi_dT_Array;
    NFL.Peak_R2_Phi_Array = Peak_R2_Phi_Array;
    NFL.int_Q_dR_Array = int_Q_dR_Array;
    NFL.R1insp = lin_R1_insp_Array;
    NFL.R1exp = lin_R1_exp_Array;
    NFL.slope_Array = slope_Array;
    NFL.work_resistance_dV = work_resistance_dV;
    NFL.work_resistance_dV_Array = work_resistance_dV_Array;
    NFL.Peff_Array = Peff_Array;
    NFL.WOB_dV = WOB_dV;
    NFL.E_Array = E_Array;
    NFL.WOBe_dV = WOBe_dV;
    NFL.WOBr_dV = WOBr_dV;
    NFL.Peff_min = Peff_min;
    save('NFL_results.mat','NFL')
end

if save_FL == 1
    FL.time = data.T';
    FL.time = data.T;
    FL.pressure = pressure;
    FL.flow = flow;
    FL.volume = volume;
    FL.linear_model = linear_model;
    FL.nonlinear_model = nonlinear_model;
    FL.RMSE_linear = RMSE_linear;
    FL.RMSE_nonlinear = RMSE_nonlinear;
    FL.resistance = resistance;
    FL.exp_point = exp_point_Array;
    FL.meanR2Phi = Ave_R2_Phi_Array;
    FL.Int_R2_Phi_dT_Array = Int_R2_Phi_dT_Array;
    FL.Peak_R2_Phi_Array = Peak_R2_Phi_Array;
    FL.int_Q_dR_Array = int_Q_dR_Array;
    FL.R1insp = lin_R1_insp_Array;
    FL.R1exp = lin_R1_exp_Array;
    FL.slope_Array = slope_Array;
    FL.work_resistance_dV = work_resistance_dV;
    FL.work_resistance_dV_Array = work_resistance_dV_Array;
    FL.Peff_Array = Peff_Array;
    FL.WOB_dV = WOB_dV;
    FL.E_Array = E_Array;
    FL.WOBe_dV = WOBe_dV;
    FL.WOBr_dV = WOBr_dV;
    FL.Peff_min = Peff_min;
    save('FL_results.mat','FL')
end

% figure('Position',figPos1)
% % c = linspace(1,10,length(Ave_R2_Phi_Array));
% % scatter(Ave_R2_Phi_Array,lin_R1_insp_Array,[],c)
% 
% sz = 40;
% scatter(lin_R1_insp_Array,res.FL.R2Phi,sz,...
%                 'MarkerEdgeColor',[0 0 1],...
%                 'MarkerFaceColor',[0 0 1],...
%                 'LineWidth',1.5); hold on
% scatter(lin_R1_exp_Array,Ave_R2_Phi_Array,sz,...
%                 'MarkerEdgeColor',[1 0 0],...
%                 'MarkerFaceColor',[1 0 0],...
%                 'LineWidth',1.5); hold on
% xlabel('R')
% ylabel('mean R_2\Phi')
% legend('R_{1,exp}','R_{2,\Phi}')
% xlim([0 50])
% ylim([0 50])
% grid on
% set(gca,'FontSize',14);



% T = table(array2tablelin_R1_insp_Array,lin_R1_exp_Array,lin_RMSE_Array,...
%     nlin_R2_phi_Array,nlin_RMSE_Array,'VariableNames',{'R1insp','R1exp','linRMSE','R2Phi','nlinRMSE'});
% disp(T)
% writetable(T,'table_1.csv')
