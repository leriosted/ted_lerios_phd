function [P_model_lin,P_model_R2,P_model_Peff,E,R1_insp,R1_exp,R2,Phi,Peff,Peff_coeff,Psi,...
    RMSE_lin,RMSE_R2,RMSE_Peff,Ave_R2_Phi,Peak_R2_Phi,Ave_Peff,Max_Peff,AUC_Peff,R1_insp_Array,R1_exp_Array,R2Array]...
    = model_function_2(t_P,P,t_Q,Q,V,P0,exp_0,exp_max)
% MODEL_FUNCTION processes each patient one at a time
%   Patient data is entered to create a model
%% Set number of splines
M = 10;

%% Interpolate to increase resolution
% t_P_long = linspace(0,t_P(end),300)';
% P = interp1(t_P,P,t_P_long,'pchip');
% t_P = t_P_long;
% t_Q_long = linspace(0,t_Q(end),300)';
% Q = interp1(t_Q,Q,t_Q_long,'pchip');
% V = interp1(t_Q,V,t_Q_long,'pchip');
% t_Q = t_Q_long;

%% Identify R for linear model:
A_lin = [Q]; b_lin = [-P];
x_lin = A_lin\b_lin;
R_lin = x_lin(1);

%% Identify R1 Insp:
insp = 1:exp_0;
lb = [0];
ub = [inf];
A_insp = [Q(insp)]; b_insp = [-P(insp)]; 
x_insp = lsqlin(A_insp,b_insp,[],[],[],[],lb,ub);
R1_insp = x_insp(1);
%% Identify R1 Exp:
exp = exp_0+1:length(t_P);
lb = [0];
ub = [inf];
A_exp = [Q(exp)]; b_exp = [-P(exp)]; 
x_exp = lsqlin(A_exp,b_exp,[],[],[],[],lb,ub);
R1_exp = x_exp(1);

%% All Model:
P_model_lin = -1*(R1_insp*Q +P0);
% Check sizes here and if they aren't the same change in an if statement...
d = 2;
kw = 0.1;

[t_spline,y_spline_Phi] = b_spline_basis_functions(M,2,(length(t_P)-(exp_0-1))/100);

Phi = zeros(length(t_P),M); % \Phi = Splines_Overall_R2




Phi(exp_0:end,:) = y_spline_Phi;
lb = [0];
ub = [inf];
A_R2 = [Q.*Phi];
b_R2 = -P;
x_R2 = lsqlin(A_R2,b_R2,[],[],[],[],lb,ub);
R2 = x_R2(1:M);

%% P_eff

[t_spline,y_spline_Peff] = b_spline_basis_functions(M,2,length((t_P))/100);
Psi = zeros(length(t_P),M); % \Psi = Splines_Overall_Peff
Psi = y_spline_Peff(1:end-1,:);
lb = [0;0];
ub = [inf;inf];

% ind = [1:exp-1,exp_max-5:exp_max+5];
% ind = find(Q==min(Q));
% ind = (ind-5):(ind+5);

aa = ceil(length(Q)*3/4);

aa =(aa):(aa+5);% aa:length(data.Q);

E= lsqnonneg(V(aa)- V(1), P(aa));

% E = lsqnonneg(V(ind)-V(1),-P(ind));

Peff_coeff = lsqlin(Psi, P-E*(V-V(1)));
Peff = Psi(:,:)*Peff_coeff;
P_model_Peff = E*(V-V(1)) - Peff;

P_model_R2 = [R1_insp*Q(insp); Phi(exp,:)*R2.*Q(exp)];
P_model_R2 = -1*P_model_R2;


%% RMSE
RMSE_R2 = sqrt(mean((P - P_model_R2).^2));
RMSE_Peff = sqrt(mean((P - P_model_Peff).^2));
RMSE_lin = sqrt(mean((P - P_model_lin).^2));

%% Averages
Ave_R2_Phi = mean(Phi(exp,:)*R2);
Peak_R2_Phi = max(Phi(exp,:)*R2);
Ave_Peff = mean(Psi*Peff_coeff);
Max_Peff = max(abs(Peff));

%% P_effort
AUC_Peff = sum(cumtrapz(Peff));

%% Arrays for plotting
R1_insp_Array = R1_insp.*ones(length(t_P),1);
R1_exp_Array = R1_exp.*ones(length(t_P),1);
R2Array = Phi*R2;


end

