close all; clear all; clc

load young_results.mat
load old_results.mat
load FL_results.mat
load NFL_results.mat

VarNames = ["Patient","Catagory","R1insp","R1exp","linear RMSE","mean R2Phi","nonlinear RMSE"];

patientNo(:,1) = 1:100;

catagory = cell(100,1);
catagory(1:20) = {'Y'};
catagory(21:40) = {'O'};
catagory(41:65) = {'N'};
catagory(66:100) = {'F'};

R1insp(1:20,1) = round(young.R1insp,1);
R1insp(21:40,1) = round(old.R1insp,1);
R1insp(41:65,1) = round(NFL.R1insp,1);
R1insp(66:100,1) = round(FL.R1insp,1);

R1exp(1:20,1) = round(young.R1exp,1);
R1exp(21:40,1) = round(old.R1exp,1);
R1exp(41:65,1) = round(NFL.R1exp,1);
R1exp(66:100,1) = round(FL.R1exp,1);

RMSE_linear(1:20,1) = round(young.RMSE_linear,1);
RMSE_linear(21:40,1) = round(old.RMSE_linear,1);
RMSE_linear(41:65,1) = round(NFL.RMSE_linear,1);
RMSE_linear(66:100,1) = round(FL.RMSE_linear,1);

meanR2Phi(1:20,1) = round(young.meanR2Phi,1);
meanR2Phi(21:40,1) = round(old.meanR2Phi,1);
meanR2Phi(41:65,1) = round(NFL.meanR2Phi,1);
meanR2Phi(66:100,1) = round(FL.meanR2Phi,1);

RMSE_nonlinear(1:20,1) = round(young.RMSE_nonlinear,1);
RMSE_nonlinear(21:40,1) = round(old.RMSE_nonlinear,1);
RMSE_nonlinear(41:65,1) = round(NFL.RMSE_nonlinear,1);
RMSE_nonlinear(66:100,1) = round(FL.RMSE_nonlinear,1);

T = table(patientNo,catagory,R1insp,R1exp,RMSE_linear,meanR2Phi,RMSE_nonlinear,'VariableNames',VarNames);
disp(T)

writetable(T,'tableResults.csv')


