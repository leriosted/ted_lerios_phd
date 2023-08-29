close all; clear all; clc

load young_results.mat
load old_results.mat
load FL_results.mat
load NFL_results.mat

for ii = 1:20
    young.WOB(ii,1) = young.WOB_dV{1, ii}(end,1);
    young.WOBr(ii,1) = young.WOBr_dV{1, ii}(end,1);
    young.WOBe(ii,1) = young.WOBe_dV{1, ii}(end,1);
end

for ii = 1:20
    old.WOB(ii,1) = old.WOB_dV{1, ii}(end,1);
    old.WOBr(ii,1) = old.WOBr_dV{1, ii}(end,1);
    old.WOBe(ii,1) = old.WOBe_dV{1, ii}(end,1);
end

for ii = 1:25
    NFL.WOB(ii,1) = NFL.WOB_dV{1, ii}(end,1);
    NFL.WOBr(ii,1) = NFL.WOBr_dV{1, ii}(end,1);
    NFL.WOBe(ii,1) = NFL.WOBe_dV{1, ii}(end,1);
end

for ii = 1:35
    FL.WOB(ii,1) = FL.WOB_dV{1, ii}(end,1);
    FL.WOBr(ii,1) = FL.WOBr_dV{1, ii}(end,1);
    FL.WOBe(ii,1) = FL.WOBe_dV{1, ii}(end,1);
end

young.R1insp_quartiles = quantile(young.R1insp,3);
old.R1insp_quartiles = quantile(old.R1insp,3);
NFL.R1insp_quartiles = quantile(NFL.R1insp,3);
FL.R1insp_quartiles = quantile(FL.R1insp,3);

young.meanR2Phi_quartiles = quantile(young.meanR2Phi,3);
old.meanR2Phi_quartiles = quantile(old.meanR2Phi,3);
NFL.meanR2Phi_quartiles = quantile(NFL.meanR2Phi,3);
FL.meanR2Phi_quartiles = quantile(FL.meanR2Phi,3);

young.E_Array_quartiles = quantile(young.E_Array,3);
old.E_Array_quartiles = quantile(old.E_Array,3);
NFL.E_Array_quartiles = quantile(NFL.E_Array,3);
FL.E_Array_quartiles = quantile(FL.E_Array,3);

young.Peff_min_quartiles = quantile(young.Peff_min,3);
old.Peff_min_quartiles = quantile(old.Peff_min,3);
NFL.Peff_min_quartiles = quantile(NFL.Peff_min,3);
FL.Peff_min_quartiles = quantile(FL.Peff_min,3);

young.WOB_Peff_quartiles = quantile(young.WOB,3);
old.WOB_Peff_quartiles = quantile(old.WOB,3);
NFL.WOB_Peff_quartiles = quantile(NFL.WOB,3);
FL.WOB_Peff_quartiles = quantile(FL.WOB,3);

young.WOBr_quartiles = quantile(young.WOBr,3);
old.WOBr_quartiles = quantile(old.WOBr,3);
NFL.WOBr_quartiles = quantile(NFL.WOBr,3);
FL.WOBr_quartiles = quantile(FL.WOBr,3);

young.WOBe_quartiles = quantile(young.WOBe,3);
old.WOBe_quartiles = quantile(old.WOBe,3);
NFL.WOBe_quartiles = quantile(NFL.WOBe,3);
FL.WOBe_quartiles = quantile(FL.WOBe,3);

VarNames = ["Patient Catagory","R1insp MED [IQR]","mean R2Phi MED [IQR]","E MED [IQR]",...
    "Peff_min MED [IQR]","WOB_Peff MED [IQR]","WOB_r MED [IQR]","WOB_e MED [IQR]","VT MED [IQR]","RMSE R1insp"];

catagory = cell(4,1);
catagory(1) = {'Y'};
catagory(2) = {'O'};
catagory(3) = {'N'};
catagory(4) = {'F'};

R1insp = cell(4,1);
R1insp(1) = {[num2str(round(median(young.R1insp),1)),' [',num2str(round(young.R1insp_quartiles(1,1),1)),'-',num2str(round(young.R1insp_quartiles(1,3),1)),']']};
R1insp(2) = {[num2str(round(median(old.R1insp),1)),' [',num2str(round(old.R1insp_quartiles(1,1),1)),'-',num2str(round(old.R1insp_quartiles(1,3),1)),']']};
R1insp(3) = {[num2str(round(median(NFL.R1insp),1)),' [',num2str(round(NFL.R1insp_quartiles(1,1),1)),'-',num2str(round(NFL.R1insp_quartiles(1,3),1)),']']};
R1insp(4) = {[num2str(round(median(FL.R1insp),1)),' [',num2str(round(FL.R1insp_quartiles(1,1),1)),'-',num2str(round(FL.R1insp_quartiles(1,3),1)),']']};

meanR2Phi = cell(4,1);
meanR2Phi(1) = {[num2str(round(median(young.meanR2Phi),1)),' [',num2str(round(young.meanR2Phi_quartiles(1,1),1)),'-',num2str(round(young.meanR2Phi_quartiles(1,3),1)),']']};
meanR2Phi(2) = {[num2str(round(median(old.meanR2Phi),1)),' [',num2str(round(old.meanR2Phi_quartiles(1,1),1)),'-',num2str(round(old.meanR2Phi_quartiles(1,3),1)),']']};
meanR2Phi(3) = {[num2str(round(median(NFL.meanR2Phi),1)),' [',num2str(round(NFL.meanR2Phi_quartiles(1,1),1)),'-',num2str(round(NFL.meanR2Phi_quartiles(1,3),1)),']']};
meanR2Phi(4) = {[num2str(round(median(FL.meanR2Phi),1)),' [',num2str(round(FL.meanR2Phi_quartiles(1,1),1)),'-',num2str(round(FL.meanR2Phi_quartiles(1,3),1)),']']};

E = cell(4,1);
E(1) = {[num2str(round(median(young.E_Array),1)),' [',num2str(round(young.E_Array_quartiles(1,1),1)),'-',num2str(round(young.E_Array_quartiles(1,3),1)),']']};
E(2) = {[num2str(round(median(old.E_Array),1)),' [',num2str(round(old.E_Array_quartiles(1,1),1)),'-',num2str(round(old.E_Array_quartiles(1,3),1)),']']};
E(3) = {[num2str(round(median(NFL.E_Array),1)),' [',num2str(round(NFL.E_Array_quartiles(1,1),1)),'-',num2str(round(NFL.E_Array_quartiles(1,3),1)),']']};
E(4) = {[num2str(round(median(FL.E_Array),1)),' [',num2str(round(FL.E_Array_quartiles(1,1),1)),'-',num2str(round(FL.E_Array_quartiles(1,3),1)),']']};

Peff_min = cell(4,1);
Peff_min(1) = {[num2str(round(median(young.Peff_min),1)),' [',num2str(round(young.Peff_min_quartiles(1,1),1)),'-',num2str(round(young.Peff_min_quartiles(1,3),1)),']']};
Peff_min(2) = {[num2str(round(median(old.Peff_min),1)),' [',num2str(round(old.Peff_min_quartiles(1,1),1)),'-',num2str(round(old.Peff_min_quartiles(1,3),1)),']']};
Peff_min(3) = {[num2str(round(median(NFL.Peff_min),1)),' [',num2str(round(NFL.Peff_min_quartiles(1,1),1)),'-',num2str(round(NFL.Peff_min_quartiles(1,3),1)),']']};
Peff_min(4) = {[num2str(round(median(FL.Peff_min),1)),' [',num2str(round(FL.Peff_min_quartiles(1,1),1)),'-',num2str(round(FL.Peff_min_quartiles(1,3),1)),']']};

WOB_Peff = cell(4,1);
WOB_Peff(1) = {[num2str(round(median(young.WOB),1)),' [',num2str(round(young.WOB_Peff_quartiles(1,1),1)),'-',num2str(round(young.WOB_Peff_quartiles(1,3),1)),']']};
WOB_Peff(2) = {[num2str(round(median(old.WOB),1)),' [',num2str(round(old.WOB_Peff_quartiles(1,1),1)),'-',num2str(round(old.WOB_Peff_quartiles(1,3),1)),']']};
WOB_Peff(3) = {[num2str(round(median(NFL.WOB),1)),' [',num2str(round(NFL.WOB_Peff_quartiles(1,1),1)),'-',num2str(round(NFL.WOB_Peff_quartiles(1,3),1)),']']};
WOB_Peff(4) = {[num2str(round(median(FL.WOB),1)),' [',num2str(round(FL.WOB_Peff_quartiles(1,1),1)),'-',num2str(round(FL.WOB_Peff_quartiles(1,3),1)),']']};

WOBr = cell(4,1);
WOBr(1) = {[num2str(round(median(young.WOBr),1)),' [',num2str(round(young.WOBr_quartiles(1,1),1)),'-',num2str(round(young.WOBr_quartiles(1,3),1)),']']};
WOBr(2) = {[num2str(round(median(old.WOBr),1)),' [',num2str(round(old.WOBr_quartiles(1,1),1)),'-',num2str(round(old.WOBr_quartiles(1,3),1)),']']};
WOBr(3) = {[num2str(round(median(NFL.WOBr),1)),' [',num2str(round(NFL.WOBr_quartiles(1,1),1)),'-',num2str(round(NFL.WOBr_quartiles(1,3),1)),']']};
WOBr(4) = {[num2str(round(median(FL.WOBr),1)),' [',num2str(round(FL.WOBr_quartiles(1,1),1)),'-',num2str(round(FL.WOBr_quartiles(1,3),1)),']']};

WOBe = cell(4,1);
WOBe(1) = {[num2str(round(median(young.WOBe),1)),' [',num2str(round(young.WOBe_quartiles(1,1),1)),'-',num2str(round(young.WOBe_quartiles(1,3),1)),']']};
WOBe(2) = {[num2str(round(median(old.WOBe),1)),' [',num2str(round(old.WOBe_quartiles(1,1),1)),'-',num2str(round(old.WOBe_quartiles(1,3),1)),']']};
WOBe(3) = {[num2str(round(median(NFL.WOBe),1)),' [',num2str(round(NFL.WOBe_quartiles(1,1),1)),'-',num2str(round(NFL.WOBe_quartiles(1,3),1)),']']};
WOBe(4) = {[num2str(round(median(FL.WOBe),1)),' [',num2str(round(FL.WOBe_quartiles(1,1),1)),'-',num2str(round(FL.WOBe_quartiles(1,3),1)),']']};

VT = cell(4,1);
VT(1) = {quantile(max(young.volume),3)};
VT(2) = {quantile(max(old.volume),3)};
VT(3) = {quantile(max(NFL.volume),3)};
VT(4) = {quantile(max(FL.volume),3)};

RMSE_R1insp = cell(4,1);
RMSE_R1insp(1) = {quantile(young.RMSE_linear,3)};
RMSE_R1insp(2) = {quantile(old.RMSE_linear,3)};
RMSE_R1insp(3) = {quantile(NFL.RMSE_linear,3)};
RMSE_R1insp(4) = {quantile(FL.RMSE_linear,3)};

T = table(catagory,R1insp,meanR2Phi,E,Peff_min,WOB_Peff,WOBr,WOBe,VT,RMSE_R1insp,'VariableNames',VarNames);
disp(T)

T2 = rows2vars(T);
disp(T2)

writetable(T,'tableIQR_MED.csv')
writetable(T2,'tableIQR_MED_2.csv')


