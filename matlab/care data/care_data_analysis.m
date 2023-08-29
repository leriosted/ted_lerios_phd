clear all; close all; clc
atWork = 0;
print_on = 1;

colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
colorGreen = '#009A17';
colorOrange = '#FF6600';
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
    figPos1 = [-1900, 150, 900, 800]; % Home
    figPos2 = [-1250, 500, 600, 600];
    figPos3 = [-600, 500, 600, 600];
    figPos4 = [-1900, 150, 600, 600];
    figPos5 = [-1250, 150, 600, 600];
    figPos6 = [-600, 150, 600, 300];
end
%% Load Data
dataCatagory = 'care_data';
files = dir(fullfile(dataCatagory,'*.mat'));
%% Process data (Calculate volume and run model)

for ii = 1:1%length(files)
    thisFileName = fullfile(dataCatagory,files(ii).name);
    fprintf('Processing "%s".\n',thisFileName);
    load(thisFileName)


    

    %% E
    for jj = 1:20
        V_ends(ii,jj) = data(jj).V(end);
        E(ii,jj) = (data(jj).P(data(jj).ind_insp(end))-data(jj).PEEP)./(max(data(jj).V));
    end
    E_med(ii,1) = median(E(ii,:));
    %% R
    
    for kk = 1:20
        P = data(kk).P(data(kk).ind_exp); 
        Q = data(kk).Q(data(kk).ind_exp); 
        V = data(kk).V(data(kk).ind_exp);
        PEEP = data(kk).PEEP;
        pp = find(P>(PEEP+1.5),1,'last');
        pp = pp:length(P)-5;
        b = P(pp)-E_med(ii)*V(pp)-PEEP;
        A = [Q(pp), Q(pp).*abs(Q(pp))];
        R(:,ii) = lsqnonneg(A,b); 
        
        figure('Position',figPos1)
        subplot(4,2,1)
        plot(data(kk).P(data(kk).ind_exp),data(kk).Q(data(kk).ind_exp),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
        xlabel('')
        ylabel('P')
        grid on
        subplot(4,2,3)
        plot(data(kk).P(data(kk).ind_exp),data(kk).Q(data(kk).ind_exp),'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
        xlabel('')
        ylabel('P')
        grid on
    end
    
end



        %% Plots


% vnames = ["Test No.","h","p"];
% 
% TT = table([1;2;3;4;5],[h1;h2;h3;h4;h5],[p1;p2;p3;p4;p5],'Variablenames',vnames);
% disp(TT)
% 
% writetable(TT,'table_ttest.csv')