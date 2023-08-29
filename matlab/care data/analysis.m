clear
close all
files = dir(fullfile(cd,'*.mat'));

for kk = 1:length(files)
    load(fullfile(files(kk).folder,files(kk).name));

    for ii = 1:10%length(data)
        %PULL OUT DATA
        T =data(ii).T;    P = data(ii).P;    Q = data(ii).Q;    V = data(ii).V;
        PEEP = data(ii).PEEP

        ind_exp = data(ii).ind_exp;
        ind_insp = data(ii).ind_insp;
        
        %BITS FOR ESTIMATING E
        DP = (P(ind_insp(end))-PEEP); DV = (max(V)-V(1)); Qie = Q(ind_insp(end));

        %JUST THE FLAT BIT OF EXPIRATION
        P = P(ind_exp); Q = Q(ind_exp); V = V(ind_exp);
        pp = find(P>(PEEP+1.5),1,'last');
        pp = pp:length(P)-5;

        %ITERATE TO GET E AND R
        R = [0,0];
        for jj = 1:10
            E = (DP-(R(1)*Qie+R(2)*abs(Qie))*Qie)/DV;
            b = P(pp)-E*V(pp)-PEEP;
            A = [Q(pp), Q(pp).*abs(Q(pp))];
            R = lsqnonneg(A,b);
        end

        Elin(ii,kk)=E;
        Rlin(ii,kk) =mean(R(1)+R(2)*abs(Q(pp)));

        figure
        subplot(3,1,1)
        plot(P(pp))
        hold on
        Pmod = A*R-b+P(pp);
        plot(Pmod)
        legend({'data','model - with Rohers'})
                grid on
        ylabel('P')

        %Fit dynamic R assuming linear E
        A = diag(Q(pp));
        b = P(pp)-E*V(pp)-PEEP;
        Rnonlin{ii,kk} = A\b;

        %Fit dynamic E assuming Rpher's R
        A = diag(V(pp));
        b = P(pp)-(R(1)+R(2)*abs(Q(pp))).*Q(pp)-PEEP;
        Enonlin{ii,kk} = A\b;

        subplot(3,1,2)
        plot((R(1)+R(2)*abs(Q(pp))))
        hold on
        plot(Rnonlin{ii,kk})
        grid on
%         ylim([-50,100])
        ylabel('R')
        legend({'Rohers R','Dynamic R'})
        pause(0.5)

                subplot(3,1,3)
        plot([1,length(pp)],Elin(ii,kk)*[1,1])
        hold on
        plot(Enonlin{ii,kk})
        grid on
        ylim([-50,100])
        ylabel('E')
        legend({'Linear E','Dynamic E (Rohrs R)'})
        pause(0.5)

    end
    figure
    boxplot([Rlin(:,kk),Elin(:,kk)])
end