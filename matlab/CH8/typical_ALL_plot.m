clear all; close all; clc
colorBlue = '#0000FF';
colorRed = '#FF0000';
colorBlack = '#000000';
colorPurple = "#FF00FF";
colorGreen = "#D95319";

pt = 0.5;
FS = 10;

load('results.mat')

figure('Position',[-1900 50 700 900])
for kk = 1:6
    subplot(4,2,kk)
    plot(A1.T{1,1,kk},A1.Paw{1,1,kk},'Color',colorBlue,'Linewidth',pt,'Linestyle','-'); hold on
    plot(A1.T{1,1,kk},A1.Pes{1,1,kk},'Color',colorBlack,'Linewidth',pt,'Linestyle','-'); hold on
    plot(A1.T{1,1,kk},A1.Peff{1,1,kk},'Color',colorRed,'Linewidth',pt,'Linestyle','--'); hold on
    plot(A2.T{1,1,kk},A2.Peff{1,1,kk},'Color',colorBlue,'Linewidth',pt,'Linestyle','--'); hold on
    plot(A3.T{1,1,kk},A3.Peff{1,1,kk},'Color',colorBlack,'Linewidth',pt,'Linestyle','--'); hold on
    str = append('Patient ',num2str(kk));
    title(str)
    grid on
    xlim([0 3])
    ylim([-20 25])
    xlabel('t (s)')
    ylabel('P (cmH_2O)')
    set(gca,'FontSize',FS);
end

leg1 = legend('P_{aw}','P_{es}','P_{eff,1}','P_{eff,2}','P_{eff,3}','FontSize',10,...
        'Orientation','horizontal');
    legend boxoff
    leg1.Position(1:2) = [.35 .24];
set(gca,'FontSize',FS);
    filename1 = append('typical_ALL');
    print('-r400','-dpng',filename1);
