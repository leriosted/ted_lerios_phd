clear all; close all; clc

load FL_00.txt

t = FL_00(:,1);             % time                              [s]
Q = FL_00(:,2);             % flow                              [L/s]
V = FL_00(:,3);             % absolute volume                   [L]
V_shift = FL_00(:,4);       % shift volume machine output       [mL]
P = FL_00(:,5);             % Alveolar Pressure from V_shift    [cmH2O]



%% Plot Setup:

figPos1 = [-2500, 1600, 900, 800]; % Home Desktop
figPos2 = [-1800, 1600, 900, 800];
figPos3 = [-1900, 400, 900, 900];
figPos4 = [-900, 400, 900, 900];

figPos1 = [-1900, 300, 1200, 800]; % Publishing size
figPos2 = [-1250, 500, 1200, 600];
figPos3 = [-600, 500, 600, 600];
figPos4 = [-1900, 150, 600, 600];
figPos5 = [-1250, 150, 600, 600];
figPos6 = [-600, 150, 600, 300];

% figPos1 = [-900, 20, 800, 700]; % Work Desktop
% figPos2 = [-900, 20, 800, 700];
% figPos3 = [-900, 20, 800, 700];
% figPos4 = [-1900, 20, 800, 700];

pt = 1.6;
FS = 11;

color1 = '#0000FF';     % Blue
color2 = '#321EE1';     
color3 = '#4D1CB7';     
color4 = '#611A96';     
color5 = '#761874';    
color6 = '#8B1654';      
color7 = '#A2142F';     % Burgandy

RGBcolor1 = [0 0 255];
RGBcolor2 = [50 30 225];
RGBcolor3 = [77 28 183];
RGBcolor4 = [97 26 150];
RGBcolor5 = [118 24 116];
RGBcolor6 = [139 22 84];
RGBcolor7 = [162 20 47];

% % Divergent Color Palette
% color1 = '#488f31';     
% color2 = '#7a9c33';     
% color3 = '#a6a73e';     
% color4 = '#f4a05c';     
% color5 = '#f18255';     
% color6 = '#ea6356';     
% color7 = '#de425b';     

%% Time Series Raw Data Plot
margin = [0.08,0.08];
fig1 = figure('Position',figPos1);
p(1) = subplot_tight(4,6,[1 2 3 4],margin);
plot(t,Q,'Color','#0072BD', 'Linewidth', pt); hold on
xlabel('t (s)')
ylabel('Q (L/s)')
grid on
set(gca,'FontSize',FS);
p(2) = subplot_tight(4,6,[7 8 9 10],margin);
plot(t,V,'Color','#0072BD', 'Linewidth', pt); hold on
xlabel('t (s)')
ylabel('V (abs) (L)')
grid on
set(gca,'FontSize',FS);
p(3) = subplot_tight(4,6,[13 14 15 16],margin);
plot(t,V_shift,'Color','#0072BD', 'Linewidth', pt); hold on
xlabel('t (s)')
ylabel('V (shift) (mL)')
grid on
set(gca,'FontSize',FS);
p(4) = subplot_tight(4,6,[19 20 21 22],margin);
plot(t,P,'Color','#0072BD', 'Linewidth', pt); hold on
xlabel('t (s)')
ylabel('P (cmH_2O)')
grid on
set(gca,'FontSize',FS);
p(5) = subplot_tight(4,6,[5 6 11 12],margin);
plot(P,V,'Color','#0072BD', 'Linewidth', pt); hold on
xlabel('P (cmH_2O)')
ylabel('V (abs) (L)')
grid on
set(gca,'FontSize',FS);
p(6) = subplot_tight(4,6,[17 18 23 24],margin);
plot(P,Q,'Color','#0072BD', 'Linewidth', pt); hold on
xlabel('P (cmH_2O)')
ylabel('Q (L/s)')
grid on
sgtitle('Flow Limited Patient 00')




fig1_savename = ('FL_Patient_00.png');
print('-r400','-dpng',fig1_savename);








