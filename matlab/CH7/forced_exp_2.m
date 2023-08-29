clear all; close all; clc

load young_results.mat

Pa = young.pressure(17,:)';
V = young.volume(17,:)';
Q = young.flow(17,:)';
T = young.time';

aa = find(Q<=0,1,"first"):length(V);
V2 = V(aa) - V(end);
Vexp = max(V)-V(end);
Texp = T(end)-T(find(Q<=0,1,"first"));
plot(-(V2-Vexp),-Q(aa))
A = pi*(Vexp)/(2*Texp);
t = T(aa) - T(find(Q<=0,1,"first"));
PassiveQ = A*sin(pi*t/Texp);
v2 = cumtrapz(t,PassiveQ);
plot(v2,PassiveQ,'--'); hold on
plot(-(V2-Vexp),-Q(aa),'r--')

load FL_results.mat

Pa = FL.pressure(17,:)';
V = FL.volume(17,:)';
Q = FL.flow(17,:)';
T = FL.time';
figure()
aa = find(Q<=0,1,"first"):length(V);
V2 = V(aa) - V(end);
Vexp = max(V)-V(end);
Texp = T(end)-T(find(Q<=0,1,"first"));
plot(-(V2-Vexp),-Q(aa)); hold on
A = pi*(Vexp)/(2*Texp);
t = T(aa) - T(find(Q<=0,1,"first"));
PassiveQ = A*sin(pi*t/Texp);
v2 = cumtrapz(t,PassiveQ);
plot(v2,PassiveQ,'--'); 
plot(-(V2-Vexp),-Q(aa),'r--')