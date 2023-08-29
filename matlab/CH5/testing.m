clear all; close all; clc

A = rand(20,1);
B = rand(20,1);
C = rand(20,1);

RMSEA = sqrt(mean((A - C).^2));
RMSEB = sqrt(mean((B - C).^2));

RMSE_both = sqrt(RMSEA^2 + RMSEB^2);
RMSE_both_check = sqrt(mean(((A) - C).^2));