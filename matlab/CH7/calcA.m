clear all; close all; clc

syms Qp A t texp tend Vp

Qp = A*sin(t/tend - texp);

Vp = int(Qp,t)

