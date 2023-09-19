%% Problema 5.3: Método de diferenças finitas — perfil de temperaturas numa resistência elétrica cilíndrica
clear all
close all
clc

u = 10^(-3); %kg/m
L = 1; %m
T = 10^3; %N

h = 1E-6;
R = 1E-3; % m
Q = 2.1E6; %W m^-3
lambda = 0.1; % W K^-1 m^-1
r = 0:h:R;
N = length(r);

A1 = diag(repmat(-2,[1 N])); % matriz com 2 na diagonal
A1(end,end) = 1; % por um 1 no ultimo elemento da matriz

A2 = diag(1 + h./(2.*r(1:(end-1))),1);
A2(1,2) = 2;
A3 = diag(1 - h./(2.*r(2:(end))),-1);
A3(end,end-1) = 0;

A = A1 + A2 + A3;

b = repmat(-h^2 * Q / lambda, [N 1]);
b(end) = 20;

T = linsolve(A,b);

plot(r,T,'.-'), xlabel('r'), ylabel('T'), grid

