%% Problema 9.2 Função de onda radial do átomo de hidrogénio
clear all
close all
clc

% S(x) = 0;
% g(x) = 2(E + 1/r - l*(l+1)/(2*r^2))
% y(x) = u(x)

% Constantes
n = 1; 
L = 0; 
E1 = (-1/2)*n^(-2);

% pré-alocações
h = 0.001;
r = 0:h:50;
Nr = length(r);
u = zeros(1,Nr);

% Condições fronteira
u(1) = 0;
u(Nr) = 0;
u(Nr-1) = h*10^(-3); % lim_u(rmax) = 0

E = -0.6;
for i = Nr-1:-1:3
    g = 2*((E + 1./r - L*(L+1)./(2*r.^2)));
    u(i-1) = (1+h^2*g(i-1)/12)^(-1)*(2*(1-5*h^2*g(i)/12)*u(i)-(1+h^2*g(i+1)/12)*u(i+1));
end
u(1) = interp1(r(2:5),u(2:5),0,'spline');
result = u(1);