%% Problema 9.1 Poço de potencial infinito a uma dimensão
clear all
close all
clc

% S(x) = 0;
% g(x) = 2(E-V(x))
% y(x) = psi(x)

a = 1;
h = 0.001;

% En = n^2 * (pi^2)/(8*(a^2))
E1 = pi^2/(8*(a^2)); % n = 1
disp(['E1 : ',num2str(E1),' Ha'])

% Pré alocações
x = -a:h:a;
Nx = length(x);

psi = zeros(1,Nx);
dpsi = zeros(1,Nx);
V = zeros(1,Nx);
g = zeros(1,Nx);

% Condições fronteira
psi(1) = 0;
psi(2) = h*10^(-3); % previsão necessária pela fórmula do método de Numerov => psi(-a+h)

E = 1.5;
for k = 2:Nx-1
    g(k) = 2*(E-0);
    psi(k+1) = (1 + h^2/12 * g(k-1))^-1 * (-(1 + h^2/12 * g(k-1))*psi(k-1) + 2*(1 - 5 * h^2/12 * g(k))*psi(k));
end
plot(x,psi)
