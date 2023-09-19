%% Problema 4.3: Oscilador de van der Pol — ode45

clear all
close all
clc

t0 = 0; %s
tf = 100; %s
v0 = [7, -2, 0.7]; %m/s
y0 = [2, -5, 0.2]; %m
eps = 0.1;

% ode45
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;

for index = 1:length(y0)
    options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
    [t_ode45,sol] = ode45(@f,[t0 tf],[y0(index) v0(index)],options,eps);

    figure(1)
    subplot(1,3,1)
    plot(t_ode45,sol(:,1)), xlabel('t'), ylabel('y')
    hold all
    subplot(1,3,2)
    plot(t_ode45,sol(:,2)), xlabel('t'), ylabel('v')
    hold all
    subplot(1,3,3)
    plot(sol(:,1),sol(:,2)), xlabel('y'), ylabel('v')
    hold all
end


function derivadas = f(t,solucao,eps)


derivadas = zeros(2,1);

% O vetor solucao tem os valores de x e v para o tempo t em que a função é chamada pela rotina ode45.
derivadas(1) = solucao(2); % = v
% Escreva acima a expressão da derivada de x em função de solucao(1) e de solucao(2).

derivadas(2) = -eps .* (solucao(1).^2 - 1) .* solucao(2)-solucao(1); % = -eps*(Y.^2-1)*V-Y
% Escreva a acima expressão da derivada de v em função de solucao(1) e de solucao(2). 
% Se achar que torna o programa mais claro, pode incluir na função as linhas,
% y = solucao(1);
% v = solucao(2);
end