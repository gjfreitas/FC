%% Problema 3.3: Oscilador harmónico simples — Runge–Kutta de passo adaptativo

clear all
close all
clc

h = 0.1;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);
t = t0:h:tf;

% solução exata
x_exact = x0*cos(w.*t); 
v_exact = -w*x0*sin(w.*t);

% ode45
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t_ode45,sol] = ode45(@f,[t0 tf],[x0 v0],options);

figure(1)
subplot(2,2,1)
plot(t,x_exact,'k.-',t_ode45,sol(:,1),'r--'), xlabel('t'), ylabel('x')
subplot(2,2,2)
plot(t,v_exact,'k.-',t_ode45,sol(:,2),'r--'), xlabel('t'), ylabel('v')

Et_exact = 0.5 * m * v_exact.^2 + 0.5 * K * x_exact.^2;
Et_ode45 = 0.5 * m * sol(:,2).^2 + 0.5 * K * sol(:,1).^2;

subplot(2,2,3)
plot(x_exact,v_exact,'k.-',sol(:,1),sol(:,2),'r--'), xlabel('x'), ylabel('v')
subplot(2,2,4)
plot(t,Et_exact,'k',t_ode45,Et_ode45,'r')

function derivadas = f(t,solucao)
K = 16;m = 1;
derivadas = zeros(2,1);

% O vetor solucao tem os valores de x e v para o tempo t em que a função é chamada pela rotina ode45.
derivadas(1) = solucao(2); % = v
% Escreva acima a expressão da derivada de x em função de solucao(1) e de solucao(2).

derivadas(2) = -K/m * solucao(1); % -K/m * x
% Escreva a acima expressão da derivada de v em função de solucao(1) e de solucao(2). 
% Se achar que torna o programa mais claro, pode incluir na função as linhas,
% x = solucao(1);
% v = solucao(2);
end