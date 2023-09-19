clc, clear all, close all


% 1) a.

% Condições iniciais:
g = 9.8;
l = 1;
b = 1;
va0 = 0;
angIni = 0.2;
tf = 30;
h = 0.01;

[t,ang,va] = rk2(h,tf,angIni,va0,b,g,l)

% Gráfico do ângulo em função do tempo.
plot(t,ang), xlabel('Tempo (s)'), ylabel('Angulo (rad)')
