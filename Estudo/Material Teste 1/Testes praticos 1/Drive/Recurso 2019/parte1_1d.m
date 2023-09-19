
clear all
close all
clc

w = 9;
g = 9.8;
R = 0.2;

ti = 0;
tf = 10;
h = 0.1;
t = ti:h:tf;

%condicoes iniciais
theta0 = acos(g/(w^2*R));
dtheta0 = 0;

%tolerancias dadas no enunciado
reltol = 3E-14;
abstol_1=1E-13;
abstol_2=1E-13;

options = odeset('RelTol', reltol, 'AbsTol', [abstol_1 abstol_2]);      %opcoes para o ode45
[t, sol] = ode45(@f, t, [theta0 dtheta0], options);              %usa a funcao f em ficheiro separado
%devolve um vetor t e um vetor sol      (tempo = t)
%a primeira coluna de sol tem todas as posições      (x = sol(:,1)
%a segunda coluna de sol tem todas as velocidades     (vx = sol(:,2)  

theta = sol(:,1);

plot(t, theta, 'y')
set(gca, 'color', 'k')
title('Método Runge-Kutta 4a ordem')
xlabel 'tempo (s)'
ylabel 'theta (rad)'


  