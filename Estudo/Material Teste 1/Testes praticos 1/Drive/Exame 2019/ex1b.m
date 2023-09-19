% Física Computacional 2018/2019  Tiago Fernandes
% Problema 3.3

clear all
close all
clc

ti = 0;
tf = 2.5;
h = 0.01;
t = ti:h:tf;            %vetor com todos os instantes

%condicoes iniciais
v0 = 30;
x0 = 0;
z0 = 1;
vx0 = v0 * cosd(60);
vz0 = v0 * sind(60);


%tolerancias dadas no enunciado
reltol = 3E-14;
abstol_1=1E-13;
abstol_2=1E-13;


options = odeset('RelTol', reltol, 'AbsTol', [abstol_1 abstol_2 abstol_1 abstol_2]);      %opcoes para o ode45
[t, sol] = ode45(@f, t, [x0 vx0 z0 vz0], options);              %usa a funcao f em ficheiro separado
%devolve um vetor t e um vetor sol      (tempo = t)
%a primeira coluna de sol tem todas as posições      (x = sol(:,1)
%a segunda coluna de sol tem todas as velocidades     (vx = sol(:,2)

plot(sol(:,1), sol(:,3), '-y')
set(gca,'Color','k')
title('Volante de badminton - ode45')
xlabel 'Abcissa (m)'
ylabel 'Altura (m)'
  