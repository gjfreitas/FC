% Física Computacional 2018/2019  Tiago Fernandes
% 1º teste 2017 - Problema 1a

clear all
close all
clc

L = 67.0;
omega = 7.292*10^-5;
phi = 48.865;           %angulo em graus
sphi = sind(phi);
g = 9.8;

ti = 0;
tf = 200;
h = 0.01;
t = ti:h:tf;
n = length(t);

x = zeros(1, n);
y = zeros(1, n);
vx = zeros(1,n);
vy = zeros(1,n);

%condicoes iniciais
x(1) = 2.00;
y(1) = 0.00;
vx1(1) = 0;
vy1(1) = 0;

%funcoes anonimas para as derivadas
fx = @(vx) vx;
fy = @(vy) vy;
fvx = @ (x, vy) 2 * omega * sphi * vy - g/L * x;
fvy = @ (y, vx) -2 * omega * sphi * vx - g/L * y;

for i = 1:n-1
    vx(i+1) = vx(i) + fvx(x(i), vy(i)) * h;
    x(i+1) = x(i) + vx(i+1) * h;
    
    vy(i+1) = vy(i) + fvy( y(i), vx(i)) * h;
    y(i+1) = y(i) + vy(i+1) * h;
end

plot(x, y, 'y')
set(gca, 'color', 'k')
xlabel 'x(m)'
ylabel 'y(m)'
title('Pêndulo de Foucalt')