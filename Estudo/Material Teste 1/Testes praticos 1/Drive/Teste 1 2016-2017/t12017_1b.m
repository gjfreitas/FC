% Física Computacional 2018/2019  Tiago Fernandes
% 1º teste 2017 - Problema 1b

clear all
close all
clc

L = 67.0;
omega = 7.292*10^-5;
phi = 48.865;           %angulo em graus
sphi = sind(phi);
g = 9.8;

ti = 0;
tf = 500;
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

fx = @(vx) vx;
fy = @(vy) vy;
fvx = @ (x, vy) 2 * omega * sphi * vy - g/L * x;
fvy = @ (y, vx) -2 * omega * sphi * vx - g/L * y;

ind_max = [];       %matriz que vai conter os indices dos maximos locais de x

for i = 1:n-1
    vx(i+1) = vx(i) + fvx(x(i), vy(i)) * h;
    x(i+1) = x(i) + vx(i+1) * h;
    
    vy(i+1) = vy(i) + fvy( y(i), vx(i)) * h;
    y(i+1) = y(i) + vy(i+1) * h;
    if i > 1
        if x(i) > x(i-1) && x(i) > x(i+1)       %se x for maior que o x anterior e o x seguinte
            ind_max = [ind_max i];
        end
    end
end

nmax = length(ind_max);     %numero de maximos locais
t_max = zeros(1, nmax);
theta_max = zeros(1, nmax);

for j = 1:length(ind_max)
    ind = ind_max(j);           %indice atual
    t_max(j) = t(ind);
    theta_max(j) = mod (atan2 (y(ind), x(ind)) , 2*pi);     %calculo de theta usando mod
end

plot(t_max, theta_max, '*y')
set(gca, 'color', 'k')
xlabel 'tempo (s)'
ylabel 'theta (rad)'
title('Máximos locais de x')
lsline
p = polyfit(t_max, theta_max, 1);

T_precessao = 2*pi / abs(p(1))      %periodo de precessão em segundos
