% Física Computacional 2018/2019  Tiago Fernandes
% 1º teste 2017 - Problema 2c

clear all
close all
clc

L = 0.100;
ro = 2.70*10^3;
sigma = 5.6703*10^-8;
A = L^2;
m = A*L / ro;
c = 0.91*10^3;

ti = 0;
tf = 3600*2;
h = 0.01;
t = ti:h:tf;
n = length(t);

T = zeros(1, n);
T(1) = 310;

fT = @(T, Tc) -sigma*A/(m*c) * (T^4 - Tc^4);
Tc = @(t) 283 + 1.0*10^-3 * t;
    
for i = 1:n-1
    
    r1T = fT(T(i), Tc(t(i)));                       %sempre igual no RK
    r2T = fT(T(i) + r1T * h/2, Tc(t(i)+0.5*h));     %T somam com rT
    r3T = fT(T(i) + r2T * h/2, Tc(t(i))+0.5*h);
    r4T = fT(T(i) + r3T * h, Tc(t(i))+h);
    T(i+1) = T(i) + 1/6*(r1T + 2*r2T + 2*r3T + r4T)*h;
    
end

plot( t, T, 'y')
set(gca, 'color', 'k')
xlabel 'tempo (s)'
ylabel 'Temperatura (K)'
title('Temperatura do bloco de alumínio')
