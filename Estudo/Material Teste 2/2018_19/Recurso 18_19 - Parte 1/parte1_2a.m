%% Exame Prático de Recurso — Parte 1 Física Computacional — 2018/2019 9 de julho de 2019
% 2 a)
clear all
close all
clc

lambda = 3.6094;

xi = 0;
xf = pi;
h = 0.01;
x = xi:h:xf;
n = length(x);

y = zeros(1,n); y(1) = 0;
vy = zeros(1,n); vy(1) = 1;

% Metódo de Euler-Cromer
fy = @(vy) vy;
fvy = @ (x, y, vy) ((x-lambda)*y+1.6*sin(x)*cos(x)*vy)/(1-0.8*sin(x)^2);
%fv = dv/dt

for i = 1:n-1
    vy(i+1) = vy(i) + fvy( x(i), y(i), vy(i)) * h;
    y(i+1) = y(i) + vy(i+1) * h;            % Euler-Cromer usa vx(i+1)
end

%grafico y(x)
plot(x,y, 'r'), xlabel('x(m)'), ylabel('y(m)')

