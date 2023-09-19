%% Exame Prático de Recurso — Parte 1 Física Computacional — 2018/2019 9 de julho de 2019
% 2 b)
clear all
close all
clc

lambda = 2.2;

xi = 0;
xf = pi;
h = 0.01;
x = xi:h:xf;
N = length(x);

y = zeros(1,N); y(1) = 0;
v = zeros(1,N); v(1) = 1;

% Shooting
nshots = 100;
B = 0; % y(L) = 0            % no metodo de shooting B é a solução conhecida para y(xfinal)
tol = 10^-12;

result = zeros(1,nshots);

guess = [39 45];

for ishot = 1:nshots
    
    y = zeros(1,N); y(1) = 0;
    v = zeros(1,N); v(1) = 1;
    
    lambda = guess(ishot);
    fy = @(vy) vy;
    fv = @ (x, y, v) ((x-lambda)*y+1.6*sin(x)*cos(x)*v)/(1-0.8*sin(x)^2);
    
    for i = 1:N-1
        v(i+1) = v(i) + fv( x(i), y(i), v(i)) * h;
        y(i+1) = y(i) + v(i+1) * h;
    end

    result(ishot) = y(end);         % result é a o ultimo y -> y(L)
    
    if ishot >= 2           % a partir da segunda iteração começamos a avaliar os results usando o metodo da secante
       m = (result(ishot)-result(ishot-1))/(guess(ishot)-guess(ishot-1));
       guess(ishot+1) = guess(ishot) + (B - result(ishot))/m;
       if abs(guess(ishot) - guess(ishot-1)) < tol
            break
       end      
    end
end

sol = guess(ishot);
disp(['lambda: ', num2str(sol)])
