%% Problema 9.2 Função de onda radial do átomo de hidrogénio
clear all
close all
clc

% S(x) = 0;
% g(x) = 2(E + 1/r - l*(l+1)/(2*r^2))
% y(x) = u(x)

% Constantes
n = 3; 
L = 1; 
E3 = (-1/2)*n^(-2);

% pré-alocações
h = 0.001;
r = 0:h:50;
Nr = length(r);

u = zeros(1,Nr);

% Condições fronteira
u(1) = 0;
u(Nr) = 0;
u(Nr-1) = h*10^(-3); % lim_u(rmax) = 0

% Resultado analitico
R_a = (4*sqrt(2))/(27*sqrt(3))*(1 - r./6).* r .* exp(-r./3);

% guesses
guess(1) = -0.6;
guess(2) = -0.7;

% Goal do shooting
B = 0; % u0 = 0
tol = 1E-12;

for i = Nr-1:-1:3
    g = 2*((guess(1) + 1./r - L*(L+1)./(2*r.^2)));
    u(i-1) = (1+h^2*g(i-1)/12)^(-1)*(2*(1-5*h^2*g(i)/12)*u(i)-(1+h^2*g(i+1)/12)*u(i+1));
end
u(1) = interp1(r(2:5),u(2:5),0,'spline');
result(1) = u(1);

while abs(guess(2)-guess(1)) > tol
    for i = Nr-1:-1:3
        g = 2*((guess(1) + 1./r - L*(L+1)./(2*r.^2)));
        u(i-1) = (1+h^2*g(i-1)/12)^(-1)*(2*(1-5*h^2*g(i)/12)*u(i)-(1+h^2*g(i+1)/12)*u(i+1));
    end
    u(1) = interp1(r(2:5),u(2:5),0,'spline');
    result(2) = u(1);
    m = (result(2)-result(1))/(guess(2)-guess(1));
    guess = [guess(2), guess(2) + (B-result(2))/m];
end
sol = guess(end-1);

C1 = trapz(r,abs(u).^2);
u_norm = u/sqrt(C1);
R(2:Nr) = u_norm(2:Nr)./r(2:Nr);
R(1) = interp1(r(2:10),R(2:10),0,'spline');

plot(r,R,'o',r,R_a)