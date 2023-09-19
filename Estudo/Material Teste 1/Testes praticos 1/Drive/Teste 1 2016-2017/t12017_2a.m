% Física Computacional 2018/2019  Tiago Fernandes
% 1º teste 2017 - Problema 2a

L = 0.100;
ro = 2.70*10^3;
sigma = 5.6703*10^-8;
A = L^2;
m = A*L / ro;
c = 0.91*10^3;
Tc = 0;

ti = 0;
tf = 1000;
h = 0.01;
t = ti:h:tf;
N = length(t);

T = zeros(1, N);

T(1) = 2000;

% metodo de crank-nicolson
for i = 1:N-1
    C = [ sigma * A/ (m * c) * h/2 , 0, 0, 1, - (T(i) - sigma * A/ (m * c) * h/2 *(T(i)^4 - 2 * Tc^4))];
    sol = roots(C);
    [valor, indice] = min(abs(sol-T(i)));
    T(i+1) = sol(indice);
end

plot( t, T, 'y')
set(gca, 'color', 'k')
xlabel 'tempo (s)'
ylabel 'Temperatura (K)'
title('Arrefecimento de um bloco de alumínio')

