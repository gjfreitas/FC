% Física Computacional 2018/2019  Tiago Fernandes
% 1º teste 2017 - Problema 2b

clear all
close all
clc

L = 0.100;
ro = 2.70*10^3;
sigma = 5.6703*10^-8;
A = L^2;
m = A*L / ro;
c = 0.91*10^3;
Tc = 0;
T0 = 2000;

ti = 0;
tf = 1000;
hs = [ 0.01 0.02 0.025 0.04 0.05];
nh = length(hs);

Tf_an = (T0^-3 + 3*sigma*A / (m*c) * tf)^(-1/3);            %temperatura final analitica

Erro = zeros(1, nh);            %erro para cada h       

for j = 1:nh
    h = hs(j);
    t = ti:h:tf;
    n = length(t);
    T = zeros(1, n);
    T(1) = T0;
    
    for i = 1:n-1
        %vetor com os coeficientes do polinomio de quarto grau
        C = [sigma*A*h / (2*m*c), 0, 0, 1, -(T(i) - sigma*A*h/(2*m*c)*(T(i)^4 - 2 * Tc^4))];
        sol = roots(C);
        [valor, indice] = min(abs(sol-T(i)));
        T(i+1) = sol(indice);
    end
    
    Erro(j) = abs(Tf_an - T(end));  %modulo da diferença entre o Tf obtido e o Tf analitico
    
end

plot( log(hs), log(Erro), '*y')
lsline
set(gca, 'color', 'k')
xlabel 'log(h) '
ylabel 'log(Erro)'
title('Arrefecimento de um bloco de alumínio')

p = polyfit(log(hs), log(Erro), 1);
ordem_do_metodo = p(1)         %ordem do metodo é o decluve da reta
