function derivadas = f(t,solucao)
w = 9;
g = 9.8;
R = 0.2;

derivadas = zeros(2,1);
% Esta linha � necess�ria e faz do vetor de sa�da um vetor coluna.
theta = solucao(1);
dtheta = solucao(2);
% O vetor solucao tem os valores de x e v para o tempo t em que a
% fun��o � chamada pela rotina ode45.
derivadas(1) = dtheta;
% Escreva acima a express�o da derivada de x em fun��o de solucao(1)
% e de solucao(2).
derivadas(2) = (w^2*cos(theta)-g/R)*sin(theta);
% Escreva a acima express�o da derivada de v em fun��o de solucao(1)
% e de solucao(2).