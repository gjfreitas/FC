function derivadas = f(t,solucao)
g = 9.8;
vlim = 6.8;
derivadas = zeros(4,1);
% Esta linha � necess�ria e faz do vetor de sa�da um vetor coluna.

% x = solucao(1);
% vx = solucao(2);
% z = solucao(3);
% vz = solucao(4);

% O vetor solucao tem os valores de x e v para o tempo t em que a
% fun��o � chamada pela rotina ode45.
derivadas(1) = solucao(2);
% Escreva acima a express�o da derivada de x em fun��o de solucao(1)
% e de solucao(2).
derivadas(2) = (- g/(vlim)^2 * solucao(2) * abs(solucao(2)));
% Escreva a acima express�o da derivada de v em fun��o de solucao(1)
% e de solucao(2).
derivadas(3) = solucao(4);
derivadas(4) = (- g/(vlim)^2 * solucao(4) * abs(solucao(4)) - g);
