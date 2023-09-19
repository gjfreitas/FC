function derivadas = f(t,solucao)
g = 9.8;
vlim = 6.8;
derivadas = zeros(4,1);
% Esta linha é necessária e faz do vetor de saída um vetor coluna.

% x = solucao(1);
% vx = solucao(2);
% z = solucao(3);
% vz = solucao(4);

% O vetor solucao tem os valores de x e v para o tempo t em que a
% função é chamada pela rotina ode45.
derivadas(1) = solucao(2);
% Escreva acima a expressão da derivada de x em função de solucao(1)
% e de solucao(2).
derivadas(2) = (- g/(vlim)^2 * solucao(2) * abs(solucao(2)));
% Escreva a acima expressão da derivada de v em função de solucao(1)
% e de solucao(2).
derivadas(3) = solucao(4);
derivadas(4) = (- g/(vlim)^2 * solucao(4) * abs(solucao(4)) - g);
