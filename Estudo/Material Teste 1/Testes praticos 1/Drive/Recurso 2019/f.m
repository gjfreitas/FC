function derivadas = f(t,solucao)
w = 9;
g = 9.8;
R = 0.2;

derivadas = zeros(2,1);
% Esta linha é necessária e faz do vetor de saída um vetor coluna.
theta = solucao(1);
dtheta = solucao(2);
% O vetor solucao tem os valores de x e v para o tempo t em que a
% função é chamada pela rotina ode45.
derivadas(1) = dtheta;
% Escreva acima a expressão da derivada de x em função de solucao(1)
% e de solucao(2).
derivadas(2) = (w^2*cos(theta)-g/R)*sin(theta);
% Escreva a acima expressão da derivada de v em função de solucao(1)
% e de solucao(2).