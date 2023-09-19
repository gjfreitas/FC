function derivadas = f(t,solucao,alfa)
m = 1.5; %Kg
K = 2; %N/m
% alfa = -0.2;

derivadas = zeros(2,1);

% O vetor solucao tem os valores de x e v para o tempo t em que a função é chamada pela rotina ode45.
derivadas(1) = solucao(2); % = v
% Escreva acima a expressão da derivada de x em função de solucao(1) e de solucao(2).

derivadas(2) = -K/m * solucao(1)*(1 + 3/2 * alfa * solucao(1));
% Escreva a acima expressão da derivada de v em função de solucao(1) e de solucao(2). 
% Se achar que torna o programa mais claro, pode incluir na função as linhas,
% x = solucao(1);
% v = solucao(2);
end