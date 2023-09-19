function derivadas = f(x,solucao,lambda)
q = 5;
% x = 0:0.01:5;
derivadas = zeros(2,1);

% O vetor solucao tem os valores de x e v para o tempo t em que a função é chamada pela rotina ode45.
derivadas(1) = solucao(2);

% Escreva acima a expressão da derivada de x em função de solucao(1) e de solucao(2).

derivadas(2) = -(lambda - 2 * q * cos(2.*x))*solucao(1);
% Escreva a acima expressão da derivada de v em função de solucao(1) e de solucao(2). 
% Se achar que torna o programa mais claro, pode incluir na função as linhas,
% y = solucao(1);
% dy/dx = solucao(2);
end