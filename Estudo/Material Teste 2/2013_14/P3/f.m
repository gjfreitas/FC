function derivadas = f(x,solucao)
const = 0.1;
Tamb = 0;
derivadas = zeros(2,1);

% O vetor solucao tem os valores de x e v para o tempo t em que a função é chamada pela rotina ode45.
derivadas(1) = solucao(2);

% Escreva acima a expressão da derivada de x em função de solucao(1) e de solucao(2).

derivadas(2) = const*(solucao(1)-Tamb);
% Escreva a acima expressão da derivada de v em função de solucao(1) e de solucao(2). 
% Se achar que torna o programa mais claro, pode incluir na função as linhas,
% y = solucao(1);
% dy/dx = solucao(2);
end