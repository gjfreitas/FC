function derivadas = f(x,solucao)

derivadas = zeros(2,1);

derivadas(1) = solucao(2);

derivadas(2) = -abs(solucao(1));

% y = solucao(1);
% dy/dx = solucao(2);
end