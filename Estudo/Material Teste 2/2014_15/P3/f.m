function derivadas = f(r,solucao,const)

derivadas = zeros(2,1);

derivadas(1) = solucao(2);

derivadas(2) = -(1/r)*solucao(2)-const;

% T = solucao(1);
% U = solucao(2);
end