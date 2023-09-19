function derivadas = f(x,solucao,const)
Tamb = 20;

derivadas = zeros(2,1);

derivadas(1) = solucao(2);

derivadas(2) = const*(solucao(1)-Tamb);

% T = solucao(1);
% dT/dx = solucao(2);
end