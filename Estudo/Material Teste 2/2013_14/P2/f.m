function derivadas = f(x,solucao,lambda)

derivadas = zeros(2,1);
y = solucao(1);
dy = solucao(2);

derivadas(1) = dy;

derivadas(2) = (-2*sech(x)^2-lambda^2)*y;

% y = solucao(1);
% dy/dx = solucao(2);
end