function derivadas = f(t,solucao)

theta = solucao(1);
Q = solucao(2);

derivadas = zeros(2,1);

derivadas(1) = Q;

derivadas(2) = -sin(theta);

end