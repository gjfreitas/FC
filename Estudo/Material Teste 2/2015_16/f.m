function derivadas = f(x,solucao,w)
L =3;
T = 5.0E4;
alfa = 5.0E-8;

derivadas = zeros(2,1);
y = solucao(1);
dy = solucao(2);

derivadas(1) = dy;

derivadas(2) = 2*alfa*T*y + alfa*w*x*(L-x);

% y = solucao(1);
% dy/dx = solucao(2);
end