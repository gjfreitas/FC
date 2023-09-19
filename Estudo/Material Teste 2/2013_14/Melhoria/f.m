function derivadas = f(t,solucao,const)

g = 9.8;
derivadas = zeros(4,1);

x = solucao(1);
vx = solucao(2);
y = solucao(3);
vy = solucao(4);


derivadas(1) = vx;
derivadas(2) = -const * sqrt(vx^2 + vy^2) * vx;

derivadas(3) = vy;
derivadas(4) = -const * sqrt(vx^2 + vy^2) * vy - g;


% Se achar que torna o programa mais claro, pode incluir na função as linhas,

end