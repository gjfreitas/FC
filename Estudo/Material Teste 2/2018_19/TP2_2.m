%% 2º Teste Prático de Avaliação Discreta Física Computacional — 2018/2019
%% 2 deve tar mal
clear all
close all
clc

lbd = 3;
h = 0.01;
x = 0:h:5;

V = -0.5*lbd*(lbd+1) * (sech(x).^2);
plot(x,V)

% 1/2(phi(i+1) - 2phi(i) + phi(i-1)) - h^2*V^2*phi(i) = -h^2*E*phi(i)
% 1/2 * phi(i+1) + (-1 - h^2*V^2)phi(i)  + 1/2 * phi(i-1) = -h^2*E * phi(i)

N = length(x);
% const = h^2*E;

% Construção da matriz A
for i = 1:N
    A1(i,1) = -1-h^2*(V(i)^2);
end
A1 = diag(A1); % transforma a matriz A1 numa matriz diagonal
A2 = diag(1/2*ones([1 N-1]),1); % sobe um posição relativamente a diagonal
A3 = diag(1/2*ones([1 N-1]),-1);% desce uma posição relativamente a diagonal
A = A1+A2+A3;

% phi = y = > Ay = const*E*y (ver ultimo slide - Aula 6 e trabalho 5.2 (2020/21))

% calcular os vetores e valores próprios
[vec,val]=eigs(A,3,'sm');
const = -h^2;

disp(['Energia do estado fundamental: ',num2str(val(1,1)/const)])
disp(['Energia do 1º estado excitado: ',num2str(val(2,2)/const)])
disp(['Energia do 2º estado excitado: ',num2str(val(3,3)/const)])

plot(x,vec(:,1)), xlabel('x'), ylabel('phi'), title('Função de onda')

% d)
B = legendre(lbd,tanh(x));

plot(x,vec(:,1),'b.-',x,-1*B(1,:),'r.-')
