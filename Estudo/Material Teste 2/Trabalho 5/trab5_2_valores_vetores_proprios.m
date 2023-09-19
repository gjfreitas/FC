%% Problema 5.2: Método das diferenças finitas — determinação da frequência de vários modos normais de vibração
clear all
close all
clc

u = 10^(-3); %kg/m
L = 1; %m
T = 10^3; %N

h = 0.01;
x = 0:h:L;
N = length(x);

% Construção da Matriz
A1 = diag(repmat(-2,1,N-2)); % matriz com -2 na diagonal
A2 = diag(ones([1 N-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(ones([1 N-3]),-1);% desce uma posição relativamente a diagonal
A = A1+A2+A3;

sol1 = eigs(A,3,'sm'); % valores proprios de A

sol = sqrt(-sol1 * T / (u*h^2));
[vec,val] = eigs(A,3,'sm');% vetores proprios de A

disp(['Omega do primeiro modo : ',num2str(sol(1)), ' rads/s'])
disp(['Omega do segundo modo : ',num2str(sol(2)), ' rads/s'])
disp(['Omega do segundo modo : ',num2str(sol(3)), ' rads/s'])

figure(1)
xx = x(2:end-1);
plot(xx,vec(:,1),'r.-',xx,vec(:,2),'b.-',xx,vec(:,3),'m.-'), xlabel('x'), ylabel('y'), title('vtores próprios numéricos de y(x)')
legend('Vetor próprio 1', 'Vetor próprio 2', 'Vetor próprio 3')