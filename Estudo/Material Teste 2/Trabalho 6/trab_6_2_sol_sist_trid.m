%% Trabalho Prático 6 Condução de calor
%% Problema 6.2 - Crank Nicolson
clear all
close all
clc

L = 50; %cm
k = 0.93; % cal/(s*cm*ºC) -> condutividade da barra
c = 0.094; %cal/(g*ºC) -> calor especifico
p = 8.9; %g/cm^3 -> densidade
a = k/(c*p);

Ti = 100; %ºC

t0 = 0; %s
tf = 500; %s
dt = 0.1; %s
t = t0:dt:tf; %s
Nt = length(t);

x0 = 0; %cm
xf = L; %cm
dx = 0.5; %cm
x = x0:dx:xf; %cm
Nx = length(x);

T = zeros(Nx, Nt);
% Condições iniciais
T(2:Nx-1,1) = Ti;
% Condições fronteira
T(1,:) = 0;
T(Nx,:) = 0;

% Implementação das matrizes com base no Slide 10 - Aula 6
eta = a *(dt/dx^2);

A1 = diag(repmat(2/eta + 2,1,Nx-2)); % (2/n + 2) na diagonal

A2 = diag(-ones([1 Nx-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(-ones([1 Nx-3]),-1);% desce uma posição relativamente a diagonal


A = A1 + A2 + A3; % matriz no slide
b = zeros(Nx-2,1);

const = 2/eta - 2;

for n = 1:Nt-1
    for i = 1: Nx-2
        b(i) = T(i,n) + const*T(i+1,n) + T(i+2,n);
    end
%     b(1) = b(1) + T(1,j+1); % pelos os slides, acrescentar o elemento na primeira linha a vermelho apesar de nao ser preciso pq é 0
%     b(Nx-2) = b(Nx-2) + T(Nx,j+1); % pelos os slides, acrescentar o elemento na ultima linha a vermelho apesar de nao ser preciso pq é 0
    % a)
    T(2:Nx-1,n+1) = linsolve(A,b);
%     % b)
%     T(2:Nx-1,n+1) = sol_sist_trid(A,b);
%     % c)
%     [L, U, P] = lu(A);
%     y = L\b;
%     T(2:Nx-1,n+1) = U\y;
end

% outra maneira
% for j = 1: Nt -1
%     b = T(1:Nx-2,j) + const * T(2:Nx-1,j) + T(3:Nx,j); % pelos slides o que esta a preto
%     b(1) = b(1) + T(1,j+1); % pelos os slides, acrescentar o elemento na primeira linha a vermelho apesar de nao ser preciso pq é 0
%     b(Nx-2) = b(Nx-2) + T(Nx,j+1); % pelos os slides, acrescentar o elemento na ultima linha a vermelho apesar de nao ser preciso pq é 0
%         % a)
%     T(2:Nx-1,n+1) = linsolve(A,b); % devolve uma matriz Nx linhas por 1 coluna
%     % b)
%     T(2:Nx-1,n+1) = sol_sist_trid(A,b);
%     % c)
%     [L, U, P] = lu(A);
%     y = L\b;
%     T(2:Nx-1,n+1) = U\y;
% 
% end



figure(1)
mesh(t,x,T), xlabel('t(s)'), ylabel('x(cm)'), zlabel('T(ºC)')
figure(2)
contourf(x,t,T'), xlabel('x(cm)'), ylabel('t(s)')  % T' transposta de T para dar igual as soluções
h = colorbar; % barra de cores
set(get(h,'label'),'string','T(ºC)'); % label da barra de cores