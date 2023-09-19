%% Trabalho Prático 6 Condução de calor
%% Problema 6.1 - Euler

% a)
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

% T(x,t+1) = T(x,n) + k/(c*p) * (dt/dx^2)[T(x-1,t) - 2*T(x,t) + T(x+1,t)]
% T é uma matriz com 0 na primeira e ultima linha e 100 nas restantes linhas da primeira coluna

T = nan(Nx, Nt);
% Condições iniciais
T(2:Nx-1,1) = Ti; % da 2 a penultima linha por Ti (= 100)
% Condições fronteira
T(1,:) = 0; % 0's na 1 linha
T(Nx,:) = 0; % 0's na ultima linha

% T = [ 0  0  0 ...
% 100 
% 100 
%  .
%  .
% 100
%  0  0  0 ...]

% Método de Euler
% Pelos slides, T(x,t+1) = T(x,t) + ((k/cp)*dt/dx^2) * (T(x-1,t) - 2T(x,t) + T(x+1,t));
for j = 1: Nt - 1 % todas as colunas
    for i = 2:Nx - 1 % desde a 2 linha a penultima
        T(i,j+1) = T(i,j)+ a*dt/dx^2*(T(i-1,j)-2*T(i,j)+T(i+1,j));
    end
end

figure(1)
mesh(t,x,T), xlabel('t(s)'), ylabel('x(cm)'), zlabel('T(ºC)')
figure(2)
contourf(x,t,T'), xlabel('x(cm)'), ylabel('t(s)')  % T' transposta de T para dar igual as soluções
h = colorbar; % barra de cores
set(get(h,'label'),'string','T(ºC)'); % label da barra de cores

% b) Mudar dx e dt
% a = k/(c*p);
% O metódo é estavel se: n = a *(dt/dx^2) <= 1/2
eta = a *(dt/dx^2) <= 1/2;  % 1 = V, 0 =F

% c)
%indice posição a 1/4 da extremidade da barra
index_x = find(abs(x - L/4) < dx/2);
figure(3)
plot(t,T(index_x,:),'r.-'), xlabel('t(s)'), ylabel('T(ºC)')
disp(['x(i) = L/4 <=> i = ',num2str(index_x)])

%indice posição para T = 50ºC na posição 1/4 da exremidade da barra
[~,index_t] = min(abs(T(index_x,:) - 50)); % calcula o indicie do tempo em que a temperatura esta mais proximo do valor pretendido
disp(['T(L/4,t) = 50ºC <=> t ~ ',num2str(t(index_t)),' s'])
