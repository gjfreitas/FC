%% Trabalho Prático 6 Condução de calor
%% Problema 6.3 - Crank Nicolson
clear all
close all
clc

% a)
L = 50; %cm
k = 0.93; % cal/(s*cm*ºC) -> condutividade da barra
c = 0.094; %cal/(g*ºC) -> calor especifico
p = 8.9; %g/cm^3 -> densidade
a = k/(c*p);

t0 = 0; %s
tf = 150; %s
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
T(:,1) = 50*sin(2*pi*x/L); % muda os valores da primeira coluna para todas as linhas
% Condições fronteira
T(1,:) = 0;
T(Nx,:) = 0;

% Implementação das matrizes com base no Slide 10 - Aula 6
n = a *(dt/dx^2);

A1 = diag(repmat(2/n + 2,1,Nx-2)); % (2/n + 2) na diagonal

A2 = diag(-ones([1 Nx-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(-ones([1 Nx-3]),-1);% desce uma posição relativamente a diagonal


A = A1 + A2 + A3; % matriz no slide
b = zeros(Nx-2,1);

const = 2/n - 2;

for n = 1:Nt-1
    for i = 1: Nx-2
        b(i) = T(i,n) + const*T(i+1,n) + T(i+2,n);
    end
    T(2:Nx-1,n+1) = linsolve(A,b);
end

figure(1)
mesh(t,x,T), xlabel('t(s)'), ylabel('x(cm)'), zlabel('T(ºC)')
figure(2)
contourf(x,t,T'), xlabel('x(cm)'), ylabel('t(s)')  % T' transposta de T para dar igual as soluções
h = colorbar; % barra de cores
set(get(h,'label'),'string','T(ºC)'); % label da barra de cores

%% b)
clear all
close all
clc

L = 50; %cm
k = 0.93; % cal/(s*cm*ºC) -> condutividade da barra
p = 8.9; %g/cm^3 -> densidade

t0 = 0; %s
tf = 250; %s
dt = 0.1; %s
t = t0:dt:tf; %s
Nt = length(t);

x0 = 0; %cm
xf = L; %cm
dx = 0.5; %cm
x = x0:dx:xf; %cm
Nx = length(x);

c = nan(Nx,1); %cal/(g*ºC) -> calor especifico
c(:,1) = 0.094; % 1 coluna com 0.094
c(floor(Nx/2):end) = 0.188; % muda desde Nx/2 até ao fim para 0.188

T = zeros(Nx, Nt);
% Condições iniciais
T(:,1) = 50*sin(2*pi*x/L); % muda os valores da primeira coluna para todas as linhas
% Condições fronteira
T(1,:) = 0;
T(Nx,:) = 0;

% Implementação das matrizes com base no Slide 10 - Aula 6
eta = k .* dt ./ (c .* p .*dx^2);

A1 = diag((2./eta(2:end-1) + 2)); % (2/n + 2) na diagonal

A2 = diag(-ones([1 Nx-3]),1); % sobe um posição relativamente a diagonal
A3 = diag(-ones([1 Nx-3]),-1);% desce uma posição relativamente a diagonal


A = A1 + A2 + A3; % matriz no slide
b = zeros(Nx-2,1);


for n = 1:Nt-1
    for i = 2: Nx-1
        b(i-1) = T(i-1,n)+( 2/eta(i) - 2)*T(i,n) + T(i+1,n);
    end
    T(2:Nx-1,n+1) = linsolve(A,b);
end

figure(1)
mesh(t,x,T), xlabel('t(s)'), ylabel('x(cm)'), zlabel('T(ºC)')
figure(2)
contourf(x,t,T'), xlabel('x(cm)'), ylabel('t(s)')  % T' transposta de T para dar igual as soluções
h = colorbar; % barra de cores
set(get(h,'label'),'string','T(ºC)'); % label da barra de cores